
#Note this was done for all 5 co-assemblies
#create contigs database
contig_list=$(cat ../'text_file_of_assembly_names.txt')
echo "$contig_list"
CON=/path/to/contigs/fasta
for a in $contig_list
  do 
    anvi-script-reformat-fasta $CON/"$a"_final.contigs.fa -o $CON/"$a".contigs-fixed.fa -l 1000 --simplify-names
    anvi-gen-contigs-database -f $CON/"$a".contigs-fixed.fa -o $CON/"$a".CONTIGS.db -n 'CONTIG DB of ZEBRAFISH SAMPLES'
    echo -e "yay! "$a" DB is generated! Ready to annotate "$a" contig database"
    echo -e "running hmms with "$a" yay!"
    anvi-run-hmms -c $CON/"$a".CONTIGS.db --num-threads 10
    echo -e "doing cogs with "$a" yay!"
    anvi-run-ncbi-cogs -c $CON/"$a".CONTIGS.db --num-threads 10 --cog-data-dir /projects/mjolnir1/data/databases/anvio/ncbi-cogs/20221104/
    echo -e "doing pfams with "$a" yay!"
    anvi-run-pfams -c $CON/"$a".CONTIGS.db --num-threads 10 --pfam-data-dir /projects/mjolnir1/data/databases/anvio/pfams/20221104/
    echo -e "doing kegg with "$a" yay!"
    anvi-run-kegg-kofams -c $CON/"$a".CONTIGS.db --num-threads 10 --kegg-data-dir /projects/mjolnir1/data/databases/anvio/kegg-kofams/20221104/
    echo -e "yay! annotations are done for "$a", please proceed to profiling"
done

#map metageomic reads back to assembly prior to profiling
FASTA_DIR=/path/to/unmapped/reads
ASSEMBLY=/path/to/assembly/fasta
WORK_DIR=/path/to/work/directory
sample_list=$(cat $FASTA_DIR/'list_of_sample_names.txt')
# Check the list for sanity
echo $sample_list
# run mapping
for stuff in $sample_list
do
minimap2 -ax sr $ASSEMBLY $FASTA_DIR/"$stuff"_unmapped_1.fastq.gz $FASTA_DIR/"$stuff"_unmapped_2.fastq.gz|samtools view -F 4 -bS > $WORK_DIR/"$stuff"_mapping.bam
done

#now generating profiles for each sample
MAP=/path/to/bam/mapping/files
CON=/path/to/contigs/database
sample_list=$(cat 'list_of_samples.txt')
for b in $sample_list
do
anvi-init-bam $MAP/"$b"_mapping.bam -o $MAP/"$b".bam
anvi-profile -i $MAP/"$b".bam -c $CON/W.CONTIGS.db -o $CON/profiles/"$b" --cluster-contigs --profile-SCVs --min-contig-length 1000 -T 10
done

#merge the profiles
PROFILES=/path/to/profiles
CONTIGS=/path/to/contigs
anvi-merge $PROFILES/*/PROFILE.db -o $PROFILES/MERGED_Profile --enforce-hierarchical-clustering -c $CONTIGS/CONTIGS.db

## Assign Taxonomy with KAIJU
anvi-get-sequences-for-gene-calls -c $CONTIGS/CONTIGS.db -o gene_calls.fa

kaiju_path='path/to/kaiju/2020-05-25_nr_euk'

kaiju -t $kaiju_path/nodes.dmp \
      -f $kaiju_path/kaiju_db_nr_euk.fmi \
      -i gene_calls.fa \
      -o gene_calls_nr.out \
      -z 10 \
      -v
#
addTaxonNames -t $kaiju_path/nodes.dmp \
              -n $kaiju_path/names.dmp \
              -i gene_calls_nr.out \
              -o gene_calls_nr.names \
              -r superkingdom,phylum,order,class,family,genus,species
#
anvi-import-taxonomy-for-genes -i gene_calls_nr.names \
                                 -c $CONTIGS/CONTIGS.db \
                                 -p kaiju \
                                 --just-do-it

#
# Use CONCOCT for binning contigs
anvi-cluster-contigs -p $PROFILES/MERGED_Profile/PROFILE.db -c $CONTIGS/CONTIGS.db -C CONCOCT --driver concoct --just-do-it

# Call MAGs
anvi-rename-bins -c $CONTIGS/CONTIGS.db \
               -p $PROFILES/PROFILE.db \
               --prefix SZ \
               --collection-to-read CONCOCT \
               --collection-to-write MAGs \
               --report-file rename.txt \
               --min-completion-for-MAG 50 \
               --max-redundancy-for-MAG 10 \
               --call-MAGs
