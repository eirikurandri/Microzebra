#FastqC for checking quality

find . -name *.fq.gz > sample_list.txt
for a in $sample_list
  do
    fastqc -o ../fastqc/fastqc_init/ "$a"
done
#then seperate directories were made for the fastqc reports for the forward and reverse reads respectively
#And the following command run to create a MultiQC report for the forward and reverse reads

multiqc forward/ -o multiqc_forward
multiqc reverse/ -o multiqc_reverse

#Aadapter and quality trimming 

sed -e 's/_1.fq.gz//g' -e 's/_2.fq.gz//g' sample_list.txt > sample_list_temp.txt
uniq sample_list_temp.txt > sample_list
sample_list=$(cat ../sample_list)
echo sample_list
#Do the quality filtering and trimming adaptors using AdapterRemoval
for a in $sample_list
  do
    AdapterRemoval --file1 "$a"_1.fq.gz --file2 "$a"_2.fq.gz --basename "$a" --output1 "$a"_filtered_1.fq.gz --output2 "$a"_filtered_2.fq.gz --adapter1 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT --adapter2 GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGATGACTATCTCGTATGCCGTCTTCTGCTTG --trimqualities --trimns --minlength 50 --threads 20 --gzip
done

#We then ran fastqc and MultiQC again as a sanity-check and for the sake of happiness (same as above just different filenames)

#We then removed duplicate sequences using sequit´s rmdup 

for a in $sample_list
do
zcat "$a"_filtered_1.fq.gz| seqkit rmdup -j 20 -s -o "$a"_1_clean.fa.gz
zcat "$a"_filtered_2.fq.gz| seqkit rmdup -j 20 -s -o "$a"_2_clean.fa.gz
done


#This was followed by re-paring the reads with bbmaps repair.sh script

for a in $sample_list
do
repair.sh in="$a"_1_clean.fq.gz in2="$a"_2_clean.fq.gz out="$a"_1_repaired.fq.gz out2="$a"_2_repaired.fq.gz outsingle="$a"_singletons.fq.gz overwrite=t
done

#Then we removed host, phix and human reads by mapping to the reference and keeping the unmapped reads in each step using minimap2
#This is the setup used for zebrafish, human and phix 
#(Just included once so its not extremely redundant, its all the same except "REF_DIR" and "Ref" are different between organimis)
REF_DIR='/path/to/reference'
Ref='name of reference file' 
#run stuff
for stuff in $sample_list
do
minimap2 -ax sr $REF_DIR/$Ref "$stuff"_1_repaired.fq.gz "$stuff"_2_repaired.fq.gz| samtools view -bS > "$stuff"_mapped_and_unmapped.bam
  samtools view -b -f 12 "$stuff"_mapped_and_unmapped.bam|samtools sort -n >  "$stuff"_sorted.bam
    bedtools bamtofastq -i "$stuff"_sorted.bam -fq "$stuff"_nohost_1.fastq -fq2 "$stuff"_nohost_2.fastq
      pigz -p 40 "$stuff"_nohost_*.fastq
done

#We also got the mean coverage accross the zebrafish genome by using samtools sort and samtools depth on the ""$stuff"_mapped_and_unmapped.bam" files
#Sort part
BAM='path/to/bam/' 
for i in $sample_list
do
samtools sort -o $BAM/sorted/"$i".sorted.bam $BAM/"$i"mapped_and_unmapped.bam
done


#depth part to output the average coverage of each sample and sample name into a file
for i in $samples
do
echo -n "$i     " >> nr.txt
samtools depth -a $BAM/sorted/"$i".sorted.bam  |  awk '{sum+=$3} END { print "Average = ",sum/NR}' >>nr.txt
done

#After getting out the unmapped reads we then created 5 co-assemblies using Megahit 
#this is just a script for generating one co-assembly
FASTA_DIR='path/to/appropriate/unmapped_fasta/files'
WORK_DIR='path/to/appropriate/work/dir'
R1=$(ls $FASTA_DIR/*1.fastq.gz | python -c 'import sys; print(",".join([x.strip() for x in sys.stdin.readlines()]))')
R2=$(ls $FASTA_DIR/*2.fastq.gz | python -c 'import sys; print(",".join([x.strip() for x in sys.stdin.readlines()]))')
# Check the FASTA list
echo $R1
echo $R2
# run megahit
megahit -1 $R1 -2 $R2 --min-contig-len 1000 -t 20 --presets meta-sensitive -o $WORK_DIR

#the assemblies were then evaluted and QC´d using Quast
assembly_list=(cat assembly_list.txt) #just a text file with my assembly names
for a in $assembly_list
  do
    quast.py -o Quast_Report/ "$a"_final.contigs.fa
done

#now we moved to the anvio part of the pipeline which can be found in the file anvio.pipeline
