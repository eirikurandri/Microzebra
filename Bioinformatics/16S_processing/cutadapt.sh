#trim reads using cutadapt:
samples=$(cat samples.txt)
for sample in $samples
do
cutadapt -e 0.15 -g ^CCTANGGGNNGCANCAG -G ^GGACTACNNGGGTATCTAAT \
--discard-untrimmed \
-o ${sample}_1a_trimmed.fq -p ${sample}_2a_trimmed.fq \
demultiplexed-${sample}_R1.fastq.gz demultiplexed-${sample}_R2.fastq.gz \
>> cutadapt_primer_trimming_stats.txt 2>&1
echo "On sample: $sample"
cutadapt -e 0.15 -g ^GGACTACNNGGGTATCTAAT -G ^CCTANGGGNNGCASCAG \
--discard-untrimmed \
-o ${sample}_1b_trimmed.fq -p ${sample}_2b_trimmed.fq \
demultiplexed-${sample}_R1.fastq.gz demultiplexed-${sample}_R2.fastq.gz \
>> cutadapt_primer_trimming_stats.txt 2>&1

#Merge both files
cat ${sample}_1a_trimmed.fq ${sample}_2b_trimmed.fq > ${sample}_1_trimmed.fq
cat ${sample}_2a_trimmed.fq ${sample}_1b_trimmed.fq > ${sample}_2_trimmed.fq
rm ${sample}_1a_trimmed.fq ${sample}_2b_trimmed.fq ${sample}_2a_trimmed.fq ${sample}_1b_trimmed.fq
done
