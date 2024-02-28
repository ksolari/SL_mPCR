samples=$(ls FastqFiles/ | cut -dR -f1 | sort | uniq)
for s in $samples;
do
  	cutadapt -g ^NNNNNNGGGTTGGTAAATTTCGTGCCAGC -G ^NNNNNNCATAGTGGGGTATCTAATCCCAGTTTG \
        -o SL_mimam_ALL_cutadapt_out/${s}R1_001.fastq -p SL_mimam_ALL_cutadapt_out/${s}R2_001.fastq --discard-untrimmed \
        FastqFiles/${s}R1_001.fastq.gz FastqFiles/${s}R2_001.fastq.gz ;
done
