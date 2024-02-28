samples=$(ls SL_mimam_ALL_cutadapt_out/ | cut -dR -f1 | sort | uniq)
for s in $samples;
do
  	pear -f SL_mimam_ALL_cutadapt_out/${s}R1_001.fastq -r SL_mimam_ALL_cutadapt_out/${s}R2_001.fastq \
        -o SL_mimam_ALL_pear_out/$s -q 26 -v 20;
done
