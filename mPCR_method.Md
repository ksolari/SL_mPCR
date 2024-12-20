# mPCR SNP Panel Protocol Summary 

__Katherine A. Solari, Shakeel Ahmad, Ellie E. Armstrong, Michael G. Campana, Hussain Ali, Shoaib Hameed, Jami Ullah, Barkat Ullah Khan, Muhammad A. Nawaz, Dmitri A. Petrov, 2024__  
Stanford University
Quaid-i-Azam University
University of California, Riverside
Smithsonian’s National Zoo and Conservation Biology Institute
Snow Leopard Foundation
Qatar University
Chan Zuckerberg BioHub

# SNP panel design
__1. Identify target SNPs__

Target SNPs will depend on the questions at hand and can be chosen using the mPCRselect pipeline [1] (https://github.com/ellieearmstrong/mPCRselect)

__2. Design primers for target SNPs__ 

It is important to identify primers that don't negatively interact when combined in one tube. One program that can be used to design primers for mulitplex PCR is NGS-PrimerPlex [2] (https://github.com/aakechin/NGS-PrimerPlex)
NGS-PrimerPlex will only be able to design primers that can be run in multiplex for a small subset of the target SNPs that you input. In this study, we input 1000 SNPs and it was able to design a primer pool that included 144 of our target loci. So, be sure to input more target loci than you want to end up with. 
Another note on NGS-PrimerPlex, the first few times you run it, it will likely error saying that it can not design primers for specific loci. Remove these loci from your input file and rerun the program. 
Example NGS-PrimerPlex submission:
`singularity exec ngs-primerplex_native.sif python3 /NGS-PrimerPlex/NGS_primerplex.py -regions TargetSNPs.txt -ref reference_genome.fasta -ad1 TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -ad2 GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -run Output_Name -minampllen 26 -optampllen 50 -maxampllen 120 -maxprimerlen 26 -minprimerlen 18 -blast -skip -primernum1 25`

# SNP panel output processing
__1. Map raw fastq files to reference genome, sort, and index__

- install bwa [3] (https://github.com/lh3/bwa)
- install samtools [4] (https://www.htslib.org/download/)
- map, sort, and index each run:
`bwa mem reference_genome.fasta sample1_R1_001.fastq.gz sample1_R2_001.fastq.gz | samtools view -bS - | samtools sort - > sample1_sorted.bam`
`samtools index sample1_sorted.bam`

__2. Calculate number of reads (optional QC step)__

`samtools flagstat sample1_sorted.bam`

__3. Filter reads with Mapping Quality of 0__

- These are low quality reads that we do not want to include in our analyses
`samtools view -h -q 1 sample1_sorted.bam -o sample1_sorted_MQ1.bam`
`samtools index sample1_sorted_MQ1.bam`

__4. Calculate number of reads after filtering based on mapping quality (optional QC step)__

`samtools flagstat sample1_sorted_MQ1.bam`

__5. Filter bam files to only include the target amplicon regions__

- Install bedtools [5] (https://bedtools.readthedocs.io/en/latest/content/installation.html)
`bedtools intersect -a sample1_sorted_MQ1.bam -b amplicons.bed > sample1_sorted_MQ1_AmpliconOnly.bam`
- The amplicon bed file should have no header and three columns: 1- Chromsome name, 2- Amplicon start coordinate, 3- Amplicon stop coordinate 
- The information necessary for the bed file can be obtained from the NGS-PrimerPlex output
- example bed file provided from our dataset - `144SNPs_1bp_FINAL.bed`

__6. Calculate the number of reads with mapping quality greater than 0 and mapping to the target amplicons (optional QC step)__

`samtools flagstat sample1_sorted_MQ1_AmpliconOnly.bam`

__7. Call genotypes for target SNPs__

- Make a text file that lists out all of the bam files that you want to call SNPs in:
e.g. - `ls *_MQ1.bam > SortedMq1BamList.txt`
- Install bcftools [6] (https://www.htslib.org/download/)
- Call SNPs using bcftools:
`bcftools mpileup -A -a AD,DP -R SNPlist.bed -f reference_genome.fasta -b SortedMq1BamList.txt  | bcftools call -m -Oz -f GQ -o output.vcf`
  - annotate (-a) with allelic depth (AD) and number of high-quality bases (DP)
  - output genotype quality (GQ) with flag -f

__8. Filter based on read depth (DP) and mapping quality (MQ)__

Install vcftools [7] (https://vcftools.github.io/)
Here, we remove genotypes that have a DP less than 6 or a MQ less than 20
`vcftools --vcf output.vcf --out output_minDP6minGQ20 --minDP  6 --minGQ 20 --recode`
Example output files from our dataset: `PKSL_n604.recode.vcf` and `SFzoo_n30.recode.vcf`

__9. Calculate the number of SNPs successfully genotyped__

`vcftools --vcf output_minDP6minGQ20.recode.vcf --out output_minDP6minGQ20_missingness --missing-indv`
This outputs the file output_minDP6minGQ20_missingness.imiss which lists the number of missing SNPs for each sample in your vcf file

__10. Remove runs with less than 20 SNPs genotyped__ 

Use the output_minDP6minGQ20_missingness.imiss file to identify which runs have less than 20 SNPs (or whatever cut off you decide to use). Make a text file of all of the runs with less than 20 SNPs (one run name per line).
`vcftools --vcf output_minDP6minGQ20.recode.vcf --out output_minDP6minGQ20_min20SNPs --remove Samples_to_remove.txt --recode`

__11. Calculate the number of shared alleles between all pairwise comparisons to identify samples from the same individual__

Install plink version 1.90b5.3 [8] (https://www.cog-genomics.org/plink/)
`plink --vcf output_minDP6minGQ20_min20SNPs.recode.vcf --out output_minDP6minGQ20_min20SNPs_IBS --allow-extra-chr --double-id --genome full`
This outputs the file output_minDP6minGQ20_min20SNPs_IBS.genome. In this file, the column labeled "DST" is the proportion of alleles identical by state. In our project, we called all run pairs with a DST value of 0.95 or higher and with 20 or more SNPs being compared to be the same individual. 
The total number of SNPs being compared in each pairwise comparison can be calculated from this output file as IBS0+IBS1+IBS2.

__12. Make PCA of SNP calls to identify population assignment (Optional)__

If you are combining your new vcf with existing data, then zip and index vcf files to be merged. Make sure that all vcf files have the same SNP loci. 
`bgzip -c output_minDP6minGQ20.recode.vcf > output_minDP6minGQ20.vcf.gz`
`bgzip -c existing.vcf > existing.vcf.gz`
`bcftools index -t output_minDP6minGQ20.vcf.gz`
`bcftools index -t existing.vcf.gz`
Merge the files and unzip
`vcf-merge output_minDP6minGQ20.vcf.gz existing.vcf.gz | bgzip -c > merged.vcf.gz`
`gunzip merged.vcf.gz` 
Note that more than two vcf files can be merged
Example output file from our dataset - `143SNP_WGSn33_mPCR_n5_56_merge.vcf`
Run PCA analysis with plink
`plink2 --make-bed --vcf merged.vcf --allow-extra-chr --set-all-var-ids @:# --out merged`
`plink2 --vcf merged.vcf --freq  --allow-extra-chr --out merged`
`plink2 --bfile merged --read-freq merged.afreq --allow-extra-chr --pca 50 --out merged`
Change the number after --pca is be the number of samples that you have in your vcf file. 


__Additional optional analyses to run:__

- Calculate allele frequency of each SNP
`vcftools --vcf output_minDP6minGQ20.recode.vcf --out output_minDP6minGQ20_freq --freq`
- Calculate the heterozygosity of each sample (this can be used to identify non-target species in your sample set)
`bcftools stats output_minDP6minGQ20.recode.vcf -S SampleList.txt > SampleList_het.txt`



## References  
1. Armstrong, E. E., Li, C., Campana, M. G., Ferrari, T., Kelley, J. L., Petrov, D. A., Solari, K. A., & Mooney, J. A. (2024). Recommendations for Population and Individual Diagnostic SNP Selection in Non-Model Species. In bioRxiv (p. 2024.07.03.601943). https://doi.org/10.1101/2024.07.03.601943
2. Kechin, A., Borobova, V., Boyarskikh, U., Khrapov, E., Subbotin, S., & Filipenko, M. (2020). NGS-PrimerPlex: High-throughput primer design for multiplex polymerase chain reactions. PLoS Computational Biology, 16(12), e1008468.
3. Li, H. (2013). Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. In arXiv [q-bio.GN]. http://arxiv.org/abs/1303.3997
4. Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., Durbin, R., & 1000 Genome Project Data Processing Subgroup. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics , 25(16), 2078–2079.
5. Quinlan, A. R., & Hall, I. M. (2010). BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics , 26(6), 841–842.
6. Li, H. A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. Bioinformatics (2011) 27(21) 2987-93
7. Danecek, P., Auton, A., Abecasis, G., Albers, C. A., Banks, E., DePristo, M. A., Handsaker, R. E., Lunter, G., Marth, G. T., Sherry, S. T., McVean, G., Durbin, R., & 1000 Genomes Project Analysis Group. (2011). The variant call format and VCFtools. Bioinformatics , 27(15), 2156–2158.
8. Purcell, S., Neale, B., Todd-Brown, K., Thomas, L., Ferreira, M. A. R., Bender, D., Maller, J., Sklar, P., de Bakker, P. I. W., Daly, M. J., & Sham, P. C. (2007). PLINK: a tool set for whole-genome association and population-based linkage analyses. American Journal of Human Genetics, 81(3), 559–575.
