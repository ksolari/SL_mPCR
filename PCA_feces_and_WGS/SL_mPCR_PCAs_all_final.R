library(ggpubr)
library(tidyverse)
library(gridExtra)
setwd("~/Documents/R")

cbPalette = c("Mongolia (South)" = "#FFBF7F","Mongolia (SW)" = "#FF7F00","Mongolia (NW)" = "#B2FF8C","Mongolia (unknown)" = "#32FF00", 
              "Russia" = "#A5EDFF", "captive" = "black",  "Kyrgyzstan" = "#19B2FF", "Tadjikistan"="#CCBFFF",  
              "Afghanistan" = "#654CFF","Pakistan" = "#FF99BF", "India" = "#E51932" , "SF_zoo_fecal" = "grey", "Pakistan_fecal" = "gold" )

cbPaletteB = c("Mongolia (South)" = "#FFBF7F","Mongolia (SW)" = "#FF7F00","Mongolia (NW)" = "#B2FF8C","Mongolia (unknown)" = "#32FF00", 
               "Russia" = "#A5EDFF", "captive" = "black",  "Kyrgyzstan" = "#19B2FF", "Tadjikistan"="#CCBFFF",  
               "Afghanistan" = "#654CFF","Pakistan" = "#FF99BF", "India" = "red3" , "SF_zoo_fecal" = "black", "Pakistan_fecal" = "#FF99BF" )


######PCA with all SNPs - n33 - WGS without relatives
####/oak/stanford/groups/dpetrov/ksolari/SL_relativesRemoved/CombinedAnnotations_wHeaderD_filtered_pass_1nonref_90maxmis_IndRename_no10xU01U08AF06KGf4_1nonref.recode.vcf

eigenvec <- read.table('CombinedAnnotations_wHeaderD_filtered_pass_1nonref_90maxmis_IndRename_no10xU01U08AF06KGf4_1nonref_pca.eigenvec',head=F)
eigenval <- read.table('CombinedAnnotations_wHeaderD_filtered_pass_1nonref_90maxmis_IndRename_no10xU01U08AF06KGf4_1nonref_pca.eigenval',head=F)
id <- read.csv('SampleInfo_SL_no02_12_20_10x_U01U08AF06KGf4_multMong.csv', head = T, sep = ',')

eigenval$V2 <- eigenval$V1/(sum(eigenval$V1))*100
eigenval$V2

eigenvec_rename <- rename(eigenvec, c("V2"="Individual")) 

eigenvec_w_sub <- merge(eigenvec_rename, id,by="Individual")

PCall <-ggplot(eigenvec_w_sub, aes(x=V3, y=V4, color=Location)) +
  scale_color_manual(name="Location", values=cbPalette) +
  geom_point(size=3, alpha = 3/4) +
  xlab("PC1 (7.69%)") + ylab("PC2 (6.53%)") +
  theme(legend.position = "right", legend.title = element_blank()) +
  theme_bw()
PCall

######PCA with 1000 (980) SNPs - n33 - WGS without relatives
####/oak/stanford/groups/dpetrov/ksolari/SL_Gencove_redownload_8_2023/GencoveSLdata/bam/completed/SL_WGS_n33_AD_DP_1000snps_1bp_rename_biallelic_noRelatives.recode.vcf 

eigenvec <- read.table('SL_WGS_n33_AD_DP_1000snps_1bp_rename_biallelic_noRelatives.eigenvec',head=F)
eigenval <- read.table('SL_WGS_n33_AD_DP_1000snps_1bp_rename_biallelic_noRelatives.eigenval',head=F)
id <- read.csv('SampleInfo_SL_no02_12_20_10x_U01U08AF06KGf4_multMong.csv', head = T, sep = ',')

eigenval$V2 <- eigenval$V1/(sum(eigenval$V1))*100
eigenval$V2

eigenvec_rename <- rename(eigenvec, c("V2"="Individual")) 

eigenvec_w_sub <- merge(eigenvec_rename, id,by="Individual")

PC1000 <-ggplot(eigenvec_w_sub, aes(x=V3, y=V4, color=Location)) +
  scale_color_manual(name="Location", values=cbPalette) +
  geom_point(size=3, alpha = 3/4) +
  xlab("PC1 (17.1%)") + ylab("PC2 (9.25%)") +
  theme(legend.position = "right", legend.title = element_blank()) +
  theme_bw()
PC1000


######PCA with 144 SNPs - n33 - WGS without relatives
####oak/stanford/groups/dpetrov/ksolari/SL_Gencove_redownload_8_2023/GencoveSLdata/bam/completed/SL_WGS_n33_AD_DP_144snps_1bp_rename_biallelic_noRelatives.recode.vcf

eigenvec <- read.table('SL_WGS_n33_AD_DP_144snps_1bp_rename_biallelic_noRelatives.eigenvec',head=F)
eigenval <- read.table('SL_WGS_n33_AD_DP_144snps_1bp_rename_biallelic_noRelatives.eigenval',head=F)
id <- read.csv('SampleInfo_SL_no02_12_20_10x_U01U08AF06KGf4_multMong.csv', head = T, sep = ',')

eigenval$V2 <- eigenval$V1/(sum(eigenval$V1))*100
eigenval$V2

eigenvec_rename <- rename(eigenvec, c("V2"="Individual")) 

eigenvec_w_sub <- merge(eigenvec_rename, id,by="Individual")

PC144 <-ggplot(eigenvec_w_sub, aes(x=V3, y=V4, color=Location)) +
  scale_color_manual(name="Location", values=cbPalette) +
  geom_point(size=3, alpha = 3/4) +
  xlab("PC1 (18.5%)") + ylab("PC2 (9.12%)") +
  theme(legend.position = "right", legend.title = element_blank()) +
  theme_bw()
PC144

######PCA with 144 SNPs - n33 - WGS without relatives + fecal - 5 SFzoo and 47 PK
#####/oak/stanford/groups/dpetrov/ksolari/SL_mPCR_AllPakistan_AllZoo/All_fastq/143SNP_WGSn33_mPCR_n5_47_merge.vcf


eigenvec <- read.table('143SNP_WGSn33_mPCR_n5_47_merge.eigenvec',head=F)
eigenval <- read.table('143SNP_WGSn33_mPCR_n5_47_merge.eigenval',head=F)
id <- read.csv('mPCR_WGS_fecal_143snps_PCA_sampleInfo_5SF47PK_33WGS.csv', head = T, sep = ',')

eigenval$V2 <- eigenval$V1/(sum(eigenval$V1))*100
eigenval$V2

eigenvec_rename <- rename(eigenvec, c("V2"="Individual")) 

eigenvec_w_sub <- merge(eigenvec_rename, id,by="Individual")

PC144wFecal <-ggplot(eigenvec_w_sub, aes(x=V3, y=V4, group=Location)) +
  geom_point(aes(shape=Location, color=Location, size=Location, alpha = 0.05))+
  scale_color_manual(name="Location", values=cbPaletteB) +
  scale_shape_manual(values=c(19, 19, 19,19, 19, 19,19, 19, 19,8, 19, 8,19))+
  scale_size_manual(values=c(3,3,3,3,3,3,3,3,3,3,3,3,3 ))+
  xlab("PC1 (16.3%)") + ylab("PC2 (5.8%)") +
  theme(legend.position = "right", legend.title = element_blank()) +
  theme_bw()
PC144wFecal


legend <- get_legend(PC144wFecal)
PCall <- PCall + theme(legend.position="none") 
PC1000 <- PC1000 + theme(legend.position="none") 
PC144 <- PC144 + theme(legend.position="none") 
PC144wFecal <- PC144wFecal + theme(legend.position="none") 
grid.arrange(PCall, PC1000, PC144, PC144wFecal, ncol=4, nrow =1, widths = c(2.5, 2.5, 2.5, 2.5))

