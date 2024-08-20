setwd("~/Documents/R")
library(ggpubr)
library(tidyverse)
library(gridExtra)
library(ggplot2)

cbPalette = c("Pakistan" = "black","Other_species" = "#E51932","Null" = "#FF99BF",
              "SF_zoo" = "#19B2FF")

mPCR_runStats <- read.csv('mPCR_OffTargetStats_updatedNulls.csv',head=T)

mPCR_runStats <- data.frame(mPCR_runStats)

het_practice <- het_practice %>%
  mutate( Location=factor(Location,levels=c("captive", "India", "Pakistan", "Afghanistan", "Tajikistan", "Kyrgyzstan", "Russia", "Mongolia (unknown)", "Mongolia (NW)", "Mongolia (SW)", "Mongolia (South)")) )

p1<-ggplot(mPCR_runStats, aes(x=group, y=Total_raw_reads_flagstat, color=group)) + 
  geom_boxplot() +
  geom_jitter(height = 0, width = 0.1, alpha = 0.5) +
  scale_color_manual(name = "group", values = cbPalette) +
  labs(y="Total_raw_reads_flagstat") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 13), 
        axis.text.y = element_text(size = 10, hjust = 1, vjust = 1), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 13),
        strip.text = element_text(size = 13),
        panel.grid.major.y = element_line(color = "grey90"),
        legend.position="none",
        plot.margin = margin(.5, .5, .5, .5, "cm")) 
p1

library(psych)
describeBy(mPCR_runStats$Total_raw_reads_flagstat, mPCR_runStats$group)

p1<-ggplot(mPCR_runStats, aes(x=group, y=Total_reads_after_MQ0_flagstat, color=group)) + 
  geom_boxplot() +
  geom_jitter(height = 0, width = 0.1, alpha = 0.5) +
  scale_color_manual(name = "group", values = cbPalette) +
  labs(y="Total_reads_after_MQ0_flagstat") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 13), 
        axis.text.y = element_text(size = 10, hjust = 1, vjust = 1), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 13),
        strip.text = element_text(size = 13),
        panel.grid.major.y = element_line(color = "grey90"),
        legend.position="none",
        plot.margin = margin(.5, .5, .5, .5, "cm")) 
p1

describeBy(mPCR_runStats$Total_reads_after_MQ0_flagstat, mPCR_runStats$group)

####Look at percent of reads read remaining after filtering for MQ0
p1<-ggplot(mPCR_runStats, aes(x=group, y=Percent_tot_wMQover0, color=group)) + 
  geom_boxplot() +
  geom_jitter(height = 0, width = 0.1, alpha = 0.5) +
  scale_color_manual(name = "group", values = cbPalette) +
  labs(y="Percent_tot_wMQover0") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 13), 
        axis.text.y = element_text(size = 10, hjust = 1, vjust = 1), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 13),
        strip.text = element_text(size = 13),
        panel.grid.major.y = element_line(color = "grey90"),
        legend.position="none",
        plot.margin = margin(.5, .5, .5, .5, "cm")) 
p1

describeBy(mPCR_runStats$Percent_tot_wMQover0, mPCR_runStats$group)

#####total number of reads after MQ0 filter and on Amplicon (+/- 50 bp) only

p1<-ggplot(mPCR_runStats, aes(x=group, y=Total_reads_after_MQ0_flagstat_onAMPonly, color=group)) + 
  geom_boxplot() +
  geom_jitter(height = 0, width = 0.1, alpha = 0.5) +
  scale_color_manual(name = "group", values = cbPalette) +
  labs(y="Total_reads_after_MQ0_flagstat_onAMPonly") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 13), 
        axis.text.y = element_text(size = 10, hjust = 1, vjust = 1), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 13),
        strip.text = element_text(size = 13),
        panel.grid.major.y = element_line(color = "grey90"),
        legend.position="none",
        plot.margin = margin(.5, .5, .5, .5, "cm")) 
p1

describeBy(mPCR_runStats$Total_reads_after_MQ0_flagstat_onAMPonly, mPCR_runStats$group)

#####% total reads after MQ0 filter that map to Amplicon (+/- 50 bp) only

p1<-ggplot(mPCR_runStats, aes(x=group, y=Percent_readsONamp_afterMQ0filter, color=group)) + 
  geom_boxplot() +
  geom_jitter(height = 0, width = 0.1, alpha = 0.5) +
  scale_color_manual(name = "group", values = cbPalette) +
  labs(y="Percent_readsONamp_afterMQ0filter") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 13), 
        axis.text.y = element_text(size = 10, hjust = 1, vjust = 1), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 13),
        strip.text = element_text(size = 13),
        panel.grid.major.y = element_line(color = "grey90"),
        legend.position="none",
        plot.margin = margin(.5, .5, .5, .5, "cm")) 
p1

describeBy(mPCR_runStats$Percent_readsONamp_afterMQ0filter, mPCR_runStats$group)





p1<-ggplot(mPCR_runStats, aes(x=group, y=OFF_AMPLICON_BASES, color=group)) + 
  geom_boxplot() +
  geom_jitter(height = 0, width = 0.1, alpha = 0.5) +
  scale_color_manual(name = "group", values = cbPalette) +
  labs(y="OFF_AMPLICON_BASES") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 13), 
        axis.text.y = element_text(size = 10, hjust = 1, vjust = 1), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 13),
        strip.text = element_text(size = 13),
        panel.grid.major.y = element_line(color = "grey90"),
        legend.position="none",
        plot.margin = margin(.5, .5, .5, .5, "cm")) 
p1

p1<-ggplot(mPCR_runStats, aes(x=group, y=ON_AMPLICON_BASES, color=group)) + 
  geom_boxplot() +
  geom_jitter(height = 0, width = 0.1, alpha = 0.5) +
  scale_color_manual(name = "group", values = cbPalette) +
  labs(y="ON_AMPLICON_BASES") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 13), 
        axis.text.y = element_text(size = 10, hjust = 1, vjust = 1), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 13),
        strip.text = element_text(size = 13),
        panel.grid.major.y = element_line(color = "grey90"),
        legend.position="none",
        plot.margin = margin(.5, .5, .5, .5, "cm")) 
p1


p1<-ggplot(mPCR_runStats, aes(x=group, y=PCT_OFF_AMPLICON, color=group)) + 
  geom_boxplot() +
  geom_jitter(height = 0, width = 0.1, alpha = 0.5) +
  scale_color_manual(name = "group", values = cbPalette) +
  labs(y="PCT_OFF_AMPLICON") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 13), 
        axis.text.y = element_text(size = 10, hjust = 1, vjust = 1), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 13),
        strip.text = element_text(size = 13),
        panel.grid.major.y = element_line(color = "grey90"),
        legend.position="none",
        plot.margin = margin(.5, .5, .5, .5, "cm")) 
p1


p1<-ggplot(mPCR_runStats, aes(x=group, y=MEAN_AMPLICON_COVERAGE, color=group)) + 
  geom_boxplot() +
  geom_jitter(height = 0, width = 0.1, alpha = 0.5) +
  scale_color_manual(name = "group", values = cbPalette) +
  labs(y="MEAN_AMPLICON_COVERAGE") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 13), 
        axis.text.y = element_text(size = 10, hjust = 1, vjust = 1), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 13),
        strip.text = element_text(size = 13),
        panel.grid.major.y = element_line(color = "grey90"),
        legend.position="none",
        plot.margin = margin(.5, .5, .5, .5, "cm")) 
p1
