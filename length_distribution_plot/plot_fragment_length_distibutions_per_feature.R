#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(ggplot2)
library(dplyr)
library(gridExtra)
library(reshape2)
library(cowplot) 

outdir <- args[1]
sample_name <- args[2]

SSU_TSS  <- read.csv(args[3], header=TRUE)
LSU_TSS  <- read.csv(args[4], header=TRUE)

SSU_TSS$fraction<-"SSU"
LSU_TSS$fraction<-"LSU"

run_combined<-rbind(SSU_TSS, LSU_TSS)
run_combined$feature_type <- factor(run_combined$feature_type, levels = c("start_of_transcript", "leader5", "start_codon", "CDS", "stop_codon", "trailer3", "end_of_transcript"))
run_combined$fraction <- factor(run_combined$fraction, levels = c("SSU","LSU"))
#un_combined <- run_combined %>% filter(X.read_length >= 17, X.read_length <= 65) %>% group_by(fraction, feature_type) %>% mutate(proportion_read_count = count/sum(count))
run_combined <- run_combined %>% filter(X.read_length >= 17) %>% group_by(fraction, feature_type) %>% mutate(proportion_read_count = count/sum(count))

run_density_plot <- ggplot(data=run_combined, aes(x=X.read_length, y=proportion_read_count, fill=fraction)) + 
  #geom_vline(xintercept = 27, colour="grey", linetype = "longdash") + 
  #geom_vline(xintercept = 40, colour="grey", linetype = "longdash") + 
  geom_vline(xintercept = 25, colour="lightblue", linetype = "dotted") +
  geom_vline(xintercept = 35, colour="lightblue", linetype = "dotted") +
  geom_density(stat="identity") +
  facet_grid(feature_type ~ fraction) + 
  scale_x_continuous(expand = c(0, 0), limits=c(15,75), breaks=c(20,30,40,50,60,70)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  xlab("Read length (nt)") + 
  ylab("Proportional counts") + 
  theme(legend.position="none")

ggsave(paste(outdir, sample_name, ".pdf", sep=""), run_density_plot, width = 10, height = 12)
ggsave(paste(outdir, sample_name, ".png", sep=""), run_density_plot, width=200, height=250, unit="mm", dpi=100)
