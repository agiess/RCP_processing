#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(ggplot2)
library(dplyr)
library(gridExtra)
library(viridis)

Palette1 <- c('skyblue4', 'orange')

outdir <- args[1]
sample_name <- args[2]

#coverage plots: Leader, CDS, Traier plots scaled to 100nt
coverage <- read.csv(args[3], header=T)

coverage$Feature  <- factor(coverage$Feature, levels=c("leader", "cds", "trailer"), labels=c("leader", "CDS", "trailer"))
coverage$Fraction <- factor(coverage$Fraction, levels=c("SSU", "LSU"), labels = c("SSU", "LSU"))

coverage.sum <- coverage %>% group_by(Gene, Fraction) %>% mutate(gene_sum = sum(Count)) #seperate
coverage.sum.43.non.zero <- filter(coverage.sum, Fraction == "SSU", gene_sum > 0) #select genes that have at least one 43S TCP-seq count 
unique.genes=unique(coverage.sum.43.non.zero$Gene)
coverage.sum.subset <- coverage.sum %>%  filter (Gene %in% unique.genes)  #subset the 43S/80S matrix to those genes with at least one 43S read

coverage.sum.subset.0_100.zscore     <- coverage.sum.subset                  %>% ungroup() %>% group_by(Gene, Fraction)             %>% mutate(windowSD=sd(Count), windowMean=mean(Count)) %>% mutate(zscore=(Count-windowMean)/windowSD)
coverage.sum.subset.0_100.zscore.sum <- coverage.sum.subset.0_100.zscore     %>% ungroup() %>% group_by(Fraction, Positon, Feature) %>% summarise(zscore_sum=sum(zscore, na.rm = T), zscore_mean=mean(zscore, na.rm = T))
coverage.sum.subset.0_100.zscore.sum <- coverage.sum.subset.0_100.zscore.sum %>% ungroup() %>% group_by(Fraction)                   %>% mutate(fraction_min=min(zscore_mean))

plot_01 <- ggplot(data=as.data.frame(coverage.sum.subset.0_100.zscore.sum), aes(x=Positon, ymax=zscore_mean, ymin=fraction_min, y=zscore_mean, colour = as.factor(Fraction))) +
  geom_ribbon(stat="identity", position = "identity", aes(fill= as.factor(Fraction), alpha=0.5)) +
  geom_line() +
  theme_bw() +   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values=Palette1) +
  scale_color_manual(values=Palette1) +
  ggtitle(paste0("0-10 Genes n=", length(unique(coverage.sum.subset$Gene))) ) +
  xlab("Scaled position in transcript") + ylab("Zscore mean over transcript") +
  theme(legend.position="none") +
  facet_grid(Fraction ~ Feature, scales = "free")

ggsave(paste(outdir, sample_name, "_ranked.pdf", sep=""), plot_01, width=200, height=150, unit="mm", dpi=100, limitsize = FALSE)
ggsave(paste(outdir, sample_name, "_ranked.png", sep=""), plot_01, width=150, height=100, unit="mm", dpi=100)
