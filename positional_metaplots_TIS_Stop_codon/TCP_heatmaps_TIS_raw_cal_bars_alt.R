#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(ggplot2)
library(dplyr)
library(gridExtra)

##
# 1 unshifted heatmaps
##

start_m5<-read.table(args[1])
start_m3<-read.table(args[2])

mono.start5.raw<-read.csv(args[3], header=FALSE, check.names=FALSE, sep=",")
mono.start3.raw<-read.csv(args[4], header=FALSE, check.names=FALSE, sep=",")

start_m5_50 <- subset(start_m5, V1<=50 & V1 >= -100 & V2 < 69 & V2 >=20);
start_m3_50 <- subset(start_m3, V1<=100 & V1 >= -50 & V2 < 69 & V2 >=20);

#hack to fix plotting problem with NaN values and the heatmap legend
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

start_m3_50[is.nan(start_m3_50)] <- NA
start_m5_50[is.nan(start_m5_50)] <- NA

mLen5<-ggplot(subset(start_m5_50, V1 <= 50 & V1 >= -100 & V2 >=20) , aes(x=V1, y=V2, fill=V3)) + geom_tile()  +
               #scale_fill_gradientn(colours=c("yellow", "lightblue", "blue", "navy"), values=c(0,0.05,0.15,0.2,0.5,1), name="sum count") +   
               scale_fill_gradientn(colours=c("yellow","lightblue","blue", "navy"), values=c(0,30,60,90,120,650), rescaler = function(x,...) x, oob = identity, name="sum count") +             
               #scale_fill_gradientn(colours=c("yellow", "lightblue", "lightblue", "blue", "blue", "navy", "navy", "navy"), name="sum count") +
               #scale_fill_gradientn(colours=c("yellow", "lightblue", "blue", "navy"),breaks=c(100,200,300,400), name="sum count") +
               xlab("Position relative to start codon") + ylab("Protected fragment length") +
               scale_x_continuous(limits = c(-101,51), expand = c(0, 0), breaks=c(-100,-90,-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50)) +
               scale_y_continuous(expand = c(0, 0)) +
               theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
               theme(text = element_text(size = 12))

mLen3<-ggplot(subset(start_m3_50, V1 <= 100 & V1 >= -50 & V2 >=20) , aes(x=V1, y=V2, fill=V3)) + geom_tile()  +
               #scale_fill_gradientn(colours=c("yellow", "lightblue", "blue", "navy"), values=c(0,0.05,0.15,0.2,0.5,1), name="sum count") +    
               scale_fill_gradientn(colours=c("yellow","lightblue","blue", "navy"), values=c(0,30,60,90,120,200), rescaler = function(x,...) x, oob = identity, name="sum count") +
               #scale_fill_gradientn(colours=c("yellow", "lightblue", "blue", "blue", "navy", "navy"), name="sum count") +
               #scale_fill_gradientn(colours=c("yellow", "lightblue", "lightblue", "blue", "blue", "navy", "navy", "navy"), name="sum count") +
               xlab("Position relative to start codon") + ylab("Protected fragment length") +
               scale_x_continuous(limits = c(-51,101), expand = c(0, 0), breaks=c(-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80,90,100)) +
               scale_y_continuous(expand = c(0, 0)) +
               theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
               theme(text = element_text(size = 12))

##
# 2 unshifted barplots
##


#scale the data by the sum of the window
#subset_mono5.raw<-filter(mono.start5.raw, V1 >=-100, V1 <= 50)
#subset_mono5.raw <- mutate(subset_mono5.raw, proportion = V2/sum(V2))

#subset_mono3.raw<-filter(mono.start3.raw, V1 >=-50, V1 <= 100)
#subset_mono3.raw <- mutate(subset_mono3.raw, proportion = V2/sum(V2))

#max_y <- max(subset_mono3.raw$proportion,subset_mono5.raw$proportion)


#calculate the bars from the values in the heatmaps (to for selected readd lengths)

start_m5_50 <- subset(start_m5, V1<=50 & V1 >= -100 & V2 < 69 & V2 >=20);
start_m3_50 <- subset(start_m3, V1<=100 & V1 >= -50 & V2 < 69 & V2 >=20);

subset_mono3.bar <- start_m3_50 %>% group_by(V1) %>% summarise(count=sum(V3))
subset_mono5.bar <- start_m5_50 %>% group_by(V1) %>% summarise(count=sum(V3))

subset_mono3.bar <- mutate(subset_mono3.bar, proportion = count/sum(count))
subset_mono5.bar <- mutate(subset_mono5.bar, proportion = count/sum(count))

max_y <- max(subset_mono3.bar$proportion,subset_mono5.bar$proportion)

#mb5 <- ggplot(data=subset_mono5.raw, aes(x=V1, y=proportion, fill=as.factor((V1 %% 3)+1))) +
mb5 <- ggplot(data=subset_mono5.bar, aes(x=V1, y=proportion, fill=as.factor((V1 %% 3)+1))) +
             geom_bar(stat="identity",width=1, colour="black", size=0.3) +
             theme_bw() +
             xlab("Position relative to start codon") +
             ylab("Proportion of reads") +
             theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
             scale_x_continuous(limits = c(-101,51), expand = c(0, 0), breaks=c(-100,-90,-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50)) +
             scale_y_continuous(limits = c(0, max_y), expand = c(0, 0)) +
             #scale_fill_manual(values=c("#db811a", "#9a9a99", "#6d6e63"), name = "Frame") +
             scale_fill_manual(values=c("#909090", "#909090", "#909090"), name = "Frame") +
             theme(text = element_text(size = 12))

#mb3 <- ggplot(data=subset_mono3.raw, aes(x=V1, y=proportion, fill=as.factor((V1 %% 3)+1))) +
mb3 <- ggplot(data=subset_mono3.bar, aes(x=V1, y=proportion, fill=as.factor((V1 %% 3)+1))) +
             geom_bar(stat="identity",width=1, colour="black", size=0.3) +
             theme_bw() +
             xlab("Position relative to start codon") +
             ylab("Proportion of reads") +
             theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
             scale_x_continuous(limits = c(-50,100), expand = c(0, 0), breaks=c(-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80,90,100)) +
             scale_y_continuous(limits = c(0, max_y), expand = c(0, 0)) +
             #scale_fill_manual(values=c("#db811a", "#9a9a99", "#6d6e63"), name = "Frame") +
             scale_fill_manual(values=c("#909090", "#909090", "#909090"), name = "Frame") +
             theme(text = element_text(size = 12))

#grid.arrange(mb5,mb3)

###
# Pretify plots
###

# Function to save legend
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Save the legend
#+++++++++++++++++++++++
legend_heatmapL5<- get_legend(mLen5)
legend_heatmapL3<- get_legend(mLen3)
legend_barchart5 <- get_legend(mb5)
legend_barchart3 <- get_legend(mb3)

mLen5 <- mLen5 + theme(legend.position="none")
mLen3 <- mLen3 + theme(legend.position="none")
mb5 <- mb5 + theme(legend.position="none")
mb3 <- mb3 + theme(legend.position="none")

##
# Output
##

lay <- rbind(c(1,1,1,1,1,1,1,1,2,3,3,3,3,3,3,3,3,4),
             c(5,5,5,5,5,5,5,5,6,7,7,7,7,7,7,7,7,8))

out1<-arrangeGrob(mb5, legend_barchart5, mb3, legend_barchart3, mLen5, legend_heatmapL5, mLen3, legend_heatmapL3, layout_matrix = lay)
ggsave(file=,args[5], out1, width=500, height=250, unit="mm", dpi=300) 
