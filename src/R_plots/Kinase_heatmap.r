library("pheatmap")
library(gplots)
library(fastcluster)
library(RColorBrewer)
library(colorRamps)
library(pheatmap)   
library(ape)

my_palette <- colorRampPalette(c("white", "black", "red"))(n = 100)
if (nrow(gm1) > 100) stop("Too many rows for heatmap, who can read?!")
fontsize_row = 10 - nrow(gm1) / 15

# whole dataset
kin <- read.table("merged.count.normalized.txt", header=TRUE,sep="\t",row.names=1)
gm <-data.matrix(kin[,1:23])
pheatmap(gm, col=my_palette,main ="Kinase heatmap - whole dataset (normalized by proteome size)", cluster_cols=F,fontsize_row=fontsize_row, border_color=NA, scale="none");

# split datasets
kin1 <- read.table("merged.count.normalized_split1.txt", header=TRUE,sep="\t",row.names=1)
gm1 <-data.matrix(kin1[,1:23])
pheatmap(gm1, col=my_palette,main ="Kinase heatmap - split 1 (normalized by proteome size)", cluster_cols=F,fontsize_row=fontsize_row, border_color=NA, scale="none");

kin2 <- read.table("merged.count.normalized_split2.txt", header=TRUE,sep="\t",row.names=1)
gm2 <-data.matrix(kin2)
pheatmap(gm2, col=my_palette,main ="Kinase heatmap - split 2 (normalized by proteome size)", cluster_cols=F,fontsize_row=fontsize_row, border_color=NA, scale="none");

kin3 <- read.table("merged.count.normalized_split3.txt", header=TRUE,sep="\t",row.names=1)
gm3 <-data.matrix(kin3[,1:23])
pheatmap(gm3, col=my_palette,main ="Kinase heatmap - split 3 (normalized by proteome size)", cluster_cols=F,fontsize_row=fontsize_row, border_color=NA);

kin4 <- read.table("merged.count.normalized_split4.txt", header=TRUE,sep="\t",row.names=1)
gm4 <-data.matrix(kin4[,1:23])
pheatmap(gm4, col=my_palette,main ="Kinase heatmap - split 4 (normalized by proteome size)", cluster_cols=F,fontsize_row=fontsize_row, border_color=NA);
