source("https://bioconductor.org/biocLite.R")
library(cummeRbund)
require(xlsx)
library(gplots)
library(ComplexHeatmap)
library(circlize)
library(dplyr)

# This script plots the heatmap of specifc genes of interest belonging to different category
path <- "E:/Chris_UM/Analysis/CoreData/15_ZhuBo_RNASeq2/cuffdiff4"
setwd(path)

# File with Gene name and category column
family <- read.table("targetGeneFamilies.txt", header = T)


####################################################################
#IDs of the samples to plot in the heatmap
sampleIds <- c("PG_WT", "PV_WT")

#This file has the pairs of sample names. This will be used to draw a second heatmap with trace for which comparison a row is significant
p1 <- c("PG_WT")
p2 <- c("PV_WT")
pairs <- data.frame(p1, p2)
pairs$name <- paste(pairs$p1, pairs$p2, sep = " vs ")

# plotTitle <- paste("Top ",topN, " regulated genes for each sample pair comparison", sep = " ")
plotTitle <- paste(pairs$name[1], "sample comparison for specific gene family list")


#read all the cufflinks data using cummeRbund package
cuff <- readCufflinks(dir = path)

#sample names
sample.names<-samples(genes(cuff))

#Get the gene set for the gene list belonging to gene family of interest
# geneSet <- getGenes(cuff, as.vector(family$Gene))

#take the top genes
gene.diff<-diffData(object = genes(cuff), x = pairs$p1[1], y = pairs$p2[1])

#filter using p_value and fold_change cutoff
filtered <- subset(x = gene.diff, status == "OK" & q_value < 0.05 & (log2_fold_change > 1 | log2_fold_change < -1))

#Find which genes from family list are present in the filtered gene list
regulated <- semi_join(x = family, y = filtered, by = c("Gene" = "gene_id"))

#select top N gene IDs
topGenes <- getGenes(object = cuff, geneIdList = regulated$Gene, sampleIdList = sampleIds)

#FPKM matrix for top N genes
expressionMat <- log2(data.matrix(repFpkmMatrix(topGenes, fullnames = TRUE)) + 0.005)

# z-score with respect to each gene
for(i in 1:nrow(expressionMat)){
  expressionMat[i,] <- (expressionMat[i,] - mean(expressionMat[i,]))/sd(expressionMat[i,])
}

#rename the row names of expressionMat with just GeneIDs as this will be required to match the geneIds from genesDF 
geneIdToName <- as.data.frame(do.call(rbind, strsplit(row.names(expressionMat), split = "|", fixed = TRUE)))
names(geneIdToName) <- c("GeneName", "GeneID")
rownames(expressionMat) <- geneIdToName$GeneName

#geneID, geneName, gene family class in one data frame
geneIdToName <- inner_join(geneIdToName, regulated, by = c("GeneID" = "Gene"))
geneCategory <- data.frame("Family" = geneIdToName$Category, row.names = geneIdToName$GeneName)

########### ***************************************
#Order geneIdToName same as expressionMat : Do this
########### ***************************************


#heatmap annotations:
haRow = rowAnnotation(df = geneCategory, name = "Gene Family", 
                      width = unit(0.5, "cm"), 
                      col = list(Family = c("EGF" = "red", "TGF" = "blue", "NOTCH_EGF" = "orange", "IGF" = "purple")), 
                      show_annotation_name = FALSE, 
                      annotation_name_gp = gpar(fontsize = 11),
                      annotation_legend_param = list(Family = list(title = "Gene Family", title_gp = gpar(fontsize = 11))))

haUp <- HeatmapAnnotation(b1 = anno_boxplot(expressionMat, which = "column", axis = TRUE), 
                          annotation_height = 10)


#plot heatmap
rpkmHeatmap <- Heatmap(expressionMat, row_title = "Genes", column_title = plotTitle,
                       col = colorRamp2(breaks = c(min(expressionMat), 0, max(expressionMat)), c("green", "black", "red"), space = "LAB"), 
                       top_annotation = haUp, top_annotation_height = unit(3, "cm"), 
                       show_row_names = TRUE, row_names_side = "left", column_names_gp = gpar(fontsize = 10), 
                       show_row_dend = FALSE, cluster_columns = FALSE, 
                       width = unit(10, "cm"), row_names_max_width = unit(15, "cm"), 
                       heatmap_legend_param = list(title = "z-score of log2(FPKM + 1)", color_bar = "continuous")) 

png(filename = "PG_WT_vs_PV_WT_geneFamily.png", width=1000, height=1200, res = 120)

rpkmHeatmap + haRow

dev.off()

