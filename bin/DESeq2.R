#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=T)

# prepare input data
countdata <- read.table(args[1], header=TRUE, row.names=1)
coldata <- read.table(args[2], header=FALSE, row.names=1)

countdata <- as.matrix(countdata[rownames(coldata)])
countdata <- round(countdata)

colnames(coldata) <- "condition"
coldata <- as.data.frame(coldata)
coldata$condition <- as.factor(coldata$condition)

# run DESeq2 analysis
library("DESeq2")

dds <- DESeqDataSetFromMatrix(countdata, coldata, design=~condition)
dds$condition <- relevel(dds$condition, ref=args[3])
dds <- DESeq(dds)
res <- results(dds, alpha=0.01)

resultsNames(dds)
resLFC <- lfcShrink(dds, coef=2, res=res, type="apeglm")
resFile <- na.omit(resLFC)
resFile <- resFile[order(resFile$padj),]

write.table(data.frame("feature"=rownames(resFile),resFile), file="output.tsv", quote=F, sep="\t", row.names=FALSE)
