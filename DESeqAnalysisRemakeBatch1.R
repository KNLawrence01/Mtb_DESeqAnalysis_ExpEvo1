install.packages("ggplot2")
install.packages("dplyr")
install.packages("tidyverse")
install.packages("ggpubr")
install.packages("stringr")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(DESeq2)
library(stringr)

### ANCESTOR: BvP ###
# read metadata file
sampleData <- read.delim("mtb_ExpEvo_RNA_metadata.txt", sep="", header = TRUE)

# reformat metadata
sampleTable <- data.frame(SampleName = sampleData$LibraryName,
                          CountsFile = sampleData$CountsFile,
                          Strain = factor(sampleData$Strain),
                          Genotype = factor(sampleData$Genotype),
                          Condition = factor(sampleData$Condition,levels = c("Planktonic","Biofilm")),
                          WetWeight = as.numeric(sampleData$WetWeight),
                          Clade = factor(sampleData$Clade),
                          SampleID = sampleData$SampleID)

# subset to only ancestral samples
A.table <- sampleTable[sampleTable$Genotype == "Ancestral",]

# DE calculations by condition (biofilm vs plantkonic)
A.DESeq <- DESeqDataSetFromHTSeqCount(sampleTable = A.table, directory = ".", design = ~ Condition)

# filter out genes with 0 counts in all samples
A.DESeq <- A.DESeq[rowSums(counts(A.DESeq)) > 1,]

# run DESeq analysis
A.DESeq <- DESeq(A.DESeq)
resultsNames(A.DESeq)

# results of biofilm vs planktonic comparison
ABvP.result <- results(A.DESeq, alpha = 0.05,lfcThreshold = 0, contrast=c("Condition","Biofilm","Planktonic"))
write.csv(x = ABvP.result,"2024_07_10_all_ABvP.csv",quote=F)

strains <- c("31","49","55","72","345","540")

for (strain in strains){
  st <- A.table[A.table$Strain == strain,]
  t <- DESeqDataSetFromHTSeqCount(sampleTable = st, directory = ".", design = ~ Condition)
  t <- t[rowSums(counts(t)) > 1,]
  t <- DESeq(t)
  r <- results(t, alpha = 0.05, lfcThreshold = 0, contrast=c("Condition","Biofilm","Planktonic"))
  ###sub <- subset(r, padj < 0.05)
  write.csv(r,paste(paste("2024_07_10",strain,"ancBvP",sep="_"),".csv",sep=""),quote=F)
}

New_ancBvP_31 <- data.frame(read.csv("2024_07_10_31_ancBvP.csv",header=T))
New_ancBvP_49 <- data.frame(read.csv("2024_07_10_49_ancBvP.csv",header=T))
New_ancBvP_55 <- data.frame(read.csv("2024_07_10_55_ancBvP.csv",header=T))
New_ancBvP_72 <- data.frame(read.csv("2024_07_10_72_ancBvP.csv",header=T))
New_ancBvP_345 <- data.frame(read.csv("2024_07_10_345_ancBvP.csv",header=T))
New_ancBvP_540 <- data.frame(read.csv("2024_07_10_540_ancBvP.csv",header=T))

### EVOLVED: BvP ###
# subset to only evolved samples
E.table <- sampleTable[sampleTable$Genotype == "Evolved",]

# DE calculations by condition (biofilm vs plantkonic)
E.DESeq <- DESeqDataSetFromHTSeqCount(sampleTable = E.table, directory = ".", design = ~ Condition)

# filter out genes with 0 counts in all samples
E.DESeq <- E.DESeq[rowSums(counts(E.DESeq)) > 1,]

# run DESeq analysis
E.DESeq <- DESeq(E.DESeq)
resultsNames(E.DESeq)

# results of biofilm vs planktonic comparison
EBvP.result <- results(E.DESeq, alpha = 0.05,lfcThreshold = 0, contrast=c("Condition","Biofilm","Planktonic"))
write.csv(x = EBvP.result,"2024_07_10_all_EBvP.csv",quote=F)

strains <- c("31","49","55","72","345","540")

for (strain in strains){
  st <- E.table[E.table$Strain == strain,]
  t <- DESeqDataSetFromHTSeqCount(sampleTable = st, directory = ".", design = ~ Condition)
  t <- t[rowSums(counts(t)) > 1,]
  t <- DESeq(t)
  r <- results(t, alpha = 0.05, lfcThreshold = 0, contrast=c("Condition","Biofilm","Planktonic"))
  ###sub <- subset(r, padj < 0.05)
  write.csv(r,paste(paste("2024_07_10",strain,"evoBvP",sep="_"),".csv",sep=""),quote=F)
}

New_evoBvP_31 <- data.frame(read.csv("2024_07_10_31_evoBvP.csv",header=T))
New_evoBvP_49 <- data.frame(read.csv("2024_07_10_49_evoBvP.csv",header=T))
New_evoBvP_55 <- data.frame(read.csv("2024_07_10_55_evoBvP.csv",header=T))
New_evoBvP_72 <- data.frame(read.csv("2024_07_10_72_evoBvP.csv",header=T))
New_evoBvP_345 <- data.frame(read.csv("2024_07_10_345_evoBvP.csv",header=T))
New_evoBvP_540 <- data.frame(read.csv("2024_07_10_540_evoBvP.csv",header=T))

### BIOFILM: AvE ###
# subset to only biofilm samples
B.table <- sampleTable[sampleTable$Condition == "Biofilm",]

# create sample table for DESeq2 input
B.DESeq <- DESeqDataSetFromHTSeqCount(sampleTable = B.table, directory = ".",design = ~ Genotype)

# remove rows with 0 counts
B.DESeq <- B.DESeq[rowSums(counts(B.DESeq)) > 1,]

# run DESeq analysis
B.DESeq <- DESeq(B.DESeq)
resultsNames(B.DESeq)

# results of biofilm: AvE
BEvA.result <- results(B.DESeq, alpha = 0.05, lfcThreshold = 0, contrast=c("Genotype","Evolved","Ancestral"))
write.csv(x = BEvA.result,"2024_07_10_all_BEvA.csv",quote=F)

# separate by individual strains (Supplementary Data 3 - BFEvA for each strain)
strains <- c("31","49","55","72","345","540")

for (strain in strains){
  st <- B.table[B.table$Strain == strain,]
  t <- DESeqDataSetFromHTSeqCount(sampleTable = st, directory = ".", design = ~ Genotype)
  t <- t[rowSums(counts(t)) > 1,]
  t <- DESeq(t)
  r <- results(t, alpha = 0.05, lfcThreshold = 0, contrast=c("Genotype","Evolved","Ancestral"))
  ###sub <- subset(r, padj < 0.05)
  write.csv(r,paste(paste("2024_07_10",strain,"BFEvA",sep="_"),".csv",sep=""),quote=F)
}

New_BFEvA_31 <- data.frame(read.csv("2024_07_10_31_BFEvA.csv",header=T))
New_BFEvA_49 <- data.frame(read.csv("2024_07_10_49_BFEvA.csv",header=T))
New_BFEvA_55 <- data.frame(read.csv("2024_07_10_55_BFEvA.csv",header=T))
New_BFEvA_72 <- data.frame(read.csv("2024_07_10_72_BFEvA.csv",header=T))
New_BFEvA_345 <- data.frame(read.csv("2024_07_10_345_BFEvA.csv",header=T))
New_BFEvA_540 <- data.frame(read.csv("2024_07_10_540_BFEvA.csv",header=T))

### PLANKTONIC: AvE ###
# subset to only biofilm samples
P.table <- sampleTable[sampleTable$Condition == "Planktonic",]

# create sample table for DESeq2 input
P.DESeq <- DESeqDataSetFromHTSeqCount(sampleTable = P.table, directory = ".",design = ~ Genotype)

# remove rows with 0 counts
P.DESeq <- P.DESeq[rowSums(counts(P.DESeq)) > 1,]

# run DESeq analysis
P.DESeq <- DESeq(P.DESeq)
resultsNames(P.DESeq)

# results of biofilm: AvE
PEvA.result <- results(P.DESeq, alpha = 0.05, lfcThreshold = 0, contrast=c("Genotype","Evolved","Ancestral"))
write.csv(x = PEvA.result,"2024_07_10_all_PEvA.csv",quote=F)

# separate by individual strains
strains <- c("31","49","55","72","345","540")

for (strain in strains){
  st <- P.table[P.table$Strain == strain,]
  t <- DESeqDataSetFromHTSeqCount(sampleTable = st, directory = ".", design = ~ Genotype)
  t <- t[rowSums(counts(t)) > 1,]
  t <- DESeq(t)
  r <- results(t, alpha = 0.05, lfcThreshold = 0, contrast=c("Genotype","Evolved","Ancestral"))
  ###sub <- subset(r, padj < 0.05)
  write.csv(r,paste(paste("2024_07_10",strain,"PEvA",sep="_"),".csv",sep=""),quote=F)
}

New_PEvA_31 <- data.frame(read.csv("2024_07_10_31_PEvA.csv",header=T))
New_PEvA_49 <- data.frame(read.csv("2024_07_10_49_PEvA.csv",header=T))
New_PEvA_55 <- data.frame(read.csv("2024_07_10_55_PEvA.csv",header=T))
New_PEvA_72 <- data.frame(read.csv("2024_07_10_72_PEvA.csv",header=T))
New_PEvA_345 <- data.frame(read.csv("2024_07_10_345_PEvA.csv",header=T))
New_PEvA_540 <- data.frame(read.csv("2024_07_10_540_PEvA.csv",header=T))
