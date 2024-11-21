# Otutab delete contamination
# setwd('E:/Shangri-La_2024/network')

# 加载包
library("readxl")
library("tidyverse")
library("dplyr")
library("vegan")
library("phyloseq")
library("ape")
library("devtools")
library("ggplot2")
library("phyloseq")
library("decontam")


rm(list = ls())

# 从所有样品(土壤样品、根系样品和空白对照样品)中，筛选出菌根实验中根系和空白对照样品
otutab <- read.table("data/EufungiMe_all_otutab.txt", header = T, row.names = 1, sep = "\t", check.names = F)
sum(otutab)
dim(otutab)
library(dplyr)
# # ?select()
# otutab1 <- select(otutab, !starts_with("FIP"))
# sum(otutab1)
# dim(otutab1)
# otutab1 <- otutab1[rowSums(otutab1) > 0, ]
# otutab1 <- otutab1[, colSums(otutab1) > 0]
# sum(otutab1)
# dim(otutab1)
# write.table(otutab1, "data/EusoFMe_all_otutab.txt", sep = "\t", col.names = NA)

otutab2 <- select(otutab, starts_with("FIP"))
sum(otutab2)
dim(otutab2)
otutab2 <- otutab2[rowSums(otutab2) > 0, ]
otutab2 <- otutab2[, colSums(otutab2) > 0]
sum(otutab2)
dim(otutab2)
write.table(otutab2, "data/FIPMe_all_otutab.txt", sep = "\t", col.names = NA)

########

rm(list = ls())

# 去掉 Control Sample(BL)中的contaminant OTU
otutabe <- read.table("data/FIPMe_all_otutab.txt", header = T, row.names = 1, sep = "\t", check.names = F)
sum(otutabe)
taxon <- read.table("data/EufungiMe_all_otus_taxonomy.txt", header = T, row.names = 1, sep = "\t", check.names = F)
sample <- read.table("data/sampledata.txt", header = T, row.names = 1, sep = "\t", check.names = F)

# 构建phyloseq对象
physeq = phyloseq(otu_table(otutab, taxa_are_rows = TRUE), sample_data(sample), tax_table(as.matrix(taxon)))
ps <- physeq

head(sample_data(ps))
head(otu_table(ps)[, 1:6])

df <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize), ]
df$Index <- seq(nrow(df))
ggplot(data = df, aes(x = Index, y = LibrarySize, color = Sample_or_Control)) + geom_point()

# 使用的第二种污染物识别方法，method = "prevalence"
sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Control Sample"
contamdf.prev <- isContaminant(ps, method = "prevalence", neg = "is.neg", threshold = 0.5)
table(contamdf.prev$contaminant)
head(which(contamdf.prev$contaminant))
taxon_contaminant <- merge(taxon, contamdf.prev, by = 0)

# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Control Sample", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "True Sample", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos = taxa_sums(ps.pa.pos), pa.neg = taxa_sums(ps.pa.neg), 
                    contaminant = contamdf.prev$contaminant)
ggplot(data = df.pa, aes(x = pa.neg, y = pa.pos, color = contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

hist(contamdf.prev$p, n = 100)

# 去掉otutabe中的contaminant OTU
taxon_sav_2 <- subset(taxon_contaminant, contaminant == "FALSE")[, 1:8]
rownames(taxon_sav_2) <- taxon_sav_2$Row.names # 重新设置行名
taxon_sav_2$Row.names <- NULL # 删除row.names列
write.table(taxon_sav_2, "data/taxonomy_decontam.txt", sep = "\t", col.names = NA)

################################################################################
rm(list = ls())

# 去除非真菌OTU (remove non-fungi OTU)
otutabe <- read.table("data/FIPMe_all_otutab.txt", header=T, row.names=1, sep='\t', check.names=F)
# sum(otutabe) 
taxon <- read.table("data/taxonomy_decontam.txt", header=T, row.names=1, sep='\t', check.names=F)
# remove non fungi OTU
taxon_fungi <- subset(taxon,Kingdom == "Fungi")
design <- read_excel("data/design.xlsx", sheet="Sheet1")

# OTU_idx <- match(rownames(taxon_fungi), rownames(otutabe))
# Sample_idx <- match(design$Sample_Fungi,colnames(otutabe))
# otutab <- otutabe[OTU_idx,Sample_idx]

# otutabe1 <- subset(otutabe, rownames(otutabe) %in% rownames(taxon))
otutabe1 <- otutabe[rownames(otutabe) %in% rownames(taxon_fungi),]
sum(otutabe1)
otutabe2 <- otutabe1[,colnames(otutabe1) %in% design$Sample_Fungi]
sum(otutabe2)
otutab3 <- otutabe2

sum(otutab3) 

write.table(otutab3, "data/otutab_decontam_fungi.txt", sep = '\t', col.names = NA)
write.table(taxon_fungi, "data/taxonomy_decontam_fingi.txt", sep = '\t', col.names = NA)

################################################################################
################################################################################
################################################################################