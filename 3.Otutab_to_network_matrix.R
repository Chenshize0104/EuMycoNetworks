# Otutab filtered to network matrix
# setwd('E:/Shangri-La_2024/network')

library(readxl)
library(tidyverse)
library(vegan)


rm(list = ls())

otutab <- read.table("data/otutab_decontam_fungi.txt", header=T, row.names=1, sep='\t', check.names=F)
# sum(otutab)
taxon <- read.table("data/taxonomy_decontam_fingi.txt", header=T, row.names=1, sep='\t', check.names=F)
design <- read_excel("data/design.xlsx", sheet="Sheet1")

# otutab_colsums <- merge(design[,1:10],colSums(otutab),by.x=1,by.y=0)
# colnames(otutab_colsums)[which(colnames(otutab_colsums) =='y')] ='Allfungi_colsums'
# View(otutab_colsums)



# remove reads less 0.1% of total reads each sample
otutab.dmin <- otutab
for (i in 1:ncol(otutab.dmin)) {
  col_sum <- sum(otutab.dmin[,i])
  threshold <- col_sum * 0.001
  otutab.dmin[otutab.dmin[,i] <= threshold, i] <- 0
}
sum(otutab.dmin) 

otutab.dmin2 <- otutab.dmin[rowSums(otutab.dmin)>0,]
otutab.dmin2 <- as.data.frame(otutab.dmin)
sum(otutab.dmin2) 

# otutab_colsums <- merge(otutab_colsums,colSums(otutab.dmin2),by.x=1,by.y=0)
# colnames(otutab_colsums)[which(colnames(otutab_colsums) =='y')] ='dmin_colsums'
# View(otutab_colsums)



# remove OTU reads less 0.01% of total reads all samples
otutab.dotu <- otutab.dmin2
sum(rowSums(otutab.dotu))*(1/100000)
otutab.dotu[which(rowSums(otutab.dotu) < sum(rowSums(otutab.dotu))*(1/100000)),] <- 0
sum(otutab.dotu)
summary(colSums(otutab.dotu))

otutab.dotu2 <- otutab.dotu[rowSums(otutab.dotu)>0,]
summary(rowSums(otutab.dotu))
sum(otutab.dotu2) 

# otutab_colsums <- merge(otutab_colsums,colSums(otutab.dotu2),by.x=1,by.y=0)
# colnames(otutab_colsums)[which(colnames(otutab_colsums) =='y')] ='dmin_rare_dotu_colsums'
# View(otutab_colsums)

################################################################################
# All Fungi
taxon.clean <- subset(taxon, rownames(taxon) %in% rownames(otutab.dotu2))
write.table(taxon.clean, "data/taxonomy_fungi_clean.txt", sep = '\t', col.names = NA)

write.table(otutab.dotu2, "data/otutab_fungi_clean.txt", sep = '\t', col.names = NA)



# Glomeromycota
taxon.clean.glom <- subset(taxon.clean,Phylum == "Glomeromycota")
write.table(taxon.clean.glom, "data/taxonomy_glom.txt", sep = '\t', col.names = NA)

otutab.glom <- subset(otutab.dotu2, rownames(otutab.dotu2) %in% rownames(taxon.clean.glom))
sum(otutab.glom)
sort(rowSums(otutab.glom))

# otutab_colsums <- merge(otutab_colsums,colSums(otutab.glom),by.x=1,by.y=0)
# colnames(otutab_colsums)[which(colnames(otutab_colsums) =='y')] ='glom_colsums'
# View(otutab_colsums)

otutab.glom.binary <- otutab.glom
otutab.glom.binary[otutab.glom.binary > 1] <- 1
sum(otutab.glom.binary) 

# otutab_colsums <- merge(otutab_colsums,colSums(otutab.glom.binary),by.x=1,by.y=0)
# colnames(otutab_colsums)[which(colnames(otutab_colsums) =='y')] ='glom_binary_colsums'
# View(otutab_colsums)

write.table(otutab.glom, "data/otutab_glom.txt", sep = '\t', col.names = NA)
write.table(otutab.glom.binary, "data/otutab_glom_binary.txt", sep = '\t', col.names = NA)

write.table(otutab_colsums, "data/otutab_glom_binary_colsums.txt", sep = '\t', col.names = NA)

################################################################################
# network matrix

# network.otutabe <- otutab.glom.binary # 需重新读入csv数据
network.otutabe <- read.table("data/otutab_glom_binary.txt", header=T, row.names=1, sep='\t', check.names=F)
# sum(network.otutabe)

taxon.clean.glom <- read.table("data/taxonomy_glom.txt", header=T, row.names=1, sep='\t', check.names=F)
design <- read_excel("data/design.xlsx", sheet="Sheet1")

network.otutabe2 <- t(network.otutabe)

# 获取矩阵的行名和列名
sample.names <- rownames(network.otutabe2)
AMFOTU.names <- colnames(network.otutabe2)

# 初始化一个空的数据框
network_df <- data.frame()
# 遍历矩阵的行和列
for (i in 1:length(sample.names)) {
  for (j in 1:length(AMFOTU.names)) {
    frequency <- network.otutabe2[i, j]
    if (frequency > 0) {
      network_df <- rbind(network_df, data.frame(Row = sample.names[i], Column = AMFOTU.names[j], Frequency = frequency))
    }
  }
}

network_df2 <- merge(network_df,taxon.clean.glom,by.x = "Column", by.y = 0)
network_df3 <- merge(network_df2,design,by.x = "Row", by.y = "Sample_Fungi")
colnames(network_df3)[which(colnames(network_df3) == "Column")] ="AMFungi_OTUID"
colnames(network_df3)[which(colnames(network_df3) == "Row")] ="Sample"
network_df4 <- network_df3[,c(1,11:24,3,2,4:10)]

network.matrix <- network_df4

write.csv(network.matrix,"data/network_matrix.csv")

################################################################################
################################################################################
################################################################################