# Beta diversity
# setwd('E:/Shangri-La_2024/network')

# 加载包
library(readxl)
library(tidyverse)
library(bipartite)
library(betapart)
library(betalink)


# Plant beta diversity
rm(list=ls())

design <- read_excel("data/design.xlsx", sheet="Sheet1")
design$Site_Plot <- paste0(design$Site,"_",design$Plot)

# 提取样地Plot分组信息
level_design <- design %>% 
  group_by(Site,SiteName,Plot,Site_Plot) %>% 
  summarize(n=n(),.groups='drop') # 使用.groups='drop'来取消分组
sum(level_design$n)

# 汇总每个Plot植物物种数，含多度信息
plantsp_count <- design %>% 
  count(Site_Plot,Plant_Species)
sum(plantsp_count$n)

# 转化为Plot，species矩阵（使用tidyr的spread()函数将计数结果转换为宽格式数据框）
plantsp_matrix <- plantsp_count %>% 
  spread(Plant_Species,n,fill=0)
sum(plantsp_matrix[,2:dim(plantsp_matrix)[2]]) # 1620

# binary，将Plot，species矩阵去除多度，转化为丰度矩阵
plantsp_richness_matrix <- plantsp_matrix %>% 
  mutate(across(2:dim(plantsp_matrix)[2],~ifelse(.>1,1,.)))
sum(plantsp_richness_matrix[,2:dim(plantsp_richness_matrix)[2]]) # 574


# 合并样地Plot分组信息和Plot的植物物种数richness
level_plantsp_richness <- left_join(level_design,plantsp_richness_matrix,by=join_by(Site_Plot))
sum(level_plantsp_richness[,6:dim(level_plantsp_richness)[2]])

# 清除临时变量，保留部分变量
rm(list=ls()[!ls() %in% c("design","level_design","plantsp_matrix","plantsp_richness_matrix","level_plantsp_richness")])

#########
# 计算beta diversity组分
sum(plantsp_matrix[,2:dim(plantsp_matrix)[2]]) # 1620
sum(plantsp_richness_matrix[,2:dim(plantsp_richness_matrix)[2]]) # 574
comun_matrix <- plantsp_richness_matrix
comun_matrix <- as.data.frame(comun_matrix)
rownames(comun_matrix) <- comun_matrix$Site_Plot
comun_matrix$Site_Plot <- NULL
sum(comun_matrix)

dist <- beta.pair(comun_matrix,index.family="jaccard")
#turnover partition
dist[[1]]
plant_turnover <- as.data.frame(as.matrix(dist[[1]]))
#nestedness partition
dist[[2]]
plant_nestedness <- as.data.frame(as.matrix(dist[[2]]))
#all beta diversity
dist[[3]]
plant_all_beta_diversity <- as.data.frame(as.matrix(dist[[3]]))
plant_beta <- rbind(plant_turnover,plant_nestedness,plant_all_beta_diversity)

write.csv(plant_beta,"beta_diversity/plantbeta_betapart.csv")

################################################################################
# AMF beta diversity
rm(list=ls())

design <- read_excel("data/design.xlsx", sheet="Sheet1")
# design <- subset(design,Plant_Species!="Euphorbia_jolkinii")
design$Site_Plot <- paste0(design$Site,"_",design$Plot)

otutabe <- read_excel("data/otutab_glom.xlsx", sheet="Sheet1")
otutabe <- as.data.frame(otutabe)
rownames(otutabe) <- otutabe[[1]]
otutabe <- otutabe[,-1]
sum(otutabe)
otutabe[otutabe > 1] <- 1
sum(otutabe)

design_otutabe <- merge(design,otutabe,by.x="Sample",by.y=0)

which(names(design_otutabe)=="OTU_172")-1 #16
design <- design_otutabe[,1:16]

# 提取样地Plot分组信息
level_design <- design %>% 
  group_by(Site,SiteName,Plot,Site_Plot) %>% 
  summarize(n=n(),.groups='drop') # 使用.groups='drop'来取消分组
sum(level_design$n)

which(names(design_otutabe)=="Site_Plot") #16
which(names(design_otutabe)=="OTU_172") #17
AMF_otu_matrix <- design_otutabe[,c(16,17:dim(design_otutabe)[2])]
sum(AMF_otu_matrix[,2:dim(AMF_otu_matrix)[2]])

# 将同一plot植物根系中的AMF_otu_matrix进行合并，得到每个plot AMF_otu_matrix物种数，含多度信息
AMF_otu_matrix <- aggregate(AMF_otu_matrix[,-1],by=list(AMF_otu_matrix[,1]),FUN=sum,na.rm=TRUE,drop=F)
sum(AMF_otu_matrix[,2:dim(AMF_otu_matrix)[2]])
colnames(AMF_otu_matrix)[which(colnames(AMF_otu_matrix)=='Group.1')]='Site_Plot'

AMF_otu_richness_matrix <- AMF_otu_matrix %>% 
  mutate(across(2:dim(AMF_otu_matrix)[2],~ifelse(.>1,1,.)))
sum(AMF_otu_richness_matrix[,2:dim(AMF_otu_richness_matrix)[2]]) # 3996

# 清除临时变量，保留部分变量
rm(list=ls()[!ls() %in% c("design","otutabe","level_design","AMF_otu_matrix","AMF_otu_richness_matrix")])

#########
# 计算beta diversity组分
sum(AMF_otu_matrix[,2:dim(AMF_otu_matrix)[2]]) # 19670
sum(AMF_otu_richness_matrix[,2:dim(AMF_otu_richness_matrix)[2]]) # 3996
comun_matrix <- AMF_otu_richness_matrix
rownames(comun_matrix) <- comun_matrix$Site_Plot
comun_matrix$Site_Plot <- NULL
sum(comun_matrix)

dist <- beta.pair(comun_matrix, index.family="jaccard")
#turnover partition
dist[[1]]
fungi_turnover <- as.data.frame(as.matrix(dist[[1]]))
#nestedness partition
dist[[2]]
fungi_nestedness <- as.data.frame(as.matrix(dist[[2]]))
#all beta diversity
dist[[3]]
AMF_all_beta_diversity <- as.data.frame(as.matrix(dist[[3]]))
AMF_beta <- rbind(fungi_turnover,fungi_nestedness,fungi_all_beta_diversity)

write.csv(AMF_beta,"beta_diversity/AMF_beta_betapart.csv")
# write.csv(AMF_beta,"beta_diversity/AMF_beta_betapart_Eup.jo_excluded.csv")

################################################################################
# Interactions beta diversity
rm(list=ls())

network_data <- read.csv("data/network_matrix_analyses.csv")
# network_data <- subset(network_data,Plant_Species!="Euphorbia_jolkinii")
head(network_data)

network_data$lower <- network_data$AMFungi_OTUID
network_data$higher <- network_data$Plant_Species
network_data$webID <- paste0(network_data$Site,network_data$Plot)
network_data$freq <- network_data$Frequency

network_matrix <- network_data[,c("lower", "higher", "webID", "freq")]
network_matrix <- as.data.frame(network_matrix)
typeof(network_matrix)

# network <- frame2webs(network_matrix,varnames=c("lower","higher","webID","freq"),type.out="list")
network <- frame2webs(network_matrix,varnames=c("lower","higher","webID","freq"),emptylist = FALSE,type.out = "array")

##物种组分不拆分
##partitioning----计算ST和OS方法：poisot, commondenom
##如果partitioning="poisot")，结果和betalink包一致
##S (Species), OS (Only Shared species links), WN (Whole Network links), ST (Species Turnover links)
(beta <- betalinkr_multi(network, partitioning="commondenom"))


##物种组分继续拆分
(betap <- betalinkr_multi(network, partitioning="commondenom", partition.st=T))

write.csv(betap,"beta_diversity/interactions_beta_betapart.csv")
# write.csv(betap,"beta_diversity/interactions_beta_betapart_Eup.jo_excluded.csv")

#################################################################################
#################################################################################
#################################################################################