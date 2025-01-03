# Community rarefaction
# setwd('E:/Shangri-La_2024/network')

# 加载包
library(readxl)
library(tidyverse)
library(bipartite)
library(psych)

rm(list=ls())

design <- read_excel("data/design.xlsx", sheet="Sheet1")
design$webID <- paste0(design$Site,design$Plot)
head(design)
# design <- subset(design,Plant_Species!="Euphorbia_jolkinii")
design <- design[,c(1:8,14)]

network <- read_excel("data/network_matrix.xlsx", sheet="Sheet1")
head(network)
# network <- subset(network,Plant_Species!="Euphorbia_jolkinii")

# network$lower <- network$AMFungi_OTUID
# network$higher <- network$Plant_Species
# network$webID_o <- paste0(network$Site,network$Plot)
# network$freq <- network$Frequency

network_matrix <- network
network_matrix <- as.data.frame(network_matrix)

rm(list=ls()[!ls() %in% c("design","network","network_matrix")])
########
set.seed(123)
# 定义函数进行随机化抽样，从0.1到0.9
rarefaction_sampling <- function(design, network_matrix, webID, num_samples = 1000) {
  design_Site <- subset(design, webID == web)
  network_Site <- subset(network_matrix, webID == web)
  sample_results <- vector("list", num_samples)
  for (i in 1:num_samples) {
    sample_data <- design_Site
    link_data <- network_Site
    # 按比例随机抽样，0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9
    sample_indices <- sample(1:nrow(sample_data), size=nrow(sample_data)*0.9, replace=FALSE)
    sub.sample <- sample_data[sample_indices, ]$Sample
    sub.network <- link_data[link_data$Sample %in% sub.sample, ]
    sample_results[[i]] <- sub.network
  }
  # 返回抽样结果
  return(sample_results)
}

# 对每个样地进行随机化抽样
rarefaction_all <- list()
for (web in unique(design$webID)) {
  rarefaction_all[[web]] <- rarefaction_sampling(design, network_matrix, web)
}

########
# 将网络数据，转化为网络矩阵
rarefaction_network_all <- list()
for (webID in unique(design$webID)) {
  for (i in 1:1000){
    rarefaction_network_all[[webID]][i] <- frame2webs(rarefaction_all[[webID]][[i]],type.out= "list")
  }
}

########
# 计算每个网络的1000个rarefaction的mean、sd
rarefactions_result <- function(rarefaction_network, index){
  raref.value <- unlist(sapply(rarefaction_network, networklevel, index=index))
  mean.raref <- mean(raref.value)
  sd.raref <- sd(raref.value)
  return(cbind(mean.raref, sd.raref))
}

# weighted connectance
Wcon_results <- list()
for (i in names(rarefaction_network_all)) {
  Wcon_results[[i]] <- rarefactions_result(rarefaction_network_all[[i]], index="weighted connectance")
}
Wcon_results
Wcon_results.all <- do.call(rbind,Wcon_results)
Wcon_results.all <- as.data.frame(Wcon_results.all)
Wcon_results.all <- cbind(webID=names(rarefaction_network_all), Wcon_results.all)
Wcon_results.all$rarefaction_index <- "weighted connectance"
Wcon_results.all



# interaction evenness (intereven="prod")
evenness_results <- list()
for (i in names(rarefaction_network_all)) {
  evenness_results[[i]] <- rarefactions_result(rarefaction_network_all[[i]], index="interaction evenness")
}
evenness_results
evenness_results.all <- do.call(rbind,evenness_results)
evenness_results.all <- as.data.frame(evenness_results.all)
evenness_results.all <- cbind(webID=names(rarefaction_network_all), evenness_results.all)
evenness_results.all$rarefaction_index <- "interaction evenness"
evenness_results.all



# weighted NODF
Wnodf_results <- list()
for (i in names(rarefaction_network_all)) {
  Wnodf_results[[i]] <- rarefactions_result(rarefaction_network_all[[i]], index="weighted NODF")
}
Wnodf_results
Wnodf_results.all <- do.call(rbind,Wnodf_results)
Wnodf_results.all <- as.data.frame(Wnodf_results.all)
Wnodf_results.all <- cbind(webID=names(rarefaction_network_all), Wnodf_results.all)
Wnodf_results.all$rarefaction_index <- "weighted NODF"
Wnodf_results.all

rarefaction_data <- bind_rows(evenness_results.all,Wnodf_results.all,Wcon_results.all)
rarefaction_data

write.csv(rarefaction_data,"results/network_analyses_rarefaction_0.9.csv")
# write.csv(rarefaction_data,"results/network_analyses_rarefaction_Eup.jo_excluded_0.9.csv")

save.image("results/network_analyses_rarefaction_0.9.RData")
# save.image("results/network_analyses_rarefaction_Eup.jo_excluded_0.9.RData")

#################################################################################
# weighted modularity，计算每个网络的1000个rarefaction的mean、sd
rarefactions_mod_result <- function(rarefaction_network){
  raref.value <- unlist(sapply(rarefaction_network,function(x)computeModules(x, method="Beckett")@likelihood))
  mean.raref <- mean(raref.value)
  sd.raref <- sd(raref.value)
  return(cbind(mean.raref, sd.raref))
}

Wmod_results <- list()
for (i in names(rarefaction_network_all)) {
  Wmod_results[[i]] <- rarefactions_mod_result(rarefaction_network_all[[i]])
}
Wmod_results
Wmod_results.all <- do.call(rbind,Wmod_results)
Wmod_results.all <- as.data.frame(Wmod_results.all)
Wmod_results.all <- cbind(webID=names(rarefaction_network_all), Wmod_results.all)
Wmod_results.all$rarefaction_index <- "weighted modularity"
Wmod_results.all


rarefaction_data2 <- bind_rows(evenness_results.all,Wnodf_results.all,Wcon_results.all,Wmod_results.all)
rarefaction_data2

write.csv(rarefaction_data2,"results/network_analyses_rarefaction2_0.9.csv")
# write.csv(rarefaction_data,"results/network_analyses_rarefaction_Eup.jo_excluded2_0.9.csv")

save.image("results/network_analyses_rarefaction2_0.9.RData")
# save.image("results/network_analyses_rarefaction_Eup.jo_excluded2_0.9.RData")

#################################################################################
#################################################################################
#################################################################################