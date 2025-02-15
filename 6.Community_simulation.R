# Community simulation
# setwd('E:/Shangri-La_2024/network')

# 加载包
library(readxl)
library(tidyverse)
library(bipartite)
library(psych)

rm(list=ls())

design <- read_excel("data/design.xlsx", sheet="Sheet1")
head(design)
# design <- subset(design,Plant_Species!="Euphorbia_jolkinii")
design <- design[,1:8]

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
# 定义函数进行随机化抽样
simulation_sampling <- function(design, network_matrix, site, num_samples = 1000) {
  design_Site <- subset(design, Site == site)
  network_Site <- subset(network_matrix, Site == site)
  sample_results <- vector("list", num_samples)
  for (i in 1:num_samples) {
    sample_data <- design_Site
    link_data <- network_Site
    # 随机抽样一半的sample
    sample_indices <- sample(1:nrow(sample_data), size=nrow(sample_data)/2, replace=FALSE)
    S1.sample <- sample_data[sample_indices, ]$Sample
    S1 <- link_data[link_data$Sample %in% S1.sample, ]
    S1$webID <- paste0(site, "_S1")
    # 将剩余的sample作为 S2
    S2.sample <- sample_data[-sample_indices, ]$Sample
    S2 <- link_data[link_data$Sample %in% S2.sample, ]
    S2$webID <- paste0(site, "_S2")
    simu_new <- rbind(S1, S2)
    sample_results[[i]] <- simu_new
  }
  # 返回抽样结果
  return(sample_results)
}

# 对每个样地进行随机化抽样
simulation_all <- list()
for (site in unique(design$Site)) {
  simulation_all[[site]] <- simulation_sampling(design, network_matrix, site)
}

########
# 将网络数据，转化为网络矩阵
simulation_network_all <- list()
for (site in unique(design$Site)) {
  simulation_network_all[[unique(simulation_all[[site]][[1]]$webID)[1]]] <- list()
  simulation_network_all[[unique(simulation_all[[site]][[1]]$webID)[2]]] <- list()
  for (i in 1:1000){
    simulation_network_all[[unique(simulation_all[[site]][[i]]$webID)[1]]][[i]] <- frame2webs(subset(simulation_all[[site]][[i]],webID==unique(simulation_all[[site]][[i]]$webID)[1]),type.out= "list")[[1]]
    
    simulation_network_all[[unique(simulation_all[[site]][[i]]$webID)[2]]][[i]] <- frame2webs(subset(simulation_all[[site]][[i]],webID==unique(simulation_all[[site]][[i]]$webID)[2]),type.out= "list")[[1]]
  }
}

########
# 计算每个网络的1000个simulation的mean、sd
simulations_result <- function(simulation_network, index){
  simul.value <- unlist(sapply(simulation_network, networklevel, index=index))
  mean.simul <- mean(simul.value)
  sd.simul <- sd(simul.value)
  return(cbind(mean.simul, sd.simul))
}

# weighted connectance
Wcon_results <- list()
for (i in names(simulation_network_all)) {
  Wcon_results[[i]] <- simulations_result(simulation_network_all[[i]], index="weighted connectance")
}
Wcon_results
Wcon_results.all <- do.call(rbind,Wcon_results)
Wcon_results.all <- as.data.frame(Wcon_results.all)
Wcon_results.all <- cbind(webID=names(simulation_network_all), Wcon_results.all)
Wcon_results.all$simulation_index <- "weighted connectance"
Wcon_results.all



# interaction evenness (intereven="prod")
evenness_results <- list()
for (i in names(simulation_network_all)) {
  evenness_results[[i]] <- simulations_result(simulation_network_all[[i]], index="interaction evenness")
}
evenness_results
evenness_results.all <- do.call(rbind,evenness_results)
evenness_results.all <- as.data.frame(evenness_results.all)
evenness_results.all <- cbind(webID=names(simulation_network_all), evenness_results.all)
evenness_results.all$simulation_index <- "interaction evenness"
evenness_results.all



# weighted NODF
Wnodf_results <- list()
for (i in names(simulation_network_all)) {
  Wnodf_results[[i]] <- simulations_result(simulation_network_all[[i]], index="weighted NODF")
}
Wnodf_results
Wnodf_results.all <- do.call(rbind,Wnodf_results)
Wnodf_results.all <- as.data.frame(Wnodf_results.all)
Wnodf_results.all <- cbind(webID=names(simulation_network_all), Wnodf_results.all)
Wnodf_results.all$simulation_index <- "weighted NODF"
Wnodf_results.all

simulation_data <- bind_rows(evenness_results.all,Wnodf_results.all,Wcon_results.all)
simulation_data

write.csv(simulation_data,"results/network_analyses_simulation.csv")
# write.csv(simulation_data,"results/network_analyses_simulation_Eup.jo_excluded.csv")

save.image("results/network_analyses_simulation.RData")
# save.image("results/network_analyses_simulation_Eup.jo_excluded.RData")

#################################################################################
# weighted modularity，计算每个网络的1000个simulation的mean、sd
simulations_mod_result <- function(simulation_network){
  simul.value <- unlist(sapply(simulation_network,function(x)computeModules(x, method="Beckett")@likelihood))
  mean.simul <- mean(simul.value)
  sd.simul <- sd(simul.value)
  return(cbind(mean.simul, sd.simul))
}

Wmod_results <- list()
for (i in names(simulation_network_all)) {
  Wmod_results[[i]] <- simulations_mod_result(simulation_network_all[[i]])
}
Wmod_results
Wmod_results.all <- do.call(rbind,Wmod_results)
Wmod_results.all <- as.data.frame(Wmod_results.all)
Wmod_results.all <- cbind(webID=names(simulation_network_all), Wmod_results.all)
Wmod_results.all$simulation_index <- "weighted modularity"
Wmod_results.all


simulation_data2 <- bind_rows(evenness_results.all,Wnodf_results.all,Wcon_results.all,Wmod_results.all)
simulation_data2

write.csv(simulation_data2,"results/network_analyses_simulation2.csv")
# write.csv(simulation_data,"results/network_analyses_simulation_Eup.jo_excluded2.csv")

save.image("results/network_analyses_simulation2.RData")
# save.image("results/network_analyses_simulation_Eup.jo_excluded2.RData")

#################################################################################
#################################################################################
#################################################################################