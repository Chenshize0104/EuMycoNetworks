# Network analyses and nullmodel
# setwd('E:/Shangri-La_2024/network')

# 加载包
library(readxl)
library(tidyverse)
library(bipartite)

rm(list=ls())

set.seed(123)
network_data <- read_excel("data/network_matrix.xlsx", sheet="Sheet1")
# network_data <- subset(network_data,Plant_Species!="Euphorbia_jolkinii")
head(network_data)

# network_data$lower <- network_data$AMFungi_OTUID
# network_data$higher <- network_data$Plant_Species
# network_data$webID <- paste0(network_data$Site, network_data$Plot)
# network_data$freq <- network_data$Frequency

network_matrix <- network_data[, c("lower", "higher", "webID", "freq")]
network_matrix <- as.data.frame(network_matrix)
typeof(network_matrix)

network <- frame2webs(network_matrix, varnames=c("lower", "higher", "webID", "freq"), type.out="list")
names(network)

##########
# 每个网络产生1000个null model
length(network)
procedure_null <- function(network, n){
  nulls <- list()
  for (i in names(network)) {
    nulls[[i]] <- nullmodel(network[[i]], N=n, method="r2dtable")
  }
  return(nulls)
}
nulls_all <- procedure_null(network, 1000)


##########
# 计算每个网络的observe与1000个null model的mean、sd和Zscore
observe_null <- function(network, nulls, index){
  obs <- networklevel(network, index=index)
  null.value <- unlist(sapply(nulls, networklevel, index=index))
  mean.null <- mean(null.value)
  sd.null <- sd(null.value)
  Z <- (obs-mean.null)/sd(null.value)
  return(cbind(obs, mean.null, sd.null, Z))
}

# weighted connectance
Wcon_results <- list()
for (i in names(network)) {
  Wcon_results[[i]] <- observe_null(network[[i]], nulls_all[[i]], index="weighted connectance")
}
Wcon_results[1]
Wcon_results
Wcon_results.all <- do.call(rbind,Wcon_results)
Wcon_results.all <- as.data.frame(Wcon_results.all)
Wcon_results.all <- cbind(webID=names(Wcon_results), Wcon_results.all)
Wcon_results.all$null_index <- "weighted connectance"
Wcon_results.all$null_model <- "r2dtable"



# interaction evenness (intereven="prod")
evenness_results <- list()
for (i in names(network)) {
  evenness_results[[i]] <- observe_null(network[[i]], nulls_all[[i]], index="interaction evenness")
}
evenness_results[1]
evenness_results
evenness_results.all <- do.call(rbind,evenness_results)
evenness_results.all <- as.data.frame(evenness_results.all)
evenness_results.all <- cbind(webID=names(evenness_results), evenness_results.all)
evenness_results.all$null_index <- "evenness_prod"
evenness_results.all$null_model <- "r2dtable"



# weighted NODF
Wnodf_results <- list()
for (i in names(network)) {
  Wnodf_results[[i]] <- observe_null(network[[i]], nulls_all[[i]], index="weighted NODF")
}
Wnodf_results[1]
Wnodf_results
Wnodf_results.all <- do.call(rbind,Wnodf_results)
Wnodf_results.all <- as.data.frame(Wnodf_results.all)
Wnodf_results.all <- cbind(webID=names(Wnodf_results), Wnodf_results.all)
Wnodf_results.all$null_index <- "weighted NODF"
Wnodf_results.all$null_model <- "r2dtable"

null_model_data <- bind_rows(evenness_results.all,Wnodf_results.all,Wcon_results.all)

write.csv(null_model_data,"network_analyses_nullmodel.csv")
# write.csv(null_model_data,"network_analyses_nullmodel_Eup.jo_excluded.csv")

#################################################################################
# weighted modularity; 计算每个网络的observe与1000个null model的mean、sd和Zscore
observe_null.mod <- function(network, nulls){
  obs <- computeModules(network)@likelihood
  null.value <- sapply(nulls,function(x)computeModules(x, method="Beckett")@likelihood)
  mean.null <- mean(null.value)
  sd.null <- sd(null.value)
  Z <- (obs-mean.null)/sd(null.value)
  return(cbind(obs, mean.null, sd.null, Z))
}

Wmod_results <- list()
for (i in names(network)) {
  Wmod_results[[i]] <- observe_null.mod(network[[i]], nulls_all[[i]])
}
Wmod_results[1]

Wmod_results

Wmod_results.all <- do.call(rbind,Wmod_results)
Wmod_results.all <- as.data.frame(Wmod_results.all)
Wmod_results.all <- cbind(webID=names(Wmod_results), Wmod_results.all)
Wmod_results.all$null_index <- "weighted modularity"
Wmod_results.all$null_model <- "r2dtable"

null_model_data2 <- bind_rows(evenness_results.all,Wnodf_results.all,Wcon_results.all,H2_results.all,generality.HL_results.all,vulnerability.LL_results.all,Wmod_results.all)

write.csv(null_model_data2,"network_analyses_nullmodel2.csv")
# write.csv(null_model_data2,"network_analyses_nullmodel_Eup.jo_excluded2.csv")
save.image("network_analyses_nullmodel.RData")
# save.image("network_analyses_nullmodel_Eup.jo_excluded.RData")

#################################################################################
#################################################################################
#################################################################################
