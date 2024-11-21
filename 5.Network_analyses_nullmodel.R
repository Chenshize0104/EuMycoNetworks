# Network analyses and nullmodel
# setwd('E:/Shangri-La_2024/network')

# 加载包
library(readxl)
library(tidyverse)
library(bipartite)

rm(list=ls())

set.seed(123)
network_data <- read.read_excel("data/network_matrix.xlsx", sheet="Sheet1")
# network_data <- subset(network_data,Plant_Species!="Euphorbia_jolkinii")
head(network_data)

network_data$lower <- network_data$AMFungi_OTUID
network_data$higher <- network_data$Plant_Species
network_data$webID <- paste0(network_data$Site, network_data$Plot)
network_data$freq <- network_data$Frequency

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
  # null.value <- apply(nulls, 3, function(x)networklevel(x, index=index))
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
# weighted modularity; 每个网络产生100个null model, 每个网络null model计算30次，取最大值
####定量模块---Beckett方法算零模型，不需要设置step---针对多个网
# 耗时很长，10个样地分开算
set.seed(123)
module_null <- function(network,n){
  result <- list()
  result1 <- list()
  result2 <- list()
  result3 <- list()
  likeliresult <- list()
  for(i in 1:length(network)){
    set.seed(123)
    nulls <- nullmodel(network[[i]], N=n, method="r2dtable")
    for (j in 1:30) {
      likeliresult[[j]] <- sapply(nulls,function(x)computeModules(x, method="Beckett")@likelihood)
    }
    result[[i]] <- t(matrix(unlist(likeliresult),length(nulls),30))
    result1[[i]] <- (apply(result[[i]],2,max))
    result2[[i]] <- (mean(result1[[i]]))
    result3[[i]] <- (sd(result1[[i]]))
  }
  module_null_res <- matrix(cbind(unlist(result2),unlist(result3)),length(network),2)
  rownames(module_null_res) <- names(network)
  colnames(module_null_res) <- c("mean","sd")
  module_null_res
}

#n为零模型个数--对于定量模块为100
Wmod <- module_null(network,100)
write.csv(Wmod, "Wmod_null.csv")
# write.csv(Wmod, "Wmod_null_Eup.jo_excluded.csv")

save.image("network_analyses_nullmodel.RData")
# save.image("network_analyses_nullmodel_Eup.jo_excluded.RData")

#################################################################################
#################################################################################
#################################################################################