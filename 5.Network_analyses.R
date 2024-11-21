# Network analyses
# setwd('E:/Shangri-La_2024/network')

# 加载包
library(readxl)
library(bipartite)
library(maxnodf)
library(tidyverse)
library(vegan)
library(ggplot2)

rm(list=ls())

set.seed(123)
network_data <- read.read_excel("data/network_matrix.xlsx", sheet="Sheet1")
# network_data <- subset(network_data,Plant_Species!="Euphorbia_jolkinii")
head(network_data)

network_data$lower <- network_data$AMFungi_OTUID
network_data$higher <- network_data$Plant_Species
network_data$webID <- paste0(network_data$Site,network_data$Plot)
network_data$freq <- network_data$Frequency

network_matrix <- network_data[,c("lower", "higher", "webID", "freq")]
network_matrix <- as.data.frame(network_matrix)
typeof(network_matrix)

network <- frame2webs(network_matrix,varnames=c("lower","higher","webID","freq"),type.out="list")
names(network)

########

# "weighted connectance", "weighted NODF"
network_index <- sapply(network, networklevel, index=c("weighted connectance", "weighted NODF"))
network_index

# "interaction evenness"
evenness_prod <- sapply(network, networklevel, index="interaction evenness", intereven="prod")
evenness_prod

#################################################################################
# weighted modularity; 定量模块观察值，每个网络重复三十次，取最大值
set.seed(123)
repeatmodule <- function(network,n){
  likeliresult <- list()
  m <- length(network)
  for (i in 1:n) {
    likeliresult[[i]] <- sapply(network,function(x)computeModules(x, method="Beckett")@likelihood)
  }
  result <- t(matrix(unlist(likeliresult),m,n))
  rownames(result) <- 1:n
  colnames(result) <- names(network)
  result
}

module_rep <- repeatmodule(network,30)
module_index <- apply(module_rep,2,max)
module_index

save.image("network_analyses.RData")
# save.image("network_analyses_Eup.jo_excluded.RData")

#################################################################################
#################################################################################
#################################################################################