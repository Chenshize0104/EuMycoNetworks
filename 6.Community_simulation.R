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

network <- read.read_excel("data/network_matrix.xlsx", sheet="Sheet1")
head(network)
# network <- subset(network,Plant_Species!="Euphorbia_jolkinii")

network$lower <- network$AMFungi_OTUID
network$higher <- network$Plant_Species
network$webID_o <- paste0(network$Site,network$Plot)
network$freq <- network$Frequency

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
  # simul.value <- apply(simulation_network, 3, function(x)networklevel(x, index=index))
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
# weighted modularity; 每个网络产生100个随机群落, 每个网络计算30次，取最大值
design <- read_excel("data/design.xlsx", sheet="Sheet1")
head(design)
# design <- subset(design,Plant_Species!="Euphorbia_jolkinii")
design <- design[,1:8]

network <- read.read_excel("data/network_matrix.xlsx", sheet="Sheet1")
head(network)
# network <- subset(network,Plant_Species!="Euphorbia_jolkinii")

network$lower <- network$AMFungi_OTUID
network$higher <- network$Plant_Species
network$webID <- paste0(network$Site,network$Plot)
network$freq <- network$Frequency

network_matrix <- network
network_matrix <- as.data.frame(network_matrix)


########
# 对样地内所有样品（包括没有真菌的样品）进行随机化抽样

# 耗时很长，10个样地分开算
unique(design$Site)

result_net1 <- list()
result_net2 <- list()
result_all <- list()
for (site in unique(design$Site)){
	design_Site <- subset(design,Site==site)
	# dim(design_Site)[1]/2
	# nrow(design_Site)/2

	network_Site <- subset(network_matrix,Site==site)
	# dim(network_Site)[1]/2
	# nrow(network_Site)/2

	set.seed(123)
	sample_data <- design_Site
	link_data <- network_Site

	# 设置抽样次数
	num_samples <- 100

	# 创建一个空列表，用于存放抽样结果
	sample_results <- vector("list", num_samples)

	# 对数据库进行抽样100次
	for (i in 1:num_samples) {
	  # 随机抽样一半的sample
	  sample_indices <- sample(1:nrow(sample_data), size = nrow(sample_data)/2, replace = FALSE)
	  s1.sample <- sample_data[sample_indices, ]$Sample
	  print(length(s1.sample))
	  s1 <- link_data[link_data$Sample %in% s1.sample, ]
	  length(unique(s1$Sample))
	  s1$webID_simulation <- paste0(s1$Site,"_s1")
	  # 将剩余的sample作为 s2
	  s2.sample  <- sample_data[-sample_indices, ]$Sample
	  length(s2.sample)
	  s2 <- link_data[link_data$Sample %in% s2.sample, ]
	  print(length(unique(s2$Sample)))
	  s2$webID_simulation <- paste0(s2$Site,"_s2")
	  simu_new <- rbind(s1,s2)
	  # 将抽样结果保存到列表中
	  sample_results[[i]] <- simu_new
	}
	# 检查抽样结果的前几个
	head(sample_results[[1]])
	network_result <- data.frame()
	for (i in 1:num_samples) {
	  link_comun <- sample_results[[i]]
	  link_comun$webID_o <- link_comun$webID
	  link_comun$webID <- link_comun$webID_simulation
	  
	  network <- frame2webs(link_comun,type.out= "list")
	  
	  #定量模块观察值，每个网络重复三十次，取最大值
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
	  Wmodul <- repeatmodule(network,30)
	  Wmodulmax <- t(t(apply(Wmodul,2,max)))
	  colnames(Wmodulmax) <- "weighted modularity"
	  network_result2 <- Wmodulmax
	  network_result_all <- rbind(network_result,network_result2)
	  network_result <- network_result_all
	}

	network_result_s1 <- network_result_all[grepl("_s1", rownames(network_result_all)), ]
	network_result_s2 <- network_result_all[grepl("_s2", rownames(network_result_all)), ]

	library(psych)
	# ?describe
	result1 <- describe(network_result_s1)
	result2 <- describe(network_result_s2)
	result_net1 <- list()
	result_net2 <- list()
	result_net1[[site]] <- result1
	result_net2[[site]] <- result2

	network_result_cha <- network_result_s1-network_result_s2
	result <- describe(network_result_cha)

	result_all <- list()
	result_all[[site]] <- result
}

write.csv(result_all,paste0("results/network_simulation_module", site, ".csv"))
# write.csv(result_all,paste0("results/network_simulation_module", site, "_Eup.jo_excluded.csv"))

result_all.all1 <- do.call(rbind,result_net1)
result_all.all2 <- do.call(rbind,result_net2)

result_all.all.all <- cbind(result_all.all1,result_all.all2)

write.csv(result_all.all.all,paste0("results/network_simulation_module", site, "all.csv"))
# write.csv(result_all.all.all,paste0("results/network_simulation_module", site, "all_Eup.jo_excluded.csv"))

save.image(paste0("results/network_simulation_module_", site, ".RData"))
# save.image(paste0("results/network_simulation_module_", site, "_Eup.jo_excluded.RData"))

#################################################################################
#################################################################################
#################################################################################