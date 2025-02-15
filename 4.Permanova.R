# Beta diversity
# setwd('E:/Shangri-La_2024/network')

# 加载包
library(readxl)
library(tidyverse)
library(vegan)
library(ggplot2)

library(zCompositions)
library(CoDaSeq)


# Plant beta diversity
rm(list=ls())

set.seed(123)
plant_comm <- read.csv("plant_comm.csv")
plant_comm <- as.data.frame(plant_comm)
head(plant_comm)
sum(plant_comm[,7:dim(plant_comm)[2]])

comm <- plant_comm[,6:dim(plant_comm)[2]]
comm <- as.data.frame(comm)
row.names(comm) <- comm$CommID
comm$CommID <- NULL

sum(comm)
rowSums(comm)
min(rowSums(comm))
max(rowSums(comm))

# group <- paste0(plant_comm$Site,plant_comm$Plot)
group <- paste0(plant_comm$Plot)
group
rownames(comm)
dis <- vegdist(comm, method="bray")

# perm_anova <- adonis2(dis~plant_comm$Plot,permutations=999)
# perm_anova

# perm_anova <- adonis2(dis~plant_comm$Site/plant_comm$Plot,permutations=999,by="terms")
# perm_anova


perm_anova <- adonis2(dis~plant_comm$Site:plant_comm$Plot,permutations=999,by="terms")
perm_anova


# nested, strata
perm_anova <- adonis2(dis~plant_comm$Site:plant_comm$Plot,permutations=999,strata=plant_comm$Site,by="terms")
perm_anova


########
dispersion <- betadisper(dis, group=group)
dispersion
permutest(dispersion)

plot(dispersion, hull=FALSE, ellipse=TRUE) ##sd ellipse
boxplot(dispersion)
anova(dispersion)
permutest(dispersion)
TukeyHSD(dispersion)

################################################################################
# AMF beta diversity
rm(list=ls())

set.seed(123)
otutabe_o <- read_excel("otutabe.xlsx", sheet="Sheet1")
otutabe_o <- as.data.frame(otutabe_o)
head(otutabe_o)
sum(otutabe_o[,7:dim(otutabe_o)[2]])

# otutabe <- otutabe_o[,c(2:4,6,7:dim(otutabe_o)[2])]
otutabe <- otutabe_o %>% group_by(Site,Plot,PlotID,CommID) %>% summarise(across(7:dim(otutabe_o)[2]-4, sum, na.rm = TRUE))
sum(otutabe[,5:dim(otutabe)[2]])

otutabe.plot <- otutabe[,4:dim(otutabe)[2]]
comm <- otutabe.plot
comm <- as.data.frame(comm)
row.names(comm) <- comm$CommID
comm$CommID <- NULL

sum(comm)
rowSums(comm)
min(rowSums(comm))
max(rowSums(comm))

# group <- paste0(otutabe$Site,otutabe$Plot)
group <- paste0(otutabe$Plot)
group
rownames(comm)
dis <- vegdist(comm, method="bray")

# perm_anova <- adonis2(dis~otutabe$Plot,permutations=999)
# perm_anova

# perm_anova <- adonis2(dis~otutabe$Site/otutabe$Plot,permutations=999,by="terms")
# perm_anova


perm_anova <- adonis2(dis~otutabe$Site:otutabe$Plot,permutations=999,by="terms")
perm_anova


# nested, strata
perm_anova <- adonis2(dis~otutabe$Site:otutabe$Plot,permutations=999,strata=otutabe$Site,by="terms")
perm_anova


########
dispersion <- betadisper(dis, group=group)
dispersion
permutest(dispersion)

plot(dispersion, hull=FALSE, ellipse=TRUE) ##sd ellipse
boxplot(dispersion)
anova(dispersion)
permutest(dispersion)
TukeyHSD(dispersion)


#################################################################################
#################################################################################
#################################################################################