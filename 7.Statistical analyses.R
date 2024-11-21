# Data statistical analyses & figure plots
# setwd('E:/Shangri-La_2024/network')

library(readxl)
library(tidyverse)
library(vegan)
library(ggplot2)
library(glmmTMB)
library(MuMIn)
library(effects)

rm(list=ls())

######################
# Fig. 3 Plant_AMF.top10

data <- read_excel("source_data.xlsx", sheet="Plant_AMF.top10")
head(data)

data$Plant_Species <- as.factor(data$Plant_Species)
data$Plant_Species_abbreviation <- as.factor(data$Plant_Species_abbreviation)
data$Plot <- as.factor(data$Plot)
data$Type <- as.factor(data$Type)
head(data)

#####
data.plant <- subset(data, Type=="Plant_relative_abundance")
data.plant.Euadjacent <- subset(data.plant, Plot=="Eu-adjacent")
data.plant.Eudistant <- subset(data.plant, Plot=="Eu-distant")

p <- ggplot(data.plant.Euadjacent,aes(Plant_Species_abbreviation,Value))+
  geom_bar(aes(fill=Plot),stat = "identity",width=0.8)+
  facet_grid(Plot~.)+
  theme_classic()+
  scale_color_manual(values=c("gray50"))+
  scale_fill_manual(values=c("gray50"))+
  scale_y_continuous(labels=scales::percent,limits=c(0,0.08))+
  ylab("Plant relative abundance")+# 设置Y轴标题
  theme(axis.title.x=element_blank(),#修改X轴标题文本
        # axis.title.x=element_text(size=8),#修改X轴标题文本
        axis.title.y=element_text(size=8,angle=90),#修改y轴标题文本
        axis.text.y=element_text(size=6.5),#修改x轴刻度标签文本
        # axis.line.x = element_line(arrow = arrow(length = unit(0.2,"lines"))),
        axis.text.x=element_text(angle=30, size=6.5, hjust=0.5),#修改y轴刻度标签文本
        panel.grid=element_blank(),#隐藏网格线
        legend.title=element_text(size=6.5),#修改图例标题大小
        legend.text=element_text(size=5.5),#修改图例文本大小
        legend.key.size=unit(0.1,"lines"),#缩小图例符号的大小
        legend.background=element_blank(),#设置图例背景为空白
        legend.position="inside",#修改图例位置
        # legend.justification=c("right","top"),#修改图例位置
        legend.position.inside=c(0.95,0.9),#修改图例位置
        legend.margin=margin(0,0,0,0),#修改图例位置
        axis.line.x=element_line(),
        strip.text = element_blank(), strip.background = element_blank())  # 取消分面标签)
p
ggsave(filename="Fig3a_plant_top10_Euadjacent.png",plot=p,width=89*1.3,height=59*1.3/2,units="mm",dpi=900)
ggsave(filename="Fig3a_plant_top10_Euadjacent.pdf",plot=p,width=89*1.3,height=59*1.3/2,units="mm")

p <- ggplot(data.plant.Eudistant,aes(Plant_Species_abbreviation,Value))+
  geom_bar(aes(fill=Plot),stat = "identity",width=0.8)+
  facet_grid(Plot~.)+
  theme_classic()+
  scale_color_manual(values=c("gray50"))+
  scale_fill_manual(values=c("gray50"))+
  scale_y_continuous(labels=scales::percent)+
  ylab("Plant relative abundance")+# 设置Y轴标题
  theme(axis.title.x=element_blank(),#修改X轴标题文本
        # axis.title.x=element_text(size=8),#修改X轴标题文本
        axis.title.y=element_text(size=8,angle=90),#修改y轴标题文本
        axis.text.y=element_text(size=6.5),#修改x轴刻度标签文本
        # axis.line.x = element_line(arrow = arrow(length = unit(0.2,"lines"))),
        axis.text.x=element_text(angle=30, size=6.5, hjust=0.5),#修改y轴刻度标签文本
        panel.grid=element_blank(),#隐藏网格线
        legend.title=element_text(size=6.5),#修改图例标题大小
        legend.text=element_text(size=5.5),#修改图例文本大小
        legend.key.size=unit(0.1,"lines"),#缩小图例符号的大小
        legend.background=element_blank(),#设置图例背景为空白
        legend.position="inside",#修改图例位置
        # legend.justification=c("right","top"),#修改图例位置
        legend.position.inside=c(0.95,0.9),#修改图例位置
        legend.margin=margin(0,0,0,0),#修改图例位置
        axis.line.x=element_line(),
        strip.text = element_blank(), strip.background = element_blank())  # 取消分面标签)
p
ggsave(filename="Fig3b_plant_top10_Eudistant.png",plot=p,width=89*1.3,height=59*1.3/2,units="mm",dpi=900)
ggsave(filename="Fig3b_plant_top10_Eudistant.pdf",plot=p,width=89*1.3,height=59*1.3/2,units="mm")

#####
data.fungi <- subset(data, Type=="AMF_richness")
data.fungi.Euadjacent <- subset(data.fungi, Plot=="Eu-adjacent")
data.fungi.Eudistant <- subset(data.fungi, Plot=="Eu-distant")

p <- ggplot(data.fungi.Euadjacent, aes(x = Plant_Species_abbreviation, y = Value)) +
  geom_bar(aes(fill=Plot),stat = "identity", position = position_dodge(width = 0.8)) +
  geom_text(aes(label=paste0("n=",Sample_count)),size=1.5,vjust=0)+
  facet_grid(Plot ~ .) +
  theme_classic()+
  scale_color_manual(values=c("gray50"))+
  scale_fill_manual(values=c("gray50"))+
  ylab("AMF richness")+# 设置Y轴标题
  theme(axis.title.x=element_blank(),#修改X轴标题文本
        # axis.title.x=element_text(size=8),#修改X轴标题文本
        axis.title.y=element_text(size=8,angle=90),#修改y轴标题文本
        axis.text.y=element_text(size=6.5),#修改x轴刻度标签文本
        # axis.line.x = element_line(arrow = arrow(length = unit(0.2,"lines"))),
        axis.text.x=element_text(angle=30, size=6.5, hjust=0.5),#修改y轴刻度标签文本
        panel.grid=element_blank(),#隐藏网格线
        legend.title=element_text(size=6.5),#修改图例标题大小
        legend.text=element_text(size=5.5),#修改图例文本大小
        legend.key.size=unit(0.1,"lines"),#缩小图例符号的大小
        legend.background=element_blank(),#设置图例背景为空白
        legend.position="inside",#修改图例位置
        # legend.justification=c("right","top"),#修改图例位置
        legend.position.inside=c(0.95,0.9),#修改图例位置
        legend.margin=margin(0,0,0,0),#修改图例位置
        axis.line.x=element_line(),
        strip.text = element_blank(), strip.background = element_blank())  # 取消分面标签)
p
ggsave(filename="Fig3c_AMF_top10_Euadjacent.png",plot=p,width=89*1.3,height=59*1.3/2,units="mm",dpi=900)
ggsave(filename="Fig3c_AMF_top10_Euadjacent.pdf",plot=p,width=89*1.3,height=59*1.3/2,units="mm")


p <- ggplot(data.fungi.Eudistant, aes(x = Plant_Species_abbreviation, y = Value)) +
  geom_bar(aes(fill=Plot),stat = "identity", position = position_dodge(width = 0.8)) +
  geom_text(aes(label=paste0("n=",Sample_count)),size=1.5,vjust=0)+
  facet_grid(Plot ~ .) +
  theme_classic()+
  scale_color_manual(values=c("gray50"))+
  scale_fill_manual(values=c("gray50"))+
  ylab("AMF richness")+# 设置Y轴标题
  theme(axis.title.x=element_blank(),#修改X轴标题文本
        # axis.title.x=element_text(size=8),#修改X轴标题文本
        axis.title.y=element_text(size=8,angle=90),#修改y轴标题文本
        axis.text.y=element_text(size=6.5),#修改x轴刻度标签文本
        # axis.line.x = element_line(arrow = arrow(length = unit(0.2,"lines"))),
        axis.text.x=element_text(angle=30, size=6.5, hjust=0.5),#修改y轴刻度标签文本
        panel.grid=element_blank(),#隐藏网格线
        legend.title=element_text(size=6.5),#修改图例标题大小
        legend.text=element_text(size=5.5),#修改图例文本大小
        legend.key.size=unit(0.1,"lines"),#缩小图例符号的大小
        legend.background=element_blank(),#设置图例背景为空白
        legend.position="inside",#修改图例位置
        # legend.justification=c("right","top"),#修改图例位置
        legend.position.inside=c(0.95,0.9),#修改图例位置
        legend.margin=margin(0,0,0,0),#修改图例位置
        axis.line.x=element_line(),
        strip.text = element_blank(), strip.background = element_blank())  # 取消分面标签)
p
ggsave(filename="Fig3d_AMF_top10_Eudistant.png",plot=p,width=89*1.3,height=59*1.3/2,units="mm",dpi=900)
ggsave(filename="Fig3d_AMF_top10_Eudistant.pdf",plot=p,width=89*1.3,height=59*1.3/2,units="mm")

################################################################################

# Fig. 4 Species_richness_interactions
rm(list=ls())

richness <- read_excel("source_data.xlsx", sheet="Species_richness_interactions")
head(richness)

richness$Site <- as.factor(richness$Site)
richness$Plot <- as.factor(richness$Plot)
richness$Plot <- as.factor(richness$Plot)
summary(richness)

# Plant_richness 比较Plot间植物richness是在统计上否有显著性差异
{
  Euadjacent.plot <- subset(richness,Plot=="Eu-adjacent")
  Eudistant.plot <- subset(richness,Plot=="Eu-distant")
  Euadjacent.plot$PlotID==Eudistant.plot$PlotID
  x <- Euadjacent.plot$Plant_richness
  y <- Eudistant.plot$Plant_richness
  t.test(x,y,paired=T)
  mean(x);mean(y)
  sd(x);sd(y)
  var(x);var(y)
}

# modlmer.TMB0 <- glmmTMB(Plant_richness ~ 1 + (1|Site), family=poisson, data=richness)
modlmer.TMB1 <- glmmTMB(Plant_richness ~ Plot + (1|Site), family=poisson, data=richness,control=glmmTMBControl(optimizer=optim, optArgs=list(method="CG")))
# (msAICc <- model.sel(modlmer.TMB0,modlmer.TMB1))

summary(modlmer.TMB1)
performance::r2(modlmer.TMB1)
allEffects(modlmer.TMB1)
plot(allEffects(modlmer.TMB1))

p <- ggplot(richness,aes(x=Plot,y=Plant_richness,color=Plot))+
  geom_jitter(shape=16,height=0.25,width=0.15,size=0.4,alpha=0.4)+
  geom_point(aes(x=Plot,y=mean(Plant_richness),color="black"),size=0.6)+
  geom_errorbar(aes(ymin=mean(Plant_richness)-sd(Plant_richness), ymax=mean(Plant_richness)+sd(Plant_richness),color="black"),width=0,linewidth=0.4)+
  scale_color_manual(values=c("#bc40bf","#bc40bf","black"))+
  scale_fill_manual(values=c("#bc40bf","#bc40bf"))+ 
  theme_classic()+
  ylab("Plant richness")+
  # scale_x_discrete(breaks=c("Eu-adjacent","Eu-distant"),labels=c("Eu-adjacent", "Eu-distant"))+
  theme(axis.title.x=element_blank(),#修改X轴标题文本
        # axis.title.x=element_text(size=8),#修改X轴标题文本
        axis.line.x = element_line(arrow = arrow(length = unit(0.2,"lines"))),
        axis.title.y=element_text(size=8,angle=90),#修改y轴标题文本
        axis.text.y=element_text(size=6.5),#修改x轴刻度标签文本
        axis.text.x=element_text(size=6.5),#修改y轴刻度标签文本
        panel.grid=element_blank(),#隐藏网格线
        legend.title=element_text(size=8),#修改图例标题大小
        legend.text=element_text(size=6.5),#修改图例文本大小
        legend.position="none")+#修改图例位置
  annotate("text",x=1.7,y=max(richness$Plant_richness),hjust=1,vjust=1,label=paste0("p > ",0.05),size=2)
p

ggsave(filename="Fig4a_plant_richness_glmm.png",plot=p,width=89/2,height=59/2,units="mm",dpi=900)
ggsave(filename="Fig4a_plant_richness_glmm.pdf",plot=p,width=89/2,height=59/2,units="mm")

#####
# Biomass of Plant
# Aboveground_biomass
{
  Euadjacent.plot <- subset(richness,Plot=="Eu-adjacent")
  Eudistant.plot <- subset(richness,Plot=="Eu-distant")
  Euadjacent.plot$PlotID==Eudistant.plot$PlotID
  x <- Euadjacent.plot$Aboveground_biomass
  y <- Eudistant.plot$Aboveground_biomass
  t.test(x,y,paired=T)
  mean(x);mean(y)
  sd(x);sd(y)
  var(x);var(y)
}
modlmer0 <- glmmTMB(Aboveground_biomass~1+(1|Site),family=Gamma(link="log"),data=richness)
modlmer1 <- glmmTMB(Aboveground_biomass~Plot+(1|Site),family=Gamma(link="log"),data=richness)
modlmer2 <- glmmTMB(Aboveground_biomass~Plot*Plant_richness+(1|Site),family=Gamma(link="log"),data=richness)
modlmer3 <- glmmTMB(Aboveground_biomass~Plot*AMF_richness+(1|Site),family=Gamma(link="log"),data=richness)
modlmer4 <- glmmTMB(Aboveground_biomass~Plot*Interactions_numbers+(1|Site),family=Gamma(link="log"),data=richness)
modlmer5 <- glmmTMB(Aboveground_biomass~Plot+Plant_richness*AMF_richness+(1|Site),family=Gamma(link="log"),data=richness)
anova(modlmer0,modlmer1,modlmer2,modlmer3,modlmer4,modlmer5)
anova(modlmer1,modlmer2,modlmer3,modlmer4,modlmer5)
anova(modlmer1,modlmer3)
summary(modlmer1)
performance::r2(modlmer1)
allEffects(modlmer1)
plot(allEffects(modlmer1))

# Belowground_biomass
{
  Euadjacent.plot <- subset(richness,Plot=="Eu-adjacent")
  Eudistant.plot <- subset(richness,Plot=="Eu-distant")
  Euadjacent.plot$PlotID==Eudistant.plot$PlotID
  x <- Euadjacent.plot$Belowground_biomass
  y <- Eudistant.plot$Belowground_biomass
  t.test(x,y,paired=T)
  mean(x);mean(y)
  sd(x);sd(y)
  var(x);var(y)
}
modlmer0 <- glmmTMB(Belowground_biomass~1+(1|Site),family=Gamma(link="log"),data=richness)
modlmer1 <- glmmTMB(Belowground_biomass~Plot+(1|Site),family=Gamma(link="log"),data=richness)
modlmer2 <- glmmTMB(Belowground_biomass~Plot*Plant_richness+(1|Site),family=Gamma(link="log"),data=richness)
modlmer3 <- glmmTMB(Belowground_biomass~Plot*AMF_richness+(1|Site),family=Gamma(link="log"),data=richness)
modlmer4 <- glmmTMB(Belowground_biomass~Plot*Interactions_numbers+(1|Site),family=Gamma(link="log"),data=richness)
modlmer5 <- glmmTMB(Belowground_biomass~Plot+Plant_richness*AMF_richness+(1|Site),family=Gamma(link="log"),data=richness)
anova(modlmer0,modlmer1,modlmer2,modlmer3,modlmer4,modlmer5)
anova(modlmer1,modlmer2,modlmer3,modlmer4,modlmer5)
anova(modlmer1,modlmer3)
summary(modlmer1)
performance::r2(modlmer1)
allEffects(modlmer1)
plot(allEffects(modlmer1))


biomass <- richness[,c("Site","Plot","PlotID","Aboveground_biomass","Belowground_biomass")]
biomass_long <- pivot_longer(biomass, cols = c(Aboveground_biomass, Belowground_biomass), names_to = "Type", values_to = "biomass")

biomass_long[which(biomass_long$Type=="Aboveground_biomass"),]$Type <- "Above"
biomass_long[which(biomass_long$Type=="Belowground_biomass"),]$Type <- "Below"

# 计算均值和标准差
biomass_long.summary <- biomass_long %>%
  group_by(Plot, Type) %>%
  summarise(biomass.mean = mean(biomass),biomass.sd = sd(biomass))

p <- ggplot(biomass_long,aes(x=Plot,y=biomass,color=Type,group=Type))+
  geom_jitter(shape=16,size=0.4,position=position_jitterdodge(jitter.height=0,jitter.width=0.25,dodge.width=0.5),alpha=0.4,show.legend=T)+
  geom_point(data=biomass_long.summary,aes(x=Plot,y=biomass.mean,group=Type),shape=16,size=0.6,color="black",position=position_dodge(0.5),alpha=0.6)+
  geom_errorbar(data=biomass_long.summary,aes(x=Plot,y=biomass.mean,ymin=biomass.mean-biomass.sd, ymax=biomass.mean+biomass.sd,group=Type),width=0,linewidth=0.2,color="black",position=position_dodge(0.5),alpha=0.6)+
  theme_classic()+
  scale_color_manual(values=c("#bc40bf","#4BA465","black"))+
  scale_fill_manual(values=c("#bc40bf","#4BA465","black"))+
  ylab("Biomass of plant")+
  # scale_x_discrete(breaks=c("Eu-adjacent","Eu-distant"),labels=c("Eu-adjacent", "Eu-distant"))+
  theme(axis.title.x=element_blank(),#修改X轴标题文本
        # axis.title.x=element_text(size=8),#修改X轴标题文本
        axis.title.y=element_text(size=8,angle=90),#修改y轴标题文本
        axis.text.y=element_text(size=6.5),#修改x轴刻度标签文本
        axis.line.x = element_line(arrow = arrow(length = unit(0.2,"lines"))),
        axis.text.x=element_text(size=6.5),#修改y轴刻度标签文本
        panel.grid=element_blank(),#隐藏网格线
        # legend.title=element_text(size=6.5),#修改图例标题大小
        legend.title=element_blank(),#修改图例标题大小
        legend.text=element_text(size=5.5),#修改图例文本大小
        legend.key.size=unit(0.1,"lines"),#缩小图例符号的大小
        legend.background=element_blank(),#设置图例背景为空白
        legend.position="inside",#修改图例位置
        # legend.justification=c("right","top"),#修改图例位置
        legend.position.inside=c(0.2,0.9),#修改图例位置
        legend.margin=margin(0,0,0,0))+#修改图例位置
  annotate("text",x=1.9,y=max(biomass_long$biomass)*0.85,hjust=1,vjust=1,label=paste0("p < ",0.001),size=2,color="black")
p

ggsave(filename="Fig4b_biomass_glmm.png",plot=p,width=89/2,height=59/2,units="mm",dpi=900)
ggsave(filename="Fig4b_biomass_glmm.pdf",plot=p,width=89/2,height=59/2,units="mm")

#####
# AMF_richness 比较Plot间AMF richness是在统计上否有显著性差异
{
  Euadjacent.plot <- subset(richness,Plot=="Eu-adjacent")
  Eudistant.plot <- subset(richness,Plot=="Eu-distant")
  Euadjacent.plot$PlotID==Eudistant.plot$PlotID
  x <- Euadjacent.plot$AMF_richness
  y <- Eudistant.plot$AMF_richness
  t.test(x,y,paired=T)
  mean(x);mean(y)
  sd(x);sd(y)
  var(x);var(y)
}

modlmer.TMB0 <- glmmTMB(AMF_richness ~ 1 + (1|Site), family=poisson, data=richness)
modlmer.TMB1 <- glmmTMB(AMF_richness ~ Plot + (1|Site), family=poisson, data=richness)
modlmer.TMB2 <- glmmTMB(AMF_richness ~ Plot + (1|Site), family=nbinom2, data=richness)
modlmer.TMB3 <- glmmTMB(AMF_richness ~ Plot + (1|Site), family=nbinom1, data=richness)
(msAICc <- model.sel(modlmer.TMB0,modlmer.TMB1,modlmer.TMB2,modlmer.TMB3))
summary(modlmer.TMB3)
performance::r2(modlmer.TMB3)
allEffects(modlmer.TMB3)
plot(allEffects(modlmer.TMB3))

{
  Euadjacent.plot <- subset(richness,Plot=="Eu-adjacent")
  Eudistant.plot <- subset(richness,Plot=="Eu-distant")
  Euadjacent.plot$PlotID==Eudistant.plot$PlotID
  x <- Euadjacent.plot$AMF_richness_Eup.jo_excluded
  y <- Eudistant.plot$AMF_richness_Eup.jo_excluded
  t.test(x,y,paired=T)
  mean(x);mean(y)
  sd(x);sd(y)
  var(x);var(y)
}
modlmer.TMB0 <- glmmTMB(AMF_richness_Eup.jo_excluded ~ 1 + (1|Site), family=poisson, data=richness)
modlmer.TMB1 <- glmmTMB(AMF_richness_Eup.jo_excluded ~ Plot + (1|Site), family=poisson, data=richness)
modlmer.TMB2 <- glmmTMB(AMF_richness_Eup.jo_excluded ~ Plot + (1|Site), family=nbinom2, data=richness)
modlmer.TMB3 <- glmmTMB(AMF_richness_Eup.jo_excluded ~ Plot + (1|Site), family=nbinom1, data=richness)
(msAICc <- model.sel(modlmer.TMB0,modlmer.TMB1,modlmer.TMB2,modlmer.TMB3))
summary(modlmer.TMB3)
performance::r2(modlmer.TMB3)
allEffects(modlmer.TMB3)
plot(allEffects(modlmer.TMB3))


richness.AMF <- richness[,c("Site","Plot","PlotID","AMF_richness","AMF_richness_Eup.jo_excluded")]
richness.AMF_long <- pivot_longer(richness.AMF, cols = c(AMF_richness, AMF_richness_Eup.jo_excluded), names_to = "Type", values_to = "AMF_richness")

richness.AMF_long[which(richness.AMF_long$Type=="AMF_richness"),]$Type <- "included"
richness.AMF_long[which(richness.AMF_long$Type=="AMF_richness_Eup.jo_excluded"),]$Type <- "excluded"
richness.AMF_long$Type <- factor(richness.AMF_long$Type,levels=c("included","excluded"))

# 计算均值和标准差
richness.AMF_long.summary <- richness.AMF_long %>%
  group_by(Plot, Type) %>%
  summarise(AMF_richness.mean = mean(AMF_richness),AMF_richness.sd = sd(AMF_richness))

p <- ggplot(richness.AMF_long,aes(x=Plot,y=AMF_richness,color=Type,group=Type))+
  geom_jitter(shape=16,size=0.4,position=position_jitterdodge(jitter.height=0,jitter.width=0.25,dodge.width=0.5),alpha=0.4,show.legend=T)+
  geom_point(data=richness.AMF_long.summary,aes(x=Plot,y=AMF_richness.mean,group=Type),shape=16,size=0.6,color="black",position=position_dodge(0.5),alpha=0.6)+
  geom_errorbar(data=richness.AMF_long.summary,aes(x=Plot,y=AMF_richness.mean,ymin=AMF_richness.mean-AMF_richness.sd, ymax=AMF_richness.mean+AMF_richness.sd,group=Type),width=0,linewidth=0.2,color="black",position=position_dodge(0.5),alpha=0.6)+
  theme_classic()+
  scale_color_manual(values=c("#bc40bf","#4BA465","black"))+
  scale_fill_manual(values=c("#bc40bf","#4BA465","black"))+
  ylab("AMF richness")+
  # scale_x_discrete(breaks=c("Eu-adjacent","Eu-distant"),labels=c("Eu-adjacent", "Eu-distant"))+
  theme(axis.title.x=element_blank(),#修改X轴标题文本
        # axis.title.x=element_text(size=8),#修改X轴标题文本
        axis.title.y=element_text(size=8,angle=90),#修改y轴标题文本
        axis.text.y=element_text(size=6.5),#修改x轴刻度标签文本
        axis.line.x = element_line(arrow = arrow(length = unit(0.2,"lines"))),
        axis.text.x=element_text(size=6.5),#修改y轴刻度标签文本
        panel.grid=element_blank(),#隐藏网格线
        # legend.title=element_text(size=6.5),#修改图例标题大小
        legend.title=element_blank(),#修改图例标题大小
        legend.text=element_text(size=5.5),#修改图例文本大小
        legend.key.size=unit(0.1,"lines"),#缩小图例符号的大小
        legend.background=element_blank(),#设置图例背景为空白
        legend.position="inside",#修改图例位置
        # legend.justification=c("right","top"),#修改图例位置
        legend.position.inside=c(0.9,0.9),#修改图例位置
        legend.margin=margin(0,0,0,0))+#修改图例位置
  annotate("text",x=1.9,y=max(richness.AMF_long$AMF_richness)*0.85,hjust=1,vjust=1,label=paste0("p > ",0.05),size=2,color="black")
p

ggsave(filename="Fig4c_AMF_richness_glmm.png",plot=p,width=89/2,height=59/2,units="mm",dpi=900)
ggsave(filename="Fig4c_AMF_richness_glmm.pdf",plot=p,width=89/2,height=59/2,units="mm")

#####
# Interaction_number 比较Plot间Interaction number是在统计上否有显著性差异
{
  Euadjacent.plot <- subset(richness,Plot=="Eu-adjacent")
  Eudistant.plot <- subset(richness,Plot=="Eu-distant")
  Euadjacent.plot$PlotID==Eudistant.plot$PlotID
  x <- Euadjacent.plot$Interaction_number
  y <- Eudistant.plot$Interaction_number
  t.test(x,y,paired=T)
  mean(x);mean(y)
  sd(x);sd(y)
  var(x);var(y)
}

modlmer.TMB0 <- glmmTMB(Interaction_number ~ 1 + (1|Site), family=poisson, data=richness)
modlmer.TMB1 <- glmmTMB(Interaction_number ~ Plot + (1|Site), family=poisson, data=richness)
modlmer.TMB2 <- glmmTMB(Interaction_number ~ Plot + (1|Site), family=nbinom2, data=richness)
modlmer.TMB3 <- glmmTMB(Interaction_number ~ Plot + (1|Site), family=nbinom1, data=richness)
(msAICc <- model.sel(modlmer.TMB0,modlmer.TMB1,modlmer.TMB2,modlmer.TMB3))
summary(modlmer.TMB3)
performance::r2(modlmer.TMB3)
allEffects(modlmer.TMB3)
plot(allEffects(modlmer.TMB3))

{
  Euadjacent.plot <- subset(richness,Plot=="Eu-adjacent")
  Eudistant.plot <- subset(richness,Plot=="Eu-distant")
  Euadjacent.plot$PlotID==Eudistant.plot$PlotID
  x <- Euadjacent.plot$Interaction_number_Eup.jo_excluded
  y <- Eudistant.plot$Interaction_number_Eup.jo_excluded
  t.test(x,y,paired=T)
  mean(x);mean(y)
  sd(x);sd(y)
  var(x);var(y)
}
modlmer.TMB0 <- glmmTMB(Interaction_number_Eup.jo_excluded ~ 1 + (1|Site), family=poisson, data=richness)
modlmer.TMB1 <- glmmTMB(Interaction_number_Eup.jo_excluded ~ Plot + (1|Site), family=poisson, data=richness)
modlmer.TMB2 <- glmmTMB(Interaction_number_Eup.jo_excluded ~ Plot + (1|Site), family=nbinom2, data=richness)
modlmer.TMB3 <- glmmTMB(Interaction_number_Eup.jo_excluded ~ Plot + (1|Site), family=nbinom1, data=richness)
(msAICc <- model.sel(modlmer.TMB0,modlmer.TMB1,modlmer.TMB2,modlmer.TMB3))
summary(modlmer.TMB3)
performance::r2(modlmer.TMB3)
allEffects(modlmer.TMB3)
plot(allEffects(modlmer.TMB3))


richness.interaction <- richness[,c("Site","Plot","PlotID","Interaction_number","Interaction_number_Eup.jo_excluded")]
richness.interaction_long <- pivot_longer(richness.interaction, cols = c(Interaction_number, Interaction_number_Eup.jo_excluded), names_to = "Type", values_to = "Interaction_number")

richness.interaction_long[which(richness.interaction_long$Type=="Interaction_number"),]$Type <- "included"
richness.interaction_long[which(richness.interaction_long$Type=="Interaction_number_Eup.jo_excluded"),]$Type <- "excluded"
richness.interaction_long$Type <- factor(richness.interaction_long$Type,levels=c("included","excluded"))

# 计算均值和标准差
richness.interaction_long.summary <- richness.interaction_long %>%
  group_by(Plot, Type) %>%
  summarise(Interaction_number.mean = mean(Interaction_number),Interaction_number.sd = sd(Interaction_number))

p <- ggplot(richness.interaction_long,aes(x=Plot,y=Interaction_number,color=Type,group=Type))+
  geom_jitter(shape=16,size=0.4,position=position_jitterdodge(jitter.height=0,jitter.width=0.25,dodge.width=0.5),alpha=0.4,show.legend=T)+
  geom_point(data=richness.interaction_long.summary,aes(x=Plot,y=Interaction_number.mean,group=Type),shape=16,size=0.6,color="black",position=position_dodge(0.5),alpha=0.6)+
  geom_errorbar(data=richness.interaction_long.summary,aes(x=Plot,y=Interaction_number.mean,ymin=Interaction_number.mean-Interaction_number.sd, ymax=Interaction_number.mean+Interaction_number.sd,group=Type),width=0,linewidth=0.2,color="black",position=position_dodge(0.5),alpha=0.6)+
  theme_classic()+
  scale_color_manual(values=c("#bc40bf","#4BA465","black"))+
  scale_fill_manual(values=c("#bc40bf","#4BA465","black"))+
  ylab("Number of interaction")+
  # scale_x_discrete(breaks=c("Eu-adjacent","Eu-distant"),labels=c("Eu-adjacent", "Eu-distant"))+
  theme(axis.title.x=element_blank(),#修改X轴标题文本
        # axis.title.x=element_text(size=8),#修改X轴标题文本
        axis.title.y=element_text(size=8,angle=90),#修改y轴标题文本
        axis.text.y=element_text(size=6.5),#修改x轴刻度标签文本
        axis.line.x = element_line(arrow = arrow(length = unit(0.2,"lines"))),
        axis.text.x=element_text(size=6.5),#修改y轴刻度标签文本
        panel.grid=element_blank(),#隐藏网格线
        # legend.title=element_text(size=6.5),#修改图例标题大小
        legend.title=element_blank(),#修改图例标题大小
        legend.text=element_text(size=5.5),#修改图例文本大小
        legend.key.size=unit(0.1,"lines"),#缩小图例符号的大小
        legend.background=element_blank(),#设置图例背景为空白
        legend.position="inside",#修改图例位置
        # legend.justification=c("right","top"),#修改图例位置
        legend.position.inside=c(0.9,0.9),#修改图例位置
        legend.margin=margin(0,0,0,0))+#修改图例位置
  annotate("text",x=2,y=max(richness.interaction_long$Interaction_number)*0.85,hjust=1,vjust=1,label=paste0("p > 0.05"," / ","p < 0.01"),size=2,color="black")
p

ggsave(filename="Fig4d_interaction_number_glmm.png",plot=p,width=89/2,height=59/2,units="mm",dpi=900)
ggsave(filename="Fig4d_interaction_number_glmm.pdf",plot=p,width=89/2,height=59/2,units="mm")

################################################################################
# Fig. 5 Beta_diversity
rm(list=ls())

Beta_diversity <- read_excel("source_data.xlsx", sheet="Beta_diversity")
head(Beta_diversity)

Beta_diversity$Site <- as.factor(Beta_diversity$Site)
Beta_diversity$Between_plot <- as.factor(Beta_diversity$Between_plot)
Beta_diversity$Beta_diversity_type <- as.factor(Beta_diversity$Beta_diversity_type)
Beta_diversity$Components <- as.factor(Beta_diversity$Components)
Beta_diversity$Components <- factor(Beta_diversity$Components,levels=c("Total","Turnover","Nestedness","βWN","βOS","βST","βST.h","βST.l","βST.lh"))
summary(Beta_diversity)

plant_beta <- subset(Beta_diversity,Beta_diversity_type=="Plant_beta_diversity")
x <- subset(plant_beta,Components=="Turnover")$Observe
y <- subset(plant_beta,Components=="Nestedness")$Observe
t.test(x,y,paired=T)

# 计算均值和标准差
plant_beta_combined <- plant_beta %>%
  group_by(Beta_diversity_type, Components) %>%
  summarise(Observe.mean = mean(Observe),Observe.sd = sd(Observe))

p <- ggplot(plant_beta,aes(x=Components,y=Observe,color=Components))+
  geom_jitter(shape=16,size=0.4,position=position_jitterdodge(jitter.height=0,jitter.width=0.25,dodge.width=0.5),alpha=0.8,show.legend=T)+
  geom_point(data=plant_beta_combined,aes(x=Components,y=Observe.mean),shape=16,size=0.6,color="black",position=position_dodge(0.5),alpha=0.8)+
  geom_errorbar(data=plant_beta_combined,aes(x=Components,y=Observe.mean,ymin=Observe.mean-Observe.sd, ymax=Observe.mean+Observe.sd),width=0,linewidth=0.2,color="black",position=position_dodge(0.5),alpha=0.6)+
  scale_color_manual(values=c("#bc40bf","#bc40bf","#bc40bf"))+
  scale_fill_manual(values=c("#bc40bf","#bc40bf","#bc40bf"))+
  theme_classic()+
  ylab("Plant beta diversity")+# 设置Y轴标题
  geom_vline(xintercept=1.5,linetype=2,linewidth=0.2,color="gray50",alpha=0.8)+
  theme(axis.title.x=element_blank(),#修改X轴标题文本
        # axis.title.x=element_text(size=8),#修改X轴标题文本
        axis.title.y=element_text(size=8,angle=90),#修改y轴标题文本
        axis.text.y=element_text(size=6.5),#修改x轴刻度标签文本
        axis.text.x=element_text(size=6.5),#修改y轴刻度标签文本
        panel.grid=element_blank(),#隐藏网格线
        legend.title=element_text(size=8),#修改图例标题大小
        legend.text=element_text(size=6.5),#修改图例文本大小
        legend.position="none")+#修改图例位置
  annotate("text",x=3,y=max(plant_beta$Observe)*0.85,hjust=1,vjust=1,label=paste0("p < ",0.001),size=2,color="black")
p
ggsave(filename="Fig5a_plant_beta_diversity.png",plot=p,width=89/2,height=59/2,units="mm",dpi=900)
ggsave(filename="Fig5a_plant_beta_diversity.pdf",plot=p,width=89/2,height=59/2,units="mm")

#####

AMF_beta <- subset(Beta_diversity,Beta_diversity_type=="AMF_beta_diversity")
x <- subset(AMF_beta,Components=="Turnover")$Observe
y <- subset(AMF_beta,Components=="Nestedness")$Observe
t.test(x,y,paired=T)

AMF_beta_long <- pivot_longer(AMF_beta, cols = c(Observe, Observe_Eup.jo_excluded), names_to = "Type", values_to = "Value")
AMF_beta_long[which(AMF_beta_long$Type=="Observe"),]$Type <- "included"
AMF_beta_long[which(AMF_beta_long$Type=="Observe_Eup.jo_excluded"),]$Type <- "excluded"
AMF_beta_long$Type <- factor(AMF_beta_long$Type,levels=c("included","excluded"))

# 计算均值和标准差
AMF_beta_long.summary <- AMF_beta_long %>%
  group_by(Components, Type) %>%
  summarise(Value.mean = mean(Value),Value.sd = sd(Value))

p <- ggplot(AMF_beta_long,aes(x=Components,y=Value,color=Type,group=Type))+
  geom_jitter(shape=16,size=0.4,position=position_jitterdodge(jitter.height=0,jitter.width=0.25,dodge.width=0.5),alpha=0.8,show.legend=T)+
  geom_point(data=AMF_beta_long.summary,aes(x=Components,y=Value.mean,group=Type),shape=16,size=0.6,color="black",position=position_dodge(0.5),alpha=0.8)+
  geom_errorbar(data=AMF_beta_long.summary,aes(x=Components,y=Value.mean,ymin=Value.mean-Value.sd, ymax=Value.mean+Value.sd,group=Type),width=0,linewidth=0.2,color="black",position=position_dodge(0.5),alpha=0.6)+
  scale_color_manual(values=c("#bc40bf","#4BA465","black"))+
  scale_fill_manual(values=c("#bc40bf","#4BA465","black"))+
  theme_classic()+
  ylab("AMF beta diversity")+# 设置Y轴标题
  geom_vline(xintercept=1.5,linetype=2,linewidth=0.2,color="gray50",alpha=0.8)+
  theme(axis.title.x=element_blank(),#修改X轴标题文本
        # axis.title.x=element_text(size=8),#修改X轴标题文本
        axis.title.y=element_text(size=8,angle=90),#修改y轴标题文本
        axis.text.y=element_text(size=6.5),#修改x轴刻度标签文本
        axis.text.x=element_text(size=6.5),#修改y轴刻度标签文本
        panel.grid=element_blank(),#隐藏网格线
        # legend.title=element_text(size=6.5),#修改图例标题大小
        legend.title=element_blank(),
        legend.text=element_text(size=5.5),#修改图例文本大小
        legend.key.size=unit(0.1,"lines"),#缩小图例符号的大小
        legend.background=element_blank(),#设置图例背景为空白
        legend.position="inside",#修改图例位置
        # legend.justification=c("right","top"),#修改图例位置
        legend.position.inside=c(0.92,0.9),#修改图例位置
        legend.margin=margin(0,0,0,0))+#修改图例位置
  annotate("text",x=3,y=max(AMF_beta_long$Value)*0.85,hjust=1,vjust=1,label=paste0("p < ",0.001),size=2,color="black")
p
ggsave(filename="Fig5b_AMF_beta_diversity.png",plot=p,width=89/2,height=59/2,units="mm",dpi=900)
ggsave(filename="Fig5b_AMF_beta_diversity.pdf",plot=p,width=89/2,height=59/2,units="mm")

#####
interaction_beta <- subset(Beta_diversity,Beta_diversity_type=="Interaction_beta_diversity")
x <- subset(interaction_beta,Components=="βOS")$Observe
y <- subset(interaction_beta,Components=="βST")$Observe
t.test(x,y,paired=T)
x <- subset(interaction_beta,Components=="βOS")$Observe_Eup.jo_excluded
y <- subset(interaction_beta,Components=="βST")$Observe_Eup.jo_excluded
t.test(x,y,paired=T)

x <- subset(interaction_beta,Components=="βST.l")$Observe
y <- subset(interaction_beta,Components=="βST.lh")$Observe
t.test(x,y,paired=T)
x <- subset(interaction_beta,Components=="βST.l")$Observe_Eup.jo_excluded
y <- subset(interaction_beta,Components=="βST.lh")$Observe_Eup.jo_excluded
t.test(x,y,paired=T)


interaction_beta_long <- pivot_longer(interaction_beta, cols = c(Observe, Observe_Eup.jo_excluded), names_to = "Type", values_to = "Value")
interaction_beta_long[which(interaction_beta_long$Type=="Observe"),]$Type <- "included"
interaction_beta_long[which(interaction_beta_long$Type=="Observe_Eup.jo_excluded"),]$Type <- "excluded"
interaction_beta_long$Type <- factor(interaction_beta_long$Type,levels=c("included","excluded"))

# 计算均值和标准差
interaction_beta_long.summary <- interaction_beta_long %>%
  group_by(Components, Type) %>%
  summarise(Value.mean = mean(Value),Value.sd = sd(Value))

p <- ggplot(interaction_beta_long,aes(x=Components,y=Value,color=Type,group=Type))+
  geom_jitter(shape=16,size=0.4,position=position_jitterdodge(jitter.height=0,jitter.width=0.25,dodge.width=0.5),alpha=0.8,show.legend=T)+
  geom_point(data=interaction_beta_long.summary,aes(x=Components,y=Value.mean,group=Type),shape=16,size=0.6,color="black",position=position_dodge(0.5),alpha=0.8)+
  geom_errorbar(data=interaction_beta_long.summary,aes(x=Components,y=Value.mean,ymin=Value.mean-Value.sd, ymax=Value.mean+Value.sd,group=Type),width=0,linewidth=0.2,color="black",position=position_dodge(0.5),alpha=0.6)+
  scale_color_manual(values=c("#bc40bf","#4BA465","black"))+
  scale_fill_manual(values=c("#bc40bf","#4BA465","black"))+
  theme_classic()+
  ylab("Interaction beta diversity")+# 设置Y轴标题
  geom_vline(xintercept=1.5,linetype=2,linewidth=0.2,color="gray50",alpha=0.8)+
  geom_vline(xintercept=3.5,linetype=2,linewidth=0.2,color="gray50",alpha=0.8)+
  theme(axis.title.x=element_blank(),#修改X轴标题文本
        # axis.title.x=element_text(size=8),#修改X轴标题文本
        axis.title.y=element_text(size=8,angle=90),#修改y轴标题文本
        axis.text.y=element_text(size=6.5),#修改x轴刻度标签文本
        axis.text.x=element_text(size=6.5),#修改y轴刻度标签文本
        panel.grid=element_blank(),#隐藏网格线
        # legend.title=element_text(size=6.5),#修改图例标题大小
        legend.title=element_blank(),
        legend.text=element_text(size=5.5),#修改图例文本大小
        legend.key.size=unit(0.1,"lines"),#缩小图例符号的大小
        legend.background=element_blank(),#设置图例背景为空白
        legend.position="inside",#修改图例位置
        # legend.justification=c("right","top"),#修改图例位置
        legend.position.inside=c(0.92,0.9),#修改图例位置
        legend.margin=margin(0,0,0,0))+#修改图例位置
  annotate("text",x=2.6,y=max(interaction_beta_long$Value)*0.85,hjust=1,vjust=1,label=paste0("p > ",0.05),size=2,color="black")+
  annotate("text",x=5.8,y=max(interaction_beta_long$Value)*0.85/4,hjust=1,vjust=1,label=paste0("p < ",0.001),size=2,color="black")
p
ggsave(filename="Fig5c_interaction_beta_diversity.png",plot=p,width=89,height=59/2,units="mm",dpi=900)
ggsave(filename="Fig5c_interaction_beta_diversity.pdf",plot=p,width=89,height=59/2,units="mm")

################################################################################
# Fig. 6 Network index
rm(list=ls())

network_index <- read_excel("source_data.xlsx", sheet="Network_index")
head(network_index)

network_index$Site <- as.factor(network_index$Site)
network_index$Plot <- as.factor(network_index$Plot)
network_index$webID <- as.factor(network_index$webID)
network_index$Index_type <- as.factor(network_index$Index_type)
summary(network_index)


WCON <- subset(network_index,Index_type=="Weighted connectance")
x <- subset(WCON,Plot=="Eu-adjacent")$Observe
y <- subset(WCON,Plot=="Eu-distant")$Observe
t.test(x,y,paired=T)

x <- subset(WCON,Plot=="Eu-adjacent")$Observe_Eup.jo_excluded
y <- subset(WCON,Plot=="Eu-distant")$Observe_Eup.jo_excluded
t.test(x,y,paired=T)

WCON_long <- pivot_longer(WCON, cols = c(Observe, Observe_Eup.jo_excluded), names_to = "Type", values_to = "Value")
WCON_long[which(WCON_long$Type=="Observe"),]$Type <- "included"
WCON_long[which(WCON_long$Type=="Observe_Eup.jo_excluded"),]$Type <- "excluded"
WCON_long$Type <- factor(WCON_long$Type,levels=c("included","excluded"))

# 计算均值和标准差
WCON_long.summary <- WCON_long %>%
  group_by(Plot, Type) %>%
  summarise(Value.mean = mean(Value),Value.sd = sd(Value))

p <- ggplot(WCON_long,aes(x=Plot,y=Value,color=Type,group=Type))+
  geom_jitter(shape=16,size=0.4,position=position_jitterdodge(jitter.height=0,jitter.width=0.25,dodge.width=0.5),alpha=0.8,show.legend=T)+
  geom_point(data=WCON_long.summary,aes(x=Plot,y=Value.mean,group=Type),shape=16,size=0.6,color="black",position=position_dodge(0.5),alpha=0.8)+
  geom_errorbar(data=WCON_long.summary,aes(x=Plot,y=Value.mean,ymin=Value.mean-Value.sd, ymax=Value.mean+Value.sd,group=Type),width=0,linewidth=0.2,color="black",position=position_dodge(0.5),alpha=0.6)+
  theme_classic()+
  scale_color_manual(values=c("#bc40bf","#4BA465","black"))+
  scale_fill_manual(values=c("#bc40bf","#4BA465","black"))+
  ylab("Weighted connectance")+# 设置Y轴标题
  # scale_x_discrete(breaks=c("Eu-adjacent","Eu-distant"),labels=c("Eu-adjacent", "Eu-distant"))+
  theme(axis.title.x=element_blank(),#修改X轴标题文本
        # axis.title.x=element_text(size=8),#修改X轴标题文本
        axis.title.y=element_text(size=8,angle=90),#修改y轴标题文本
        axis.text.y=element_text(size=6.5),#修改x轴刻度标签文本
        axis.line.x = element_line(arrow = arrow(length = unit(0.2,"lines"))),
        axis.text.x=element_text(size=6.5),#修改y轴刻度标签文本
        panel.grid=element_blank(),#隐藏网格线
        # legend.title=element_text(size=6.5),#修改图例标题大小
        legend.title=element_blank(),#修改图例标题大小
        legend.text=element_text(size=5.5),#修改图例文本大小
        legend.key.size=unit(0.1,"lines"),#缩小图例符号的大小
        legend.background=element_blank(),#设置图例背景为空白
        legend.position="inside",#修改图例位置
        # legend.justification=c("right","top"),#修改图例位置
        legend.position.inside=c(0.9,0.9),#修改图例位置
        legend.margin=margin(0,0,0,0))+#修改图例位置
  annotate("text",x=1.8,y=max(WCON_long$Value),hjust=1,vjust=1,label=paste0("p < ",0.01),size=2,color="black")
p

ggsave(filename="Fig6a_network_index_Weighted_connectance.png",plot=p,width=89/2,height=59/2,units="mm",dpi=900)
ggsave(filename="Fig6a_network_index_Weighted_connectance.pdf",plot=p,width=89/2,height=59/2,units="mm")

#####
Evenness <- subset(network_index,Index_type=="Interaction evenness")
x <- subset(Evenness,Plot=="Eu-adjacent")$Observe
y <- subset(Evenness,Plot=="Eu-distant")$Observe
t.test(x,y,paired=T)

x <- subset(Evenness,Plot=="Eu-adjacent")$Observe_Eup.jo_excluded
y <- subset(Evenness,Plot=="Eu-distant")$Observe_Eup.jo_excluded
t.test(x,y,paired=T)

Evenness_long <- pivot_longer(Evenness, cols = c(Observe, Observe_Eup.jo_excluded), names_to = "Type", values_to = "Value")
Evenness_long[which(Evenness_long$Type=="Observe"),]$Type <- "included"
Evenness_long[which(Evenness_long$Type=="Observe_Eup.jo_excluded"),]$Type <- "excluded"
Evenness_long$Type <- factor(Evenness_long$Type,levels=c("included","excluded"))

# 计算均值和标准差
Evenness_long.summary <- Evenness_long %>%
  group_by(Plot, Type) %>%
  summarise(Value.mean = mean(Value),Value.sd = sd(Value))

p <- ggplot(Evenness_long,aes(x=Plot,y=Value,color=Type,group=Type))+
  geom_jitter(shape=16,size=0.4,position=position_jitterdodge(jitter.height=0,jitter.width=0.25,dodge.width=0.5),alpha=0.8,show.legend=T)+
  geom_point(data=Evenness_long.summary,aes(x=Plot,y=Value.mean,group=Type),shape=16,size=0.6,color="black",position=position_dodge(0.5),alpha=0.8)+
  geom_errorbar(data=Evenness_long.summary,aes(x=Plot,y=Value.mean,ymin=Value.mean-Value.sd, ymax=Value.mean+Value.sd,group=Type),width=0,linewidth=0.2,color="black",position=position_dodge(0.5),alpha=0.6)+
  theme_classic()+
  scale_color_manual(values=c("#bc40bf","#4BA465","black"))+
  scale_fill_manual(values=c("#bc40bf","#4BA465","black"))+
  ylab("Interaction evenness")+# 设置Y轴标题
  # scale_x_discrete(breaks=c("Eu-adjacent","Eu-distant"),labels=c("Eu-adjacent", "Eu-distant"))+
  theme(axis.title.x=element_blank(),#修改X轴标题文本
        # axis.title.x=element_text(size=8),#修改X轴标题文本
        axis.title.y=element_text(size=8,angle=90),#修改y轴标题文本
        axis.text.y=element_text(size=6.5),#修改x轴刻度标签文本
        axis.line.x = element_line(arrow = arrow(length = unit(0.2,"lines"))),
        axis.text.x=element_text(size=6.5),#修改y轴刻度标签文本
        panel.grid=element_blank(),#隐藏网格线
        # legend.title=element_text(size=6.5),#修改图例标题大小
        legend.title=element_blank(),#修改图例标题大小
        legend.text=element_text(size=5.5),#修改图例文本大小
        legend.key.size=unit(0.1,"lines"),#缩小图例符号的大小
        legend.background=element_blank(),#设置图例背景为空白
        legend.position="inside",#修改图例位置
        # legend.justification=c("right","top"),#修改图例位置
        legend.position.inside=c(0.9,0.9),#修改图例位置
        legend.margin=margin(0,0,0,0))+#修改图例位置
  annotate("text",x=1.8,y=max(Evenness_long$Value),hjust=1,vjust=1,label=paste0("p > ",0.05),size=2,color="black")
p

ggsave(filename="Fig6b_network_index_evenness.png",plot=p,width=89/2,height=59/2,units="mm",dpi=900)
ggsave(filename="Fig6b_network_index_evenness.pdf",plot=p,width=89/2,height=59/2,units="mm")

#####
WNODF <- subset(network_index,Index_type=="Weighted NODF")
x <- subset(WNODF,Plot=="Eu-adjacent")$Observe
y <- subset(WNODF,Plot=="Eu-distant")$Observe
t.test(x,y,paired=T)

x <- subset(WNODF,Plot=="Eu-adjacent")$Observe_Eup.jo_excluded
y <- subset(WNODF,Plot=="Eu-distant")$Observe_Eup.jo_excluded
t.test(x,y,paired=T)

WNODF_long <- pivot_longer(WNODF, cols = c(Observe, Observe_Eup.jo_excluded), names_to = "Type", values_to = "Value")
WNODF_long[which(WNODF_long$Type=="Observe"),]$Type <- "included"
WNODF_long[which(WNODF_long$Type=="Observe_Eup.jo_excluded"),]$Type <- "excluded"
WNODF_long$Type <- factor(WNODF_long$Type,levels=c("included","excluded"))

# 计算均值和标准差
WNODF_long.summary <- WNODF_long %>%
  group_by(Plot, Type) %>%
  summarise(Value.mean = mean(Value),Value.sd = sd(Value))

p <- ggplot(WNODF_long,aes(x=Plot,y=Value,color=Type,group=Type))+
  geom_jitter(shape=16,size=0.4,position=position_jitterdodge(jitter.height=0,jitter.width=0.25,dodge.width=0.5),alpha=0.8,show.legend=T)+
  geom_point(data=WNODF_long.summary,aes(x=Plot,y=Value.mean,group=Type),shape=16,size=0.6,color="black",position=position_dodge(0.5),alpha=0.8)+
  geom_errorbar(data=WNODF_long.summary,aes(x=Plot,y=Value.mean,ymin=Value.mean-Value.sd, ymax=Value.mean+Value.sd,group=Type),width=0,linewidth=0.2,color="black",position=position_dodge(0.5),alpha=0.6)+
  theme_classic()+
  scale_color_manual(values=c("#bc40bf","#4BA465","black"))+
  scale_fill_manual(values=c("#bc40bf","#4BA465","black"))+
  ylab("Weighted NODF")+# 设置Y轴标题
  # scale_x_discrete(breaks=c("Eu-adjacent","Eu-distant"),labels=c("Eu-adjacent", "Eu-distant"))+
  theme(axis.title.x=element_blank(),#修改X轴标题文本
        # axis.title.x=element_text(size=8),#修改X轴标题文本
        axis.title.y=element_text(size=8,angle=90),#修改y轴标题文本
        axis.text.y=element_text(size=6.5),#修改x轴刻度标签文本
        axis.line.x = element_line(arrow = arrow(length = unit(0.2,"lines"))),
        axis.text.x=element_text(size=6.5),#修改y轴刻度标签文本
        panel.grid=element_blank(),#隐藏网格线
        # legend.title=element_text(size=6.5),#修改图例标题大小
        legend.title=element_blank(),#修改图例标题大小
        legend.text=element_text(size=5.5),#修改图例文本大小
        legend.key.size=unit(0.1,"lines"),#缩小图例符号的大小
        legend.background=element_blank(),#设置图例背景为空白
        legend.position="inside",#修改图例位置
        # legend.justification=c("right","top"),#修改图例位置
        legend.position.inside=c(0.9,0.9),#修改图例位置
        legend.margin=margin(0,0,0,0))+#修改图例位置
  annotate("text",x=1.8,y=max(WNODF_long$Value),hjust=1,vjust=1,label=paste0("p < ",0.01),size=2,color="black")
p

ggsave(filename="Fig6c_network_index_Weighted_NODF.png",plot=p,width=89/2,height=59/2,units="mm",dpi=900)
ggsave(filename="Fig6c_network_index_Weighted_NODF.pdf",plot=p,width=89/2,height=59/2,units="mm")

#####
WMOD <- subset(network_index,Index_type=="Weighted modularity")
x <- subset(WMOD,Plot=="Eu-adjacent")$Observe
y <- subset(WMOD,Plot=="Eu-distant")$Observe
t.test(x,y,paired=T)

x <- subset(WMOD,Plot=="Eu-adjacent")$Observe_Eup.jo_excluded
y <- subset(WMOD,Plot=="Eu-distant")$Observe_Eup.jo_excluded
t.test(x,y,paired=T)

WMOD_long <- pivot_longer(WMOD, cols = c(Observe, Observe_Eup.jo_excluded), names_to = "Type", values_to = "Value")
WMOD_long[which(WMOD_long$Type=="Observe"),]$Type <- "included"
WMOD_long[which(WMOD_long$Type=="Observe_Eup.jo_excluded"),]$Type <- "excluded"
WMOD_long$Type <- factor(WMOD_long$Type,levels=c("included","excluded"))

# 计算均值和标准差
WMOD_long.summary <- WMOD_long %>%
  group_by(Plot, Type) %>%
  summarise(Value.mean = mean(Value),Value.sd = sd(Value))

p <- ggplot(WMOD_long,aes(x=Plot,y=Value,color=Type,group=Type))+
  geom_jitter(shape=16,size=0.4,position=position_jitterdodge(jitter.height=0,jitter.width=0.25,dodge.width=0.5),alpha=0.8,show.legend=T)+
  geom_point(data=WMOD_long.summary,aes(x=Plot,y=Value.mean,group=Type),shape=16,size=0.6,color="black",position=position_dodge(0.5),alpha=0.8)+
  geom_errorbar(data=WMOD_long.summary,aes(x=Plot,y=Value.mean,ymin=Value.mean-Value.sd, ymax=Value.mean+Value.sd,group=Type),width=0,linewidth=0.2,color="black",position=position_dodge(0.5),alpha=0.6)+
  theme_classic()+
  scale_color_manual(values=c("#bc40bf","#4BA465","black"))+
  scale_fill_manual(values=c("#bc40bf","#4BA465","black"))+
  ylab("Weighted modularity")+# 设置Y轴标题
  # scale_x_discrete(breaks=c("Eu-adjacent","Eu-distant"),labels=c("Eu-adjacent", "Eu-distant"))+
  theme(axis.title.x=element_blank(),#修改X轴标题文本
        # axis.title.x=element_text(size=8),#修改X轴标题文本
        axis.title.y=element_text(size=8,angle=90),#修改y轴标题文本
        axis.text.y=element_text(size=6.5),#修改x轴刻度标签文本
        axis.line.x = element_line(arrow = arrow(length = unit(0.2,"lines"))),
        axis.text.x=element_text(size=6.5),#修改y轴刻度标签文本
        panel.grid=element_blank(),#隐藏网格线
        # legend.title=element_text(size=6.5),#修改图例标题大小
        legend.title=element_blank(),#修改图例标题大小
        legend.text=element_text(size=5.5),#修改图例文本大小
        legend.key.size=unit(0.1,"lines"),#缩小图例符号的大小
        legend.background=element_blank(),#设置图例背景为空白
        legend.position="inside",#修改图例位置
        # legend.justification=c("right","top"),#修改图例位置
        legend.position.inside=c(0.9,0.9),#修改图例位置
        legend.margin=margin(0,0,0,0))+#修改图例位置
  annotate("text",x=1.8,y=max(WMOD_long$Value),hjust=1,vjust=1,label=paste0("p < ",0.05),size=2,color="black")
p

ggsave(filename="Fig6d_network_index_Weighted_modularity.png",plot=p,width=89/2,height=59/2,units="mm",dpi=900)
ggsave(filename="Fig6d_network_index_Weighted_modularity.pdf",plot=p,width=89/2,height=59/2,units="mm")

#####
# simulation
rm(list=ls())

network_index_simulation <- read_excel("source_data.xlsx", sheet="Network_index_Simulation")
head(network_index_simulation)

network_index_simulation$Site <- as.factor(network_index_simulation$Site)
network_index_simulation$Between_plot <- as.factor(network_index_simulation$Between_plot)
network_index_simulation$Index_type <- as.factor(network_index_simulation$Index_type)
network_index_simulation$Delta_type <- as.factor(network_index_simulation$Delta_type)
summary(network_index_simulation)


WCON <- subset(network_index_simulation,Index_type=="Weighted connectance")
x <- subset(WCON,Delta_type=="delta_Observed")$Delta_value
y <- subset(WCON,Delta_type=="delta_Simulation")$Delta_value
t.test(x,y,paired=T)

x <- subset(WCON,Delta_type=="delta_Observed")$Delta_value_Eup.jo_excluded
y <- subset(WCON,Delta_type=="delta_Simulation")$Delta_value_Eup.jo_excluded
t.test(x,y,paired=T)

WCON_long <- pivot_longer(WCON, cols = c(Delta_value, Delta_value_Eup.jo_excluded), names_to = "Type", values_to = "Value")
WCON_long[which(WCON_long$Type=="Delta_value"),]$Type <- "included"
WCON_long[which(WCON_long$Type=="Delta_value_Eup.jo_excluded"),]$Type <- "excluded"
WCON_long$Type <- factor(WCON_long$Type,levels=c("included","excluded"))

# 计算均值和标准差
WCON_long.summary <- WCON_long %>%
  group_by(Delta_type, Type) %>%
  summarise(Value.mean = mean(Value),Value.sd = sd(Value))

p <- ggplot(WCON_long,aes(x=Delta_type,y=Value,color=Type,group=Type))+
  geom_jitter(shape=16,size=0.4,position=position_jitterdodge(jitter.height=0,jitter.width=0.25,dodge.width=0.5),alpha=0.8,show.legend=T)+
  geom_point(data=WCON_long.summary,aes(x=Delta_type,y=Value.mean,group=Type),shape=16,size=0.6,color="black",position=position_dodge(0.5),alpha=0.8)+
  geom_errorbar(data=WCON_long.summary,aes(x=Delta_type,y=Value.mean,ymin=Value.mean-Value.sd, ymax=Value.mean+Value.sd,group=Type),width=0,linewidth=0.2,color="black",position=position_dodge(0.5),alpha=0.6)+
  theme_classic()+
  scale_color_manual(values=c("#bc40bf","#4BA465","black"))+
  scale_fill_manual(values=c("#bc40bf","#4BA465","black"))+
  ylab("Weighted connectance")+# 设置Y轴标题
  scale_x_discrete(breaks=c("delta_Observed","delta_Simulation"),labels=c("ΔObserved", "ΔSimulation"))+
  theme(axis.title.x=element_blank(),#修改X轴标题文本
        # axis.title.x=element_text(size=8),#修改X轴标题文本
        axis.title.y=element_text(size=8,angle=90),#修改y轴标题文本
        axis.text.y=element_text(size=6.5),#修改x轴刻度标签文本
        axis.line.x = element_line(arrow = arrow(length = unit(0.2,"lines"))),
        axis.text.x=element_text(size=6.5),#修改y轴刻度标签文本
        panel.grid=element_blank(),#隐藏网格线
        # legend.title=element_text(size=6.5),#修改图例标题大小
        legend.title=element_blank(),#修改图例标题大小
        legend.text=element_text(size=5.5),#修改图例文本大小
        legend.key.size=unit(0.1,"lines"),#缩小图例符号的大小
        legend.background=element_blank(),#设置图例背景为空白
        legend.position="inside",#修改图例位置
        # legend.justification=c("right","top"),#修改图例位置
        legend.position.inside=c(0.9,0.9),#修改图例位置
        legend.margin=margin(0,0,0,0))+#修改图例位置
  annotate("text",x=1.8,y=max(WCON_long$Value),hjust=1,vjust=1,label=paste0("p < ",0.01),size=2,color="black")
p

ggsave(filename="Fig8a_network_index_Weighted_connectance_simulation.png",plot=p,width=89/2,height=59/2,units="mm",dpi=900)
ggsave(filename="Fig8a_network_index_Weighted_connectance_simulation.pdf",plot=p,width=89/2,height=59/2,units="mm")

#####
Evenness <- subset(network_index_simulation,Index_type=="Interaction evenness")
x <- subset(Evenness,Delta_type=="delta_Observed")$Delta_value
y <- subset(Evenness,Delta_type=="delta_Simulation")$Delta_value
t.test(x,y,paired=T)

x <- subset(Evenness,Delta_type=="delta_Observed")$Delta_value_Eup.jo_excluded
y <- subset(Evenness,Delta_type=="delta_Simulation")$Delta_value_Eup.jo_excluded
t.test(x,y,paired=T)

Evenness_long <- pivot_longer(Evenness, cols = c(Delta_value, Delta_value_Eup.jo_excluded), names_to = "Type", values_to = "Value")
Evenness_long[which(Evenness_long$Type=="Delta_value"),]$Type <- "included"
Evenness_long[which(Evenness_long$Type=="Delta_value_Eup.jo_excluded"),]$Type <- "excluded"
Evenness_long$Type <- factor(Evenness_long$Type,levels=c("included","excluded"))

# 计算均值和标准差
Evenness_long.summary <- Evenness_long %>%
  group_by(Delta_type, Type) %>%
  summarise(Value.mean = mean(Value),Value.sd = sd(Value))

p <- ggplot(Evenness_long,aes(x=Delta_type,y=Value,color=Type,group=Type))+
  geom_jitter(shape=16,size=0.4,position=position_jitterdodge(jitter.height=0,jitter.width=0.25,dodge.width=0.5),alpha=0.8,show.legend=T)+
  geom_point(data=Evenness_long.summary,aes(x=Delta_type,y=Value.mean,group=Type),shape=16,size=0.6,color="black",position=position_dodge(0.5),alpha=0.8)+
  geom_errorbar(data=Evenness_long.summary,aes(x=Delta_type,y=Value.mean,ymin=Value.mean-Value.sd, ymax=Value.mean+Value.sd,group=Type),width=0,linewidth=0.2,color="black",position=position_dodge(0.5),alpha=0.6)+
  theme_classic()+
  scale_color_manual(values=c("#bc40bf","#4BA465","black"))+
  scale_fill_manual(values=c("#bc40bf","#4BA465","black"))+
  ylab("Interaction evenness")+# 设置Y轴标题
  scale_x_discrete(breaks=c("delta_Observed","delta_Simulation"),labels=c("ΔObserved", "ΔSimulation"))+
  theme(axis.title.x=element_blank(),#修改X轴标题文本
        # axis.title.x=element_text(size=8),#修改X轴标题文本
        axis.title.y=element_text(size=8,angle=90),#修改y轴标题文本
        axis.text.y=element_text(size=6.5),#修改x轴刻度标签文本
        axis.line.x = element_line(arrow = arrow(length = unit(0.2,"lines"))),
        axis.text.x=element_text(size=6.5),#修改y轴刻度标签文本
        panel.grid=element_blank(),#隐藏网格线
        # legend.title=element_text(size=6.5),#修改图例标题大小
        legend.title=element_blank(),#修改图例标题大小
        legend.text=element_text(size=5.5),#修改图例文本大小
        legend.key.size=unit(0.1,"lines"),#缩小图例符号的大小
        legend.background=element_blank(),#设置图例背景为空白
        legend.position="inside",#修改图例位置
        # legend.justification=c("right","top"),#修改图例位置
        legend.position.inside=c(0.9,0.9),#修改图例位置
        legend.margin=margin(0,0,0,0))+#修改图例位置
  annotate("text",x=1.8,y=max(Evenness_long$Value),hjust=1,vjust=1,label=paste0("p > ",0.05),size=2,color="black")
p

ggsave(filename="Fig8b_network_index_evenness_simulation.png",plot=p,width=89/2,height=59/2,units="mm",dpi=900)
ggsave(filename="Fig8b_network_index_evenness_simulation.pdf",plot=p,width=89/2,height=59/2,units="mm")

#####
WNODF <- subset(network_index_simulation,Index_type=="Weighted NODF")
x <- subset(WNODF,Delta_type=="delta_Observed")$Delta_value
y <- subset(WNODF,Delta_type=="delta_Simulation")$Delta_value
t.test(x,y,paired=T)

x <- subset(WNODF,Delta_type=="delta_Observed")$Delta_value_Eup.jo_excluded
y <- subset(WNODF,Delta_type=="delta_Simulation")$Delta_value_Eup.jo_excluded
t.test(x,y,paired=T)

WNODF_long <- pivot_longer(WNODF, cols = c(Delta_value, Delta_value_Eup.jo_excluded), names_to = "Type", values_to = "Value")
WNODF_long[which(WNODF_long$Type=="Delta_value"),]$Type <- "included"
WNODF_long[which(WNODF_long$Type=="Delta_value_Eup.jo_excluded"),]$Type <- "excluded"
WNODF_long$Type <- factor(WNODF_long$Type,levels=c("included","excluded"))

# 计算均值和标准差
WNODF_long.summary <- WNODF_long %>%
  group_by(Delta_type, Type) %>%
  summarise(Value.mean = mean(Value),Value.sd = sd(Value))

p <- ggplot(WNODF_long,aes(x=Delta_type,y=Value,color=Type,group=Type))+
  geom_jitter(shape=16,size=0.4,position=position_jitterdodge(jitter.height=0,jitter.width=0.25,dodge.width=0.5),alpha=0.8,show.legend=T)+
  geom_point(data=WNODF_long.summary,aes(x=Delta_type,y=Value.mean,group=Type),shape=16,size=0.6,color="black",position=position_dodge(0.5),alpha=0.8)+
  geom_errorbar(data=WNODF_long.summary,aes(x=Delta_type,y=Value.mean,ymin=Value.mean-Value.sd, ymax=Value.mean+Value.sd,group=Type),width=0,linewidth=0.2,color="black",position=position_dodge(0.5),alpha=0.6)+
  theme_classic()+
  scale_color_manual(values=c("#bc40bf","#4BA465","black"))+
  scale_fill_manual(values=c("#bc40bf","#4BA465","black"))+
  ylab("Weighted NODF")+# 设置Y轴标题
  scale_x_discrete(breaks=c("delta_Observed","delta_Simulation"),labels=c("ΔObserved", "ΔSimulation"))+
  theme(axis.title.x=element_blank(),#修改X轴标题文本
        # axis.title.x=element_text(size=8),#修改X轴标题文本
        axis.title.y=element_text(size=8,angle=90),#修改y轴标题文本
        axis.text.y=element_text(size=6.5),#修改x轴刻度标签文本
        axis.line.x = element_line(arrow = arrow(length = unit(0.2,"lines"))),
        axis.text.x=element_text(size=6.5),#修改y轴刻度标签文本
        panel.grid=element_blank(),#隐藏网格线
        # legend.title=element_text(size=6.5),#修改图例标题大小
        legend.title=element_blank(),#修改图例标题大小
        legend.text=element_text(size=5.5),#修改图例文本大小
        legend.key.size=unit(0.1,"lines"),#缩小图例符号的大小
        legend.background=element_blank(),#设置图例背景为空白
        legend.position="inside",#修改图例位置
        # legend.justification=c("right","top"),#修改图例位置
        legend.position.inside=c(0.9,0.9),#修改图例位置
        legend.margin=margin(0,0,0,0))+#修改图例位置
  annotate("text",x=1.8,y=max(WNODF_long$Value),hjust=1,vjust=1,label=paste0("p < ",0.01),size=2,color="black")
p

ggsave(filename="Fig8c_network_index_Weighted_NODF_simulation.png",plot=p,width=89/2,height=59/2,units="mm",dpi=900)
ggsave(filename="Fig8c_network_index_Weighted_NODF_simulation.pdf",plot=p,width=89/2,height=59/2,units="mm")

#####
WMOD <- subset(network_index_simulation,Index_type=="Weighted modularity")
x <- subset(WMOD,Delta_type=="delta_Observed")$Delta_value
y <- subset(WMOD,Delta_type=="delta_Simulation")$Delta_value
t.test(x,y,paired=T)

x <- subset(WMOD,Delta_type=="delta_Observed")$Delta_value_Eup.jo_excluded
y <- subset(WMOD,Delta_type=="delta_Simulation")$Delta_value_Eup.jo_excluded
t.test(x,y,paired=T)

WMOD_long <- pivot_longer(WMOD, cols = c(Delta_value, Delta_value_Eup.jo_excluded), names_to = "Type", values_to = "Value")
WMOD_long[which(WMOD_long$Type=="Delta_value"),]$Type <- "included"
WMOD_long[which(WMOD_long$Type=="Delta_value_Eup.jo_excluded"),]$Type <- "excluded"
WMOD_long$Type <- factor(WMOD_long$Type,levels=c("included","excluded"))

# 计算均值和标准差
WMOD_long.summary <- WMOD_long %>%
  group_by(Delta_type, Type) %>%
  summarise(Value.mean = mean(Value),Value.sd = sd(Value))

p <- ggplot(WMOD_long,aes(x=Delta_type,y=Value,color=Type,group=Type))+
  geom_jitter(shape=16,size=0.4,position=position_jitterdodge(jitter.height=0,jitter.width=0.25,dodge.width=0.5),alpha=0.8,show.legend=T)+
  geom_point(data=WMOD_long.summary,aes(x=Delta_type,y=Value.mean,group=Type),shape=16,size=0.6,color="black",position=position_dodge(0.5),alpha=0.8)+
  geom_errorbar(data=WMOD_long.summary,aes(x=Delta_type,y=Value.mean,ymin=Value.mean-Value.sd, ymax=Value.mean+Value.sd,group=Type),width=0,linewidth=0.2,color="black",position=position_dodge(0.5),alpha=0.6)+
  theme_classic()+
  scale_color_manual(values=c("#bc40bf","#4BA465","black"))+
  scale_fill_manual(values=c("#bc40bf","#4BA465","black"))+
  ylab("Weighted modularity")+# 设置Y轴标题
  scale_x_discrete(breaks=c("delta_Observed","delta_Simulation"),labels=c("ΔObserved", "ΔSimulation"))+
  theme(axis.title.x=element_blank(),#修改X轴标题文本
        # axis.title.x=element_text(size=8),#修改X轴标题文本
        axis.title.y=element_text(size=8,angle=90),#修改y轴标题文本
        axis.text.y=element_text(size=6.5),#修改x轴刻度标签文本
        axis.line.x = element_line(arrow = arrow(length = unit(0.2,"lines"))),
        axis.text.x=element_text(size=6.5),#修改y轴刻度标签文本
        panel.grid=element_blank(),#隐藏网格线
        # legend.title=element_text(size=6.5),#修改图例标题大小
        legend.title=element_blank(),#修改图例标题大小
        legend.text=element_text(size=5.5),#修改图例文本大小
        legend.key.size=unit(0.1,"lines"),#缩小图例符号的大小
        legend.background=element_blank(),#设置图例背景为空白
        legend.position="inside",#修改图例位置
        # legend.justification=c("right","top"),#修改图例位置
        legend.position.inside=c(0.9,0.9),#修改图例位置
        legend.margin=margin(0,0,0,0))+#修改图例位置
  annotate("text",x=1.8,y=max(WMOD_long$Value),hjust=1,vjust=1,label=paste0("p < ",0.05),size=2,color="black")
p

ggsave(filename="Fig8d_network_index_Weighted_modularity_simulation.png",plot=p,width=89/2,height=59/2,units="mm",dpi=900)
ggsave(filename="Fig8d_network_index_Weighted_modularity_simulation.pdf",plot=p,width=89/2,height=59/2,units="mm")

#################################################################################
#################################################################################
#################################################################################