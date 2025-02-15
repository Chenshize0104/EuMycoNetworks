# Data statistical analyses & figure plots
# setwd('E:/Shangri-La_2024/network')

library(readxl)
library(tidyverse)
library(vegan)
library(ggplot2)
library(glmmTMB)
library(MuMIn)
library(DHARMa)
library(effects)
library(gridExtra)
library(cowplot)
library(car)

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

p1 <- ggplot(data.plant.Euadjacent,aes(Plant_Species_abbreviation,Value))+
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
        axis.text.y=element_text(size=6.5),#修改y轴刻度标签文本
        # axis.line.x = element_line(arrow = arrow(length = unit(0.2,"lines"))),
        axis.text.x=element_blank(),#修改x轴刻度标签文本
        # axis.text.x=element_text(angle=30, size=6.5, hjust=0.5),#修改x轴刻度标签文本
        panel.grid=element_blank(),#隐藏网格线
        legend.title=element_blank(),#修改图例标题大小
        # legend.title=element_text(size=6.5),#修改图例标题大小
        legend.text=element_text(size=5.5),#修改图例文本大小
        legend.key.size=unit(0.1,"lines"),#缩小图例符号的大小
        legend.background=element_blank(),#设置图例背景为空白
        legend.position="inside",#修改图例位置
        # legend.justification=c("right","top"),#修改图例位置
        legend.position.inside=c(0.95,0.95),#修改图例位置
        legend.margin=margin(0,0,0,0),#修改图例位置
        axis.line.x=element_line(),
        strip.text = element_blank(), strip.background = element_blank())  # 取消分面标签)
p1
ggsave(filename="Fig3a_plant_top10_Euadjacent.png",plot=p1,width=89*1.3,height=59*1.3/2,units="mm",dpi=900)
ggsave(filename="Fig3a_plant_top10_Euadjacent.pdf",plot=p1,width=89*1.3,height=59*1.3/2,units="mm")

p2 <- ggplot(data.plant.Eudistant,aes(Plant_Species_abbreviation,Value))+
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
        axis.text.y=element_text(size=6.5),#修改y轴刻度标签文本
        # axis.line.x = element_line(arrow = arrow(length = unit(0.2,"lines"))),
        axis.text.x=element_text(angle=30, size=6.5, hjust=0.5),#修改x轴刻度标签文本
        panel.grid=element_blank(),#隐藏网格线
        legend.title=element_blank(),#修改图例标题大小
        # legend.title=element_text(size=6.5),#修改图例标题大小
        legend.text=element_text(size=5.5),#修改图例文本大小
        legend.key.size=unit(0.1,"lines"),#缩小图例符号的大小
        legend.background=element_blank(),#设置图例背景为空白
        legend.position="inside",#修改图例位置
        # legend.justification=c("right","top"),#修改图例位置
        legend.position.inside=c(0.95,0.95),#修改图例位置
        legend.margin=margin(0,0,0,0),#修改图例位置
        axis.line.x=element_line(),
        strip.text = element_blank(), strip.background = element_blank())  # 取消分面标签)
p2
ggsave(filename="Fig3b_plant_top10_Eudistant.png",plot=p2,width=89*1.3,height=59*1.3/2,units="mm",dpi=900)
ggsave(filename="Fig3b_plant_top10_Eudistant.pdf",plot=p2,width=89*1.3,height=59*1.3/2,units="mm")

#####
data.fungi <- subset(data, Type=="AMF_richness")
data.fungi$AMF_richness_Estimator <- as.numeric(data.fungi$AMF_richness_Estimator)
data.fungi.Euadjacent <- subset(data.fungi, Plot=="Eu-adjacent")
data.fungi.Eudistant <- subset(data.fungi, Plot=="Eu-distant")

p3 <- ggplot(data.fungi.Euadjacent, aes(x = Plant_Species_abbreviation, y = Value)) +
  geom_bar(aes(fill=Plot),stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(y=Value,ymin=Value*0.95,ymax=AMF_richness_Estimator),width=0.8,linewidth=0.4,color="gray50",alpha=0.8)+
  geom_text(aes(label=paste0("n=",Sample_count)),size=1.5,vjust=0)+
  facet_grid(Plot ~ .) +
  theme_classic()+
  scale_color_manual(values=c("gray50"))+
  scale_fill_manual(values=c("gray50"))+
  ylab("AMF richness")+# 设置Y轴标题
  theme(axis.title.x=element_blank(),#修改X轴标题文本
        # axis.title.x=element_text(size=8),#修改X轴标题文本
        axis.title.y=element_text(size=8,angle=90),#修改y轴标题文本
        axis.text.y=element_text(size=6.5),#修改y轴刻度标签文本
        # axis.line.x = element_line(arrow = arrow(length = unit(0.2,"lines"))),
        axis.text.x=element_blank(),#修改x轴刻度标签文本
        # axis.text.x=element_text(angle=30, size=6.5, hjust=0.5),#修改x轴刻度标签文本
        panel.grid=element_blank(),#隐藏网格线
        legend.title=element_blank(),#修改图例标题大小
        # legend.title=element_text(size=6.5),#修改图例标题大小
        legend.text=element_text(size=5.5),#修改图例文本大小
        legend.key.size=unit(0.1,"lines"),#缩小图例符号的大小
        legend.background=element_blank(),#设置图例背景为空白
        legend.position="inside",#修改图例位置
        # legend.justification=c("right","top"),#修改图例位置
        legend.position.inside=c(0.95,0.95),#修改图例位置
        legend.margin=margin(0,0,0,0),#修改图例位置
        axis.line.x=element_line(),
        strip.text = element_blank(), strip.background = element_blank())  # 取消分面标签)
p3
ggsave(filename="Fig3c_AMF_top10_Euadjacent.png",plot=p3,width=89*1.3,height=59*1.3/2,units="mm",dpi=900)
ggsave(filename="Fig3c_AMF_top10_Euadjacent.pdf",plot=p3,width=89*1.3,height=59*1.3/2,units="mm")


p4 <- ggplot(data.fungi.Eudistant, aes(x = Plant_Species_abbreviation, y = Value)) +
  geom_bar(aes(fill=Plot),stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(y=Value,ymin=Value*0.95,ymax=AMF_richness_Estimator),width=0.8,linewidth=0.4,color="gray50",alpha=0.8)+
  geom_text(aes(label=paste0("n=",Sample_count)),size=1.5,vjust=0)+
  facet_grid(Plot ~ .) +
  theme_classic()+
  scale_color_manual(values=c("gray50"))+
  scale_fill_manual(values=c("gray50"))+
  ylab("AMF richness")+# 设置Y轴标题
  theme(axis.title.x=element_blank(),#修改X轴标题文本
        # axis.title.x=element_text(size=8),#修改X轴标题文本
        axis.title.y=element_text(size=8,angle=90),#修改y轴标题文本
        axis.text.y=element_text(size=6.5),#修改y轴刻度标签文本
        # axis.line.x = element_line(arrow = arrow(length = unit(0.2,"lines"))),
        axis.text.x=element_text(angle=30, size=6.5, hjust=0.5),#修改x轴刻度标签文本
        panel.grid=element_blank(),#隐藏网格线
        legend.title=element_blank(),#修改图例标题大小
        # legend.title=element_text(size=6.5),#修改图例标题大小
        legend.text=element_text(size=5.5),#修改图例文本大小
        legend.key.size=unit(0.1,"lines"),#缩小图例符号的大小
        legend.background=element_blank(),#设置图例背景为空白
        legend.position="inside",#修改图例位置
        # legend.justification=c("right","top"),#修改图例位置
        legend.position.inside=c(0.95,0.95),#修改图例位置
        legend.margin=margin(0,0,0,0),#修改图例位置
        axis.line.x=element_line(),
        strip.text = element_blank(), strip.background = element_blank())  # 取消分面标签)
p4
ggsave(filename="Fig3d_AMF_top10_Eudistant.png",plot=p4,width=89*1.3,height=59*1.3/2,units="mm",dpi=900)
ggsave(filename="Fig3d_AMF_top10_Eudistant.pdf",plot=p4,width=89*1.3,height=59*1.3/2,units="mm")

?plot_grid
p_all <- plot_grid(p1,p2,p3,p4,labels=c("a","","b",""),label_size=10,label_x=0,label_y=1,ncol=1,nrow=4)
p_all
ggsave(filename="Fig3.png",plot=p_all,width=89*1.2,height=59*2.4,units="mm",dpi=900)
ggsave(filename="Fig3.pdf",plot=p_all,width=89*1.2,height=59*2.4,units="mm")

################################################################################

# Fig. 4 Species_richness_interactions
rm(list=ls())

richness <- read_excel("source_data.xlsx", sheet="Species_richness_interactions")
head(richness)

richness$Site <- as.factor(richness$Site)
richness$Plot <- as.factor(richness$Plot)
richness$Plot <- as.factor(richness$Plot)
summary(richness)

# 如果Shapiro-Wilk检验的p值小于0.05，说明数据不符合正态分布。以下结果大于0.05，符合正态分布。
shapiro.test(richness$Plant_richness)

# Levene检验的p值如果大于0.05，说明各组之间方差齐性，可以进行方差分析，否则变换数据或选择非参数检验。
leveneTest(Plant_richness~Plot, data = richness)

# 方差分析
anova_result <- aov(Plant_richness~Plot, data = richness)
summary(anova_result)


# 非参数检验（Kruskal-Wallis检验）
kruskal_test <- kruskal.test(Plant_richness~Plot, data = richness)
# 输出检验结果
kruskal_test

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
  shapiro.test(x)
  shapiro.test(y)
}

# library(lme4)
# mod1 <- glm(Plant_richness ~ Plot, family=poisson, data=richness)
# mod2 <- glmer(Plant_richness ~ Plot + (1|Site), family=poisson, data=richness)
# (msAICc <- model.sel(mod1,mod2))
# summary(mod1)
# performance::r2(mod1)
# allEffects(mod1)
# plot(allEffects(mod1))

modlmer1 <- glmmTMB(Plant_richness ~ Plot, family=poisson, data=richness)
modlmer2 <- glmmTMB(Plant_richness ~ Plot + (1|Site), family=poisson, data=richness)
(msAICc <- model.sel(modlmer1,modlmer2))

# modlmer3 <- glmmTMB(Plant_richness ~ Plot + (1|Site), family=nbinom2, data=richness)
# modlmer4 <- glmmTMB(Plant_richness ~ Plot + (1|Site), family=nbinom1, data=richness)
# (msAICc <- model.sel(modlmer1,modlmer2,modlmer3,modlmer4))

plot(simulateResiduals(modlmer1))

summary(modlmer1)
performance::r2(modlmer1)
allEffects(modlmer1)
plot(allEffects(modlmer1))

Anova(modlmer1,type=2)

Anova(modlmer2,type=2)


##
modlmer1 <- glmmTMB(log(Plant_richness) ~ Plot, data=richness)
modlmer2 <- glmmTMB(log(Plant_richness) ~ Plot + (1|Site), data=richness)
(msAICc <- model.sel(modlmer1,modlmer2))

plot(simulateResiduals(modlmer2))

summary(modlmer2)
performance::r2(modlmer2)
allEffects(modlmer2)
plot(allEffects(modlmer2))

p1 <- ggplot(richness,aes(x=Plot,y=Plant_richness,color=Plot))+
  geom_boxplot(width=0.15,linewidth=0.25,outliers=F,alpha=1,show.legend=T,position=position_dodge(0.5))+
  stat_boxplot(geom="errorbar",width=0.15,linewidth=0.25,position=position_dodge(0.5))+
  geom_jitter(shape=16,size=0.4,position=position_jitterdodge(jitter.height=0.25,jitter.width=0.45,dodge.width=0.5),alpha=0.4,show.legend=T)+
  # geom_point(aes(x=Plot,y=mean(Plant_richness)),shape=16,size=0.8,color="black",position=position_dodge(0.5),alpha=1)+
  # geom_errorbar(aes(x=Plot,y=mean(Plant_richness),ymin=mean(Plant_richness)-sd(Plant_richness), ymax=mean(Plant_richness)+sd(Plant_richness)),width=0,linewidth=0.4,color="black",position=position_dodge(0.5),alpha=1)+
  # scale_color_manual(values=c("#bc40bf","#bc40bf","black"))+
  # scale_fill_manual(values=c("#bc40bf","#bc40bf"))+ 
  scale_color_manual(values=c("gray40","gray40","gray40"))+
  scale_fill_manual(values=c("gray40","gray40"))+ 
  theme_classic()+
  ylab("Plant richness")+
  # scale_x_discrete(breaks=c("Eu-adjacent","Eu-distant"),labels=c("Eu-adjacent", "Eu-distant"))+
  # scale_y_continuous(labels = function(x) sprintf("%03d", x))+
  theme(axis.title.x=element_blank(),#修改X轴标题文本
        # axis.title.x=element_text(size=8),#修改X轴标题文本
        # axis.line.x = element_line(arrow = arrow(length = unit(0.2,"lines"))),
        axis.title.y=element_text(size=8,angle=90),#修改y轴标题文本
        axis.text.y=element_text(size=6.5),#修改y轴刻度标签文本
        axis.text.x=element_blank(),#修改x轴刻度标签文本
        # axis.text.x=element_text(size=6.5),#修改x轴刻度标签文本
        panel.grid=element_blank(),#隐藏网格线
        legend.title=element_text(size=8),#修改图例标题大小
        legend.text=element_text(size=6.5),#修改图例文本大小
        legend.position="none")+#修改图例位置
  # annotate("text",x=1.7,y=max(richness$Plant_richness),hjust=1,vjust=1,label=paste0("p > ",0.05),size=2)+
  annotate("text",x=1.5,y=(max(richness$Plant_richness)+(max(richness$Plant_richness)-min(richness$Plant_richness))*0.1),hjust=1,vjust=1,label=paste0("NS"),size=2)
p1

ggsave(filename="Fig4a_plant_richness_glmm.png",plot=p1,width=89/2,height=59/2,units="mm",dpi=900)
ggsave(filename="Fig4a_plant_richness_glmm.pdf",plot=p1,width=89/2,height=59/2,units="mm")

#####
# Biomass of Plant
# Aboveground_biomass
{
  Euadjacent.plot <- subset(richness,Plot=="Eu-adjacent")
  Eudistant.plot <- subset(richness,Plot=="Eu-distant")
  Euadjacent.plot$PlotID==Eudistant.plot$PlotID
  x <- Euadjacent.plot$`Aboveground_biomass (kg/m2)`
  y <- Eudistant.plot$`Aboveground_biomass (kg/m2)`
  t.test(x,y,paired=T)
  mean(x);mean(y)
  sd(x);sd(y)
  var(x);var(y)
  shapiro.test(x)
  shapiro.test(y)
}

leveneTest(`Aboveground_biomass (kg/m2)`~Plot, data = richness)

# mod1 <- glm(`Aboveground_biomass (kg/m2)` ~ Plot, family=Gamma(link="log"), data=richness)
# mod2 <- glmer(`Aboveground_biomass (kg/m2)` ~ Plot + (1|Site), family=Gamma(link="log"), data=richness)
# (msAICc <- model.sel(mod1,mod2))
# summary(mod2)
# performance::r2(mod2)
# allEffects(mod2)
# plot(allEffects(mod2))

modlmer1 <- glmmTMB(`Aboveground_biomass (kg/m2)`~Plot,family=Gamma(link="log"),data=richness)
modlmer2 <- glmmTMB(`Aboveground_biomass (kg/m2)`~Plot+(1|Site),family=Gamma(link="log"),data=richness)
(msAICc <- model.sel(modlmer1,modlmer2))

modlmer3 <- glmmTMB(`Aboveground_biomass (kg/m2)`~Plot*Plant_richness+(1|Site),family=Gamma(link="log"),data=richness)
modlmer4 <- glmmTMB(`Aboveground_biomass (kg/m2)`~Plot*AMF_richness+(1|Site),family=Gamma(link="log"),data=richness)
modlmer5 <- glmmTMB(`Aboveground_biomass (kg/m2)`~Plot*Interaction_number+(1|Site),family=Gamma(link="log"),data=richness)
(msAICc <- model.sel(modlmer2,modlmer3,modlmer4,modlmer5))

plot(simulateResiduals(modlmer2))

summary(modlmer2)
performance::r2(modlmer2)
allEffects(modlmer2)
plot(allEffects(modlmer2))

Anova(modlmer2,type=2)


# Belowground_biomass
{
  Euadjacent.plot <- subset(richness,Plot=="Eu-adjacent")
  Eudistant.plot <- subset(richness,Plot=="Eu-distant")
  Euadjacent.plot$PlotID==Eudistant.plot$PlotID
  x <- Euadjacent.plot$`Belowground_biomass (kg/m2)`
  y <- Eudistant.plot$`Belowground_biomass (kg/m2)`
  t.test(x,y,paired=T)
  mean(x);mean(y)
  sd(x);sd(y)
  var(x);var(y)
  shapiro.test(x)
  shapiro.test(y)
}

leveneTest(`Belowground_biomass (kg/m2)`~Plot, data = richness)

modlmer1 <- glmmTMB(`Belowground_biomass (kg/m2)`~Plot,family=Gamma(link="log"),data=richness)
modlmer2 <- glmmTMB(`Belowground_biomass (kg/m2)`~Plot+(1|Site),family=Gamma(link="log"),data=richness)
(msAICc <- model.sel(modlmer1,modlmer2))

modlmer3 <- glmmTMB(`Belowground_biomass (kg/m2)`~Plot*Plant_richness+(1|Site),family=Gamma(link="log"),data=richness)
modlmer4 <- glmmTMB(`Belowground_biomass (kg/m2)`~Plot*AMF_richness+(1|Site),family=Gamma(link="log"),data=richness)
modlmer5 <- glmmTMB(`Belowground_biomass (kg/m2)`~Plot*Interaction_number+(1|Site),family=Gamma(link="log"),data=richness)
(msAICc <- model.sel(modlmer2,modlmer3,modlmer4,modlmer5))

plot(simulateResiduals(modlmer2))

summary(modlmer2)
performance::r2(modlmer2)
allEffects(modlmer2)
plot(allEffects(modlmer2))

Anova(modlmer2,type=2)

biomass <- richness[,c("Site","Plot","PlotID","Aboveground_biomass (kg/m2)","Belowground_biomass (kg/m2)")]
colnames(biomass) <- c("Site","Plot","PlotID","Aboveground_biomass","Belowground_biomass")
biomass_long <- pivot_longer(biomass, cols = c(Aboveground_biomass, Belowground_biomass), names_to="Type", values_to="biomass")

biomass_long[which(biomass_long$Type=="Aboveground_biomass"),]$Type <- "Above"
biomass_long[which(biomass_long$Type=="Belowground_biomass"),]$Type <- "Below"
biomass_long$Plot <- as.factor(biomass_long$Plot)
biomass_long$Type <- as.factor(biomass_long$Type)

# 计算均值和标准差
biomass_long.summary <- biomass_long %>%
  group_by(Plot, Type) %>%
  summarise(biomass.mean = mean(biomass),biomass.sd = sd(biomass))

p2 <- ggplot(biomass_long,aes(x=Plot,y=biomass,color=Type))+
  geom_boxplot(width=0.25,linewidth=0.25,outliers=F,alpha=1,show.legend=T,position=position_dodge(0.5))+
  stat_boxplot(geom="errorbar",width=0.25,linewidth=0.25,position=position_dodge(0.5))+
  geom_jitter(shape=16,size=0.4,position=position_jitterdodge(jitter.height=0,jitter.width=0.45,dodge.width=0.5),alpha=0.4,show.legend=T)+
  theme_classic()+
  scale_color_manual(values=c("#bc40bf","#4BA465","black"))+
  scale_fill_manual(values=c("#bc40bf","#4BA465","black"))+
  ylab("Biomass of plant")+
  # scale_x_discrete(breaks=c("Eu-adjacent","Eu-distant"),labels=c("Eu-adjacent", "Eu-distant"))+
  scale_y_continuous(labels = function(x) sprintf("%03d", x))+
  theme(axis.title.x=element_blank(),#修改X轴标题文本
        # axis.title.x=element_text(size=8),#修改X轴标题文本
        axis.title.y=element_text(size=8,angle=90),#修改y轴标题文本
        axis.text.y=element_text(size=6.5),#修改y轴刻度标签文本
        # axis.line.x = element_line(arrow = arrow(length = unit(0.2,"lines"))),
        axis.text.x=element_blank(),#修改x轴刻度标签文本
        # axis.text.x=element_text(size=6.5),#修改x轴刻度标签文本
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
  annotate("text",x=0.8,y=max(biomass_long$biomass)*0.95,hjust=1,vjust=1,label=expression(paste("kg/","m"^2)),size=2,color="black")+
  # annotate("text",x=1.9,y=max(biomass_long$biomass)*0.85,hjust=1,vjust=1,label=paste0("p < ",0.001),size=2,color="black")+
  annotate("text",x=1.5,y=max(biomass_long$biomass)*0.9,hjust=1,vjust=1,label=paste0("***"),size=2,color="black")+
  annotate("text",x=1.5+0.2,y=max(biomass_long$biomass)*1,hjust=1,vjust=1,label=paste0("***"),size=2,color="black")
p2

ggsave(filename="Fig4b_biomass_glmm.png",plot=p2,width=89/2,height=59/2,units="mm",dpi=900)
ggsave(filename="Fig4b_biomass_glmm.pdf",plot=p2,width=89/2,height=59/2,units="mm")

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
  shapiro.test(x)
  shapiro.test(y)
}

leveneTest(AMF_richness~Plot, data = richness)

modlmer1 <- glmmTMB(AMF_richness ~ Plot, family=poisson, data=richness)
modlmer2 <- glmmTMB(AMF_richness ~ Plot + (1|Site), family=poisson, data=richness)
(msAICc <- model.sel(modlmer1,modlmer2))

modlmer3 <- glmmTMB(AMF_richness ~ Plot + (1|Site), family=nbinom2, data=richness)
modlmer4 <- glmmTMB(AMF_richness ~ Plot + (1|Site), family=nbinom1, data=richness)
(msAICc <- model.sel(modlmer2,modlmer3,modlmer4))

plot(simulateResiduals(modlmer4))

summary(modlmer4)
performance::r2(modlmer4)
allEffects(modlmer4)
plot(allEffects(modlmer4))

Anova(modlmer4,type=2)

# {
#   Euadjacent.plot <- subset(richness,Plot=="Eu-adjacent")
#   Eudistant.plot <- subset(richness,Plot=="Eu-distant")
#   Euadjacent.plot$PlotID==Eudistant.plot$PlotID
#   x <- Euadjacent.plot$AMF_richness_Euphorbia_excluded
#   y <- Eudistant.plot$AMF_richness_Euphorbia_excluded
#   t.test(x,y,paired=T)
#   mean(x);mean(y)
#   sd(x);sd(y)
#   var(x);var(y)
# }
# 
# modlmer1 <- glmmTMB(AMF_richness_Euphorbia_excluded ~ Plot, family=poisson, data=richness)
# modlmer2 <- glmmTMB(AMF_richness_Euphorbia_excluded ~ Plot + (1|Site), family=poisson, data=richness)
# (msAICc <- model.sel(modlmer1,modlmer2))
# 
# modlmer3 <- glmmTMB(AMF_richness_Euphorbia_excluded ~ Plot + (1|Site), family=nbinom2, data=richness)
# modlmer4 <- glmmTMB(AMF_richness_Euphorbia_excluded ~ Plot + (1|Site), family=nbinom1, data=richness)
# (msAICc <- model.sel(modlmer2,modlmer3,modlmer4))
# 
# plot(simulateResiduals(modlmer4))
# 
# summary(modlmer4)
# performance::r2(modlmer4)
# allEffects(modlmer4)
# plot(allEffects(modlmer4))
# 
# richness.AMF <- richness[,c("Site","Plot","PlotID","AMF_richness","AMF_richness_Euphorbia_excluded")]
# richness.AMF_long <- pivot_longer(richness.AMF, cols = c(AMF_richness, AMF_richness_Euphorbia_excluded), names_to = "Type", values_to = "AMF_richness")
# 
# richness.AMF_long[which(richness.AMF_long$Type=="AMF_richness"),]$Type <- "included"
# richness.AMF_long[which(richness.AMF_long$Type=="AMF_richness_Euphorbia_excluded"),]$Type <- "excluded"
# richness.AMF_long$Type <- factor(richness.AMF_long$Type,levels=c("included","excluded"))
# 
# # 计算均值和标准差
# richness.AMF_long.summary <- richness.AMF_long %>%
#   group_by(Plot, Type) %>%
#   summarise(AMF_richness.mean = mean(AMF_richness),AMF_richness.sd = sd(AMF_richness))
# 
# p3 <- ggplot(richness.AMF_long,aes(x=Plot,y=AMF_richness,,color=Type))+
#   geom_boxplot(width=0.25,linewidth=0.25,outliers=F,alpha=1,show.legend=T,position=position_dodge(0.5))+
#   stat_boxplot(geom="errorbar",width=0.25,linewidth=0.25,position=position_dodge(0.5))+
#   geom_jitter(shape=16,size=0.4,position=position_jitterdodge(jitter.height=0,jitter.width=0.4,dodge.width=0.5),alpha=0.4,show.legend=T)+
#   theme_classic()+
#   scale_color_manual(values=c("#bc40bf","#4BA465","black"))+
#   scale_fill_manual(values=c("#bc40bf","#4BA465","black"))+
#   ylab("AMF richness")+
#   # scale_x_discrete(breaks=c("Eu-adjacent","Eu-distant"),labels=c("Eu-adjacent", "Eu-distant"))+
#   theme(axis.title.x=element_blank(),#修改X轴标题文本
#         # axis.title.x=element_text(size=8),#修改X轴标题文本
#         axis.title.y=element_text(size=8,angle=90),#修改y轴标题文本
#         axis.text.y=element_text(size=6.5),#修改y轴刻度标签文本
#         # axis.line.x = element_line(arrow = arrow(length = unit(0.2,"lines"))),
#         axis.text.x=element_text(size=6.5),#修改x轴刻度标签文本
#         panel.grid=element_blank(),#隐藏网格线
#         # legend.title=element_text(size=6.5),#修改图例标题大小
#         legend.title=element_blank(),#修改图例标题大小
#         legend.text=element_text(size=5.5),#修改图例文本大小
#         legend.key.size=unit(0.1,"lines"),#缩小图例符号的大小
#         legend.background=element_blank(),#设置图例背景为空白
#         legend.position="inside",#修改图例位置
#         # legend.justification=c("right","top"),#修改图例位置
#         legend.position.inside=c(0.9,0.9),#修改图例位置
#         legend.margin=margin(0,0,0,0))+#修改图例位置
#   # annotate("text",x=1.9,y=max(richness.AMF_long$AMF_richness)*0.85,hjust=1,vjust=1,label=paste0("p > ",0.05),size=2,color="black")+
#   annotate("text",x=1.5,y=(max(richness.AMF_long$AMF_richness)+(max(richness.AMF_long$AMF_richness)-min(richness.AMF_long$AMF_richness))*0.2),hjust=1,vjust=1,label=paste0("NS"),size=2,color="black")+
#   annotate("text",x=1.5+0.2,y=(max(richness.AMF_long$AMF_richness)+(max(richness.AMF_long$AMF_richness)-min(richness.AMF_long$AMF_richness))*0.1),hjust=1,vjust=1,label=paste0("NS"),size=2,color="black")
# p3

p3 <- ggplot(richness,aes(x=Plot,y=AMF_richness,color=Plot))+
  geom_boxplot(width=0.15,linewidth=0.25,outliers=F,alpha=1,show.legend=T,position=position_dodge(0.5))+
  stat_boxplot(geom="errorbar",width=0.15,linewidth=0.25,position=position_dodge(0.5))+
  geom_jitter(shape=16,size=0.4,position=position_jitterdodge(jitter.height=0.25,jitter.width=0.45,dodge.width=0.5),alpha=0.4,show.legend=T)+
  # geom_point(aes(x=Plot,y=mean(Plant_richness)),shape=16,size=0.8,color="black",position=position_dodge(0.5),alpha=1)+
  # geom_errorbar(aes(x=Plot,y=mean(Plant_richness),ymin=mean(Plant_richness)-sd(Plant_richness), ymax=mean(Plant_richness)+sd(Plant_richness)),width=0,linewidth=0.4,color="black",position=position_dodge(0.5),alpha=1)+
  # scale_color_manual(values=c("#bc40bf","#bc40bf","black"))+
  # scale_fill_manual(values=c("#bc40bf","#bc40bf"))+ 
  scale_color_manual(values=c("gray40","gray40","gray40"))+
  scale_fill_manual(values=c("gray40","gray40"))+ 
  theme_classic()+
  ylab("AMF richness")+
  # scale_x_discrete(breaks=c("Eu-adjacent","Eu-distant"),labels=c("Eu-adjacent", "Eu-distant"))+
  # scale_y_continuous(labels = function(x) sprintf("%03d", x))+
  theme(axis.title.x=element_blank(),#修改X轴标题文本
        # axis.title.x=element_text(size=8),#修改X轴标题文本
        # axis.line.x = element_line(arrow = arrow(length = unit(0.2,"lines"))),
        axis.title.y=element_text(size=8,angle=90),#修改y轴标题文本
        axis.text.y=element_text(size=6.5),#修改y轴刻度标签文本
        # axis.text.x=element_blank(),#修改x轴刻度标签文本
        axis.text.x=element_text(size=6.5),#修改x轴刻度标签文本
        panel.grid=element_blank(),#隐藏网格线
        legend.title=element_text(size=8),#修改图例标题大小
        legend.text=element_text(size=6.5),#修改图例文本大小
        legend.position="none")+#修改图例位置
  # annotate("text",x=1.7,y=max(richness$AMF_richness),hjust=1,vjust=1,label=paste0("p > ",0.05),size=2)+
  annotate("text",x=1.5,y=(max(richness$AMF_richness)+(max(richness$AMF_richness)-min(richness$AMF_richness))*0.1),hjust=1,vjust=1,label=paste0("NS"),size=2)
p3

ggsave(filename="Fig4c_AMF_richness_glmm.png",plot=p3,width=89/2,height=59/2,units="mm",dpi=900)
ggsave(filename="Fig4c_AMF_richness_glmm.pdf",plot=p3,width=89/2,height=59/2,units="mm")

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
  shapiro.test(x)
  shapiro.test(y)
}

leveneTest(Interaction_number~Plot, data = richness)

modlmer1 <- glmmTMB(Interaction_number ~ Plot, family=poisson, data=richness)
modlmer2 <- glmmTMB(Interaction_number ~ Plot + (1|Site), family=poisson, data=richness)
(msAICc <- model.sel(modlmer1,modlmer2))

modlmer3 <- glmmTMB(Interaction_number ~ Plot + (1|Site), family=nbinom2, data=richness)
modlmer4 <- glmmTMB(Interaction_number ~ Plot + (1|Site), family=nbinom1, data=richness)
(msAICc <- model.sel(modlmer2,modlmer3,modlmer4))

plot(simulateResiduals(modlmer4))

summary(modlmer4)
performance::r2(modlmer4)
allEffects(modlmer4)
plot(allEffects(modlmer4))

Anova(modlmer4,type=2)

# {
#   Euadjacent.plot <- subset(richness,Plot=="Eu-adjacent")
#   Eudistant.plot <- subset(richness,Plot=="Eu-distant")
#   Euadjacent.plot$PlotID==Eudistant.plot$PlotID
#   x <- Euadjacent.plot$Interaction_number_Euphorbia_excluded
#   y <- Eudistant.plot$Interaction_number_Euphorbia_excluded
#   t.test(x,y,paired=T)
#   mean(x);mean(y)
#   sd(x);sd(y)
#   var(x);var(y)
# }
# 
# modlmer1 <- glmmTMB(Interaction_number_Euphorbia_excluded ~ Plot, family=poisson, data=richness)
# modlmer2 <- glmmTMB(Interaction_number_Euphorbia_excluded ~ Plot + (1|Site), family=poisson, data=richness)
# (msAICc <- model.sel(modlmer1,modlmer2))
# 
# modlmer3 <- glmmTMB(Interaction_number_Euphorbia_excluded ~ Plot + (1|Site), family=nbinom2, data=richness)
# modlmer4 <- glmmTMB(Interaction_number_Euphorbia_excluded ~ Plot + (1|Site), family=nbinom1, data=richness)
# (msAICc <- model.sel(modlmer2,modlmer3,modlmer4))
# 
# plot(simulateResiduals(modlmer4))
# 
# summary(modlmer4)
# performance::r2(modlmer4)
# allEffects(modlmer4)
# plot(allEffects(modlmer4))
# 
# richness.interaction <- richness[,c("Site","Plot","PlotID","Interaction_number","Interaction_number_Euphorbia_excluded")]
# richness.interaction_long <- pivot_longer(richness.interaction, cols = c(Interaction_number, Interaction_number_Euphorbia_excluded), names_to = "Type", values_to = "Interaction_number")
# 
# richness.interaction_long[which(richness.interaction_long$Type=="Interaction_number"),]$Type <- "included"
# richness.interaction_long[which(richness.interaction_long$Type=="Interaction_number_Euphorbia_excluded"),]$Type <- "excluded"
# richness.interaction_long$Type <- factor(richness.interaction_long$Type,levels=c("included","excluded"))
# 
# # 计算均值和标准差
# richness.interaction_long.summary <- richness.interaction_long %>%
#   group_by(Plot, Type) %>%
#   summarise(Interaction_number.mean = mean(Interaction_number),Interaction_number.sd = sd(Interaction_number))
# 
# p4 <- ggplot(richness.interaction_long,aes(x=Plot,y=Interaction_number,color=Type))+
#   geom_boxplot(width=0.25,linewidth=0.25,outliers=F,alpha=1,show.legend=T,position=position_dodge(0.5))+
#   stat_boxplot(geom="errorbar",width=0.25,linewidth=0.25,position=position_dodge(0.5))+
#   geom_jitter(shape=16,size=0.4,position=position_jitterdodge(jitter.height=0,jitter.width=0.4,dodge.width=0.5),alpha=0.4,show.legend=T)+
#   theme_classic()+
#   scale_color_manual(values=c("#bc40bf","#4BA465","black"))+
#   scale_fill_manual(values=c("#bc40bf","#4BA465","black"))+
#   ylab("Number of interaction")+
#   # scale_x_discrete(breaks=c("Eu-adjacent","Eu-distant"),labels=c("Eu-adjacent", "Eu-distant"))+
#   theme(axis.title.x=element_blank(),#修改X轴标题文本
#         # axis.title.x=element_text(size=8),#修改X轴标题文本
#         axis.title.y=element_text(size=8,angle=90),#修改y轴标题文本
#         axis.text.y=element_text(size=6.5),#修改y轴刻度标签文本
#         # axis.line.x = element_line(arrow = arrow(length = unit(0.2,"lines"))),
#         axis.text.x=element_text(size=6.5),#修改x轴刻度标签文本
#         panel.grid=element_blank(),#隐藏网格线
#         # legend.title=element_text(size=6.5),#修改图例标题大小
#         legend.title=element_blank(),#修改图例标题大小
#         legend.text=element_text(size=5.5),#修改图例文本大小
#         legend.key.size=unit(0.1,"lines"),#缩小图例符号的大小
#         legend.background=element_blank(),#设置图例背景为空白
#         legend.position="inside",#修改图例位置
#         # legend.justification=c("right","top"),#修改图例位置
#         legend.position.inside=c(0.9,0.9),#修改图例位置
#         legend.margin=margin(0,0,0,0))+#修改图例位置
#   # annotate("text",x=2,y=max(richness.interaction_long$Interaction_number)*0.85,hjust=1,vjust=1,label=paste0("p > 0.05"," / ","p < 0.01"),size=2,color="black")+
#   annotate("text",x=1.5,y=(max(richness.interaction_long$Interaction_number)+(max(richness.interaction_long$Interaction_number)-min(richness.interaction_long$Interaction_number))*0.1),hjust=1,vjust=1,label=paste0("NS"),size=2,color="black")+
#   annotate("text",x=1.5+0.2,y=(max(richness.interaction_long$Interaction_number)),hjust=1,vjust=1,label=paste0("**"),size=2,color="black")
# p4

p4 <- ggplot(richness,aes(x=Plot,y=Interaction_number,color=Plot))+
  geom_boxplot(width=0.15,linewidth=0.25,outliers=F,alpha=1,show.legend=T,position=position_dodge(0.5))+
  stat_boxplot(geom="errorbar",width=0.15,linewidth=0.25,position=position_dodge(0.5))+
  geom_jitter(shape=16,size=0.4,position=position_jitterdodge(jitter.height=0.25,jitter.width=0.45,dodge.width=0.5),alpha=0.4,show.legend=T)+
  # geom_point(aes(x=Plot,y=mean(Plant_richness)),shape=16,size=0.8,color="black",position=position_dodge(0.5),alpha=1)+
  # geom_errorbar(aes(x=Plot,y=mean(Plant_richness),ymin=mean(Plant_richness)-sd(Plant_richness), ymax=mean(Plant_richness)+sd(Plant_richness)),width=0,linewidth=0.4,color="black",position=position_dodge(0.5),alpha=1)+
  # scale_color_manual(values=c("#bc40bf","#bc40bf","black"))+
  # scale_fill_manual(values=c("#bc40bf","#bc40bf"))+ 
  scale_color_manual(values=c("gray40","gray40","gray40"))+
  scale_fill_manual(values=c("gray40","gray40"))+ 
  theme_classic()+
  ylab("Number of interaction")+
  # scale_x_discrete(breaks=c("Eu-adjacent","Eu-distant"),labels=c("Eu-adjacent", "Eu-distant"))+
  # scale_y_continuous(labels = function(x) sprintf("%03d", x))+
  theme(axis.title.x=element_blank(),#修改X轴标题文本
        # axis.title.x=element_text(size=8),#修改X轴标题文本
        # axis.line.x = element_line(arrow = arrow(length = unit(0.2,"lines"))),
        axis.title.y=element_text(size=8,angle=90),#修改y轴标题文本
        axis.text.y=element_text(size=6.5),#修改y轴刻度标签文本
        # axis.text.x=element_blank(),#修改x轴刻度标签文本
        axis.text.x=element_text(size=6.5),#修改x轴刻度标签文本
        panel.grid=element_blank(),#隐藏网格线
        legend.title=element_text(size=8),#修改图例标题大小
        legend.text=element_text(size=6.5),#修改图例文本大小
        legend.position="none")+#修改图例位置
  # annotate("text",x=1.7,y=max(richness$Interaction_number),hjust=1,vjust=1,label=paste0("p > ",0.05),size=2)+
  annotate("text",x=1.5,y=(max(richness$Interaction_number)+(max(richness$Interaction_number)-min(richness$Interaction_number))*0.1),hjust=1,vjust=1,label=paste0("NS"),size=2)
p4

ggsave(filename="Fig4d_interaction_number_glmm.png",plot=p4,width=89/2,height=59/2,units="mm",dpi=900)
ggsave(filename="Fig4d_interaction_number_glmm.pdf",plot=p4,width=89/2,height=59/2,units="mm")

p_all <- plot_grid(p1,p2,p3,p4,labels=c("a","b","c","d"),label_size=8,label_x=0.25,label_y=1,ncol=2,nrow=2)
p_all
ggsave(filename="Fig4_0127.png",plot=p_all,width=89,height=59,units="mm",dpi=900)
ggsave(filename="Fig4_0127.pdf",plot=p_all,width=89,height=59,units="mm")

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
shapiro.test(x)
shapiro.test(y)
levene_test <- leveneTest(Observe ~ Components, data = plant_beta)
levene_test
t.test(x,y,paired=T)

d <- x-y
# 进行Welch's t检验
result_welch <- t.test(d, mu = 0, alternative = "two.sided", var.equal = FALSE)
result_welch

# 计算均值和标准差
plant_beta_combined <- plant_beta %>%
  group_by(Beta_diversity_type, Components) %>%
  summarise(Observe.mean = mean(Observe),Observe.sd = sd(Observe))

p1 <- ggplot(plant_beta,aes(x=Components,y=Observe,color=Components))+
  geom_boxplot(width=0.225,linewidth=0.25,outliers=F,alpha=1,show.legend=T,position=position_dodge(0.5))+
  stat_boxplot(geom="errorbar",width=0.225,linewidth=0.25,position=position_dodge(0.5))+
  geom_jitter(shape=16,size=0.6,position=position_jitterdodge(jitter.height=0,jitter.width=0.8,dodge.width=0.5),alpha=0.6,show.legend=T)+
  scale_color_manual(values=c("gray40","gray40","gray40"))+
  scale_fill_manual(values=c("gray40","gray40","gray40"))+
  theme_classic()+
  ylab("Plant beta diversity")+# 设置Y轴标题
  geom_vline(xintercept=1.5,linetype=2,linewidth=0.2,color="gray50",alpha=0.8)+
  theme(axis.title.x=element_blank(),#修改X轴标题文本
        # axis.title.x=element_text(size=8),#修改X轴标题文本
        axis.title.y=element_text(size=8,angle=90),#修改y轴标题文本
        axis.text.y=element_text(size=6.5),#修改y轴刻度标签文本
        axis.text.x=element_text(size=6.5),#修改x轴刻度标签文本
        panel.grid=element_blank(),#隐藏网格线
        legend.title=element_text(size=8),#修改图例标题大小
        legend.text=element_text(size=6.5),#修改图例文本大小
        legend.position="none")+#修改图例位置
  # annotate("text",x=3,y=max(plant_beta$Observe)*0.85,hjust=1,vjust=1,label=paste0("p < ",0.001),size=2,color="black")+
  annotate("text",x=2.5,y=(max(plant_beta$Observe)+(max(plant_beta$Observe)-min(plant_beta$Observe))*0.1),hjust=1,vjust=1,label=paste0("***"),size=2,color="black")
p1
ggsave(filename="Fig5a_plant_beta_diversity.png",plot=p1,width=89/2,height=59/2,units="mm",dpi=900)
ggsave(filename="Fig5a_plant_beta_diversity.pdf",plot=p1,width=89/2,height=59/2,units="mm")

#####

AMF_beta <- subset(Beta_diversity,Beta_diversity_type=="AMF_beta_diversity")
x <- subset(AMF_beta,Components=="Turnover")$Observe_Euphorbia_excluded
y <- subset(AMF_beta,Components=="Nestedness")$Observe_Euphorbia_excluded
shapiro.test(x)
shapiro.test(y)
levene_test <- leveneTest(Observe ~ Components, data = AMF_beta)
levene_test
t.test(x,y,paired=T)


# AMF_beta_long <- pivot_longer(AMF_beta, cols = c(Observe, Observe_Euphorbia_excluded), names_to = "Type", values_to = "Value")
# AMF_beta_long[which(AMF_beta_long$Type=="Observe"),]$Type <- "included"
# AMF_beta_long[which(AMF_beta_long$Type=="Observe_Euphorbia_excluded"),]$Type <- "excluded"
# AMF_beta_long$Type <- factor(AMF_beta_long$Type,levels=c("included","excluded"))
# 
# # 计算均值和标准差
# AMF_beta_long.summary <- AMF_beta_long %>%
#   group_by(Components, Type) %>%
#   summarise(Value.mean = mean(Value),Value.sd = sd(Value))
# 
# p2 <- ggplot(AMF_beta_long,aes(x=Components,y=Value,color=Type))+
#   geom_boxplot(width=0.4,linewidth=0.25,outliers=F,alpha=1,show.legend=T,position=position_dodge(0.5))+
#   stat_boxplot(geom="errorbar",width=0.4,linewidth=0.25,position=position_dodge(0.5))+
#   geom_jitter(shape=16,size=0.6,position=position_jitterdodge(jitter.height=0,jitter.width=0.45,dodge.width=0.5),alpha=0.6,show.legend=T)+
#   scale_color_manual(values=c("#bc40bf","#4BA465","black"))+
#   scale_fill_manual(values=c("#bc40bf","#4BA465","black"))+
#   theme_classic()+
#   ylab("AMF beta diversity")+# 设置Y轴标题
#   geom_vline(xintercept=1.5,linetype=2,linewidth=0.2,color="gray50",alpha=0.8)+
#   theme(axis.title.x=element_blank(),#修改X轴标题文本
#         # axis.title.x=element_text(size=8),#修改X轴标题文本
#         axis.title.y=element_text(size=8,angle=90),#修改y轴标题文本
#         axis.text.y=element_text(size=6.5),#修改y轴刻度标签文本
#         axis.text.x=element_text(size=6.5),#修改x轴刻度标签文本
#         panel.grid=element_blank(),#隐藏网格线
#         # legend.title=element_text(size=6.5),#修改图例标题大小
#         legend.title=element_blank(),
#         legend.text=element_text(size=5.5),#修改图例文本大小
#         legend.key.size=unit(0.1,"lines"),#缩小图例符号的大小
#         legend.background=element_blank(),#设置图例背景为空白
#         legend.position="inside",#修改图例位置
#         # legend.justification=c("right","top"),#修改图例位置
#         legend.position.inside=c(0.9,0.9),#修改图例位置
#         legend.margin=margin(0,0,0,0))+#修改图例位置
#   # annotate("text",x=3,y=max(AMF_beta_long$Value)*0.85,hjust=1,vjust=1,label=paste0("p < ",0.001),size=2,color="black")+
#   annotate("text",x=2.5,y=(max(AMF_beta_long$Value)+(max(AMF_beta_long$Value)-min(AMF_beta_long$Value))*0.1),hjust=1,vjust=1,label=paste0("***"),size=2,color="black")+
#   annotate("text",x=2.5+0.2,y=(max(AMF_beta_long$Value)),hjust=1,vjust=1,label=paste0("***"),size=2,color="black")
# p2

p2 <- ggplot(AMF_beta,aes(x=Components,y=Observe,color=Components))+
  geom_boxplot(width=0.225,linewidth=0.25,outliers=F,alpha=1,show.legend=T,position=position_dodge(0.5))+
  stat_boxplot(geom="errorbar",width=0.225,linewidth=0.25,position=position_dodge(0.5))+
  geom_jitter(shape=16,size=0.6,position=position_jitterdodge(jitter.height=0,jitter.width=0.8,dodge.width=0.5),alpha=0.6,show.legend=T)+
  scale_color_manual(values=c("gray40","gray40","gray40"))+
  scale_fill_manual(values=c("gray40","gray40","gray40"))+
  theme_classic()+
  ylab("AMF beta diversity")+# 设置Y轴标题
  geom_vline(xintercept=1.5,linetype=2,linewidth=0.2,color="gray50",alpha=0.8)+
  theme(axis.title.x=element_blank(),#修改X轴标题文本
        # axis.title.x=element_text(size=8),#修改X轴标题文本
        axis.title.y=element_text(size=8,angle=90),#修改y轴标题文本
        axis.text.y=element_text(size=6.5),#修改y轴刻度标签文本
        axis.text.x=element_text(size=6.5),#修改x轴刻度标签文本
        panel.grid=element_blank(),#隐藏网格线
        legend.title=element_text(size=8),#修改图例标题大小
        legend.text=element_text(size=6.5),#修改图例文本大小
        legend.position="none")+#修改图例位置
  # annotate("text",x=3,y=max(plant_beta$Observe)*0.85,hjust=1,vjust=1,label=paste0("p < ",0.001),size=2,color="black")+
  annotate("text",x=2.5,y=(max(plant_beta$Observe)+(max(plant_beta$Observe)-min(plant_beta$Observe))*0.1),hjust=1,vjust=1,label=paste0("***"),size=2,color="black")
p2

ggsave(filename="Fig5b_AMF_beta_diversity.png",plot=p2,width=89/2,height=59/2,units="mm",dpi=900)
ggsave(filename="Fig5b_AMF_beta_diversity.pdf",plot=p2,width=89/2,height=59/2,units="mm")

#####
interaction_beta <- subset(Beta_diversity,Beta_diversity_type=="Interaction_beta_diversity")
x <- subset(interaction_beta,Components=="βOS")$Observe
y <- subset(interaction_beta,Components=="βST")$Observe
shapiro.test(x)
shapiro.test(y)
levene_test <- leveneTest(Observe ~ Components, data = interaction_beta)
levene_test
t.test(x,y,paired=T)

d <- x-y
# 进行Welch's t检验
result_welch <- t.test(d, mu = 0, alternative = "two.sided", var.equal = FALSE)
result_welch


x <- subset(interaction_beta,Components=="βOS")$Observe_Euphorbia_excluded
y <- subset(interaction_beta,Components=="βST")$Observe_Euphorbia_excluded
shapiro.test(x)
shapiro.test(y)
levene_test <- leveneTest(Observe ~ Components, data = interaction_beta)
levene_test
t.test(x,y,paired=T)

d <- x-y
# 进行Welch's t检验
result_welch <- t.test(d, mu = 0, alternative = "two.sided", var.equal = FALSE)
result_welch


x <- subset(interaction_beta,Components=="βST.h")$Observe
y <- subset(interaction_beta,Components=="βST.l")$Observe
z <- subset(interaction_beta,Components=="βST.lh")$Observe
shapiro.test(x)
shapiro.test(y)
shapiro.test(z)
levene_test <- leveneTest(Observe ~ Components, data = interaction_beta)
levene_test
t.test(x,y,paired=T)
aov(x,y,z)

d <- x-y
# 进行Welch's t检验
result_welch <- t.test(d, mu = 0, alternative = "two.sided", var.equal = FALSE)
result_welch

interaction_beta2 <- subset(interaction_beta,Components %in% c("βST.h","βST.l","βST.lh"))[,c(1,4,5)]
interaction_beta2 <- as.data.frame(interaction_beta2)
# # 进行Friedman检验
friedman_result <- friedman.test(Observe~Components|Site,data=interaction_beta2)
friedman_result <- friedman.test(Observe~Site|Components,data=interaction_beta2)

# 使用aov函数
anova_result <- aov(Observe~Components, data = interaction_beta2)
summary(anova_result)

# x:"βST.h";y:"βST.l";z:"βST.lh"
wilcox.test(x, y, paired = TRUE);t.test(x,y,paired=T)
# V = 55, p-value = 0.001953;t = 5.6994, df = 9, p-value = 0.0002945
wilcox.test(x, z, paired = TRUE);t.test(x,z,paired=T)
# V = 55, p-value = 0.001953;t = 9.175, df = 9, p-value = 7.292e-06
wilcox.test(y, z, paired = TRUE);t.test(y,z,paired=T)
# V = 51, p-value = 0.01367;t = 3.5067, df = 9, p-value = 0.006653





x <- subset(interaction_beta,Components=="βST.h")$Observe_Euphorbia_excluded
y <- subset(interaction_beta,Components=="βST.lh")$Observe_Euphorbia_excluded
z <- subset(interaction_beta,Components=="βST.lh")$Observe_Euphorbia_excluded
shapiro.test(x)
shapiro.test(y)
shapiro.test(z)
levene_test <- leveneTest(Observe ~ Components, data = interaction_beta)
levene_test
t.test(x,y,paired=T)

d <- x-y
# 进行Welch's t检验
result_welch <- t.test(d, mu = 0, alternative = "two.sided", var.equal = FALSE)
result_welch

# interaction_beta_long <- pivot_longer(interaction_beta, cols = c(Observe, Observe_Euphorbia_excluded), names_to = "Type", values_to = "Value")
# interaction_beta_long[which(interaction_beta_long$Type=="Observe"),]$Type <- "included"
# interaction_beta_long[which(interaction_beta_long$Type=="Observe_Euphorbia_excluded"),]$Type <- "excluded"
# interaction_beta_long$Type <- factor(interaction_beta_long$Type,levels=c("included","excluded"))
# 
# # 计算均值和标准差
# interaction_beta_long.summary <- interaction_beta_long %>%
#   group_by(Components, Type) %>%
#   summarise(Value.mean = mean(Value),Value.sd = sd(Value))
# 
# p3 <- ggplot(interaction_beta_long,aes(x=Components,y=Value,color=Type))+
#   geom_boxplot(width=0.4,linewidth=0.25,outliers=F,alpha=1,show.legend=T,position=position_dodge(0.5))+
#   stat_boxplot(geom="errorbar",width=0.4,linewidth=0.25,position=position_dodge(0.5))+
#   geom_jitter(shape=16,size=0.6,position=position_jitterdodge(jitter.height=0,jitter.width=0.45,dodge.width=0.5),alpha=0.6,show.legend=T)+
#   scale_color_manual(values=c("#bc40bf","#4BA465","black"))+
#   scale_fill_manual(values=c("#bc40bf","#4BA465","black"))+
#   theme_classic()+
#   ylab("Interaction beta diversity")+# 设置Y轴标题
#   geom_vline(xintercept=1.5,linetype=2,linewidth=0.2,color="gray50",alpha=0.8)+
#   geom_vline(xintercept=3.5,linetype=2,linewidth=0.2,color="gray50",alpha=0.8)+
#   theme(axis.title.x=element_blank(),#修改X轴标题文本
#         # axis.title.x=element_text(size=8),#修改X轴标题文本
#         axis.title.y=element_text(size=8,angle=90),#修改y轴标题文本
#         axis.text.y=element_text(size=6.5),#修改y轴刻度标签文本
#         axis.text.x=element_text(size=6.5),#修改x轴刻度标签文本
#         panel.grid=element_blank(),#隐藏网格线
#         # legend.title=element_text(size=6.5),#修改图例标题大小
#         legend.title=element_blank(),
#         legend.text=element_text(size=5.5),#修改图例文本大小
#         legend.key.size=unit(0.1,"lines"),#缩小图例符号的大小
#         legend.background=element_blank(),#设置图例背景为空白
#         legend.position="inside",#修改图例位置
#         # legend.justification=c("right","top"),#修改图例位置
#         legend.position.inside=c(0.92,0.9),#修改图例位置
#         legend.margin=margin(0,0,0,0))+#修改图例位置
#   # annotate("text",x=2.6,y=max(interaction_beta_long$Value)*0.85,hjust=1,vjust=1,label=paste0("p > ",0.05),size=2,color="black")+
#   # annotate("text",x=5.3,y=max(interaction_beta_long$Value)*0.85/1.2,hjust=1,vjust=1,label=paste0("p < ",0.001),size=2,color="black")+
#   # annotate("text",x=4.8,y=max(interaction_beta_long$Value)*0.85/1.5,hjust=1,vjust=1,label=paste0("p < ",0.001),size=2,color="black")+
#   # annotate("text",x=5.8,y=max(interaction_beta_long$Value)*0.85/4,hjust=1,vjust=1,label=paste0("p < ",0.001),size=2,color="black")
#   annotate("text",x=2.5,y=max(interaction_beta_long$Value)*1,hjust=1,vjust=1,label=paste0("NS"),size=2,color="black")+
#   annotate("text",x=2.5+0.2,y=max(interaction_beta_long$Value)*0.9,hjust=1,vjust=1,label=paste0("NS"),size=2,color="black")+
#   annotate("text",x=5,y=max(interaction_beta_long$Value)*1,hjust=1,vjust=1,label=paste0("***"),size=2,color="black")+
#   annotate("text",x=5+0.2,y=max(interaction_beta_long$Value)*0.9,hjust=1,vjust=1,label=paste0("***"),size=2,color="black")+
#   annotate("text",x=4.5,y=max(interaction_beta_long$Value)*0.8,hjust=1,vjust=1,label=paste0("***"),size=2,color="black")+
#   annotate("text",x=4.5+0.2,y=max(interaction_beta_long$Value)*0.7,hjust=1,vjust=1,label=paste0("***"),size=2,color="black")+
#   annotate("text",x=5.5,y=max(interaction_beta_long$Value)*0.4,hjust=1,vjust=1,label=paste0("***"),size=2,color="black")+
#   annotate("text",x=5.5+0.2,y=max(interaction_beta_long$Value)*0.3,hjust=1,vjust=1,label=paste0("***"),size=2,color="black")
# p3
# ggsave(filename="Fig5c_interaction_beta_diversity.png",plot=p3,width=89,height=59/2,units="mm",dpi=900)
# ggsave(filename="Fig5c_interaction_beta_diversity.pdf",plot=p3,width=89,height=59/2,units="mm")

p3 <- ggplot(interaction_beta,aes(x=Components,y=Observe,color=Components))+
  geom_boxplot(width=0.225,linewidth=0.25,outliers=F,alpha=1,show.legend=T,position=position_dodge(0.5))+
  stat_boxplot(geom="errorbar",width=0.225,linewidth=0.25,position=position_dodge(0.5))+
  geom_jitter(shape=16,size=0.6,position=position_jitterdodge(jitter.height=0,jitter.width=0.8,dodge.width=0.5),alpha=0.6,show.legend=T)+
  scale_color_manual(values=c("gray40","gray40","gray40","gray40","gray40","gray40","black"))+
  scale_fill_manual(values=c("gray40","gray40","gray40","gray40","gray40","gray40","black"))+
  theme_classic()+
  ylab("Interaction beta diversity")+# 设置Y轴标题
  geom_vline(xintercept=1.5,linetype=2,linewidth=0.2,color="gray50",alpha=0.8)+
  geom_vline(xintercept=3.5,linetype=2,linewidth=0.2,color="gray50",alpha=0.8)+
  theme(axis.title.x=element_blank(),#修改X轴标题文本
        # axis.title.x=element_text(size=8),#修改X轴标题文本
        axis.title.y=element_text(size=8,angle=90),#修改y轴标题文本
        axis.text.y=element_text(size=6.5),#修改y轴刻度标签文本
        axis.text.x=element_text(size=6.5),#修改x轴刻度标签文本
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
  annotate("text",x=2.5,y=max(interaction_beta$Observe)*0.9,hjust=1,vjust=1,label=paste0("NS"),size=2,color="black")+
  annotate("text",x=5,y=max(interaction_beta$Observe)*0.82,hjust=1,vjust=1,label=paste0("***"),size=2,color="black")+
  annotate("text",x=4.5,y=max(interaction_beta$Observe)*0.7,hjust=1,vjust=1,label=paste0("***"),size=2,color="black")+
  annotate("text",x=5.5,y=max(interaction_beta$Observe)*0.3,hjust=1,vjust=1,label=paste0("***"),size=2,color="black")
p3

p4 <- plot_grid(p1,p2,labels=c("a","b"),label_size=8,label_x=0.25,label_y=1,ncol=2,nrow=1)
p4
p_all <- plot_grid(p4,p3,labels=c("","c"),label_size=8,label_x=0.125,label_y=1,ncol=1,nrow=2)
p_all
ggsave(filename="Fig5_0127.png",plot=p_all,width=89,height=59,units="mm",dpi=900)
ggsave(filename="Fig5_0127.pdf",plot=p_all,width=89,height=59,units="mm")

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
shapiro.test(x)
shapiro.test(y)
levene_test <- leveneTest(Observe ~ Plot, data = WCON)
levene_test
t.test(x,y,paired=T)

# x <- subset(WCON,Plot=="Eu-adjacent")$Observe_Euphorbia_excluded
# y <- subset(WCON,Plot=="Eu-distant")$Observe_Euphorbia_excluded
# t.test(x,y,paired=T)
# 
# WCON_long <- pivot_longer(WCON, cols = c(Observe, Observe_Euphorbia_excluded), names_to = "Type", values_to = "Value")
# WCON_long[which(WCON_long$Type=="Observe"),]$Type <- "included"
# WCON_long[which(WCON_long$Type=="Observe_Euphorbia_excluded"),]$Type <- "excluded"
# WCON_long$Type <- factor(WCON_long$Type,levels=c("included","excluded"))
# 
# # 计算均值和标准差
# WCON_long.summary <- WCON_long %>%
#   group_by(Plot, Type) %>%
#   summarise(Value.mean = mean(Value),Value.sd = sd(Value))

# p1 <- ggplot(WCON_long,aes(x=Plot,y=Value,color=Type))+
#   geom_boxplot(width=0.4,linewidth=0.25,outliers=F,alpha=1,show.legend=T,position=position_dodge(0.5))+
#   stat_boxplot(geom="errorbar",width=0.4,linewidth=0.25,position=position_dodge(0.5))+
#   geom_jitter(shape=16,size=0.6,position=position_jitterdodge(jitter.height=0,jitter.width=0.45,dodge.width=0.5),alpha=0.6,show.legend=T)+
#   theme_classic()+
#   scale_color_manual(values=c("#bc40bf","#4BA465","black"))+
#   scale_fill_manual(values=c("#bc40bf","#4BA465","black"))+
#   ylab("Weighted connectance")+# 设置Y轴标题
#   # scale_x_discrete(breaks=c("Eu-adjacent","Eu-distant"),labels=c("Eu-adjacent", "Eu-distant"))+
#   theme(axis.title.x=element_blank(),#修改X轴标题文本
#         # axis.title.x=element_text(size=8),#修改X轴标题文本
#         axis.title.y=element_text(size=8,angle=90),#修改y轴标题文本
#         axis.text.y=element_text(size=6.5),#修改y轴刻度标签文本
#         # axis.line.x = element_line(arrow = arrow(length = unit(0.2,"lines"))),
#         axis.text.x=element_blank(),#修改x轴刻度标签文本
#         # axis.text.x=element_text(size=6.5),#修改x轴刻度标签文本
#         panel.grid=element_blank(),#隐藏网格线
#         # legend.title=element_text(size=6.5),#修改图例标题大小
#         legend.title=element_blank(),#修改图例标题大小
#         legend.text=element_text(size=5.5),#修改图例文本大小
#         legend.key.size=unit(0.1,"lines"),#缩小图例符号的大小
#         legend.background=element_blank(),#设置图例背景为空白
#         legend.position="inside",#修改图例位置
#         # legend.justification=c("right","top"),#修改图例位置
#         legend.position.inside=c(0.9,0.9),#修改图例位置
#         legend.margin=margin(0,0,0,0))+#修改图例位置
#   # annotate("text",x=1.8,y=max(WCON_long$Value),hjust=1,vjust=1,label=paste0("p < ",0.01),size=2,color="black")+
#   annotate("text",x=1.5,y=(max(WCON_long$Value)+(max(WCON_long$Value)-min(WCON_long$Value))*0.2),hjust=1,vjust=1,label=paste0("**"),size=2,color="black")+
#   annotate("text",x=1.5+0.2,y=(max(WCON_long$Value)+(max(WCON_long$Value)-min(WCON_long$Value))*0.1),hjust=1,vjust=1,label=paste0("**"),size=2,color="black")
# p1

p1 <- ggplot(WCON,aes(x=Plot,y=Observe,color=Index_type))+
  # geom_boxplot(width=0.225,linewidth=0.25,outliers=F,alpha=1,show.legend=T,position=position_dodge(0.5))+
  # stat_boxplot(geom="errorbar",width=0.225,linewidth=0.25,position=position_dodge(0.5))+
  # geom_jitter(shape=16,size=0.6,position=position_jitterdodge(jitter.height=0,jitter.width=0.8,dodge.width=0.5),alpha=0.4,show.legend=T)+
  geom_boxplot(width=0.225,linewidth=0.25,outliers=F,alpha=1,show.legend=T,position=position_dodge(0.5))+
  stat_boxplot(geom="errorbar",width=0.225,linewidth=0.25,position=position_dodge(0.5))+
  geom_jitter(shape=16,size=0.6,position=position_jitterdodge(jitter.height=0,jitter.width=0.8,dodge.width=0.5),alpha=0.6,show.legend=T)+
  scale_color_manual(values=c("gray40","gray40","black"))+
  scale_fill_manual(values=c("gray40","gray40","black"))+
  theme_classic()+
  # scale_color_manual(values=c("gray40","gray40","black"))+
  # scale_fill_manual(values=c("gray40","gray40","black"))+
  ylab("Weighted connectance")+# 设置Y轴标题
  # scale_x_discrete(breaks=c("Eu-adjacent","Eu-distant"),labels=c("Eu-adjacent", "Eu-distant"))+
  theme(axis.title.x=element_blank(),#修改X轴标题文本
        # axis.title.x=element_text(size=8),#修改X轴标题文本
        axis.title.y=element_text(size=8,angle=90),#修改y轴标题文本
        axis.text.y=element_text(size=6.5),#修改y轴刻度标签文本
        # axis.line.x = element_line(arrow = arrow(length = unit(0.2,"lines"))),
        axis.text.x=element_blank(),#修改x轴刻度标签文本
        # axis.text.x=element_text(size=6.5),#修改x轴刻度标签文本
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
  # annotate("text",x=1.8,y=max(WCON$Observe),hjust=1,vjust=1,label=paste0("p < ",0.01),size=2,color="black")+
  annotate("text",x=1.5,y=(max(WCON$Observe)+(max(WCON$Observe)-min(WCON$Observe))*0.1),hjust=1,vjust=1,label=paste0("**"),size=2,color="black")
p1

ggsave(filename="Fig6a_network_index_Weighted_connectance.png",plot=p1,width=89/2,height=59/2,units="mm",dpi=900)
ggsave(filename="Fig6a_network_index_Weighted_connectance.pdf",plot=p1,width=89/2,height=59/2,units="mm")

#####
Evenness <- subset(network_index,Index_type=="Interaction evenness")
x <- subset(Evenness,Plot=="Eu-adjacent")$Observe
y <- subset(Evenness,Plot=="Eu-distant")$Observe
shapiro.test(x)
shapiro.test(y)
levene_test <- leveneTest(Observe ~ Plot, data = Evenness)
levene_test
t.test(x,y,paired=T)

# x <- subset(Evenness,Plot=="Eu-adjacent")$Observe_Euphorbia_excluded
# y <- subset(Evenness,Plot=="Eu-distant")$Observe_Euphorbia_excluded
# t.test(x,y,paired=T)
# 
# Evenness_long <- pivot_longer(Evenness, cols = c(Observe, Observe_Euphorbia_excluded), names_to = "Type", values_to = "Value")
# Evenness_long[which(Evenness_long$Type=="Observe"),]$Type <- "included"
# Evenness_long[which(Evenness_long$Type=="Observe_Euphorbia_excluded"),]$Type <- "excluded"
# Evenness_long$Type <- factor(Evenness_long$Type,levels=c("included","excluded"))
# 
# # 计算均值和标准差
# Evenness_long.summary <- Evenness_long %>%
#   group_by(Plot, Type) %>%
#   summarise(Value.mean = mean(Value),Value.sd = sd(Value))

# p2 <- ggplot(Evenness_long,aes(x=Plot,y=Value,color=Type))+
#   geom_boxplot(width=0.4,linewidth=0.25,outliers=F,alpha=1,show.legend=T,position=position_dodge(0.5))+
#   stat_boxplot(geom="errorbar",width=0.4,linewidth=0.25,position=position_dodge(0.5))+
#   geom_jitter(shape=16,size=0.6,position=position_jitterdodge(jitter.height=0,jitter.width=0.45,dodge.width=0.5),alpha=0.6,show.legend=T)+
#   theme_classic()+
#   scale_color_manual(values=c("#bc40bf","#4BA465","black"))+
#   scale_fill_manual(values=c("#bc40bf","#4BA465","black"))+
#   ylab("Interaction evenness")+# 设置Y轴标题
#   # scale_x_discrete(breaks=c("Eu-adjacent","Eu-distant"),labels=c("Eu-adjacent", "Eu-distant"))+
#   theme(axis.title.x=element_blank(),#修改X轴标题文本
#         # axis.title.x=element_text(size=8),#修改X轴标题文本
#         axis.title.y=element_text(size=8,angle=90),#修改y轴标题文本
#         axis.text.y=element_text(size=6.5),#修改y轴刻度标签文本
#         # axis.line.x = element_line(arrow = arrow(length = unit(0.2,"lines"))),
#         axis.text.x=element_blank(),#修改x轴刻度标签文本
#         # axis.text.x=element_text(size=6.5),#修改x轴刻度标签文本
#         panel.grid=element_blank(),#隐藏网格线
#         # legend.title=element_text(size=6.5),#修改图例标题大小
#         legend.title=element_blank(),#修改图例标题大小
#         legend.text=element_text(size=5.5),#修改图例文本大小
#         legend.key.size=unit(0.1,"lines"),#缩小图例符号的大小
#         legend.background=element_blank(),#设置图例背景为空白
#         legend.position="inside",#修改图例位置
#         # legend.justification=c("right","top"),#修改图例位置
#         legend.position.inside=c(0.9,0.9),#修改图例位置
#         legend.margin=margin(0,0,0,0))+#修改图例位置
#   # annotate("text",x=1.8,y=max(Evenness_long$Value),hjust=1,vjust=1,label=paste0("p > ",0.05),size=2,color="black")+
#   annotate("text",x=1.5,y=(max(Evenness_long$Value)+(max(Evenness_long$Value)-min(Evenness_long$Value))*0.2),hjust=1,vjust=1,label=paste0("NS"),size=2,color="black")+
#   annotate("text",x=1.5+0.2,y=(max(Evenness_long$Value)+(max(Evenness_long$Value)-min(Evenness_long$Value))*0.1),hjust=1,vjust=1,label=paste0("NS"),size=2,color="black")
# p2

p2 <- ggplot(Evenness,aes(x=Plot,y=Observe,color=Index_type))+
  # geom_boxplot(width=0.225,linewidth=0.25,outliers=F,alpha=1,show.legend=T,position=position_dodge(0.5))+
  # stat_boxplot(geom="errorbar",width=0.225,linewidth=0.25,position=position_dodge(0.5))+
  # geom_jitter(shape=16,size=0.6,position=position_jitterdodge(jitter.height=0,jitter.width=0.8,dodge.width=0.5),alpha=0.4,show.legend=T)+
  geom_boxplot(width=0.225,linewidth=0.25,outliers=F,alpha=1,show.legend=T,position=position_dodge(0.5))+
  stat_boxplot(geom="errorbar",width=0.225,linewidth=0.25,position=position_dodge(0.5))+
  geom_jitter(shape=16,size=0.6,position=position_jitterdodge(jitter.height=0,jitter.width=0.8,dodge.width=0.5),alpha=0.6,show.legend=T)+
  scale_color_manual(values=c("gray40","gray40","black"))+
  scale_fill_manual(values=c("gray40","gray40","black"))+
  theme_classic()+
  # scale_color_manual(values=c("gray40","gray40","black"))+
  # scale_fill_manual(values=c("gray40","gray40","black"))+
  ylab("Interaction evenness")+# 设置Y轴标题
  # scale_x_discrete(breaks=c("Eu-adjacent","Eu-distant"),labels=c("Eu-adjacent", "Eu-distant"))+
  theme(axis.title.x=element_blank(),#修改X轴标题文本
        # axis.title.x=element_text(size=8),#修改X轴标题文本
        axis.title.y=element_text(size=8,angle=90),#修改y轴标题文本
        axis.text.y=element_text(size=6.5),#修改y轴刻度标签文本
        # axis.line.x = element_line(arrow = arrow(length = unit(0.2,"lines"))),
        axis.text.x=element_blank(),#修改x轴刻度标签文本
        # axis.text.x=element_text(size=6.5),#修改x轴刻度标签文本
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
  # annotate("text",x=1.8,y=max(Evenness$Observe),hjust=1,vjust=1,label=paste0("p > ",0.05),size=2,color="black")+
  annotate("text",x=1.5,y=(max(Evenness$Observe)+(max(Evenness$Observe)-min(Evenness$Observe))*0.2),hjust=1,vjust=1,label=paste0("NS"),size=2,color="black")
p2

ggsave(filename="Fig6b_network_index_evenness.png",plot=p2,width=89/2,height=59/2,units="mm",dpi=900)
ggsave(filename="Fig6b_network_index_evenness.pdf",plot=p2,width=89/2,height=59/2,units="mm")

#####
WNODF <- subset(network_index,Index_type=="Weighted NODF")
x <- subset(WNODF,Plot=="Eu-adjacent")$Observe
y <- subset(WNODF,Plot=="Eu-distant")$Observe
shapiro.test(x)
shapiro.test(y)
levene_test <- leveneTest(Observe ~ Plot, data = WNODF)
levene_test
t.test(x,y,paired=T)

# x <- subset(WNODF,Plot=="Eu-adjacent")$Observe_Euphorbia_excluded
# y <- subset(WNODF,Plot=="Eu-distant")$Observe_Euphorbia_excluded
# t.test(x,y,paired=T)
# 
# WNODF_long <- pivot_longer(WNODF, cols = c(Observe, Observe_Euphorbia_excluded), names_to = "Type", values_to = "Value")
# WNODF_long[which(WNODF_long$Type=="Observe"),]$Type <- "included"
# WNODF_long[which(WNODF_long$Type=="Observe_Euphorbia_excluded"),]$Type <- "excluded"
# WNODF_long$Type <- factor(WNODF_long$Type,levels=c("included","excluded"))
# 
# # 计算均值和标准差
# WNODF_long.summary <- WNODF_long %>%
#   group_by(Plot, Type) %>%
#   summarise(Value.mean = mean(Value),Value.sd = sd(Value))

# p3 <- ggplot(WNODF_long,aes(x=Plot,y=Value,color=Type))+
#   geom_boxplot(width=0.4,linewidth=0.25,outliers=F,alpha=1,show.legend=T,position=position_dodge(0.5))+
#   stat_boxplot(geom="errorbar",width=0.4,linewidth=0.25,position=position_dodge(0.5))+
#   geom_jitter(shape=16,size=0.6,position=position_jitterdodge(jitter.height=0,jitter.width=0.45,dodge.width=0.5),alpha=0.6,show.legend=T)+
#   theme_classic()+
#   scale_color_manual(values=c("#bc40bf","#4BA465","black"))+
#   scale_fill_manual(values=c("#bc40bf","#4BA465","black"))+
#   ylab("Weighted NODF")+# 设置Y轴标题
#   # scale_x_discrete(breaks=c("Eu-adjacent","Eu-distant"),labels=c("Eu-adjacent", "Eu-distant"))+
#   scale_y_continuous(labels = function(x) sprintf("%03d", x))+
#   theme(axis.title.x=element_blank(),#修改X轴标题文本
#         # axis.title.x=element_text(size=8),#修改X轴标题文本
#         axis.title.y=element_text(size=8,angle=90),#修改y轴标题文本
#         axis.text.y=element_text(size=6.5),#修改y轴刻度标签文本
#         # axis.line.x = element_line(arrow = arrow(length = unit(0.2,"lines"))),
#         axis.text.x=element_text(size=6.5),#修改x轴刻度标签文本
#         panel.grid=element_blank(),#隐藏网格线
#         # legend.title=element_text(size=6.5),#修改图例标题大小
#         legend.title=element_blank(),#修改图例标题大小
#         legend.text=element_text(size=5.5),#修改图例文本大小
#         legend.key.size=unit(0.1,"lines"),#缩小图例符号的大小
#         legend.background=element_blank(),#设置图例背景为空白
#         legend.position="inside",#修改图例位置
#         # legend.justification=c("right","top"),#修改图例位置
#         legend.position.inside=c(0.9,0.9),#修改图例位置
#         legend.margin=margin(0,0,0,0))+#修改图例位置
#   # annotate("text",x=1.8,y=max(WNODF_long$Value),hjust=1,vjust=1,label=paste0("p < ",0.01),size=2,color="black")+
#   annotate("text",x=1.5,y=(max(WNODF_long$Value)+(max(WNODF_long$Value)-min(WNODF_long$Value))*0.2),hjust=1,vjust=1,label=paste0("**"),size=2,color="black")+
#   annotate("text",x=1.5+0.2,y=(max(WNODF_long$Value)+(max(WNODF_long$Value)-min(WNODF_long$Value))*0.1),hjust=1,vjust=1,label=paste0("**"),size=2,color="black")
# p3

p3 <- ggplot(WNODF,aes(x=Plot,y=Observe,color=Index_type))+
  # geom_boxplot(width=0.225,linewidth=0.25,outliers=F,alpha=1,show.legend=T,position=position_dodge(0.5))+
  # stat_boxplot(geom="errorbar",width=0.225,linewidth=0.25,position=position_dodge(0.5))+
  # geom_jitter(shape=16,size=0.6,position=position_jitterdodge(jitter.height=0,jitter.width=0.8,dodge.width=0.5),alpha=0.4,show.legend=T)+
  geom_boxplot(width=0.225,linewidth=0.25,outliers=F,alpha=1,show.legend=T,position=position_dodge(0.5))+
  stat_boxplot(geom="errorbar",width=0.225,linewidth=0.25,position=position_dodge(0.5))+
  geom_jitter(shape=16,size=0.6,position=position_jitterdodge(jitter.height=0,jitter.width=0.8,dodge.width=0.5),alpha=0.6,show.legend=T)+
  scale_color_manual(values=c("gray40","gray40","black"))+
  scale_fill_manual(values=c("gray40","gray40","black"))+
  theme_classic()+
  # scale_color_manual(values=c("gray40","gray40","black"))+
  # scale_fill_manual(values=c("gray40","gray40","black"))+
  ylab("Weighted NODF")+# 设置Y轴标题
  # scale_x_discrete(breaks=c("Eu-adjacent","Eu-distant"),labels=c("Eu-adjacent", "Eu-distant"))+
  scale_y_continuous(labels = function(x) sprintf("%03d", x))+
  theme(axis.title.x=element_blank(),#修改X轴标题文本
        # axis.title.x=element_text(size=8),#修改X轴标题文本
        axis.title.y=element_text(size=8,angle=90),#修改y轴标题文本
        axis.text.y=element_text(size=6.5),#修改y轴刻度标签文本
        # axis.line.x = element_line(arrow = arrow(length = unit(0.2,"lines"))),
        axis.text.x=element_text(size=6.5),#修改x轴刻度标签文本
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
  # annotate("text",x=1.8,y=max(WNODF$Observe),hjust=1,vjust=1,label=paste0("p < ",0.01),size=2,color="black")+
  annotate("text",x=1.5,y=(max(WNODF$Observe)+(max(WNODF$Observe)-min(WNODF$Observe))*0.2),hjust=1,vjust=1,label=paste0("**"),size=2,color="black")
p3

ggsave(filename="Fig6c_network_index_Weighted_NODF.png",plot=p3,width=89/2,height=59/2,units="mm",dpi=900)
ggsave(filename="Fig6c_network_index_Weighted_NODF.pdf",plot=p3,width=89/2,height=59/2,units="mm")

#####
WMOD <- subset(network_index,Index_type=="Weighted modularity")
x <- subset(WMOD,Plot=="Eu-adjacent")$Observe
y <- subset(WMOD,Plot=="Eu-distant")$Observe
shapiro.test(x)
shapiro.test(y)
levene_test <- leveneTest(Observe ~ Plot, data = WMOD)
levene_test
t.test(x,y,paired=T)

# x <- subset(WMOD,Plot=="Eu-adjacent")$Observe_Euphorbia_excluded
# y <- subset(WMOD,Plot=="Eu-distant")$Observe_Euphorbia_excluded
# t.test(x,y,paired=T)
# 
# WMOD_long <- pivot_longer(WMOD, cols = c(Observe, Observe_Euphorbia_excluded), names_to = "Type", values_to = "Value")
# WMOD_long[which(WMOD_long$Type=="Observe"),]$Type <- "included"
# WMOD_long[which(WMOD_long$Type=="Observe_Euphorbia_excluded"),]$Type <- "excluded"
# WMOD_long$Type <- factor(WMOD_long$Type,levels=c("included","excluded"))
# 
# # 计算均值和标准差
# WMOD_long.summary <- WMOD_long %>%
#   group_by(Plot, Type) %>%
#   summarise(Value.mean = mean(Value),Value.sd = sd(Value))

# p4 <- ggplot(WMOD_long,aes(x=Plot,y=Value,color=Type))+
#   geom_boxplot(width=0.4,linewidth=0.25,outliers=F,alpha=1,show.legend=T,position=position_dodge(0.5))+
#   stat_boxplot(geom="errorbar",width=0.4,linewidth=0.25,position=position_dodge(0.5))+
#   geom_jitter(shape=16,size=0.6,position=position_jitterdodge(jitter.height=0,jitter.width=0.45,dodge.width=0.5),alpha=0.6,show.legend=T)+
#   theme_classic()+
#   scale_color_manual(values=c("#bc40bf","#4BA465","black"))+
#   scale_fill_manual(values=c("#bc40bf","#4BA465","black"))+
#   ylab("Weighted modularity")+# 设置Y轴标题
#   # scale_x_discrete(breaks=c("Eu-adjacent","Eu-distant"),labels=c("Eu-adjacent", "Eu-distant"))+
#   theme(axis.title.x=element_blank(),#修改X轴标题文本
#         # axis.title.x=element_text(size=8),#修改X轴标题文本
#         axis.title.y=element_text(size=8,angle=90),#修改y轴标题文本
#         axis.text.y=element_text(size=6.5),#修改y轴刻度标签文本
#         # axis.line.x = element_line(arrow = arrow(length = unit(0.2,"lines"))),
#         axis.text.x=element_text(size=6.5),#修改x轴刻度标签文本
#         panel.grid=element_blank(),#隐藏网格线
#         # legend.title=element_text(size=6.5),#修改图例标题大小
#         legend.title=element_blank(),#修改图例标题大小
#         legend.text=element_text(size=5.5),#修改图例文本大小
#         legend.key.size=unit(0.1,"lines"),#缩小图例符号的大小
#         legend.background=element_blank(),#设置图例背景为空白
#         legend.position="inside",#修改图例位置
#         # legend.justification=c("right","top"),#修改图例位置
#         legend.position.inside=c(0.9,0.9),#修改图例位置
#         legend.margin=margin(0,0,0,0))+#修改图例位置
#   # annotate("text",x=1.8,y=max(WMOD_long$Value),hjust=1,vjust=1,label=paste0("p < ",0.05),size=2,color="black")+
#   annotate("text",x=1.5,y=(max(WMOD_long$Value)+(max(WMOD_long$Value)-min(WMOD_long$Value))*0.2),hjust=1,vjust=1,label=paste0("*"),size=2,color="black")+
#   annotate("text",x=1.5+0.2,y=(max(WMOD_long$Value)+(max(WMOD_long$Value)-min(WMOD_long$Value))*0.1),hjust=1,vjust=1,label=paste0("*"),size=2,color="black")
# p4

p4 <- ggplot(WMOD,aes(x=Plot,y=Observe,color=Index_type))+
  # geom_boxplot(width=0.225,linewidth=0.25,outliers=F,alpha=1,show.legend=T,position=position_dodge(0.5))+
  # stat_boxplot(geom="errorbar",width=0.225,linewidth=0.25,position=position_dodge(0.5))+
  # geom_jitter(shape=16,size=0.6,position=position_jitterdodge(jitter.height=0,jitter.width=0.8,dodge.width=0.5),alpha=0.4,show.legend=T)+
  geom_boxplot(width=0.225,linewidth=0.25,outliers=F,alpha=1,show.legend=T,position=position_dodge(0.5))+
  stat_boxplot(geom="errorbar",width=0.225,linewidth=0.25,position=position_dodge(0.5))+
  geom_jitter(shape=16,size=0.6,position=position_jitterdodge(jitter.height=0,jitter.width=0.8,dodge.width=0.5),alpha=0.6,show.legend=T)+
  scale_color_manual(values=c("gray40","gray40","black"))+
  scale_fill_manual(values=c("gray40","gray40","black"))+
  theme_classic()+
  # scale_color_manual(values=c("gray40","gray40","black"))+
  # scale_fill_manual(values=c("gray40","gray40","black"))+
  ylab("Weighted modularity")+# 设置Y轴标题
  # scale_x_discrete(breaks=c("Eu-adjacent","Eu-distant"),labels=c("Eu-adjacent", "Eu-distant"))+
  theme(axis.title.x=element_blank(),#修改X轴标题文本
        # axis.title.x=element_text(size=8),#修改X轴标题文本
        axis.title.y=element_text(size=8,angle=90),#修改y轴标题文本
        axis.text.y=element_text(size=6.5),#修改y轴刻度标签文本
        # axis.line.x = element_line(arrow = arrow(length = unit(0.2,"lines"))),
        axis.text.x=element_text(size=6.5),#修改x轴刻度标签文本
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
  # annotate("text",x=1.8,y=max(WMOD$Observe),hjust=1,vjust=1,label=paste0("p < ",0.05),size=2,color="black")+
  annotate("text",x=1.5,y=(max(WMOD$Observe)+(max(WMOD$Observe)-min(WMOD$Observe))*0.2),hjust=1,vjust=1,label=paste0("*"),size=2,color="black")
p4

ggsave(filename="Fig6d_network_index_Weighted_modularity.png",plot=p4,width=89/2,height=59/2,units="mm",dpi=900)
ggsave(filename="Fig6d_network_index_Weighted_modularity.pdf",plot=p4,width=89/2,height=59/2,units="mm")

p_all <- plot_grid(p1,p2,p3,p4,labels=c("a","b","c","d"),label_size=8,label_x=0.28,label_y=1,ncol=2,nrow=2)
p_all
ggsave(filename="Fig6_0209.png",plot=p_all,width=89,height=59,units="mm",dpi=900)
ggsave(filename="Fig6_0209.pdf",plot=p_all,width=89,height=59,units="mm")

###############
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

x <- subset(WCON,Delta_type=="delta_Observed")$Delta_value_Euphorbia_excluded
y <- subset(WCON,Delta_type=="delta_Simulation")$Delta_value_Euphorbia_excluded
t.test(x,y,paired=T)

WCON_long <- pivot_longer(WCON, cols = c(Delta_value, Delta_value_Euphorbia_excluded), names_to = "Type", values_to = "Value")
WCON_long[which(WCON_long$Type=="Delta_value"),]$Type <- "included"
WCON_long[which(WCON_long$Type=="Delta_value_Euphorbia_excluded"),]$Type <- "excluded"
WCON_long$Type <- factor(WCON_long$Type,levels=c("included","excluded"))

# 计算均值和标准差
WCON_long.summary <- WCON_long %>%
  group_by(Delta_type, Type) %>%
  summarise(Value.mean = mean(Value),Value.sd = sd(Value))

p1 <- ggplot(WCON_long,aes(x=Delta_type,y=Value,color=Type))+
  geom_boxplot(width=0.4,linewidth=0.25,outliers=F,alpha=1,show.legend=T,position=position_dodge(0.5))+
  stat_boxplot(geom="errorbar",width=0.4,linewidth=0.25,position=position_dodge(0.5))+
  geom_jitter(shape=16,size=0.6,position=position_jitterdodge(jitter.height=0,jitter.width=0.45,dodge.width=0.5),alpha=0.6,show.legend=T)+
  theme_classic()+
  scale_color_manual(values=c("#bc40bf","#4BA465","black"))+
  scale_fill_manual(values=c("#bc40bf","#4BA465","black"))+
  ylab("ΔWeighted connectance")+# 设置Y轴标题
  scale_x_discrete(breaks=c("delta_Observed","delta_Simulation"),labels=c("Observed", "Simulation"))+
  theme(axis.title.x=element_blank(),#修改X轴标题文本
        # axis.title.x=element_text(size=8),#修改X轴标题文本
        axis.title.y=element_text(size=8,angle=90),#修改y轴标题文本
        axis.text.y=element_text(size=6.5),#修改y轴刻度标签文本
        # axis.line.x = element_line(arrow = arrow(length = unit(0.2,"lines"))),
        axis.text.x=element_blank(),#修改x轴刻度标签文本
        # axis.text.x=element_text(size=6.5),#修改x轴刻度标签文本
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
  # annotate("text",x=1.8,y=max(WCON_long$Value),hjust=1,vjust=1,label=paste0("p < ",0.01),size=2,color="black")+
  annotate("text",x=1.5,y=(max(WCON_long$Value)+(max(WCON_long$Value)-min(WCON_long$Value))*0.2),hjust=1,vjust=1,label=paste0("**"),size=2,color="black")+
  annotate("text",x=1.5+0.2,y=(max(WCON_long$Value)+(max(WCON_long$Value)-min(WCON_long$Value))*0.1),hjust=1,vjust=1,label=paste0("**"),size=2,color="black")
p1

ggsave(filename="FigS1a_network_index_Weighted_connectance_simulation.png",plot=p1,width=89/2,height=59/2,units="mm",dpi=900)
ggsave(filename="FigS1a_network_index_Weighted_connectance_simulation.pdf",plot=p1,width=89/2,height=59/2,units="mm")

#####
Evenness <- subset(network_index_simulation,Index_type=="Interaction evenness")
x <- subset(Evenness,Delta_type=="delta_Observed")$Delta_value
y <- subset(Evenness,Delta_type=="delta_Simulation")$Delta_value
t.test(x,y,paired=T)

x <- subset(Evenness,Delta_type=="delta_Observed")$Delta_value_Euphorbia_excluded
y <- subset(Evenness,Delta_type=="delta_Simulation")$Delta_value_Euphorbia_excluded
t.test(x,y,paired=T)

Evenness_long <- pivot_longer(Evenness, cols = c(Delta_value, Delta_value_Euphorbia_excluded), names_to = "Type", values_to = "Value")
Evenness_long[which(Evenness_long$Type=="Delta_value"),]$Type <- "included"
Evenness_long[which(Evenness_long$Type=="Delta_value_Euphorbia_excluded"),]$Type <- "excluded"
Evenness_long$Type <- factor(Evenness_long$Type,levels=c("included","excluded"))

# 计算均值和标准差
Evenness_long.summary <- Evenness_long %>%
  group_by(Delta_type, Type) %>%
  summarise(Value.mean = mean(Value),Value.sd = sd(Value))

p2 <- ggplot(Evenness_long,aes(x=Delta_type,y=Value,color=Type))+
  geom_boxplot(width=0.4,linewidth=0.25,outliers=F,alpha=1,show.legend=T,position=position_dodge(0.5))+
  stat_boxplot(geom="errorbar",width=0.4,linewidth=0.25,position=position_dodge(0.5))+
  geom_jitter(shape=16,size=0.6,position=position_jitterdodge(jitter.height=0,jitter.width=0.45,dodge.width=0.5),alpha=0.6,show.legend=T)+
  theme_classic()+
  scale_color_manual(values=c("#bc40bf","#4BA465","black"))+
  scale_fill_manual(values=c("#bc40bf","#4BA465","black"))+
  ylab("ΔInteraction evenness")+# 设置Y轴标题
  scale_x_discrete(breaks=c("delta_Observed","delta_Simulation"),labels=c("Observed", "Simulation"))+
  theme(axis.title.x=element_blank(),#修改X轴标题文本
        # axis.title.x=element_text(size=8),#修改X轴标题文本
        axis.title.y=element_text(size=8,angle=90),#修改y轴标题文本
        axis.text.y=element_text(size=6.5),#修改y轴刻度标签文本
        # axis.line.x = element_line(arrow = arrow(length = unit(0.2,"lines"))),
        axis.text.x=element_blank(),#修改x轴刻度标签文本
        # axis.text.x=element_text(size=6.5),#修改x轴刻度标签文本
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
  # annotate("text",x=1.8,y=max(Evenness_long$Value),hjust=1,vjust=1,label=paste0("p > ",0.05),size=2,color="black")+
  annotate("text",x=1.5,y=(max(Evenness_long$Value)+(max(Evenness_long$Value)-min(Evenness_long$Value))*0.2),hjust=1,vjust=1,label=paste0("NS"),size=2,color="black")+
  annotate("text",x=1.5+0.2,y=(max(Evenness_long$Value)+(max(Evenness_long$Value)-min(Evenness_long$Value))*0.1),hjust=1,vjust=1,label=paste0("NS"),size=2,color="black")
p2

ggsave(filename="FigS1b_network_index_evenness_simulation.png",plot=p2,width=89/2,height=59/2,units="mm",dpi=900)
ggsave(filename="FigS1b_network_index_evenness_simulation.pdf",plot=p2,width=89/2,height=59/2,units="mm")

#####
WNODF <- subset(network_index_simulation,Index_type=="Weighted NODF")
x <- subset(WNODF,Delta_type=="delta_Observed")$Delta_value
y <- subset(WNODF,Delta_type=="delta_Simulation")$Delta_value
t.test(x,y,paired=T)

x <- subset(WNODF,Delta_type=="delta_Observed")$Delta_value_Euphorbia_excluded
y <- subset(WNODF,Delta_type=="delta_Simulation")$Delta_value_Euphorbia_excluded
t.test(x,y,paired=T)

WNODF_long <- pivot_longer(WNODF, cols = c(Delta_value, Delta_value_Euphorbia_excluded), names_to = "Type", values_to = "Value")
WNODF_long[which(WNODF_long$Type=="Delta_value"),]$Type <- "included"
WNODF_long[which(WNODF_long$Type=="Delta_value_Euphorbia_excluded"),]$Type <- "excluded"
WNODF_long$Type <- factor(WNODF_long$Type,levels=c("included","excluded"))

# 计算均值和标准差
WNODF_long.summary <- WNODF_long %>%
  group_by(Delta_type, Type) %>%
  summarise(Value.mean = mean(Value),Value.sd = sd(Value))

p3 <- ggplot(WNODF_long,aes(x=Delta_type,y=Value,color=Type))+
  geom_boxplot(width=0.4,linewidth=0.25,outliers=F,alpha=1,show.legend=T,position=position_dodge(0.5))+
  stat_boxplot(geom="errorbar",width=0.4,linewidth=0.25,position=position_dodge(0.5))+
  geom_jitter(shape=16,size=0.6,position=position_jitterdodge(jitter.height=0,jitter.width=0.45,dodge.width=0.5),alpha=0.6,show.legend=T)+
  theme_classic()+
  scale_color_manual(values=c("#bc40bf","#4BA465","black"))+
  scale_fill_manual(values=c("#bc40bf","#4BA465","black"))+
  ylab("ΔWeighted NODF")+# 设置Y轴标题
  scale_x_discrete(breaks=c("delta_Observed","delta_Simulation"),labels=c("Observed", "Simulation"))+
  # scale_y_continuous(labels = function(x) sprintf("%03d", x))+
  theme(axis.title.x=element_blank(),#修改X轴标题文本
        # axis.title.x=element_text(size=8),#修改X轴标题文本
        axis.title.y=element_text(size=8,angle=90),#修改y轴标题文本
        axis.text.y=element_text(size=6.5),#修改y轴刻度标签文本
        # axis.line.x = element_line(arrow = arrow(length = unit(0.2,"lines"))),
        axis.text.x=element_text(size=6.5),#修改x轴刻度标签文本
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
  # annotate("text",x=1.8,y=max(WNODF_long$Value),hjust=1,vjust=1,label=paste0("p < ",0.01),size=2,color="black")+
  annotate("text",x=1.5,y=(max(WNODF_long$Value)+(max(WNODF_long$Value)-min(WNODF_long$Value))*0.2),hjust=1,vjust=1,label=paste0("**"),size=2,color="black")+
  annotate("text",x=1.5+0.2,y=(max(WNODF_long$Value)+(max(WNODF_long$Value)-min(WNODF_long$Value))*0.1),hjust=1,vjust=1,label=paste0("**"),size=2,color="black")
p3

ggsave(filename="FigS1c_network_index_Weighted_NODF_simulation.png",plot=p3,width=89/2,height=59/2,units="mm",dpi=900)
ggsave(filename="FigS1c_network_index_Weighted_NODF_simulation.pdf",plot=p3,width=89/2,height=59/2,units="mm")

#####
WMOD <- subset(network_index_simulation,Index_type=="Weighted modularity")
x <- subset(WMOD,Delta_type=="delta_Observed")$Delta_value
y <- subset(WMOD,Delta_type=="delta_Simulation")$Delta_value
t.test(x,y,paired=T)

x <- subset(WMOD,Delta_type=="delta_Observed")$Delta_value_Euphorbia_excluded
y <- subset(WMOD,Delta_type=="delta_Simulation")$Delta_value_Euphorbia_excluded
t.test(x,y,paired=T)

WMOD_long <- pivot_longer(WMOD, cols = c(Delta_value, Delta_value_Euphorbia_excluded), names_to = "Type", values_to = "Value")
WMOD_long[which(WMOD_long$Type=="Delta_value"),]$Type <- "included"
WMOD_long[which(WMOD_long$Type=="Delta_value_Euphorbia_excluded"),]$Type <- "excluded"
WMOD_long$Type <- factor(WMOD_long$Type,levels=c("included","excluded"))

# 计算均值和标准差
WMOD_long.summary <- WMOD_long %>%
  group_by(Delta_type, Type) %>%
  summarise(Value.mean = mean(Value),Value.sd = sd(Value))

p4 <- ggplot(WMOD_long,aes(x=Delta_type,y=Value,color=Type))+
  geom_boxplot(width=0.4,linewidth=0.25,outliers=F,alpha=1,show.legend=T,position=position_dodge(0.5))+
  stat_boxplot(geom="errorbar",width=0.4,linewidth=0.25,position=position_dodge(0.5))+
  geom_jitter(shape=16,size=0.6,position=position_jitterdodge(jitter.height=0,jitter.width=0.45,dodge.width=0.5),alpha=0.6,show.legend=T)+
  theme_classic()+
  scale_color_manual(values=c("#bc40bf","#4BA465","black"))+
  scale_fill_manual(values=c("#bc40bf","#4BA465","black"))+
  ylab("ΔWeighted modularity")+# 设置Y轴标题
  scale_x_discrete(breaks=c("delta_Observed","delta_Simulation"),labels=c("Observed", "Simulation"))+
  theme(axis.title.x=element_blank(),#修改X轴标题文本
        # axis.title.x=element_text(size=8),#修改X轴标题文本
        axis.title.y=element_text(size=8,angle=90),#修改y轴标题文本
        axis.text.y=element_text(size=6.5),#修改y轴刻度标签文本
        # axis.line.x = element_line(arrow = arrow(length = unit(0.2,"lines"))),
        axis.text.x=element_text(size=6.5),#修改x轴刻度标签文本
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
  # annotate("text",x=1.8,y=max(WMOD_long$Value),hjust=1,vjust=1,label=paste0("p < ",0.05),size=2,color="black")+
  annotate("text",x=1.5,y=(max(WMOD_long$Value)+(max(WMOD_long$Value)-min(WMOD_long$Value))*0.2),hjust=1,vjust=1,label=paste0("*"),size=2,color="black")+
  annotate("text",x=1.5+0.2,y=max(max(WMOD_long$Value)+(max(WMOD_long$Value)-min(WMOD_long$Value))*0.1),hjust=1,vjust=1,label=paste0("*"),size=2,color="black")
p4

ggsave(filename="FigS1d_network_index_Weighted_modularity_simulation.png",plot=p4,width=89/2,height=59/2,units="mm",dpi=900)
ggsave(filename="FigS1d_network_index_Weighted_modularity_simulation.pdf",plot=p4,width=89/2,height=59/2,units="mm")

p_all <- plot_grid(p1,p2,p3,p4,labels=c("a","b","c","d"),label_size=8,label_x=0.28,label_y=1,ncol=2,nrow=2)
p_all
ggsave(filename="FigS1_1210.png",plot=p_all,width=89,height=59,units="mm",dpi=900)
ggsave(filename="FigS1_1210.pdf",plot=p_all,width=89,height=59,units="mm")

#################################################################################

# Reference:
# White TJ, Bruns T, Lee S, Taylor J. Amplification and direct sequencing of fungal ribosomal RNA genes for phylogenetics. In: PCR Protocols: A guide to methods and applications (eds Innis MA, Gelfand DH, Sninsky JJ, White TJ). Academic Press, Inc (1990).
# Chen S, et al. Validation of the ITS2 Region as a Novel DNA Barcode for Identifying Medicinal Plant Species. Plos One 5, e8613 (2010).
# Ihrmark K, et al. New primers to amplify the fungal ITS2 region – evaluation by 454-sequencing of artificial and natural communities. Fems Microbiol Ecol 82, 666-677 (2012).
# Kohout P, et al. Comparison of commonly used primer sets for evaluating arbuscular mycorrhizal fungal communities: Is there a universal solution? Soil Biology and Biochemistry 68, 482-493 (2014).
# Lekberg Y, et al. More bang for the buck? Can arbuscular mycorrhizal fungal communities be characterized adequately alongside other fungi using general fungal primers? New Phytol 220, 971-976 (2018).
# Fadrosh DW, et al. An improved dual-indexing approach for multiplexed 16S rRNA gene sequencing on the Illumina MiSeq platform. Microbiome 2, 6 (2014).
# Kitson JJ, Hahn C, Sands RJ, Straw NA, Evans DM, Lunt DH. Detecting host–parasitoid interactions in an invasive Lepidopteran using nested tagging DNA metabarcoding. Mol Ecol 28, 471-483 (2019).
# Korbie DJ, Mattick JS. Touchdown PCR for increased specificity and sensitivity in PCR amplification. Nat Protoc 3, 1452-1456 (2008).
# Stevens JL, Jackson RL, Olson JB. Slowing PCR ramp speed reduces chimera formation from environmental samples. J Microbiol Meth 93, 203-205 (2013).
# Illumina. Effects of Index Misassignment on Multiplexing and Downstream Analysis [White paper]. https://doi.org/https://www.illumina.com/content/dam/illumina-marketing/documents/products/whitepapers/index-hopping-white-paper-770-2017-004.pdf (2017).
# Andrews S. FastQC: a quality control tool for high throughput sequence data. https://github.com/s-andrews/FastQC (2010).
# Aronesty E. Comparison of Sequencing Utility Programs. The Open Bioinformatics Journal 7, 1-8 (2013).
# Rognes T, Flouri T, Nichols B, Quince C, Mahe F. VSEARCH: a versatile open source tool for metagenomics. Peerj 4, e2584 (2016).
# Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics 34, i884-i890 (2018).
# Martin M. Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet journal 17, 10-12 (2011).
# Edgar RC, Flyvbjerg H. Error filtering, pair assembly and error correction for next-generation sequencing reads. Bioinformatics 31, 3476-3482 (2015).
# 1Rivers A, Weber K, Gardner T, Liu S, Armstrong S. ITSxpress: Software to rapidly trim internally transcribed spacer sequences with quality scores for marker gene analysis. https://f1000research.com/articles/7-1418/v1 (2018).
# McGinnis S, Madden TL. BLAST: at the core of a powerful and diverse set of sequence analysis tools. Nucleic Acids Res 32, W20-W25 (2004).
# Abarenkov K, et al. UNITE UCHIME reference dataset. Version 9.0.). UNITE Community (2022).
# Tedersoo L, et al. Best practices in metabarcoding of fungi: From experimental design to results. Mol Ecol 31, 2769-2795 (2022).
# Abarenkov K, et al. UNITE USEARCH/UTAX release for eukaryotes. Version 9.0.). UNITE Community (2023).
# Davis NM, Proctor DM, Holmes SP, Relman DA, Callahan BJ. Simple statistical identification and removal of contaminant sequences in marker-gene and metagenomics data. Microbiome 6, 226 (2018).
# Nguyen NH, Smith D, Peay K, Kennedy P. Parsing ecological signal from noise in next generation amplicon sequencing. New Phytol 205, 1389-1393 (2015).
# Cobian GM, Egan CP, Amend AS. Plant–microbe specificity varies as a function of elevation. The ISME Journal 13, 2778-2788 (2019).
# Montesinos-Navarro A, Segarra-Moragues JG, Valiente-Banuet A, Verdú M. The network structure of plant–arbuscular mycorrhizal fungi. New Phytol 194, 536-547 (2012).
# Toju H, Sato H, Yamamoto S, Tanabe AS. Structural diversity across arbuscular mycorrhizal, ectomycorrhizal, and endophytic plant–fungus networks. Bmc Plant Biol 18, 292 (2018).
# R Core Team. R: A language and environment for statistical computing. R Foundation for Statistical Computing. https://www.R-project.org/ (2013).
# Csardi G, Nepusz T. The igraph software package for complex network research. InterJournal Complex Systems, 1965 (2006).
# Tyner S, C., Briatte F, Hofmann H. Network Visualization with ggplot2. The R Journal,  (2017).
# Wickham H, Chang W, Wickham MH. Package 'ggplot2'. Create elegant data visualisations using the grammar of graphics Version 2, 1-189 (2016).
# Hsieh T, Ma K, Chao A. iNEXT: an R package for rarefaction and extrapolation of species diversity (H ill numbers). Methods Ecol Evol 7, 1451-1456 (2016).
# Magnusson A, et al. Package 'glmmTMB'. https://cran.r-project.org/package=glmmTMB (2017).
# Bartoń K. Package 'MuMIn'. https://CRAN.R-project.org/package=MuMIn (2015).
# Nelder JA, Wedderburn RWM. Generalized Linear Models. Royal Statistical Society Journal Series A: General 135, 370-384 (1972).
# Bolker BM, et al. Generalized linear mixed models: a practical guide for ecology and evolution. Trends Ecol Evol 24, 127-135 (2009).
# Hartig F. DHARMa: residual diagnostics for hierarchical (multi-level/mixed) regression models. https://CRAN.R-project.org/package=DHARMa (2018).
# Baselga A, Orme CDL. betapart: an R package for the study of beta diversity. Methods Ecol Evol 3, 808-812 (2012).
# Baselga A. Partitioning the turnover and nestedness components of beta diversity. Global Ecology & Biogeography 19, 134-143 (2010).
# Poisot T, Canard E, Mouillot D, Mouquet N, Gravel D. The dissimilarity of species interaction networks. Ecol Lett 15, 1353-1361 (2012).
# Dormann CF, Gruber B, Fründ J. Introducing the bipartite package: analysing ecological networks. interaction 1,  (2008).
# Bersier L-F, Banašek-Richter C, Cattin M-F. QUANTITATIVE DESCRIPTORS OF FOOD-WEB MATRICES. Ecology 83, 2394-2407 (2002).
# Tylianakis JM, Tscharntke T, Lewis OT. Habitat modification alters the structure of tropical host–parasitoid food webs. Nature 445, 202-205 (2007).
# Almeida-Neto M, Ulrich W. A straightforward computational approach for measuring nestedness using quantitative matrices. Environ Modell Softw 26, 173-178 (2011).
# Dormann CF, Strauss R. A method for detecting modules in quantitative bipartite networks. Methods Ecol Evol 5, 90-98 (2014).
# Beckett SJ. Improved community detection in weighted bipartite networks. Roy Soc Open Sci 3, 140536 (2016).
# Bastolla U, Fortuna MA, Pascual-Garcia A, Ferrera A, Luque B, Bascompte J. The architecture of mutualistic networks minimizes competition and increases biodiversity. Nature 458, 1018-U1091 (2009).
# Patefield WM. Algorithm AS 159: an efficient method of generating random R× C tables with given row and column totals. Journal of the Royal Statistical Society Series C (Applied Statistics) 30, 91-97 (1981).

#################################################################################
#################################################################################
#################################################################################