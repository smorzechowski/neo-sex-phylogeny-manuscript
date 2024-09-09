# Landscape graphs for repeat divergences
# Sophie MacRae Orzechowski
# Updated March 2024


setwd("~/PhD research/Neo sex chromosome/repeatmasker")

library(tidyr)
library(dplyr)
library(ggplot2)
library(scales)
library(patchwork)
show_col(hue_pal()(20))

#data<- read.table("cyan_flye_NP_assembly_incl_addedZandW_v2_rm1kb.fasta.align.landscape.Div.Rclass.tab",header=F)
data<- read.table("Ecyan_HiC_v3.0_neoW.fa.align.landscape.Div.Rclass.tab",header=F)
data <- data[-1,]
Wlength<-83413875


data_long <- data %>% pivot_longer(-V1, names_to = "Category", values_to = "count")
data_long$Divergence <- as.numeric(gsub("V","",data_long$Category))-1
data_long$count <- as.numeric(as.character(data_long$count))
data_long$percent <- (data_long$count/Wlength)
data_long$V1 <- as.character(data_long$V1)
colnames(data_long)[1] <-'Class'



#Categories<-c("V2","V3","V4","V5","V6","V7","V8","V9","V10","V11","V12","V13","V14","V15","V16","V17","V18","V19","V20","V21","V22","V23","V24",
#"V25","V26","V27","V28","V29","V30","V31","V32","V33","V34","V35","V36","V37","V38","V39","V40","V41","V42","V43","V44","V45","V46","V47","V48","V49","V50")
# scale_x_discrete(limits = Categories,breaks=c('V2','V10','V20','V30','V40','V50'),
#                  labels=c('0','10','20','30','40','50'))+


neoW <- ggplot(data_long,aes(Divergence,percent,fill=Class))+
  geom_bar(stat="identity")+
  scale_y_continuous(labels=scales::percent)+
  #xlab('% Divergence')+
  #ylab('% of chromosome')+
  xlab('')+
  ylab('')+
  theme_bw()+
  theme(axis.title = element_text(size=13),
        axis.text = element_text(size=11),
        legend.position="none")

############### Neo-Z


data<- read.table("Ecyan_HiC_v3.0_neoZ.fa.align.landscape.Div.Rclass.tab",header=F)
data <- data[-1,]
Zlength<-129165907

data_long <- data %>% pivot_longer(-V1, names_to = "Category", values_to = "count")
data_long$Divergence <- as.numeric(gsub("V","",data_long$Category))-1
data_long$count <- as.numeric(as.character(data_long$count))
data_long$percent <- (data_long$count/Zlength)
data_long$V1 <- as.character(data_long$V1)
colnames(data_long)[1] <-'Class'

ggplot(data_long,aes(Divergence,percent,fill=Class))+
  geom_bar(stat="identity")+
  scale_y_continuous(labels=scales::percent,limits=c(0,0.025))+
  xlab('% Divergence')+
  ylab('% of chromosome')+
  theme_bw()+
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=12))

################### Autosomes


data<- read.table("Ecyan_HiC_v3.0_sans_neoZW.fa.align.landscape.Div.Rclass.tab",header=F)
data <- data[-1,]
Alength<-918782050


data_long <- data %>% pivot_longer(-V1, names_to = "Category", values_to = "count")
data_long$Divergence <- as.numeric(gsub("V","",data_long$Category))-1
data_long$count <- as.numeric(as.character(data_long$count))
data_long$percent <- (data_long$count/Alength)
data_long$V1 <- as.character(data_long$V1)
colnames(data_long)[1] <-'Class'


Autos <- ggplot(data_long,aes(Divergence,percent,fill=Class))+
  geom_bar(stat="identity")+
  scale_y_continuous(labels=scales::percent,limits=c(0,0.05))+
  #xlab('% Divergence')+
  #ylab('% of chromosome')+
  xlab('')+
  ylab('')+
  xlim(0,40)+
  theme_bw()+
  theme(axis.title = element_text(size=13),
        axis.text = element_text(size=16),
        legend.position="none")
Autos
### Ancestral W

data<- read.table("Ecyan_HiC_v3.0_ancW.fa.align.landscape.Div.Rclass.tab",header=F)
data <- data[-1,]
AncW<-21329568
  
  
data_long <- data %>% pivot_longer(-V1, names_to = "Category", values_to = "count")
data_long$Divergence <- as.numeric(gsub("V","",data_long$Category))-1
data_long$count <- as.numeric(as.character(data_long$count))
data_long$percent <- (data_long$count/AncW)
data_long$V1 <- as.character(data_long$V1)
colnames(data_long)[1] <-'Class'
  
ancW <- ggplot(data_long,aes(Divergence,percent,fill=Class))+
  geom_bar(stat="identity")+
  scale_y_continuous(labels=scales::percent)+
  #xlab('% Divergence')+
  #ylab('% of chromosome')+
  xlab('')+
  ylab('')+
  xlim(0,40)+
  theme_bw()+
  theme(axis.title = element_text(size=13),
        axis.text = element_text(size=16),
        legend.position="none")


### Added W
data<- read.table("Ecyan_HiC_v3.0_addedW.fa.align.landscape.Div.Rclass.tab",header=F)
data <- data[-1,]
AddW<-44581082


data_long <- data %>% pivot_longer(-V1, names_to = "Category", values_to = "count")
data_long$Divergence <- as.numeric(gsub("V","",data_long$Category))-1
data_long$count <- as.numeric(as.character(data_long$count))
data_long$percent <- (data_long$count/AddW)
data_long$V1 <- as.character(data_long$V1)
colnames(data_long)[1] <-'Class'

custom_colors<-c("#F8766D",'#B79F00','#39B600','#619CFF','#F564E3')

addW <- ggplot(data_long,aes(Divergence,percent,fill=Class))+
  geom_bar(stat="identity")+
  scale_y_continuous(labels=scales::percent,limits=c(0,0.05))+
  #xlab('% Divergence')+
  #ylab('% of chromosome')+
  xlab('')+
  ylab('')+
  xlim(0,40)+
  theme_bw()+
  theme(axis.title = element_text(size=13),
        axis.text = element_text(size=16),
        legend.position="none")+
  scale_fill_manual(values=custom_colors)+
  scale_color_manual(values=custom_colors)

ancW / addW

### new PAR
data<- read.table("Ecyan_HiC_v3.0_newPAR.fa.align.landscape.Div.Rclass.tab",header=F)
data <- data[-1,]
newPAR<-17503221


data_long <- data %>% pivot_longer(-V1, names_to = "Category", values_to = "count")
data_long$Divergence <- as.numeric(gsub("V","",data_long$Category))-1
data_long$count <- as.numeric(as.character(data_long$count))
data_long$percent <- (data_long$count/newPAR)
data_long$V1 <- as.character(data_long$V1)
colnames(data_long)[1] <-'Class'

newPAR <- ggplot(data_long,aes(Divergence,percent,fill=Class))+
  geom_bar(stat="identity")+
  scale_y_continuous(labels=scales::percent,limits=c(0,0.05))+
  #xlab('% Divergence')+
  #ylab('% of chromosome')+
  xlab('')+
  ylab('')+
  xlim(0,40)+
  theme_bw()+
  theme(axis.title = element_text(size=13),
        axis.text = element_text(size=16),
        legend.position="none")

### Ancestral Z
data<- read.table("Ecyan_HiC_v3.0_ancZ.fa.align.landscape.Div.Rclass.tab",header=F)
data <- data[-1,]
ancZ<-74640337


data_long <- data %>% pivot_longer(-V1, names_to = "Category", values_to = "count")
data_long$Divergence <- as.numeric(gsub("V","",data_long$Category))-1
data_long$count <- as.numeric(as.character(data_long$count))
data_long$percent <- (data_long$count/ancZ)
data_long$V1 <- as.character(data_long$V1)
colnames(data_long)[1] <-'Class'

ancZ <- ggplot(data_long,aes(Divergence,percent,fill=Class))+
  geom_bar(stat="identity")+
  scale_y_continuous(labels=scales::percent,limits=c(0,0.05))+
  #xlab('% Divergence')+
  #ylab('% of chromosome')+
  xlab('')+
  ylab('')+
  xlim(0,40)+
  theme_bw()+
  theme(axis.title = element_text(size=13),
        axis.text = element_text(size=16),
        legend.position="none")
  
### Added Z

data<- read.table("Ecyan_HiC_v3.0_addedZ.fa.align.landscape.Div.Rclass.tab",header=F)
data <- data[-1,]
addedZ<-37022348


data_long <- data %>% pivot_longer(-V1, names_to = "Category", values_to = "count")
data_long$Divergence <- as.numeric(gsub("V","",data_long$Category))-1
data_long$count <- as.numeric(as.character(data_long$count))
data_long$percent <- (data_long$count/addedZ)
data_long$V1 <- as.character(data_long$V1)
colnames(data_long)[1] <-'Class'

addZ <- ggplot(data_long,aes(Divergence,percent,fill=Class))+
  geom_bar(stat="identity")+
  scale_y_continuous(labels=scales::percent,limits=c(0,0.05))+
  #xlab('% Divergence')+
  #ylab('% of chromosome')+
  xlab('')+
  ylab('')+
  xlim(0,40)+
  theme_bw()+
  theme(axis.title = element_text(size=13),
        axis.text = element_text(size=16),
        legend.position="none")


(ancW + addW + newPAR) / (ancZ + addZ + Autos)
