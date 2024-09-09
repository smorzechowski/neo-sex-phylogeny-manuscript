# Neo-sex chromosomes repeat content
# Sophia C. M. Orzechowski
# 17 April 2022


library(ggplot2)
#install.packages("viridis")
library(viridis)
setwd("~/PhD research/Neo sex chromosome")

data <- read.csv("neo_sex_repeats.csv",header=T,stringsAsFactors = FALSE)
data <- read.csv("neo_sex_repeats_updated2024.csv",header=T,stringsAsFactors = FALSE)
data$percent <- round(data$Masked_length/data$Total_length,4)

data$Elements<- factor(data$Elements,levels=c("DNA","LINEs","LTRs","Low_complexity","Simple_repeats","Unclassified"))

ggplot(data,aes(reorder(factor(Region),Order),Masked_length/Total_length,fill=Elements))+
  geom_bar(stat="identity",position='stack')+
 # scale_color_viridis(discrete=TRUE)+
 #scale_fill_viridis(discrete=TRUE)+
  xlab("")+
  #xlab("Region")+
  ylim(0,1)+
  #ylab("Proportion of TEs per region")+
  ylab("")+
  theme_classic()+
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=15),
        axis.text.x =element_text(size=15),
        legend.text = element_text(size=15),
        legend.title = element_text(size=20))
 

  