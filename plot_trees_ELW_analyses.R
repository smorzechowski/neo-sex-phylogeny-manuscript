# Plot Reference Trees for Expected Likelihood Weight Analyses
# For Figure 6 and Supplementary Figure 2
# Sophie MacRae Orzechowski
# February 2024

library(dplyr)
library(tidyr)
library(ggplot2)
library(ape)
library(ggridges)
library(hrbrthemes)
#library(ggtree)

setwd("~/PhD research/Neo sex chromosome/Recombination suppression")

trees <- read.tree(file="RaXML/original reference tree set/reference_topologies_shortened.txt")
trees

# rename tips with Scientific Names
for(i in c(1:length(trees))) {
  trees[[i]]$tip.label[trees[[i]]$tip.label=="Gpicta"] <- "G_picta_Z"
  trees[[i]]$tip.label[trees[[i]]$tip.label=="Malbogularis"] <- "M_albogularis_Z"
  trees[[i]]$tip.label[trees[[i]]$tip.label=="HeHo"] <- "L_melanops_cassidix"
  trees[[i]]$tip.label[trees[[i]]$tip.label=="Mmel"] <- "M_melanocephala"
  trees[[i]]$tip.label[trees[[i]]$tip.label=="Striat"] <- "P_striatus"
  trees[[i]]$tip.label[trees[[i]]$tip.label=="Chry"] <- "A_chrysorrhoa"
  trees[[i]]$tip.label[trees[[i]]$tip.label=="Nleucotis_Z"] <- "N_leucotis_Z"
  trees[[i]]$tip.label[trees[[i]]$tip.label=="Ecyanotis_Z"] <- "E_cyanotis_Z"
  trees[[i]]$tip.label[trees[[i]]$tip.label=="Pcitreogularis_Z"] <- "P_citreogularis_Z"
  trees[[i]]$tip.label[trees[[i]]$tip.label=="Nleucotis_W"] <- "N_leucotis_W"
  trees[[i]]$tip.label[trees[[i]]$tip.label=="Ecyanotis_W"] <- "E_cyanotis_W"
  trees[[i]]$tip.label[trees[[i]]$tip.label=="Pcitreogularis_W"] <- "P_citreogularis_W"}

# rename tips with Common Names
for(i in c(1:length(trees))) {
  trees[[i]]$tip.label[trees[[i]]$tip.label=="Gpicta"] <- "Painted_Honeyeater_Z"
  trees[[i]]$tip.label[trees[[i]]$tip.label=="Malbogularis"] <- "White-naped_Honeyeater_Z"
  trees[[i]]$tip.label[trees[[i]]$tip.label=="HeHo"] <- "Helmeted_Honeyeater"
  trees[[i]]$tip.label[trees[[i]]$tip.label=="Mmel"] <- "Noisy_Miner"
  trees[[i]]$tip.label[trees[[i]]$tip.label=="Striat"] <- "Striated_Pardalote"
  trees[[i]]$tip.label[trees[[i]]$tip.label=="Chry"] <- "Yellow-rumped_Thornbill"
  trees[[i]]$tip.label[trees[[i]]$tip.label=="Nleucotis_Z"] <- "White-eared_Honeyeater_Z"
  trees[[i]]$tip.label[trees[[i]]$tip.label=="Ecyanotis_Z"] <- "Blue-faced_Honeyeater_Z"
  trees[[i]]$tip.label[trees[[i]]$tip.label=="Pcitreogularis_Z"] <- "Little_Friarbird_Z"
  trees[[i]]$tip.label[trees[[i]]$tip.label=="Nleucotis_W"] <- "White-eared_Honeyeater_W"
  trees[[i]]$tip.label[trees[[i]]$tip.label=="Ecyanotis_W"] <- "Blue-faced_Honeyeater_W"
  trees[[i]]$tip.label[trees[[i]]$tip.label=="Pcitreogularis_W"] <- "Little_Friarbird_W"}

# rename tips with Common Names and added Z-W
for(i in c(1:length(trees))) {
  trees[[i]]$tip.label[trees[[i]]$tip.label=="Gpicta"] <- "Painted_Honeyeater_added_Z"
  trees[[i]]$tip.label[trees[[i]]$tip.label=="Malbogularis"] <- "White-naped_Honeyeater_added_Z"
  trees[[i]]$tip.label[trees[[i]]$tip.label=="HeHo"] <- "Helmeted_Honeyeater"
  trees[[i]]$tip.label[trees[[i]]$tip.label=="Mmel"] <- "Noisy_Miner"
  trees[[i]]$tip.label[trees[[i]]$tip.label=="Striat"] <- "Striated_Pardalote"
  trees[[i]]$tip.label[trees[[i]]$tip.label=="Chry"] <- "Yellow-rumped_Thornbill"
  trees[[i]]$tip.label[trees[[i]]$tip.label=="Nleucotis_Z"] <- "White-eared_Honeyeater_added_Z"
  trees[[i]]$tip.label[trees[[i]]$tip.label=="Ecyanotis_Z"] <- "Blue-faced_Honeyeater_added_Z"
  trees[[i]]$tip.label[trees[[i]]$tip.label=="Pcitreogularis_Z"] <- "Little_Friarbird_added_Z"
  trees[[i]]$tip.label[trees[[i]]$tip.label=="Nleucotis_W"] <- "White-eared_Honeyeater_added_W"
  trees[[i]]$tip.label[trees[[i]]$tip.label=="Ecyanotis_W"] <- "Blue-faced_Honeyeater_added_W"
  trees[[i]]$tip.label[trees[[i]]$tip.label=="Pcitreogularis_W"] <- "Little_Friarbird_added_W"}


#From the list of tip labels, copy the name of the outgroup species

# Scientific Names
#trees <- root(trees,c("A_chrysorrhoa","P_striatus"),resolve.root=TRUE)
#is.rooted(trees)

# Common Names
trees <- root(trees,c("Yellow-rumped_Thornbill","Striated_Pardalote"),resolve.root=TRUE)
is.rooted(trees)

#trees[[9]]$tip.label

layout(matrix(1:10,2,5,byrow =TRUE))
#layout(matrix(1:4,1,4,byrow =TRUE))
par(mar=c(2,0.3,2,0.3))
par(mfrow=c(1,1))
plot(ladderize(trees[[1]],right=FALSE), tip.color=c("blue","blue","blue","blue","black","red",
                                                    "red","black","black","black",
                                                    "black","red"),main="Tree 1: no suppression.; monophyletic",cex=1.1);

plot(ladderize(trees[[2]],right=FALSE), tip.color=c("blue","blue","black","red","red","black",
                                                    "black","blue","blue","black",
                                                    "black","red"),main="Tree 2: no suppression.; paraphyletic",cex=1.1);
plot(ladderize(trees[[3]],right=FALSE), tip.color=c("black","black","blue","blue","red","red",
                                                    "black","black","black","black",
                                                    "black","red"),main="Tree 3: partial suppression.; monophyletic",cex=1.1);
plot(ladderize(trees[[4]],right=FALSE), tip.color=c("black","black","black","black","black","red",
                                                    "red","blue","blue","black",
                                                    "black","red"),main="Tree 4: partial suppression.; paraphyletic",cex=1.1);
plot(ladderize(trees[[5]],right=FALSE), tip.color=c("black","black","blue","blue","red","red",
                                                    "red","black","black","black",
                                                    "black","black"),main="Tree 5: full suppression.; monophyletic",cex=1.1);
plot(ladderize(trees[[6]],right=FALSE), tip.color=c("black","black","black","black","black","blue",
                                                    "blue","red","red","black",
                                                    "black","red"),main="Tree 6: partial suppression.; paraphyletic",cex=1.1);
plot(ladderize(trees[[7]],right=FALSE), tip.color=c("black","black","red","red","red","blue",
                                                    "blue","black","black","black",
                                                    "black","black"),main="Tree 7: full suppression.; monophyletic",cex=1.1);
plot(ladderize(trees[[8]],right=FALSE), tip.color=c("black","black","black","black","black","red",
                                                    "red","blue","blue","red",
                                                    "black","black"),main="Tree 8: partial suppression x2.; paraphyletic",cex=1.1);

#png("tree9_plot.png", width = 7, height = 7, units = "in", res = 600)
plot(ladderize(trees[[9]],right=FALSE), tip.color=c("black","black","red","red","red","black",
                                                    "black","black","blue","blue",
                                                    "black","black"),main="Tree 9: full suppression.; paraphyletic",cex=1.1);
#dev.off()
#trees[[9]]$tip.label

#############################################################################################
# Create 2 x 3 plot for Figure 6 

png("trees_ELW_figure.png", width = 12, height = 6, units = "in", res = 600,bg='transparent',type="cairo")
layout(matrix(1:6,2,3,byrow=FALSE))
par(mar=c(2,0.3,2,0.3))

#T1Z <- drop.tip(trees[[1]],c('N_leucotis_W','E_cyanotis_W','P_citreogularis_W'))
#T1Z$tip.label

T1Z <- drop.tip(trees[[1]],c('White-eared_Honeyeater_W','Blue-faced_Honeyeater_W','Little_Friarbird_W'))
T1Z$tip.label
T1Z$tip.label[T1Z$tip.label=="Helmeted_Honeyeater"] <- "Helmeted_Honeyeater_Z"
T1Z$tip.label[T1Z$tip.label=="Noisy_Miner"] <- "Noisy_Miner_Z"
T1Z$tip.label[T1Z$tip.label=="Striated_Pardalote"] <- "Striated_Pardalote_Z"
T1Z$tip.label[T1Z$tip.label=="Yellow-rumped_Thornbill"] <- "Yellow-rumped_Thornbill_Z"
#T1Z$tip.label[T1Z$tip.label=="White-naped_Honeyeater_added_Z"] <- "White-naped_Honeyeater_ancestral_Z"
#T1Z$tip.label[T1Z$tip.label=="Painted_Honeyeater_added_Z"] <- "Painted_Honeyeater_ancestral_Z"
#T1Z$tip.label[T1Z$tip.label=="White-eared_Honeyeater_added_Z"] <- "White-eared_Honeyeater_ancestral_Z"
#T1Z$tip.label[T1Z$tip.label=="Blue-faced_Honeyeater_added_Z"] <- "Blue-faced_Honeyeater_ancestral_Z"
#T1Z$tip.label[T1Z$tip.label=="Little_Friarbird_added_Z"] <- "Little_Friarbird_ancestral_Z"


plot(ladderize(T1Z,right=FALSE),tip.color=c("blue","blue","blue","blue",
                                            "black","black","black","black",
                                            "black"),main="Tree 1Z: mtDNA tree",cex=1.4,font=2,cex.main=2);

#T2Z <- drop.tip(trees[[2]],c('N_leucotis_W','E_cyanotis_W','P_citreogularis_W'))
#T2Z$tip.label

T2Z <- drop.tip(trees[[2]],c('White-eared_Honeyeater_W','Blue-faced_Honeyeater_W','Little_Friarbird_W'))
T2Z$tip.label
T2Z$tip.label[T2Z$tip.label=="Helmeted_Honeyeater"] <- "Helmeted_Honeyeater_Z"
T2Z$tip.label[T2Z$tip.label=="Noisy_Miner"] <- "Noisy_Miner_Z"
T2Z$tip.label[T2Z$tip.label=="Striated_Pardalote"] <- "Striated_Pardalote_Z"
T2Z$tip.label[T2Z$tip.label=="Yellow-rumped_Thornbill"] <- "Yellow-rumped_Thornbill_Z"
#T2Z$tip.label[T2Z$tip.label=="White-naped_Honeyeater_added_Z"] <- "White-naped_Honeyeater_ancestral_Z"
#T2Z$tip.label[T2Z$tip.label=="Painted_Honeyeater_added_Z"] <- "Painted_Honeyeater_ancestral_Z"
#T2Z$tip.label[T2Z$tip.label=="White-eared_Honeyeater_added_Z"] <- "White-eared_Honeyeater_ancestral_Z"
#T2Z$tip.label[T2Z$tip.label=="Blue-faced_Honeyeater_added_Z"] <- "Blue-faced_Honeyeater_ancestral_Z"
#T2Z$tip.label[T2Z$tip.label=="Little_Friarbird_added_Z"] <- "Little_Friarbird_ancestral_Z"

plot(ladderize(T2Z,right=FALSE),tip.color=c("blue","blue","black","black",
                                            "black","blue","blue","black",
                                            "black"),main="Tree 2Z: species tree",cex=1.4,font=2,cex.main=2);

plot(ladderize(trees[[5]],right=FALSE), tip.color=c("blue","blue","blue","blue","red","red",
                                                    "red","black","black","black",
                                                    "black","black"),main="Tree 5: clade-wide RS",cex=1.4,font=2,cex.main=2);

plot(ladderize(trees[[9]],right=FALSE), tip.color=c("blue","blue","red","red","red","black",
                                                    "black","black","blue","blue",
                                                    "black","black"),main="Tree 9: clade-wide RS",cex=1.4,font=2,cex.main=2);

plot(ladderize(trees[[1]],right=FALSE), tip.color=c("blue","blue","blue","blue","black","red",
                                                    "red","black","black","black",
                                                    "black","red"),main="Tree 1: no RS",cex=1.4,font=2,cex.main=2);

trees[[1]]$tip.label

plot(ladderize(trees[[2]],right=FALSE), tip.color=c("blue","blue","black","red","red","black",
                                                    "black","blue","blue","black",
                                                    "black","red"),main="Tree 2: no RS",cex=1.4,font=2,cex.main=2);


dev.off()


