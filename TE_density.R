# Visualizing TE densities across chromosomes
# Sophie MacRae Orzechowski
# 27 March 2024
# Updated 21 October 2025

setwd("~/PhD research/Neo sex chromosome/repeatmasker")

library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(scales)
library(patchwork)

#################################################################################################################
###### Blue-faced Honeyeater plots ######

data<- read.csv("Ecyan_HiC_v3.0_neoW.fa.out.csv",header=T)
data$class <- as.factor(data$class)
data$begin<-as.numeric(data$begin)
data$class_family<-paste(data$class,data$family,sep="/")

Helitrons<- sum(data$end[data$family=="Helitron"]-data$begin[data$family=="Helitron"])

# get levels in the order used by ggplot
class_levels <- levels(factor(data$class))
n <- length(class_levels)

# build the default ggplot hue palette and then replace the RC entry
default_cols <- hue_pal()(n)               # default ggplot colours
names(default_cols) <- class_levels
default_cols["RC"] <- "red"                # replace only RC

#### Neo-W
Ecyan_W<- ggplot(data,aes(x=begin/1000000, y=after_stat(count),group=class,color=class,fill=class))+
  geom_density(kernel="gaussian",trim=TRUE,linewidth=1,linetype=1,show.legend=TRUE,alpha=0.2)+
  guides(colour = guide_legend(override.aes = list(alpha = 1,linetype=1)))+
  ylab("Probability Density Count Scaled")+
  xlab("neo-W coordinates (Mb)")+
  ggtitle("A.")+
  theme_bw()+
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=12),
        plot.title=element_text(size=16,face="bold"))+
  scale_x_continuous(breaks=c(0,10,20,30,40,50,60,70,80,90))+
  geom_vline(aes(xintercept=27.286927),linetype=2,color='red')+
  scale_color_manual(values = default_cols) +
  scale_fill_manual(values = default_cols)

Ecyan_W 



#### Neo-Z

data<- read.csv("Ecyan_HiC_v3.0_neoZ.fa.out.csv",header=T)
data$class <- as.factor(data$class)
data$begin<-as.numeric(data$begin)
data_sort <- data[order(data$begin,data$repeat.),]

Helitrons<- sum(data$end[data$family=="Helitron"]-data$begin[data$family=="Helitron"])

# get levels in the order used by ggplot
class_levels <- levels(factor(data$class))
n <- length(class_levels)

# build the default ggplot hue palette and then replace the RC entry
default_cols <- hue_pal()(n)               # default ggplot colours
names(default_cols) <- class_levels
default_cols["RC"] <- "red"                # replace only RC

Ecyan_Z <- ggplot(data,aes(x=begin/1000000, y=after_stat(count),group=class,color=class,fill=class))+
  geom_density(kernel="gaussian",trim=TRUE,linewidth=1,linetype=1,show.legend=TRUE,alpha=0.2)+
  guides(colour = guide_legend(override.aes = list(alpha = 1,linetype=1)))+
  ylab("Probability Density Count Scaled")+
  xlab("neo-Z coordinates (Mb)")+
  ggtitle("B.")+
  theme_bw()+
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=12),
        plot.title=element_text(size=16,face="bold"))+
  scale_x_continuous(breaks=c(0,10,20,30,40,50,60,70,80,90,100,110,120,130))+
  geom_vline(aes(xintercept=74.640337),linetype=2,color='red')+
  ylim(0,700)+
  scale_color_manual(values = default_cols) +
  scale_fill_manual(values = default_cols)

Ecyan_Z

Ecyan_combined <- Ecyan_W / Ecyan_Z

# save plot
ggsave("C:/Users/sophi/Documents/PhD research/Manuscripts/Meliphagoidea neo-sex phylogeny/Finalized figures for revised manuscript/Ecyan_TE_landscape_neoW_neoZ.png", 
       plot = Ecyan_combined, dpi = 300, width = 8, height = 9, units = "in")

#################################################################################################################
###### Noisy Miner plots #######

### Mmel Chr1A + 5
data<- read.csv("Mmel_HiC_v1.0_Chr1A5.fa.out.csv",header=T)
data$class <- as.factor(data$class)
data$begin<-as.numeric(data$begin)
data_sort <- data[order(data$begin,data$repeat.),]
Helitrons<- sum(data$end[data$family=="Helitron"]-data$begin[data$family=="Helitron"])

mel5_1a <- ggplot(data,aes(x=begin/1000000, y=after_stat(count),group=class,color=class,fill=class))+
  geom_density(kernel="gaussian",trim=TRUE,linewidth=1,linetype=1,show.legend=TRUE,alpha=0.2)+
  guides(colour = guide_legend(override.aes = list(alpha = 1,linetype=1)))+
  ylab("Probability Density Count Scaled")+
 # ylab("")+
  xlab("Noisy Miner Chr1A+5 coordinates (Mb)")+
  ggtitle("A.")+
  geom_vline(aes(xintercept=62.451976),linetype=2,color='blue')+
  theme_bw()+
  theme(axis.title = element_text(size=16),
        plot.title=element_text(size=16,face="bold"))+
  scale_x_continuous(breaks=c(0,10,20,30,40,50,60,70,80,90,100,110,120))+
  ylim(0,600)+
  scale_color_manual(values = default_cols) +
  scale_fill_manual(values = default_cols)

mel5_1a

### Mmel ChrZ
data<- read.csv("Mmel_HiC_v1.0_Z.fa.out.csv",header=T)
data$class <- as.factor(data$class)
data$begin<-as.numeric(data$begin)
data_sort <- data[order(data$begin,data$repeat.),]

Helitrons<- sum(data$end[data$family=="Helitron"]-data$begin[data$family=="Helitron"])

ggplot(data,aes(x=begin/1000000, y=after_stat(count),group=class,color=class,fill=class))+
  geom_density(kernel="gaussian",trim=TRUE,linewidth=1,linetype=1,show.legend=TRUE,alpha=0.2)+
  guides(colour = guide_legend(override.aes = list(alpha = 1,linetype=1)))+
  ylab("Probability Density Count Scaled")+
  xlab("Noisy Miner Z coordinates (Mb)")+
  theme_bw()+
  theme(axis.title = element_text(size=16),
        plot.title=element_text(size=16,face="bold"))+
  scale_x_continuous(breaks=c(0,10,20,30,40,50,60,70,80))+
  ylim(0,600)

#### Mmel ChrZ REDO - PAR was missing
data<- read.csv("Mmel_HiC_v3.0_Zchrom.fa.out.csv",header=T)
data$class <- as.factor(data$class)
data$begin<-as.numeric(data$begin)
data_sort <- data[order(data$begin,data$repeat.),]

Helitrons<- sum(data$end[data$family=="Helitron"]-data$begin[data$family=="Helitron"])
#max(data$begin)

# flip the start coordinates of this chromosome for comparability with the other species
data$begin_rev <- 74220008 - data$begin



melZ <- ggplot(data,aes(x=begin_rev/1000000, y=after_stat(count),group=class,color=class,fill=class))+
  geom_density(kernel="gaussian",trim=TRUE,linewidth=1,linetype=1,show.legend=TRUE,alpha=0.2)+
  guides(colour = guide_legend(override.aes = list(alpha = 1,linetype=1)))+
  ylab("Probability Density Count Scaled")+
  xlab("Noisy Miner Z coordinates (Mb)")+
  theme_bw()+
  ggtitle("B.")+
  theme(axis.title = element_text(size=16),
        plot.title=element_text(size=16,face="bold"))+
  scale_x_continuous(breaks=c(0,10,20,30,40,50,60,70,80))+
  ylim(0,600)+
  scale_color_manual(values = default_cols) +
  scale_fill_manual(values = default_cols)+
  geom_vline(aes(xintercept=74.22),linetype=2,color='red')

mel_combined <- mel5_1a / melZ
 # plot_layout(guides="collect")


mel_combined
# save plot
ggsave("C:/Users/sophi/Documents/PhD research/Manuscripts/Meliphagoidea neo-sex phylogeny/Finalized figures for revised manuscript/mel_TE_landscape_Z_1a_5.png", 
       plot = mel_combined, dpi = 300, width = 8, height = 9, units = "in")


#################################################################################################################
###### Painted Honeyeater plot ######

# Gpicta scaffolded to Ecyan sans neo-W genome

data<- read.csv("Gpicta_scaffolded_to_Ecyan_sans_neo_W_scaffold_2.fasta.out.csv",header=T)
data$class <- as.factor(data$class)
data$begin<-as.numeric(data$begin)
data$class_family<-paste(data$class,data$family,sep="/")

Helitrons<- sum(data$end[data$family=="Helitron"]-data$begin[data$family=="Helitron"])

picta <- ggplot(data,aes(x=begin/1000000, y=after_stat(count),group=class,color=class,fill=class))+
  geom_density(kernel="gaussian",trim=TRUE,linewidth=1,linetype=1,show.legend=TRUE,alpha=0.2)+
  guides(colour = guide_legend(override.aes = list(alpha = 1,linetype=1)))+
  ylab("Probability Density \n Count Scaled")+
  xlab("neo-Z coordinates (Mb)")+
  theme_bw()+
  ggtitle("A. Painted Honeyeater")+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=12),
        legend.position = "none",
        plot.title=element_text(size=16,face="bold"))+
  scale_x_continuous(breaks=c(0,10,20,30,40,50,60,70,80,90,100,110,120,130))+
  geom_vline(aes(xintercept=74.640337),linetype=2,color='red')+
  scale_color_manual(values = default_cols) +
  scale_fill_manual(values = default_cols)+
  ylim(0,600)


picta

#################################################################################################################
###### Yellow-rumped Thornbill plot #######

# Chry scaffolded to Ecyan sans neo-W genome
data<- read.csv("Chry_scaffolded_to_Ecyan_sans_neoW_scaffold_2.fasta.out.csv",header=T)
data$class <- as.factor(data$class)
data$begin<-as.numeric(data$begin)
data$class_family<-paste(data$class,data$family,sep="/")

Helitrons<- sum(data$end[data$family=="Helitron"]-data$begin[data$family=="Helitron"])

#Ecyan end of ancestral Z: 74.640337
#Striat homology end of Z: ~73.5
#PAR homology
#scaffold_2_RagTag     135919842  73518127   73518243   +  scaffold_2    129165907  74608895   74609011   116      116       60
data<- data[data$begin<74000000,]

chry <- ggplot(data,aes(x=begin/1000000, y=after_stat(count),group=class,color=class,fill=class))+
  geom_density(kernel="gaussian",trim=TRUE,linewidth=1,linetype=1,show.legend=TRUE,alpha=0.2)+
  guides(colour = guide_legend(override.aes = list(alpha = 1,linetype=1)))+
  ylab("Probability Density \n Count Scaled")+
  #ylab("")+
  ylim(0,600)+
  xlab("Z coordinates (Mb)")+
  ggtitle("B. Yellow-rumped Thornbill")+
  theme_bw()+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=12),
        plot.title = element_text(size=16,face="bold"))+
  scale_x_continuous(breaks=c(0,10,20,30,40,50,60,70,80,90,100,110,120,130))+
  #geom_vline(aes(xintercept=74.640337),linetype=2,color='red')
  geom_vline(aes(xintercept=73.5),linetype=2,color='red')+
  scale_color_manual(values = default_cols) +
  scale_fill_manual(values = default_cols)

chry

#################################################################################################################
###### Striated Pardalote plot #######

# Striat scaffolded to Ecyan sans neo-W genome
#Ecyan end of ancestral Z: 74.640337
#Striat homology end of Z: ~75.9
# PAR homology
# scaffold_2_RagTag    135708702  75853070   75900517   -  scaffold_2    129165907  74534429   74579927   4032     51110    60


data<- read.csv("Striat_scaffolded_to_Ecyan_sans_neoW_scaffold_2.fasta.out.csv",header=T)
data$class <- as.factor(data$class)
data$begin<-as.numeric(data$begin)
data$class_family<-paste(data$class,data$family,sep="/")

Helitrons<- sum(data$end[data$family=="Helitron"]-data$begin[data$family=="Helitron"])

data<- data[data$begin<76000000,]

striat <- ggplot(data,aes(x=begin/1000000, y=after_stat(count),group=class,color=class,fill=class))+
  geom_density(kernel="gaussian",trim=TRUE,linewidth=1,linetype=1,show.legend=TRUE,alpha=0.2)+
  guides(colour = guide_legend(override.aes = list(alpha = 1,linetype=1)))+
  ylab("Probability Density \n Count Scaled")+
  xlab("Z coordinates (Mb)")+
  ggtitle("C. Striated Pardalote")+
  theme_bw()+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=12),
        legend.position="none",
        plot.title = element_text(size=16,face="bold"))+
  scale_x_continuous(breaks=c(0,10,20,30,40,50,60,70,80))+
  # scale_x_continuous(breaks=c(0,10,20,30,40,50,60,70,80,90,100,110,120,130))+
  #geom_vline(aes(xintercept=74.640337),linetype=2,color='red')
  geom_vline(aes(xintercept=75.9),linetype=2,color='red')+
  ylim(0,600)+
  scale_color_manual(values = default_cols) +
  scale_fill_manual(values = default_cols)

striat

plots_combined <- picta / chry / striat +
  plot_layout(guides="collect")

plots_combined

ggsave("C:/Users/sophi/Documents/PhD research/Manuscripts/Meliphagoidea neo-sex phylogeny/Finalized figures for revised manuscript/picta_outgroups_TE_landscape_Z.png", 
       plot = plots_combined, dpi = 300, width = 8, height = 10, units = "in")



#################################################################################################################
###### Helitron totals and genome distribution plots #######

# Look at helitrons, total bp masked over all genomes
cyan<- read.csv('Ecyan_Helitrons.csv',header=T)
cyan$Species='Ecyan'
cyan$helsum=cyan$end-cyan$begin
HelCyan<- sum(cyan$end-cyan$begin)

cyanchrom <-read.table('cyan_out_scaffolds_final.chrom.sizes',header=F)
cyanchrom$start <- 1
colnames(cyanchrom) <- c("query","chrom_end","chrom_start")

cyanchrom_filt <- cyanchrom[cyanchrom$chrom_end>400000,]
cyanchrom_filt <- cyanchrom_filt[cyanchrom_filt$query!="scaffold_26" & cyanchrom_filt$query!="scaffold_16"
                                 & cyanchrom_filt$query!="scaffold_14"
                                 & cyanchrom_filt$query!="scaffold_13"
                                 & cyanchrom_filt$query!="scaffold_9",]

cyanchrom_filt$query<-factor(cyanchrom_filt$query,levels=c("scaffold_1","scaffold_2","scaffold_3","scaffold_4",
                                                         "scaffold_5","scaffold_6","scaffold_7","scaffold_8",
                                                         "scaffold_10","scaffold_11","scaffold_12",
                                                         "scaffold_15",
                                                         "scaffold_17","scaffold_18","scaffold_19","scaffold_20",
                                                         "scaffold_21","scaffold_22","scaffold_23","scaffold_24",
                                                         "scaffold_25","scaffold_27","scaffold_28",
                                                         "scaffold_29","scaffold_30","scaffold_31","scaffold_32",
                                                         "scaffold_33","scaffold_34","scaffold_35","scaffold_36",
                                                         "scaffold_37","scaffold_38","scaffold_39","scaffold_40"))



cyan_combined <- left_join(cyan,cyanchrom)
cyan_combined <- cyan_combined[cyan_combined$chrom_end>400000,]

cyan_combined$end_diff <- cyan_combined$chrom_end - cyan_combined$begin

cyan_combined$query<-factor(cyan_combined$query,levels=c("scaffold_1","scaffold_2","scaffold_3","scaffold_4",
                                                         "scaffold_5","scaffold_6","scaffold_7","scaffold_8",
                                                         "scaffold_9","scaffold_10","scaffold_11","scaffold_12",
                                                         "scaffold_13","scaffold_14","scaffold_15","scaffold_16",
                                                         "scaffold_17","scaffold_18","scaffold_19","scaffold_20",
                                                         "scaffold_21","scaffold_22","scaffold_23","scaffold_24",
                                                         "scaffold_25","scaffold_26","scaffold_27","scaffold_28",
                                                         "scaffold_29","scaffold_30","scaffold_31","scaffold_32",
                                                         "scaffold_33","scaffold_34","scaffold_35","scaffold_36",
                                                         "scaffold_37","scaffold_38","scaffold_39","scaffold_40"))

#cyan_combined$query <- str_replace_all(cyan_combined$query,"scaffold_","Chr")
#cyanchrom_filt$query <- str_replace_all(cyanchrom_filt$query,"scaffold_","Chr")

ggplot(na.omit(cyan_combined),aes(query,begin/1000000))+
  geom_point()+
  geom_point(data=cyanchrom_filt,aes(x=query,y=chrom_end/1000000,color='red'))+
  coord_flip()+
  ylab('Distribution of putative Helitrons across chromosomes (Mb)')+
  xlab('')+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=14))
  

heho<- read.csv('HeHo_helitrons.csv',header=T)
heho$Species='HeHo'
heho$helsum=heho$end-heho$begin

HelHeHo<- sum(heho$end-heho$begin)

mel<- read.csv('Mmel_helitrons.csv',header=T)
mel$Species='Mmel'
mel$helsum=mel$end-mel$begin

Helmel<- sum(mel$end-mel$begin)

chry<- read.csv('Chry_helitrons.csv',header=T)
chry$Species='Chry'
chry$helsum=chry$end-chry$begin
Helchry<- sum(chry$end-chry$begin)

striat<- read.csv('Striat_helitrons.csv',header=T)
striat$Species='Striat'
striat$helsum=striat$end-striat$begin
Helstriat<- sum(striat$end-striat$begin)

malurus<- read.csv('Malurus_helitrons.csv',header=T)
malurus$Species='Malurus'
malurus$helsum=malurus$end-malurus$begin
Helmalurus<- sum(malurus$end-malurus$begin)

hola<- read.csv('HoLa_helitrons.csv',header=T)
hola$Species='HoLa'
hola$helsum=hola$end-hola$begin
Helhola<- sum(hola$end-hola$begin)



heldata <- rbind(cyan,mel,chry,striat,malurus,heho)
heldata$Species <- factor(heldata$Species,levels=c("Ecyan","HeHo","Mmel","Striat","Chry","Malurus"))
ggplot(heldata,aes(Species,helsum/1000000))+
  coord_flip()+
  geom_bar(stat='identity')+
  scale_y_continuous()+
  ylim(0,4)+
  ylab('Sum total helitron elements (Mb)')+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=14))


