# Expected Likelihood Weight Analyses
# Sophie MacRae Orzechowski
# February 2024

library(dplyr)
library(tidyr)
library(ggplot2)
library(ape)
library(ggridges)
library(hrbrthemes)
library(stringr)
library(ggpattern)
#remotes::install_github("coolbutuseless/ggpattern")

setwd("~/PhD research/Neo sex chromosome/Recombination suppression")

# gene position data
origpositiondata<-read.csv('Ecyan_gene_status_longest_transcript.csv',header=T,stringsAsFactors = FALSE,strip.white = TRUE)
newpositiondata<-read.table('Ecyan_HiC_v2.0_neoZ_position_data.txt',header=T)
newpositiondata$genepos <- as.numeric(row.names(newpositiondata))
origpositiondata_filt<-origpositiondata[,c('gene','tr_length','gene_id')]
positiondata_filt<-left_join(origpositiondata_filt,newpositiondata,join_by(gene==transcriptID))

# new PAR data
newPAR<- read.csv('NewPAR.iqtree.logtable.reports.formatted.csv',header=T,stringsAsFactors = FALSE,strip.white = TRUE)
newPAR_append<- left_join(newPAR,positiondata_filt,join_by(geneID==gene_id))
newPAR_append_filt <- newPAR_append[newPAR$sig.ELW=="+" & newPAR_append$tr_length>1000,]


ggplot(newPAR_append_filt,aes(genepos,Tree))+
  geom_raster(aes(fill=c.ELW))+
  scale_fill_continuous(high = "#660033", low = "#f2eaf2")

# ancestral PAR data
ancPAR <- read.csv('ancestralPAR.iqtree.logtable.reports.formatted.csv',header=T,stringsAsFactors = FALSE,strip.white = TRUE)
ancPAR_append<- left_join(ancPAR,positiondata_filt,join_by(geneID==gene_id))
ancPAR_append_filt <- ancPAR_append[ancPAR_append$sig.ELW=="+" & ancPAR_append$tr_length>1000,]

# newZ data
newZ<- read.csv('NewZ.iqtree.logtable.reports.formatted.csv',header=T,stringsAsFactors = FALSE,strip.white = TRUE)
newZ_append<- left_join(newZ,positiondata_filt,join_by(geneID==gene_id))
newZ_append_filt <- newZ_append[newZ_append$geneID!='ATP8B1'& newZ_append$sig.ELW=="+"& newZ_append$tr_length>1000,]


ggplot(newZ_append_filt,aes(genepos,Tree))+
  geom_raster(aes(fill=c.ELW))+
  scale_fill_continuous(high = "#660033", low = "#f2eaf2")

# oldZ data
oldZ<-read.csv('OldZ.iqtree.logtable.reports.combined.csv',header=T,stringsAsFactors = FALSE,strip.white = TRUE)
oldZ$geneID <- str_replace_all(oldZ$geneID,'OldZ_','')
oldZ$geneID <- str_replace_all(oldZ$geneID,'_rna-[X,N][N,M]_[0-9]+','')
oldZ$geneID <- str_replace_all(oldZ$geneID,'_ENSGALT[0-9]+','')
oldZ_append<- left_join(oldZ,positiondata_filt,join_by(geneID==gene_id))
oldZ_append_filt <- oldZ_append[oldZ$sig.ELW=="+"& oldZ_append$tr_length>1000,]
oldZ_append_filt$Tree <- paste(oldZ_append_filt$Tree,'Z',sep="")

# top trees with max ELW
newZ_toptree <- newZ_append_filt %>% group_by(geneID) %>% filter(c.ELW == max(c.ELW))
newPAR_toptree <- newPAR_append_filt %>% group_by(geneID) %>% filter(c.ELW == max(c.ELW))
ancPAR_toptree <- ancPAR_append_filt %>% group_by(geneID) %>% filter(c.ELW == max(c.ELW))
oldZ_toptree <- oldZ_append_filt %>% group_by(geneID) %>% filter(c.ELW == max(c.ELW))

ggplot(newZ_toptree,aes(genepos,Tree))+
  geom_tile(aes(fill=c.ELW))+
  scale_fill_continuous(high = "#660033", low = "#f2eaf2")+
  scale_y_continuous(breaks=c(1,2,3,4,5,6,7,8,9))
  
# combine datasets
combined_append_filt <- rbind(newZ_append_filt,newPAR_append_filt,oldZ_append_filt,ancPAR_append_filt)  
combined_newZ_ancPAR<- rbind(newZ_append_filt,ancPAR_append_filt)
combined_newZ_ancPAR_newPAR<- rbind(newZ_append_filt,ancPAR_append_filt,newPAR_append_filt)
combined_toptree <- as.data.frame(combined_append_filt %>% group_by(geneID) %>% filter(c.ELW == max(c.ELW)))

# Make tree factor variable
combined_toptree$Tree<-as.factor(combined_toptree$Tree)
oldZ_toptree$Tree<-as.factor(oldZ_toptree$Tree)
combined_append_filt$Tree<-as.factor(combined_append_filt$Tree)
oldZ_append_filt$Tree<-as.factor(oldZ_append_filt$Tree)
ancPAR_append_filt$Tree<-as.factor(ancPAR_append_filt$Tree)
newZ_append_filt$Tree<-as.factor(newZ_append_filt$Tree)


ggplot(combined_append_filt,aes(genepos,Tree))+
  geom_raster(aes(fill=c.ELW))+
  scale_fill_continuous(high = "#660033", low = "white")+
  scale_y_continuous(breaks=c(1,2,3,4,5,6,7,8,9))

ggplot()+
  geom_tile(data=combined_toptree,aes(genepos,Tree,fill=c.ELW))+
  scale_fill_continuous(high = "#660033", low = "white")+
  scale_y_continuous(breaks=c(1,2,3,4,5,6,7,8,9))+
#  geom_rect(aes(xmin=564,xmax=800,ymin=-Inf,ymax=Inf),alpha=0.24,fill="lightpink")+
#  geom_rect(aes(xmin=800,xmax=960,ymin=-Inf,ymax=Inf),alpha=0.24,fill="lightblue")+
#  geom_rect(aes(xmin=960,xmax=1075,ymin=-Inf,ymax=Inf),alpha=0.3,fill="lightgray")+
  xlab('Added Z position')+
  theme_classic()+
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=12))


table(combined_toptree$Tree)
table(newZ_toptree$Tree[newZ_toptree$genepos<800])
table(newZ_toptree$Tree[newZ_toptree$genepos>800])
table(newPAR_toptree$Tree)

ggplot(newZ_append_filt,aes(Tree))+
  geom_histogram(bins=9)+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9))


ggplot(combined_toptree,aes(genepos,c.ELW,fill=as.factor(Tree),color=as.factor(Tree)))+
  geom_point()+
  geom_smooth()

ggplot(combined_toptree,aes(x=genepos, y=after_stat(count),group=as.factor(Tree),color=as.factor(Tree)))+
  geom_density()+
  xlim(564,1072)


# Just the OldZ top tree genepos
ggplot(oldZ_toptree,aes(x=genepos, y=after_stat(count),group=Tree,color=Tree,fill=Tree))+
  geom_density(kernel="gaussian",trim=TRUE,linewidth=1,linetype=1,show.legend=TRUE,alpha=0.2)+
  guides(colour = guide_legend(override.aes = list(alpha = 1,linetype=1)))+
  ylab("Max ELW Tree Topology Probability Density")+
  xlab("Gene Position Across Added-Z")

# Just the OldZ top tree startpos
ggplot(oldZ_toptree,aes(x=startpos/1000000, y=after_stat(count),group=Tree,color=Tree,fill=Tree))+
  geom_density(kernel="gaussian",trim=TRUE,linewidth=1,linetype=1,show.legend=TRUE,alpha=0.2)+
  guides(colour = guide_legend(override.aes = list(alpha = 1,linetype=1)))+
  ylab("Max ELW Tree Topology Probability Density")+
  xlab("Gene Position Across Added-Z") 
  
ggplot(oldZ_append_filt,aes(x=genepos, y=after_stat(count),group=Tree,color=Tree,fill=Tree))+
  geom_density(kernel="gaussian",trim=TRUE,linewidth=1,linetype=1,show.legend=TRUE,alpha=0.2)+
  guides(colour = guide_legend(override.aes = list(alpha = 1,linetype=1)))+
  ylab("Max ELW Tree Topology Probability Density")+
  xlab("Gene Position Across Added-Z")


ggplot(combined_toptree,aes(x=genepos, y=after_stat(count),group=Tree,color=Tree,fill=Tree))+
  geom_density(kernel="gaussian",trim=TRUE,linewidth=1,linetype=1,bounds=c(1,1675),show.legend=TRUE,alpha=0.2)+
  guides(colour = guide_legend(override.aes = list(alpha = 1,linetype=1)))+
  ylab("Max ELW Tree Topology Probability Density")+
  xlab("Gene Position Across Added-Z")+
  theme_bw()+
  theme(axis.title = element_text(size=16))
#  geom_vline(xintercept=960,linetype=3)+
#  geom_vline(xintercept=800,linetype=3)+
#  annotate("text",x=1000,y=0.45,label="New PAR")+
#  annotate("text",x=600,y=0.45,label="Inversion")

ggplot(combined_toptree,aes(x=startpos/1000000, y=after_stat(count),group=Tree,color=Tree,fill=Tree))+
  geom_density(kernel="gaussian",trim=TRUE,linewidth=1,linetype=1,show.legend=TRUE,alpha=0.2)+
  guides(colour = guide_legend(override.aes = list(alpha = 1,linetype=1)))+
  ylab("Max ELW Tree Topology Probability Density")+
  xlab("Gene Position Across Added-Z")+
  theme_bw()+
  theme(axis.title = element_text(size=16))
#  geom_vline(xintercept=960,linetype=3)+
#  geom_vline(xintercept=800,linetype=3)+
#  annotate("text",x=1000,y=0.45,label="New PAR")+
#  annotate("text",x=600,y=0.45,label="Inversion")

# Separate the oldZ from the NewZ region
combined_toptree$Tree<-factor(combined_toptree$Tree,levels=c('1Z','1','2Z','2','3','4','5','6','7','8','9'))

#https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible

safe_colorblind_palette <- c("#40CEEE", "#999970", "#CC6677", "#117733", "#882255", "#888888", 
                             "#332288", "#DDCC77", "#44AA99", "#661100","#AA4499")

combined_toptree$mtDNA<-'yes'
combined_toptree$mtDNA[combined_toptree$Tree=='9']<-'no'
combined_toptree$mtDNA[combined_toptree$Tree=='2Z']<-'no'
combined_toptree$mtDNA[combined_toptree$Tree=='2']<-'no'
combined_toptree$mtDNA[combined_toptree$Tree=='4']<-'no'
combined_toptree$mtDNA[combined_toptree$Tree=='6']<-'no'
combined_toptree$mtDNA[combined_toptree$Tree=='8']<-'no'
combined_toptree$mtDNA[combined_toptree$Tree=='9']<-'no'

ggplot(combined_toptree[combined_toptree$genepos>833,],aes(x=startpos/1000000, y=after_stat(count),group=Tree,color=Tree,fill=Tree,pattern=mtDNA))+
  geom_density_pattern(data=combined_toptree[combined_toptree$genepos<=833,],aes(x=startpos/1000000, 
                                                                                 y=after_stat(count),
                                                                                 group=Tree,
                                                                                 color=Tree,
                                                                                 fill=Tree,
                                                                                 pattern=mtDNA),
                       kernel="cosine",
                       trim=TRUE,
                       linewidth=0.1,
                       linetype=1,
                       show.legend=TRUE,
                       alpha=0.3,
                       pattern_fill="black",
                       pattern_density=0.1,
                       pattern_spacing=0.025,
                       pattern_angle = 45)+
  geom_density(kernel="cosine",trim=TRUE,linewidth=0.1,linetype=1,show.legend=TRUE,alpha=0.3)+
  geom_density_pattern(position = 'identity',
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.025,
                   pattern_key_scale_factor = 1,show.legend = FALSE,
                   alpha=0.3,trim=TRUE,kernel="cosine") +
  scale_pattern_manual(values = c(yes = "stripe", no = "none")) +
  guides(colour = guide_legend(override.aes = list(alpha = 1,linetype=1)))+
  geom_rug(data=combined_toptree,aes(x=startpos/1000000, y=0))+
  ylab("Max ELW Tree Topology Probability Density")+
  xlab("Neo-Z coordinates (Mb)")+
  theme_classic()+
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=12))+
  scale_x_continuous(breaks=c(0,10,20,30,40,50,60,70,80,90,100,110,120,130))+
  scale_fill_manual(values=safe_colorblind_palette)+
  scale_color_manual(values=safe_colorblind_palette)+
  geom_vline(aes(xintercept=129.17-17.4),linetype='dotted',color='black')+
  geom_vline(aes(xintercept=74.2),linetype=2,color='red')


ggplot(combined_append_filt,aes(x=genepos, y=after_stat(count),group=Tree,color=Tree,fill=Tree))+
  geom_density(kernel="gaussian",trim=TRUE,linewidth=1,linetype=1,bounds=c(1,1072),show.legend=TRUE,alpha=0.2)+
  guides(colour = guide_legend(override.aes = list(alpha = 1,linetype=1)))+
  ylab("Tree Topology Probability Density")+
  xlab("Gene Position Across Added-Z")+
  theme_bw()+
  theme(axis.title = element_text(size=16))+
  #geom_vline(xintercept=960,linetype=3)+
  #geom_vline(xintercept=800,linetype=3)+
  #annotate("text",x=1000,y=0.45,label="New PAR")+
  #annotate("text",x=600,y=0.45,label="Inversion")

# https://datawookie.dev/blog/2022/10/scaling-density-plots/
# https://stackoverflow.com/questions/17506053/making-line-legends-for-geom-density-in-ggplot2-in-r

ggplot(combined_toptree) +
  geom_density_ridges(
    aes(x = genepos, y = after_stat(count), group=as.factor(Tree), fill = as.factor(Tree), height = after_stat(count)),
    stat="density",
    scale = 1.5,
    alpha = 0.5)

#  facet_wrap(~Tree)

# try cumulative density plot
ggplot(combined_toptree, aes(genepos,group=Tree)) +
  stat_ecdf(geom = "step")

### Zoom into the ancestral PAR region

ggplot(ancPAR_append_filt,aes(genepos,c.ELW,fill=Tree))+
  geom_bar(position='dodge',stat='identity')+
  ylab('Expected Likelihood Weight (ELW)')+
  xlab('Gene position in the ancestral PAR')

PAR_safe_colorblind_palette <- c("#882255", "#888888","#332288", "#DDCC77", "#44AA99", "#661100","#AA4499")

#"#40CEEE", "#999970", "#CC6677", "#117733"

agg.data <- aggregate(cbind(startpos,c.ELW) ~ geneID, data = ancPAR_append_filt, max)
agg.data$Tree <- as.factor(c("5","7","9","3","9","5"))
agg.data$c.ELW<- agg.data$c.ELW+0.03
#svg(filename = "AncestralPAR_ELW.svg",width = 10, height = 20)
png(filename = "AncestralPAR_ELW.png",width = 2400, height = 1200,res=600)
## The margin definition
#par(mar = c(bottom, left, top, right)
par(mar = c(1,1,1,1))

ggplot(ancPAR_append_filt,aes(startpos/1000000,c.ELW,fill=Tree))+
  geom_bar(position='dodge',stat='identity')+
  ylab('Expected Likelihood Weight (ELW)')+
  xlab('Ancestral PAR (Mb)')+
  scale_x_continuous()+
  geom_text(data=agg.data,aes(label=geneID))+
  theme_bw()+
  ylim(0,1)+
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=12))+
  scale_fill_manual(values=PAR_safe_colorblind_palette)+
  scale_color_manual(values=PAR_safe_colorblind_palette)


## Saving the file
dev.off()


  

#breaks=c(70,80,90,100,110,120,130

ggplot(newZ_append_filt,aes(genepos,c.ELW,fill=Tree))+
  geom_bar(position='dodge',stat='identity')+
  ylab('Expected Likelihood Weight (ELW)')+
  xlab('Gene position in the added Z')


combined_newZ_ancPAR$Tree<-as.factor(combined_newZ_ancPAR$Tree)
combined_newZ_ancPAR_newPAR$Tree<-as.factor(combined_newZ_ancPAR_newPAR$Tree)
ggplot(combined_newZ_ancPAR_newPAR,aes(startpos/1000000,c.ELW,fill=Tree))+
  geom_bar(position='dodge',stat='identity',width=1)+
  ylab('Expected Likelihood Weight (ELW)')+
  xlab('Gene position in the added Z')

combined_append_filt$Tree<-factor(combined_append_filt$Tree,levels=c('1Z','2Z','1','2','3','4','5','6','7','8','9'))

ggplot(combined_append_filt,aes(startpos/1000000,c.ELW,fill=Tree))+
  geom_bar(position='dodge',stat='identity',width=0.8)+
  ylab('Expected Likelihood Weight (ELW)')+
  xlab('Neo-Z coordinates (Mb)')+
  scale_fill_manual(values=safe_colorblind_palette)+
  scale_color_manual(values=safe_colorblind_palette)+
  geom_vline(aes(xintercept=129.6-17.4),linetype='dotted',color='black',size=0.8)+
  geom_vline(aes(xintercept=74.2),linetype=2,color='red',size=.8)+
  scale_x_continuous(breaks=c(0,10,20,30,40,50,60,70,80,90,100,110,120,130))+
  theme_classic()+
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=12))
  
safe_colorblind_palette <- c("#40CEEE", "#999970", "#CC6677", "#117733", "#882255", "#888888", 
                             "#332288", "#DDCC77", "#44AA99", "#661100","#AA4499")



ggplot(oldZ_append_filt,aes(genepos,c.ELW,fill=Tree))+
  geom_bar(position='dodge',stat='identity',width=3)+
  ylab('Expected Likelihood Weight (ELW)')+
  xlab('Gene position in the ancestral Z')


