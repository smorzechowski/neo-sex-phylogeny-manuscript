# R script to visualize Expected Likelihood Weight Analyses
# For creating Figures 6 and 7 in the finalized manuscript
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

# ancestral PAR data
ancPAR <- read.csv('ancestralPAR.iqtree.logtable.reports.formatted.csv',header=T,stringsAsFactors = FALSE,strip.white = TRUE)
ancPAR_append<- left_join(ancPAR,positiondata_filt,join_by(geneID==gene_id))
ancPAR_append_filt <- ancPAR_append[ancPAR_append$sig.ELW=="+" & ancPAR_append$tr_length>1000,]

# newZ data
newZ<- read.csv('NewZ.iqtree.logtable.reports.formatted.csv',header=T,stringsAsFactors = FALSE,strip.white = TRUE)
newZ_append<- left_join(newZ,positiondata_filt,join_by(geneID==gene_id))
newZ_append_filt <- newZ_append[newZ_append$geneID!='ATP8B1'& newZ_append$sig.ELW=="+"& newZ_append$tr_length>1000,]

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

# Alignment lengths
alignlength <- read.csv("alignment_lengths.csv",header=T)

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

# get estimates of total number of coding regions in added and ancestral region
newZ_append_filt_order <- newZ_append_filt[order(newZ_append_filt$geneID),]
newZ_append_filt_order_nondup <- newZ_append_filt_order[!duplicated(newZ_append_filt_order$geneID),]

newPAR_append_filt_order <- newPAR_append_filt[order(newPAR_append_filt$geneID),]
newPAR_append_filt_order_nondup <- newPAR_append_filt_order[!duplicated(newPAR_append_filt_order$geneID),]
total_added_regions <- 36 + 128

oldZ_append_filt_order <- oldZ_append_filt[order(oldZ_append_filt$geneID),]
oldZ_append_filt_order_nondup <- oldZ_append_filt_order[!duplicated(oldZ_append_filt_order$geneID),]

ancPAR_append_filt_order <- ancPAR_append_filt[order(ancPAR_append_filt$geneID),]
ancPAR_append_filt_order_nondup <- ancPAR_append_filt_order[!duplicated(ancPAR_append_filt_order$geneID),]

# get mean, max, min for transcript length
combined_append_filt_order <- combined_append_filt[order(combined_append_filt$geneID),]
combined_append_filt_order_nondup <- combined_append_filt_order[!duplicated(combined_append_filt_order$geneID),]
all_total <- 232+6+164

mean_length <- mean(combined_append_filt_order_nondup$tr_length)
median_length <- median(combined_append_filt_order_nondup$tr_length)
max_length <- max(combined_append_filt_order_nondup$tr_length)
min_length <- min(combined_append_filt_order_nondup$tr_length)

# Look at total number of trees with max ELW and Calculate alignment length!
newZ_toptree_length <- left_join(newZ_toptree,alignlength)

median(newZ_toptree_length$Alignment_length[newZ_toptree_length$c.ELW<0.5])
median(newZ_toptree_length$Alignment_length[newZ_toptree_length$c.ELW>0.5])
median(newZ_toptree_length$c.ELW[newZ_toptree_length$Tree==9 & newZ_toptree_length$c.ELW>0.5])

# How many loci show a pattern of clade-wide recombination suppression?
table(combined_toptree$Tree)
total_clade_wide <- (73+31+12)/(73+31+15+12+13+4+3)

table(newZ_toptree$Tree)
newZ_clade_wide <- (71+29+11)/(4+4+29+13+11+15+71)

# Separate the oldZ from the NewZ region
combined_toptree$Tree<-factor(combined_toptree$Tree,levels=c('1Z','1','2Z','2','3','4','5','6','7','8','9'))

#https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible

safe_colorblind_palette <- c("#40CEEE", "#999970", "#CC6677", "#117733", "#882255", "#888888", 
                             "#332288", "#DDCC77", "#44AA99", "#661100","#AA4499")

# Categorical variable: mtDNA tree yes/no
combined_toptree$mtDNA<-'yes'
combined_toptree$mtDNA[combined_toptree$Tree=='9']<-'no'
combined_toptree$mtDNA[combined_toptree$Tree=='2Z']<-'no'
combined_toptree$mtDNA[combined_toptree$Tree=='2']<-'no'
combined_toptree$mtDNA[combined_toptree$Tree=='4']<-'no'
combined_toptree$mtDNA[combined_toptree$Tree=='6']<-'no'
combined_toptree$mtDNA[combined_toptree$Tree=='8']<-'no'
combined_toptree$mtDNA[combined_toptree$Tree=='9']<-'no'

combined_toptree_nondup <- combined_toptree[!duplicated(combined_toptree$geneID),]
write.csv(combined_toptree_nondup,"C:/Users/sophi/Documents/PhD research/Manuscripts/Meliphagoidea neo-sex phylogeny/list_of_loci_length.csv")

#########################################################################################
# Maximum expected likelihood weights for Figure 6.
# For completeness including the loci with low phylogenetic signal where multiple trees have the same maximum ELW.

p <- ggplot(combined_toptree[combined_toptree$genepos>833,],aes(x=startpos/1000000, y=after_stat(count),group=Tree,color=Tree,fill=Tree,pattern=mtDNA))+
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
  #guides(colour = guide_legend(override.aes = list(alpha = 1,linetype=1,linewidth=0)))+
  guides(
    colour = guide_legend(
      override.aes = list(alpha = 1, linetype = 0, linewidth = 0)
    ),
    fill = guide_legend(
      override.aes = list(alpha = 1, linetype = 0, linewidth = 0, pattern = "none")
    )
    #    pattern = guide_legend(
    #      override.aes = list(pattern = "none", pattern_density = 0, pattern_spacing = 0, pattern_fill = NA, alpha = 1)
    #    )
  )+
  geom_rug(data=combined_toptree,aes(x=startpos/1000000, y=0))+
  ylab("Max ELW Tree Topology Probability Density")+
  xlab("Neo-Z coordinates (Mb)")+
  theme_classic()+
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=12),
        legend.position="right")+
  scale_x_continuous(breaks=c(0,10,20,30,40,50,60,70,80,90,100,110,120,130))+
  scale_fill_manual(values=safe_colorblind_palette)+
  scale_color_manual(values=safe_colorblind_palette)+
  geom_vline(aes(xintercept=129.17-17.4),linetype='dotted',color='black')+
  geom_vline(aes(xintercept=74.2),linetype=2,color='red')

p

ggsave("C:/Users/sophi/Documents/PhD research/Manuscripts/Meliphagoidea neo-sex phylogeny/Finalized figures for revised manuscript/ELW_recomb_suppr_figure_legend_all_loci.png", 
       plot = p, dpi = 600, width = 8, height = 5, units = "in")


#########################################################################################
# All expected likelihood weights for Figure 7.

combined_newZ_ancPAR$Tree<-as.factor(combined_newZ_ancPAR$Tree)
combined_newZ_ancPAR_newPAR$Tree<-as.factor(combined_newZ_ancPAR_newPAR$Tree)


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



#########################################################################################
### Zoom into the ancestral PAR region for Figure 7.


PAR_safe_colorblind_palette <- c("#882255", "#888888","#332288", "#DDCC77", "#44AA99", "#661100","#AA4499")


agg.data <- aggregate(cbind(startpos,c.ELW) ~ geneID, data = ancPAR_append_filt, max)
agg.data$Tree <- as.factor(c("5","7","9","3","9","5"))
agg.data$c.ELW<- agg.data$c.ELW+0.03

## The margin definition
#par(mar = c(bottom, left, top, right)
par(mar = c(1,1,1,1))

p <- ggplot(ancPAR_append_filt,aes(startpos/1000000,c.ELW,fill=Tree))+
  geom_bar(position='dodge',stat='identity')+
  ylab('Expected Likelihood Weight (ELW)')+
  xlab('Ancestral PAR (Mb)')+
  scale_x_continuous()+
  geom_text(data=agg.data,aes(label=geneID), fontface = "italic")+
  theme_bw()+
  ylim(0,1)+
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=12))+
  scale_fill_manual(values=PAR_safe_colorblind_palette)+
  scale_color_manual(values=PAR_safe_colorblind_palette)

p
ggsave("C:/Users/sophi/Documents/PhD research/Manuscripts/Meliphagoidea neo-sex phylogeny/AncestralPAR_ELW_italic_geneID.png", 
       plot = p, dpi = 450, width = 7, height = 4, units = "in")


