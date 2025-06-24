# Visualize sex differences in heterozygosity on added-Z/W
# Output from FindZX pipeline
# For creating Figure 4 of the finalized manuscript
# Sophia C. M. Orzechowski 
# 11 October 2022
# Updated March 2024

# Load libraries
library(ggplot2)
library(scales)
library(patchwork)
show_col(hue_pal()(10))

####1####
# Set working directory
#setwd("~/PhD research/Neo sex chromosome/findZX/cyan_synteny/output/synteny/T_guttatus/tables")
setwd("~/PhD research/Neo sex chromosome/findZX/2022-12-21/results/Entomyzon/output/no_synteny/tables")

# Load data
#data <- read.table('diffHeterozygosity.100000bp.out',header=TRUE,sep="\t",stringsAsFactors = FALSE)
#data <- read.table('diffHeterozygosity.50000bp.out',header=TRUE,sep="\t",stringsAsFactors = FALSE)
data <- read.table('diffHeterozygosity.25000bp.out',header=TRUE,sep="\t",stringsAsFactors = FALSE)
str(data)
#cyan_subset <- data[data$chr=="NC_044217.2",]
cyan_subset <- data[data$chr=="ONT_addedZ_V2_contig_4_segment0",]
cyan_subset$Species <- "Entomyzon"
cyan_subset$Species <- "Blue-faced Honeyeater"

ggplot(cyan_subset,aes(start,diff))+
  geom_point()+
  scale_x_continuous(labels=comma)

#write.table(cyan_subset,file="Cyanotis_NC_044217.2_diffHet.100kb.out")
#write.table(data_subset,file="Entomyzon_ONT_addedZ_V2_contig_4_segment0_diffHet.100kb.out")



####2####
# Set working directory
#setwd("~/PhD research/Neo sex chromosome/findZX/Nesoptilotis/Nesoptilotis/output/synteny/T_guttatus/tables")
setwd("~/PhD research/Neo sex chromosome/findZX/2022-12-21/results/Nesoptilotis/output/no_synteny/tables")


# Load data
#data <- read.table('diffHeterozygosity.100000bp.out',header=TRUE,sep="\t",stringsAsFactors = FALSE)
#data <- read.table('diffHeterozygosity.50000bp.out',header=TRUE,sep="\t",stringsAsFactors = FALSE)
data <- read.table('diffHeterozygosity.25000bp.out',header=TRUE,sep="\t",stringsAsFactors = FALSE)
str(data)
#nleu_subset <- data[data$chr=="NC_044217.2",]
nleu_subset <- data[data$chr=="ONT_addedZ_V2_contig_4_segment0",]

nleu_subset$Species <- "Nesoptilotis"
nleu_subset$Species <- "White-eared Honeyeater"
#write.table(data_subset,file="Neosptilotis_NC_044217.2_diffHet.100kb.out")



ggplot(nleu_subset,aes(start,diff))+
  geom_point()+
  scale_x_continuous(labels=comma)

####3####

# Set working directory
setwd("~/PhD research/Neo sex chromosome/findZX/Philemon/Philemon/output/synteny/T_guttatus/tables")
setwd("~/PhD research/Neo sex chromosome/findZX/2022-12-21/results/Philemon/output/no_synteny/tables")

# Load data
#data <- read.table('diffHeterozygosity.100000bp.out',header=TRUE,sep="\t",stringsAsFactors = FALSE)
#data <- read.table('diffHeterozygosity.50000bp.out',header=TRUE,sep="\t",stringsAsFactors = FALSE)
data <- read.table('diffHeterozygosity.25000bp.out',header=TRUE,sep="\t",stringsAsFactors = FALSE)

str(data)
phil_subset <- data[data$chr=="NC_044217.2",]

phil_subset <- data[data$chr=="ONT_addedZ_V2_contig_4_segment0",]

phil_subset$Species <- "Philemon"
phil_subset$Species <- "Little Friarbird"
#write.table(phil_subset,file="Philemon_NC_044217.2_diffHet.100kb.out")

ggplot(phil_subset,aes(start,diff))+
  geom_point()+
  scale_x_continuous(labels=comma)

####4####
setwd("~/PhD research/Neo sex chromosome/findZX/2022-12-21/results/Conopophila/output/no_synteny/tables")
#data <- read.table('diffHeterozygosity.100000bp.out',header=TRUE,sep="\t",stringsAsFactors = FALSE)
#data <- read.table('diffHeterozygosity.50000bp.out',header=TRUE,sep="\t",stringsAsFactors = FALSE)
data <- read.table('diffHeterozygosity.25000bp.out',header=TRUE,sep="\t",stringsAsFactors = FALSE)

cruf_subset <- data[data$chr=="ONT_addedZ_V2_contig_4_segment0",]
cruf_subset$Species <- "Conopophila"
cruf_subset$Species <- "Rufous-throated Honeyeater"

# A graph with all three species
#all_species <- rbind(cyan_subset,phil_subset,nleu_subset,cruf_subset)
all_species <- rbind(cyan_subset,phil_subset,nleu_subset)


# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# The palette that matches the sex differences graph
#sexdiffpalette <- c("#A3A500","#E76BF3","#D89000","#00B0F6")
sexdiffpalette <- c("#A3A500","#E76BF3","#00B0F6")

# Estimated the size of the ancestral Z from JuiceBox -- good enough I think for this??
#all_species$NeoZ_coords <- all_species$start+74200000
all_species$NeoZ_coords <- all_species$start+74650000
reversed_axis<-c(130,120,110,100,90,80)

all_species_filt <- all_species[all_species$NeoZ_coords<129170000,]

p2 <- ggplot(all_species,aes(NeoZ_coords/1000000,diff,fill=Species,color=Species))+
  geom_point()+
  #geom_smooth()+
#xlim(0,55)+  
 scale_x_reverse(breaks = seq(74,130,10),labels=reversed_axis)+
  #geom_vline(aes(xintercept=17.400000+74.2),linetype='dotted',color='black',size=1.1)+
  geom_vline(aes(xintercept=17.400000+74.65),linetype='dotted',color='black',size=1.1)+
  geom_vline(aes(xintercept=129.17),linetype=2,color='red',size=1.1)+
 # geom_text(label="44,058,037 bp",x=50058037,y=1.5)+
 # annotate(geom="text", x=50058037, y=1.5, label="44,058,037 bp",color="black")+
  xlab("E. cyanotis added-Z Coordinates") +
  xlab("Neo-Z Coordinates (Mbp)") +
  ylab("Δ (Female - Male) Heterozygosity")+
  ggtitle("B.")+
 theme_classic()+
  theme(axis.title = element_text(size=17),
        axis.text = element_text(size=16),
        legend.text = element_text(size=12),
        legend.title = element_blank(),
        title=element_text(size=16))+
  scale_fill_manual(values=sexdiffpalette)+
  scale_color_manual(values=sexdiffpalette)

p2

# Now run neo_sex_het_coverage_differences_figure.R to generate plot p

final_plot <- p/p2
final_plot

ggsave("C:/Users/sophi/Documents/PhD research/Manuscripts/Meliphagoidea neo-sex phylogeny/Finalized figures for revised manuscript/Figure3_Heterozygosity.png", 
       plot = final_plot, dpi = 600, width = 8, height = 10, units = "in")



