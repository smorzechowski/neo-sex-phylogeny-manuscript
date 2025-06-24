# Visualize mean sex differences in heterozygosity and read depth on added-Z/W
# Output from FindZX pipeline
# For creating Figure 4 of the finalized manuscript
# Sophia C. M. Orzechowski 
# February 2023
# Updated March 2024

# Load libraries
library(ggplot2)
library(scales)
show_col(hue_pal()(10))
####1####
# Set working directory
setwd("~/PhD research/Neo sex chromosome/findZX/2022-12-21/results/Entomyzon/output/no_synteny/tables")

# Load data
cyan <- read.table('sexDifferences_mean_SD.100000bp.window.tsv',header=TRUE,sep="\t",stringsAsFactors = FALSE)
str(data)

cyan$Species <- "Entomyzon"
cyan$Species <- "Blue-faced Honeyeater"

ggplot(cyan,aes( Mean.genome.coverage.difference..unfiltered.,Mean.heterozygosity.difference))+
  geom_point()+
  scale_x_continuous(labels=comma)



####2####
# Set working directory
setwd("~/PhD research/Neo sex chromosome/findZX/2022-12-21/results/Nesoptilotis/output/no_synteny/tables")


# Load data
nleu <- read.table('sexDifferences_mean_SD.100000bp.window.tsv',header=TRUE,sep="\t",stringsAsFactors = FALSE)

str(nleu)


nleu$Species <- "Nesoptilotis"
nleu$Species <- "White-eared Honeyeater"

ggplot(nleu,aes( Mean.genome.coverage.difference..unfiltered.,Mean.heterozygosity.difference))+
  geom_point()+
  scale_x_continuous(labels=comma)

####3####

# Set working directory
setwd("~/PhD research/Neo sex chromosome/findZX/2022-12-21/results/Philemon/output/no_synteny/tables")

# Load data
phil <- read.table('sexDifferences_mean_SD.100000bp.window.tsv',header=TRUE,sep="\t",stringsAsFactors = FALSE)


phil$Species <- "Philemon"
phil$Species <- "Little Friarbird"

ggplot(phil,aes( Mean.genome.coverage.difference..unfiltered.,Mean.heterozygosity.difference))+
  geom_point()+
  scale_x_continuous(labels=comma)

####4####
setwd("~/PhD research/Neo sex chromosome/findZX/2022-12-21/results/Conopophila/output/no_synteny/tables")
cruf <- read.table('sexDifferences_mean_SD.100000bp.window.tsv',header=TRUE,sep="\t",stringsAsFactors = FALSE)



cruf$Species <- "Conopophila"
cruf$Species <- "Rufous-throated Honeyeater"

####5####
setwd("~/PhD research/Neo sex chromosome/findZX/2022-12-21/results/Ptiloprora/output/no_synteny/tables")
ptilo <- read.table('sexDifferences_mean_SD.100000bp.window.tsv',header=TRUE,sep="\t",stringsAsFactors = FALSE)



ptilo$Species <- "Ptiloprora"
ptilo$Species <- "Rufous-backed Honeyeater"

####6####
setwd("~/PhD research/Neo sex chromosome/findZX/2022-12-21/results/Melilestes/output/no_synteny/tables")
melilest <- read.table('sexDifferences_mean_SD.100000bp.window.tsv',header=TRUE,sep="\t",stringsAsFactors = FALSE)


melilest$Species <- "Melilestes"
melilest$Species <- "Long-billed Honeyeater"
####7####
setwd("~/PhD research/Neo sex chromosome/findZX/2022-12-21/results/Melipotes/output/no_synteny/tables")
melipotes <- read.table('sexDifferences_mean_SD.100000bp.window.tsv',header=TRUE,sep="\t",stringsAsFactors = FALSE)


melipotes$Species <- "Melipotes"
melipotes$Species <- "Smoky Honeyeater"
####8####
setwd("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/findZX/striat_synteny/output/synteny/T_guttatus/tables")


striathet <- read.table('diffHeterozygosity.100000bp.out',header=TRUE,sep="\t",stringsAsFactors = FALSE)
striathet_subset <- striathet[striathet$chr=="NC_044217.2",]
mean(striathet_subset$diff)
sd(striathet_subset$diff)

striatcov <- read.table('diffGenomeCoverage.mismatch.unfiltered.100000bp.out',header=TRUE,sep="\t",stringsAsFactors = FALSE)
striatcov_subset <- striatcov[striatcov$chr=="NC_044217.2",]
mean(striatcov_subset$diff)
sd(striatcov_subset$diff)


####9####
setwd("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/findZX/chry_dataset/output/synteny/T_guttatus/tables")


chryhet <- read.table('diffHeterozygosity.100000bp.out',header=TRUE,sep="\t",stringsAsFactors = FALSE)
chryhet_subset <- chryhet[chryhet$chr=="NC_044217.2",]
mean(chryhet_subset$diff)
sd(chryhet_subset$diff)

chrycov <- read.table('diffGenomeCoverage.mismatch.unfiltered.100000bp.out',header=TRUE,sep="\t",stringsAsFactors = FALSE)
chrycov_subset <- chrycov[chrycov$chr=="NC_044217.2",]
mean(chrycov_subset$diff)
sd(chrycov_subset$diff)

####10####
setwd("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/findZX/mel_synteny_cluster/output/synteny/T_guttatus/tables")


melhet <- read.table('diffHeterozygosity.100000bp.out',header=TRUE,sep="\t",stringsAsFactors = FALSE)
melhet_subset <- melhet[melhet$chr=="NC_044217.2",]
mean(melhet_subset$diff)
sd(melhet_subset$diff)

melcov <- read.table('diffGenomeCoverage.mismatch.unfiltered.100000bp.out',header=TRUE,sep="\t",stringsAsFactors = FALSE)
melcov_subset <- melcov[melcov$chr=="NC_044217.2",]
mean(melcov_subset$diff)
sd(melcov_subset$diff)



# A graph with all species
all_species <- rbind(cyan,phil,nleu,cruf,ptilo,melilest,melipotes)

#write.csv(all_species,"mean_sex_differences.csv")

# Added the last three species manually
# Added common names manually
setwd("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/findZX")
data <- read.csv("mean_sex_differences.csv",header=T,stringsAsFactors = FALSE) 


# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

custom_palette<-c("#A3A500","#E76BF3","#00BF7D","#39B600","#FF62BC","#D89000","#00BFC4","#9590FF","#00B0F6","#F8766D")
data <- data[data$Chromosomes=="ONT_addedZ_V2_contig_4_segment0" | data$Chromosomes=="NC_044217.2",]

p <- ggplot(data,aes(Mean.genome.coverage.difference..unfiltered.,Mean.heterozygosity.difference,
                 ymin=Mean.heterozygosity.difference-SD.heterozygosity.difference,
                 ymax=Mean.heterozygosity.difference+SD.heterozygosity.difference,
                 xmin=Mean.genome.coverage.difference..unfiltered.-SD.genome.coverage.difference..unfiltered.,
                 xmax=Mean.genome.coverage.difference..unfiltered.+SD.genome.coverage.difference..unfiltered., 
                 fill=CommonName,color=CommonName))+
  geom_point(size=5)+
  geom_errorbar(size=1)+
  geom_errorbarh(size=1)+
  geom_vline(aes(xintercept=0),linetype=2,color='black')+
  geom_hline(aes(yintercept=0),linetype=2,color='black')+
  xlab("Δ (Female - Male) Genome Coverage") +
  ylab("Δ (Female - Male) Heterozygosity")+
  ggtitle("A.")+
  theme_bw()+
  theme(axis.title = element_text(size=17),
        axis.text = element_text(size=16),
        legend.text = element_text(size=12),
        legend.title = element_blank(),
        title = element_text(size=16))+
  xlim(-0.5,.5)+
  scale_fill_manual(values=custom_palette)+
  scale_color_manual(values=custom_palette)

p

