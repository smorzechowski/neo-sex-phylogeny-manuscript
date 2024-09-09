# Get lists of genes found only on Z (ancestral Z), shared between Z and W (added),
# or found only on W. This is using information from the TOGA annotation, as opposed to
# just coverage from bamstat

# Sophie Orzechowski
# March 2024

library(dplyr)
library(tidyr)
library(stringr)

setwd("~/PhD research/Neo sex chromosome/Recombination suppression")

cyanF<-read.table("Ecyan_204_F_bamstat_output_bed.txt",header=TRUE,sep='\t')
cyanM<-read.table("Ecyan_185_M_bamstat_output_bed.txt",header=TRUE,sep='\t')

Zspec<-read.table("neoZ_specific_transcripts_lengths_genecol.txt",header=FALSE,sep='\t')
Wspec<-read.table("neoW_specific_transcripts_lengths_genecol.txt",header=FALSE,sep='\t')
ZWshared<-read.table("neoZ_neoW_shared_transcripts_lengths_genecol.txt",header=FALSE,sep='\t')


# Pare down gene names
Zspec$gene_id <- str_replace_all(Zspec$V1,'_rna-[X,N][N,M]_[0-9]+\\.[[:digit:]]+\\.[[:digit:]]+','')
Zspec$gene_id <- str_replace_all(Zspec$gene_id,'_ENSGALT[0-9]+.[[:digit:]].[[:digit:]]+','')
Zspec$gene_id <- str_replace(Zspec$gene_id,'\\.[[:digit:]]+\\.[[:digit:]]+','')

ZWshared$gene_id <- str_replace_all(ZWshared$V1,'_rna-[X,N][N,M]_[0-9]+\\.[[:digit:]]+\\.[[:digit:]]+','')
ZWshared$gene_id <- str_replace_all(ZWshared$gene_id,'_ENSGALT[0-9]+.[[:digit:]].[[:digit:]]+','')
ZWshared$gene_id <- str_replace(ZWshared$gene_id,'\\.[[:digit:]]+\\.[[:digit:]]+','')


# Get the longest isoform for Z specific genes!
Zspec_li <- Zspec %>% group_by(gene_id)%>%
                    slice_max(V3,n=1,with_ties=FALSE)


# Get the longest isoform for ZW shared genes!
ZWshared_li <- ZWshared %>% group_by(gene_id)%>%
  slice_max(V3,n=1,with_ties=FALSE)

# Scaffold_2 boundaries - extracted from ragtag scaffolding onto HiC genome
# New PAR: >=
# added-Z RS <= 54.5 Mb
# Old Z > 54.5 Mb

cyanF$boundary <- "OldZ"
cyanF$boundary[cyanF$end>75035781] <- "NewZ"
cyanF$boundary[cyanF$end>112623882] <- "NewPAR"

cyanM$boundary <- "OldZ"
cyanM$boundary[cyanF$end>75035781] <- "NewZ"
cyanM$boundary[cyanF$end>112623882] <- "NewPAR"


cyanF$CDS_length <- cyanF$end - cyanF$start
cyanM$CDS_length <- cyanM$end - cyanM$start 


filt_cyanF <- cyanF[,c('sample','avgcov_0','gene','boundary','CDS_length','start')]
filt_cyanM <- cyanM[,c('sample','avgcov_0','gene','boundary','CDS_length','start')]

covtranF <- filt_cyanF %>% group_by(gene,boundary) %>% summarise(cyanF_mean = mean(avgcov_0),tr_length = sum(CDS_length),startpos=min(start))
covtranM <- filt_cyanM %>% group_by(gene,boundary) %>% summarise(cyanM_mean = mean(avgcov_0),tr_length = sum(CDS_length),startpos=min(start))
covtranF$cyanF_mean_norm <- covtranF$cyanF_mean*1.094

coverage_comp <- left_join(covtranF,covtranM)
coverage_comp$status ='present'
coverage_comp$threshold_low <- coverage_comp$cyanM_mean*0.7
coverage_comp$status[coverage_comp$cyanF_mean_norm<coverage_comp$threshold_low] <- 'W_lost'
coverage_comp$threshold_high <- coverage_comp$cyanM_mean*1.3
coverage_comp$status[coverage_comp$cyanF_mean_norm>coverage_comp$threshold_high] <- 'W_dup'

coverage_comp$status[coverage_comp$cyanM_mean<1] <- 'M_missing'
coverage_comp$status[coverage_comp$cyanM_mean==0 & coverage_comp$cyanF_mean==0] <- 'MF_missing'

# extract the gene ID, first match found for alphanumeric string with hyphen in it
# https://stackoverflow.com/questions/336210/regular-expression-for-alphanumeric-and-underscores
# need to add hyphen to beginning if pattern has a hyphen somewhere in it :)
# https://www.quora.com/How-do-I-specify-a-hyphen-in-a-regex-character-set

coverage_comp$gene_id <- str_extract(coverage_comp$gene,'^[-a-zA-Z0-9]*')



coverage_comp_filt <- coverage_comp[coverage_comp$gene_id!='tama'&
                                      coverage_comp$gene_id!='rna',]


coverage_comp_sort <- coverage_comp_filt[order(coverage_comp_filt$gene_id,-coverage_comp_filt$tr_length),]
coverage_comp_sort$duplicated <- duplicated(coverage_comp_sort$gene_id)

coverage_comp_sort_nondup <- coverage_comp_sort[coverage_comp_sort$duplicated==FALSE,]

table(coverage_comp$status,coverage_comp$boundary)
table(coverage_comp_sort_nondup$status,coverage_comp_sort_nondup$boundary)

# Join the Z specific data
Zspecific_coverage_comp <- left_join(coverage_comp_sort_nondup,Zspec_li)

# The category totals in this won't change, FYI! 
table(Zspecific_coverage_comp$status,Zspecific_coverage_comp$boundary)

# Join the Z W shared data

ZWshared_coverage_comp <- left_join(coverage_comp_sort_nondup,ZWshared_li)

write.table(coverage_comp,'Ecyan_gene_status.txt',row.names=FALSE,sep="\t",quote=FALSE)
write.table(coverage_comp_sort_nondup,'Ecyan_gene_status_longest_transcript.txt',row.names=FALSE,sep="\t",quote=FALSE)

# The longest isoforms of genes deemed 'present' based on the above coverage parsing:

write.table(coverage_comp_sort_nondup$gene_id[coverage_comp_sort_nondup$status=='present'],'Ecyan_gene_id_present_longest_transcript.txt',row.names=FALSE,sep="\t",quote=FALSE)
write.table(coverage_comp_sort_nondup$gene[coverage_comp_sort_nondup$status=='present'],'Ecyan_transcript_id_present_longest_transcript.txt',row.names=FALSE,sep="\t",quote=FALSE)


# The longest isoforms of genes present on the Old Z based on the new TOGA annotation of Ecyan HiC including both neo-Z and neo-W
write.table(Zspecific_coverage_comp$gene[Zspecific_coverage_comp$boundary=='OldZ' & !is.na(Zspecific_coverage_comp$V3)],'Ecyan_transcript_id_OldZ_specific_longest_transcript.txt',row.names=FALSE,sep="\t",quote=FALSE)


# The longest isoforms of genes in each region that are present
write.table(coverage_comp_sort_nondup$gene[coverage_comp_sort_nondup$status=='present'&
                                           coverage_comp_sort_nondup$boundary=='NewPAR'],
            'NewPAR_Ecyan_transcript_id_present_longest_transcript.txt',row.names=FALSE,sep="\t",quote=FALSE)

write.table(coverage_comp_sort_nondup$gene[coverage_comp_sort_nondup$status=='present'&
                                             coverage_comp_sort_nondup$boundary=='NewZ'],
            'NewZ_Ecyan_transcript_id_present_longest_transcript.txt',row.names=FALSE,sep="\t",quote=FALSE)

write.table(coverage_comp_sort_nondup$gene[coverage_comp_sort_nondup$status=='present'&
                                             coverage_comp_sort_nondup$boundary=='OldZ'],
            'OldZ_Ecyan_transcript_id_present_longest_transcript.txt',row.names=FALSE,sep="\t",quote=FALSE)
