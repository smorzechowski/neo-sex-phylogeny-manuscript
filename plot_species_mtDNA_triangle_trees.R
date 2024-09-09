# Plot species tree and mtDNA with triangles showing proportional size of clades
# Sophie MacRae Orzechowski
# March 2024
#



library(phytools)
library(geiger)

setwd("~/PhD research/Neo sex chromosome/UCE data/genera59_mafft_min75percent_gblocks0.65_clean_nexus_raxml")


### ** Examples

## first create our backbone tree with
## random subtree diversities
tree<-phytools:::lambdaTree(pbtree(n=10),lambda=0.5)
## create a translation table
## leaving a couple of single-taxon clades for fun
tip.label<-sample(tree$tip.label,8)
clade.label<-LETTERS[1:8]
N<-ceiling(runif(n=8,min=1,max=20))
## set crown node depth to 1/2 the maximum depth
depth<-sapply(tip.label,function(x,y)
  0.5*y$edge.length[which(tree$edge[,2]==
                            which(y$tip.label==x))],y=tree)
trans<-data.frame(tip.label,clade.label,N,depth)
rownames(trans)<-NULL
## here's what trans looks like
print(trans)
par(mar=c(1,1,1,1))
layout(matrix(1,1,1,byrow=FALSE))
obj<-phylo.toBackbone(tree,trans)
## plot
plot(obj)

par(mar=c(5.1,4.1,4.1,2.1))

#https://github.com/liamrevell/phytools/issues/83
#http://blog.phytools.org/2017/01/plotting-terminal-edges-of-tree.html
#https://stackoverflow.com/questions/54785864/exporting-a-high-res-node-labelled-phylogenetic-tree
#https://www.youtube.com/watch?v=6HLOdDtBFfA
#https://4va.github.io/biodatasci/r-ggtree.html

##########################################################################################


# NOW create for Andersen et al. 2019 species tree
tree<-readNexus('andersen_2019_backbone_tree.nex')
plot(tree)
tree_cladogram=ape::compute.brlen(tree)

absent<- c("taxon_1","taxon_2","taxon_3","taxon_4","taxon_5","taxon_9","taxon_10",
           "taxon_11", "taxon_12", "taxon_13","taxon_14")
present<-c("taxon_7","taxon_6","taxon_8")

tree_cladogram<-paintBranches(tree_cladogram,edge=sapply(present,match,tree_cladogram$tip.label),
                  state="present",anc.state="absent")

plot(tree_cladogram)

#clade.label<-LETTERS[1:14]
#clade.label<-rep(NA,14)
clade.label<-c('H','H','H',NA,'G','D','E','C','B','A',NA,'Myza',NA,NA)

# starting from clade H to base of tree including pardalotes and thornbills
N<-c(21,12,20,24,10,44,22,10,19,13,2,2,1,1)
tip.label<-tree_cladogram$tip.label

## set crown node depth to be consistent across all species and with the mtDNA tree!
depth<-rep(c(0.04545455),times=14)

#depth<-sapply(tip.label,function(x,y)
#  0.5*y$edge.length[which(tree_cladogram$edge[,2]==
#                            which(y$tip.label==x))],y=tree_cladogram)

trans<-data.frame(tip.label,clade.label,N,depth)
rownames(trans)<-NULL

## here's what trans looks like
print(trans)
obj<-phylo.toBackbone(tree_cladogram,trans)

## plot
plot(obj,col=c('black','black','black','black','red','red','red',
               'black','black','black','black','black','black','black'))

####
## Setting up the png output file
## you might want to change the width and height here!
svg(filename = "AndersenSpeciesTree_v2.svg",width = 10, height = 20)
#png(filename = "SpeciesTree_test.png",width = 800, height = 800,res=1000)
## The margin definition
#par(mar = c(bottom, left, top, right)
par(mar = c(1,1,1,1))

plot(obj,col=c('black','black','black','black','red','red','red',
               'black','black','black','black','black','black','black'))

## Saving the file
dev.off()

##########################################################################################

### NOW create for Andersen et al. 2019 mtDNA tree Fig S6 Bayesian Maximum Credibility Tree 
tree<-readNexus('mtDNA_andersen_2019_backbone_tree_v2.nex')
plot(tree)
tree_cladogram=ape::compute.brlen(tree)


absent<-c("taxon_5","taxon_6","taxon_7","taxon_8","taxon_9","taxon_10","taxon_11", "taxon_12","taxon_13","taxon_14")
present<- c("taxon_1","taxon_2","taxon_3","taxon_4")

tree_cladogram<-paintBranches(tree_cladogram,edge=sapply(present,match,tree_cladogram$tip.label),
                              state="present",anc.state="absent")

plot(tree_cladogram)

#clade.label<-LETTERS[1:14]
#clade.label<-rep(NA,14)
clade.label<-c('F','G','D','E','Spine','H','H','H','B','C','A','Myza',NA,NA)

# starting from clade F to base of tree including pardalotes and thornbills
N<-c(24,10,44,22,2,21,12,20,19,10,13,2,1,1)
tip.label<-tree_cladogram$tip.label

## set crown node depth to be consistent across all species and with the species tree!
depth<-rep(c(0.04545455),times=14)
#depth<-sapply(tip.label,function(x,y)
#  0.5*y$edge.length[which(tree_cladogram$edge[,2]==
#                            which(y$tip.label==x))],y=tree_cladogram)

trans<-data.frame(tip.label,clade.label,N,depth)
rownames(trans)<-NULL
## here's what trans looks like
print(trans)
obj<-phylo.toBackbone(tree_cladogram,trans)


## plot
plot(obj,col=c('red','red','red','red','black','black','black',
               'black','black','black','black','black','black','black'))

####
## Setting up the png output file
## you might want to change the width and height here!
svg(filename = "Andersen_mtDNATree_v2.svg",width = 10, height = 20)
#png(filename = "SpeciesTree_test.png",width = 800, height = 800,res=1000)
## The margin definition
#par(mar = c(bottom, left, top, right)
par(mar = c(1,1,1,1))

plot(obj,col=c('red','red','red','red','black','black','black',
               'black','black','black','black','black','black','black'))

## Saving the file
dev.off()



