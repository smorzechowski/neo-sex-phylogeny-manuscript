# Root reference trees with ape and save as PDF
# Sophie MacRae Orzechowski
# February 2024

library(ape)
library(phylotools)
setwd("~/PhD research/RaXML")

trees <- read.tree(file="OldZ.trees")
trees

trees[[1]]$tip.label
#From the list of tip labels, copy the name of the outgroup species

trees <- root(trees,"Chry",resolve.root=TRUE)
is.rooted(trees)



#If you want to see what your rooted gene trees look like (and to visually confirm that your trees have been rooted with the correct outgroup)

pdf("OldZ_rooted_trees.pdf")
for(i in 1:length(trees)){
  plot(trees[[i]])
}
dev.off()

ultra_trees <- trees
pdf("OldZ_rooted_ultra_trees.pdf")
for(i in 1:length(trees)){
  ultra_tree=trees[[i]]
  ultra_tree$edge.length<-NULL
  plot(ultra_tree)
}
dev.off()

#finally, let’s write out the rooted trees file

write.tree(trees,file="OldZ_rooted.trees")
