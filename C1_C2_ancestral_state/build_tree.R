library(phytools)
library(ape)
tree <- read.tree("gtdbtk.bac120.user_msa_selected.rmdup_root.nwk")
tree <- root(tree, "GCA_002174125.1_ASM217412v1_genomic", resolve.root = TRUE)
x <- read.csv("1251genomes_metadata.CSV", row.names = 1)
x <- as.matrix(x)[,1]
mtree<-make.simmap(tree,x,model="ER")

pd<-summary(mtree,plot=FALSE)
pd

cols <- read.csv("cols.csv",row.names =1)
cols <- as.matrix(cols)[,1]
plot(mtree,cols,type="fan",fsize=0.8,ftype="off")

target_group <- "C1"
group_file <- "C1vsC2_metadata.tab"
group_metadata <- read.table(group_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)

target_tips <- group_metadata[group_metadata[, 2] == target_group, 1]
target_tips <- intersect(target_tips, tree$tip.label)

length(target_tips)

target_mtree <- keep.tip(mtree, target_tips)
plot(target_mtree,cols,type="fan",fsize=0.8,ftype="off")
summary(target_mtree, plot = FALSE)
