library(phytools)
library(ape)
tree <- read.tree("gtdbtk.bac120.user_msa.rmdup.nwk")
tree <- root(tree, "GCA_002174125.1_ASM217412v1_genomic", resolve.root = TRUE)
x <- read.csv("Acinetobactin.tab", row.names = 1) # aadA.tab
x <- as.matrix(x)[,1]
mtree<-make.simmap(tree,x,model="ER")
pd<-summary(mtree,plot=FALSE)
pd
cols <- read.csv("cols.csv",row.names =1)
cols <- as.matrix(cols)[,1]
plot(mtree,cols,type="fan",fsize=0.8,ftype="off")