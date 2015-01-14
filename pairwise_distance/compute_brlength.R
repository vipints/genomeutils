##
## Compute the branch length for a Phylogenetic tree
## Usage:
## cat ${DIR}/compute_brlength.R | $R_PATH --slave --args $infile.tre $outfile.tre 

## install package 
#install.packages(ape)

### get arguments 1: INFILE, 2: OUTFILE 
args <- commandArgs()

INFILE<-args[4]
OUTFILE<-args[5]

## load the library 
library(ape)

##Load the phylo tree in tree format
myTree<-read.tree(INFILE)

#print the uploaded tree 
myTree

## Calculate the branch length
myTreeBL<-compute.brlen(myTree, method = "Grafen", power = 1)

## print the branch length included tree 
myTreeBL

## save the file in tree format
write.tree(myTreeBL, file=OUTFILE)

## Done
