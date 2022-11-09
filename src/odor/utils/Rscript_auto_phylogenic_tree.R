#### AUTOMATED PROTOCOL FOR LIGANDS POCKETS DESCRIPTION ####
#"C:\Program Files\R\R-3.6.3\bin\Rscript.exe"
## Collect arguments
args <- commandArgs(TRUE)
## Default setting when no arguments passed
if(length(args) < 1) {
  args <- c("--help")
}
## Help section
if("--help" %in% args) {
  cat("
      The R Script
      
      Arguments:
      [TO MODIFY]
      --or               - string of the or
      --output                - output path and name file
      --tree             - path to the tree structure
      
      Example:
      Rscript Rscript_auto_phylogenic_tree.R --or=\"Or5k1;Or4p20\" --tree=\"../Data/phylogenic_tree/data_HORDE/phylo_tree_PhyML_Or10ad1.tree\" --o=\".\" \n\n")
  q(save="no")
}

## LIBRARY
require(showtext)
library(ggtree)
library(ggrepel)
library(dplyr)
library(jsonlite)
library(ggimage)
library(ggtree)
library(ggtext)
library(stringi)
## PARSE ARGUMENTS
## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1
args = argsL

if(is.null(argsL$tree)) {
  args$tree = "../Data/phylogenic_tree/data_HORDE/phylo_tree_PhyML_Or10ad1.tree"
}
if(is.null(argsL$output)) {
  args$output = "./img.svg"
}
print("-/-")
print(args$output)
print(args)
fileName <- args$tree
file_txt = readChar(fileName, file.info(fileName)$size)
file_txt = gsub("[\r\n]", "", file_txt)
myTree <- ape::read.tree(text=file_txt)
options(ignore.negative.edge=TRUE)


list_or_family = stri_extract_first_regex(myTree$tip.label, "[0-9]+")
cls = NULL
for (i in sort(as.numeric(unique(list_or_family)))) {
  cls[[as.character(i)]] = myTree$tip.label[which(list_or_family == as.character(i))]
}

mTree <- groupOTU(myTree, cls)

vec_or = unlist(strsplit(args$or,";"))

image = ggtree(mTree,aes(color=group), layout="circular", size = 0.1, branch.length = 'none') + 
  geom_tiplab(data = ~ subset(., !(label %in% c(vec_or))),aes( label=label), 
              geom = "text", linetype='dashed', linesize=.1, size = 1)  +
  geom_label_repel(data = ~ subset(., label %in% c(vec_or)),
                   aes( label=label), 
                   box.padding = 1,
                   show.legend = FALSE, max.overlaps = 1000000000)
ggsave(file=args$output, plot=image, width=20, height=16)


