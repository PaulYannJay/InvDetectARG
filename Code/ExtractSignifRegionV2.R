library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
start_time <- Sys.time()
args = commandArgs(trailingOnly=TRUE)
NodeSortFL2=read.table(args[1], header=T, stringsAsFactors = F) #Genetic map
NodeSortFL=read.table(args[2], header=T, stringsAsFactors = F) #Genetic map
# NodeSortFL2=read.table("~/Projects/InvDetectARG/Output/1KGP.chr1.dated.test1.NodeStat.summarisedSignif.txt", header=T, stringsAsFactors = F) #Genetic map
# NodeSortFL=read.table("~/Projects/InvDetectARG/Output/1KGP.chr1.dated.test1.NodeStat.summarisedSignif_IndivNode.txt", header=T, stringsAsFactors = F) #Genetic map
#  
 InterestSub=NodeSortFL2[NodeSortFL2$VarSpanP<0.05,]
 InterestingNode= NodeSortFL[(NodeSortFL$Left %in% InterestSub$Left & NodeSortFL$EmpP<0.05),]
 InterestingNode = InterestingNode %>% group_by(IDID, NodeMeanId, NodeHarmoMeanId) %>% mutate(Left=min(Left), Right=max(Right)) %>%
                      distinct(IDID, NodeMeanId, NodeHarmoMeanId, .keep_all=TRUE) # merge 
write.table(InterestSub, paste(args[1], ".top", sep=""), quote = F, col.names = T)
write.table(InterestingNode, paste(args[2], ".top", sep=""), quote = F, col.names = T)

