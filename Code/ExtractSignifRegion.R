library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
start_time <- Sys.time()
args = commandArgs(trailingOnly=TRUE)
NodeSortFL2=read.table(args[1], header=T, stringsAsFactors = F) #Genetic map
NodeSortFL=read.table(args[2], header=T, stringsAsFactors = F) #Genetic map

 
 InterestSub=NodeSortFL2[NodeSortFL2$VarSpanP<0.05,]
 IntestingNode= NodeSortFL[(NodeSortFL$Left %in% InterestSub$Left & NodeSortFL$EmpP<0.05),]
 IntestingNode = IntestingNode %>% group_by(SampleSet,IDID, NSample) %>% mutate(FirstPos=min(Left), LastPos=max(Right)) %>%
                      distinct(SampleSet,IDID, NSample, .keep_all=TRUE) # merge 
write.table(InterestSub, paste(args[1], ".top", sep=""), quote = F, col.names = T)
write.table(IntestingNode, paste(args[2], ".top", sep=""), quote = F, col.names = T)

