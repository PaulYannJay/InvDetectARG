library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
start_time <- Sys.time()
args = commandArgs(trailingOnly=TRUE)
HumanRecomb=read.table(args[1], header=T, stringsAsFactors = F) #Genetic map

Recomb_bins=20 #Number of recombination bins
HumanRecomb=HumanRecomb %>% mutate(MyQuantileBins = cut(COMBINED_rate.cM.Mb., 
                                                        breaks = unique(quantile(COMBINED_rate.cM.Mb.,probs=seq.int(0,1, by=1/Recomb_bins))), 
                                                        labels= quantile(COMBINED_rate.cM.Mb.,probs=seq.int((1/Recomb_bins)/2,1-(1/Recomb_bins)/2, by=1/Recomb_bins)), # To get the media of each recombination interval
                                                        include.lowest=TRUE)) #Assign each window to a given bin of recombination rate

HumanRecomb=HumanRecomb %>% mutate(Start=lag(position)  %>% coalesce(0)) ## Genetic map indicate the recombination rate for fragment defined by the focal position and the previous position. We grep here the previous position

Node=read.table(args[2], stringsAsFactors = F) #Node table
colnames(Node)=c("Left", "Right", "ID", "Age","BranchLenght","SampleSet", "NSample")
end_time <- Sys.time()
TimeSpend=end_time - start_time
 print(paste0("Reading Node files= ",TimeSpend," mins"))

Node$IndexPos=findInterval(Node$Left, HumanRecomb$Start) #For each node, find the recombination interval in which it falls

Node$RecomInt=HumanRecomb$MyQuantileBins[Node$IndexPos] #For each node, find the recombination rates of the interval in which it falls


Node$ParentAge=Node$Age+Node$BranchLenght #Define the parent Age
Node$Size=Node$Right-Node$Left #Size of the interval span by the tree where this node was found

Node=Node %>% arrange(SampleSet, Left) %>% mutate(Follow=ifelse((lag(Right)==Left & lag(NSample)==NSample & lag(SampleSet)==SampleSet), 0, 1) %>% coalesce(0)) #Indicate with a 1 when we switch to another node, to another sample set or to a non-consecutive position
Node$RecomInt=as.numeric(as.character(Node$RecomInt))
Time_bins=100
Node=Node %>%  mutate(TimeBin = cut(Age,
                                                breaks = unique(quantile(Age,probs=seq.int(0,1, by=1/Time_bins))),
                                                labels = quantile(Age,probs=seq.int((1/Time_bins)/2,1-(1/Time_bins)/2, by=1/Time_bins)),
                                                include.lowest=TRUE)) #Define the time bin of each nodes
Node$TimeBin=as.numeric(as.character(Node$TimeBin))
end_time <- Sys.time()
TimeSpend=end_time - start_time
 print(paste0("Binning Node Time= ",TimeSpend," mins"))

 NodeSortFL =Node %>% group_by(SampleSet, NSample) %>% mutate(IDID=cumsum(Follow) + 1) %>% ##Create a group ID based on sample set and position. If a node is associated with the sample sample set than previously but with a non-consecutive position, it creare a new group
   group_by(SampleSet,IDID, NSample) %>% #Group by this ID
   mutate(SumSpan=sum(Size), 
             FirstPos=min(Left), LastPos=max(Right), 
             ParentAge=mean(ParentAge), Age=mean(Age),
             RecombRate=median(RecomInt),
             TimeBin=median(TimeBin)) %>%   #Sum the position
   mutate(TimeSurface=SumSpan*Age) #%>% 
   #group_by(Left, Right) %>% summarise(VarSpan=var(SumSpan),SumTimeSurface=sum(TimeSurface), VarTimeSurface=var(TimeSurface), )

 NodeSortFL=NodeSortFL %>% group_by(TimeBin, RecombRate) %>%  mutate (EmpP=dense_rank(desc(SumSpan))/length(unique(SumSpan))) %>% ungroup()
 NodeSortFL2=NodeSortFL %>% group_by(Left, Right) %>% summarise(VarSpan=var(SumSpan),
                                                                SumTimeSurface=sum(TimeSurface),
                                                                VarTimeSurface=var(TimeSurface),
                                                                MeanP=mean(EmpP),
                                                                RecomInt=mean(RecomInt), #Useless mean calculation, all value are the same. This is just an ugly trick
                                                                TopQuantP=quantile(EmpP,0.1))
 
 NodeSortFL2=NodeSortFL2 %>% group_by(RecomInt) %>% mutate (VarSpanP=dense_rank(desc(VarSpan))/length(unique(VarSpan)), #For each summary variable, calculate its "empirical P"
                                                            n=n(),
                                                            SumTimeSurfaceP=dense_rank(desc(SumTimeSurface))/length(unique(SumTimeSurface)),
                                                            VarTimeSurfaceP=dense_rank(desc(VarTimeSurface))/length(unique(VarTimeSurface)),
                                                            ) %>% ungroup()

write.table(NodeSortFL2, paste(args[2], ".summarisedSignif.txt", sep=""), quote = F, col.names = T)
write.table(NodeSortFL, paste(args[2], ".summarisedSignif_IndivNode.txt", sep=""), quote = F, col.names = T)
end_time <- Sys.time()

TimeSpend=end_time - start_time
 print(paste0("TotalTime= ",TimeSpend," mins"))
