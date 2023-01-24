import numpy as np
import sys
import os
import cyvcf2
import tsinfer
import getopt
import tskit
from tqdm import tqdm

def main(argv):
    global TreeFile
    global NSample 
    global Positions 
    global Prefix 
    NSample="All"
    try:
        #opts, args = getopt.getopt(argv,"hg:o:w:k:t:k:",["ifile=","ofile="])
        opts, args = getopt.getopt(argv,"hi:o:n:p:")
    except getopt.GetoptError:
        print('Bug1 ! Usage: ComputeTrees.py -i <file.trees> -o <OutputPrefix> -n NumberOfSample -p <FileWithPositions>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('Usage: LoadTs_OutEdgeStat.py -i <file.trees> -o <OutputPrefix>')
            sys.exit()
        elif opt in ("-p"):
            Positions = arg
        elif opt in ("-n"):
            NSample = arg
        elif opt in ("-i"):
            TreeFile = arg
        elif opt in ("-o"):
            Prefix = arg
    print('Input file is ', TreeFile)
    print('Output prefix is ', Prefix)

main(sys.argv[1:])

ts = tskit.load(TreeFile)

print("Tree sequence: Loaded !")
if NSample=="All":
    ts_subset = ts
else:
    subset = range(0, int(NSample))
    ts_subset = ts.simplify(subset, filter_sites=False)

print("Tree sequence: simplified !")
#Node_time=ts_subset.nodes_times #Need new version of tskit...
ID=0
NodeTim={}
for value in ts_subset.tables.nodes.time:
    NodeTim[ID] = value
    ID=ID+1

print("Node data: compiled !")

nbtree=ts_subset.num_trees
TreeList=ts_subset.trees(sample_lists=True) #Iterator of trees

TreeRoot=set([])
for tree in TreeList: #Store all roots in a set
    TreeNode=tree.roots
    for i in TreeNode:
        TreeRoot.add(i)

OldNodeNotRoot = {key:value for (key, value) in NodeTim.items() if (value >= 1000 and key not in TreeRoot)} #Only keep node that are pretty old and not roots of trees

textfileNode=open(Prefix + ".NodeStat", "w")

TreeList=ts_subset.trees(sample_lists=True) #When iterating over *.trees(), it clear the list. So we reload it
TreeID=1
for tree in TreeList:
    TreeInterval=tree.interval
    Left=str(TreeInterval[0])
    Right=str(TreeInterval[1])
    TreeNode=tree.nodes()
    for node in TreeNode:
        if node in OldNodeNotRoot:
            BranchLength=tree.branch_length(node)
            NodeId=sum(tree.samples(node))/len(list(tree.samples(node)))
            Line=Left + "\t" + Right + "\t" + str(node) + "\t" + str(OldNodeNotRoot[node]) + "\t" + str(BranchLength) + "\t" +  str(NodeId) + "\t" + str(len(list(tree.samples(node)))) + "\n"
            textfileNode.write(Line) 
    print("Tree:" + str(TreeID) + "/" + str(nbtree))
    TreeID=TreeID+1
    #print(TreeNode)

print("Finished !")
