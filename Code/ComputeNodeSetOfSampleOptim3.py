import numba as nb
import numpy as np
import sys
import os
import cyvcf2
import tsinfer
import getopt
import tskit
import statistics
import timeit
from tqdm import tqdm

def main(argv):
    global TreeFile
    global Prefix 
    NSample="All"
    try:
        opts, args = getopt.getopt(argv,"hi:o:")
    except getopt.GetoptError:
        print('Bug1 ! Usage: ComputeNodeSetOfSample.py -i <file.trees> -o <OutputPrefix>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('Usage: ComputeNodeSetOfSample.py -i <file.trees> -o <OutputPrefix>')
            sys.exit()
        elif opt in ("-i"):
            TreeFile = arg
        elif opt in ("-o"):
            Prefix = arg
    print('Input file is ', TreeFile)
    print('Output prefix is ', Prefix)

main(sys.argv[1:])
textfileNode=open(Prefix + ".NodeStat", "w") #open the output file

ts = tskit.load(TreeFile)
#samples = np.zeros(ts.num_samples, dtype=np.int32)
samples = ts.num_samples
nodes = ts.num_nodes
#nodes = np.zeros((ts.num_trees,ts.num_nodes), dtype=np.int32)
#for u in range(ts.num_samples):
#    samples[u] = u
print("sisi")
@nb.jit(nopython=True)
#@nb.jit(nopython=True, parallel=True)
def GetSampleSum(parent):
    sumNodes = np.zeros(len(parent)+1, dtype=np.int32)
    harmoSumNodes = np.zeros(len(parent)+1, dtype=np.float64)
    nSample = np.zeros(len(parent)+1, dtype=np.int32)
    for u in range(samples):
       # sumNodes[u]=u
        v = u
        while v < len(parent):
            sumNodes[parent[v]]=sumNodes[parent[v]]+u
            harmoSumNodes[parent[v]]=harmoSumNodes[parent[v]]+1/(u+1)
            nSample[parent[v]]=nSample[parent[v]]+1
            v=parent[v]
    sampleMean=sumNodes/nSample
    sampleMean=sampleMean[~np.isnan(sampleMean)]
    sampleHarmoMean=nSample/harmoSumNodes
    sampleHarmoMean=sampleHarmoMean[~np.isnan(sampleHarmoMean)]
    return sampleMean, sampleHarmoMean
    
print("sisi")
#def GetSampleSum(index):
#    for u in samples:
#        sumNodes[index, u]=u
#        v = u
#        while v != tskit.NULL:
#            sumNodes[index, parent[v]]=sumNodes[index, parent[v]]+u
#            v=parent[v]
        
#tree=ts.first()
TreeList=ts.trees(sample_lists=True) #When iterating over *.trees(), it clear the list. So we reload it
#nbtree=ts_sub.num_trees
#TreeID=1
with open("testHarmoMean.txt", 'a') as f:
    for tree in TreeList:
        #parent = np.zeros(ts.num_nodes, dtype=np.int32)
        parent=tree.parent_array
        parent=np.delete(parent, np.argwhere(parent==-1))
        uniq, uniqInd = np.unique(parent, return_inverse=True)
        starttime = timeit.default_timer()
        sampleMean, sampleHarmoMean=GetSampleSum(uniqInd+samples)
        np.savetxt(f, np.rot90([uniq,sampleMean,sampleHarmoMean, np.repeat(tree.interval[0], len(uniq)), np.repeat(tree.interval[1], len(uniq))]), delimiter=" ")
        print("The time difference is :", timeit.default_timer() - starttime)
        #sampleMean, sampleHormoMean=GetSampleSum(uniqInd)
        #print(sampleMean)
        #print(uniq)
        #print(uniqInd+samples)
        #print(uniqInd)
        #print(len(uniqInd)-samples)
       # print(sumNodes)
        print(tree.index)
       # print(nSample[tree.index, 7214])
       # print(parent)
       # #kprint(nSample)

#print("Tree sequence: Loaded !")
#nodes.tofile("test.txt", sep=" ")
#sampleMean=sumNodes/nSample
#sampleHarmoMean=nSample/harmoSumNodes
#np.savetxt("test.txt", sumNodes.astype(int), fmt='%i', delimiter=" ")
#np.savetxt("testSample.txt", nSample.astype(int), fmt='%i', delimiter=" ")
#np.savetxt("testParent.txt", parent.astype(int), fmt='%i', delimiter=" ")
#np.savetxt("testMean.txt", sampleMean, delimiter=" ")
#np.savetxt("testHarmoMean.txt", sampleHarmoMean, delimiter=" ")
#np.savetxt("testHarmo.txt", harmoSumNodes, delimiter=" ")
#ID=0
#ID=0
#NodeTim={}
#for value in ts.tables.nodes.time: #Store in a dictionary the age of nodes.
#    NodeTim[ID] = value
#    ID=ID+1
#
#print("Node data: compiled !")
#
#TreeList=ts.trees(sample_lists=True) #Iterator of trees
#TreeRoot=set([])
#for tree in TreeList: #Store all roots in a set
#    TreeNode=tree.roots
#    for i in TreeNode:
#        TreeRoot.add(i)
#
#OldNodeNotRoot = {key:value for (key, value) in NodeTim.items() if (value >= 100 and key not in TreeRoot)} #Only keep nodes that are pretty old (e.g. older than 1000 generation here, to reduce computation time) and not roots of trees
#NodeList=list(OldNodeNotRoot.keys())
#ts_sub=ts.subset(NodeList, remove_unreferenced=False)
#print(ts_sub.num_nodes)
#
##TreeList=ts_sub.trees(sample_lists=True) #When iterating over *.trees(), it clear the list. So we reload it
#TreeList=ts.trees(sample_lists=True) #When iterating over *.trees(), it clear the list. So we reload it
#nbtree=ts_sub.num_trees
#TreeID=1
#for tree in TreeList:
#    TreeInterval=tree.interval
#    Left=str(TreeInterval[0])
#    print(Left)
#    Right=str(TreeInterval[1])
#    TreeNode=tree.nodes()
#    print(tree.num_edges)
#    for node in TreeNode:
#        #if node in OldNodeNotRoot:
#            NodeIdMean=statistics.fmean(tree.samples(node)) #Mean of sample numerical ID
#            print(tree.parent(node))
#            print(NodeIdMean)
#            NodeIdHarmoMean=statistics.harmonic_mean(tree.samples(node)) #Mean of sample numerical ID
#            Line=Left + "\t" + Right + "\t" + str(node) + "\t" + str(OldNodeNotRoot[node]) + "\t" +  str(NodeIdMean) + "\t" +  str(NodeIdHarmoMean) + "\t" + str(len(list(tree.samples(node)))) + "\n" #Output Line
#            textfileNode.write(Line)  #Write output
#    print("Tree:" + str(TreeID) + "/" + str(nbtree))
#    TreeID=TreeID+1
#
#print("Done !")
