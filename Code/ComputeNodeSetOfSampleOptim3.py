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
#textfileNode=open(Prefix + ".NodeStat", "w") #open the output file

ts = tskit.load(TreeFile)
times=ts.nodes_time #np array of node age
samples = ts.num_samples #Number of samples
nodes = ts.num_nodes #Number of nodes

@nb.jit(nopython=True) #Use numba for fast calculation
#@nb.jit(nopython=True, parallel=True)
def GetSampleSum(parent): #Function to identify node using two statistic: the mean of the ID of the samples coalscing to the node and the harmonic mean of the same vector (with slight modification: see below). This function iterate through the tree and visit each node.
    sumNodes = np.zeros(len(parent)+1, dtype=np.int32) #Sum of the ID of the sample that coalesce to the nodes
    harmoSumNodes = np.zeros(len(parent)+1, dtype=np.float64) #Harmonic sum
    nSample = np.zeros(len(parent)+1, dtype=np.int32)#number of sample that coalesce to the nodes
    maxval=parent.max()#the ID of the oldest parent (the root)
    for u in range(samples):#for all samples. u is the ID of the sample
        v = u # v is the ID of the focal nodes 
        while v < maxval: #while the focal nodes is not the root
            sumNodes[parent[v]]=sumNodes[parent[v]]+u #Increment the sum of the parent of the focal node by the ID of the focal sample
            harmoSumNodes[parent[v]]=harmoSumNodes[parent[v]]+1/(u+1) #Same but with the harmonic value. Here, I add "+1" to avoid dividing by 0 for the first sample. So it is not exactly the harmonic mean but it has similar properties 
            nSample[parent[v]]=nSample[parent[v]]+1 #Same with the number of sample
            v=parent[v] #move to the parent of the focal node and loop
    sampleMean=sumNodes/nSample #Calculate the mean
    sampleMean=sampleMean[~np.isnan(sampleMean)] #the sample node (0 => node ID => samples) have nSample=0, so this give "NaN". remove them.
    sampleHarmoMean=nSample/harmoSumNodes #calculate the harmonic mean
    sampleHarmoMean=sampleHarmoMean[~np.isnan(sampleHarmoMean)] #Same
    return sampleMean, sampleHarmoMean
    
TreeList=ts.trees(sample_lists=True) #List of Tree
nbtree=ts.num_trees #Number of Tree

@nb.jit(nopython=True)
def getParentTime(parentNod): #Function to get the age of the focal node
    time=np.zeros(len(parentNod), dtype=np.float64)
    i=0
    for u in parentNod: #for all focal node
        if (u!=-1): #If the focal node is not the virtual root
            time[i]=times[u] #Grep its age in the node age array
        else:
            time[i]=-1 #else return "-1"
        i=i+1
    return time
        
def getParentArray(nodes): #Function to get the parent of the focal nodes #I do not remember why, but it seem to not work if I replace the "tree.parent(u)" by a code looking at the ts.parent_array (as done before; which should be more efficient as it allow using numba). Maybe try again latter.
    parent=np.zeros(len(nodes), dtype=np.int32)
    i=0
    for u in nodes:
        parentu=tree.parent(u)
        parent[i]=parentu
        i=i+1
    return parent[:-1]

#with open("testHarmoMean.txt", 'w') as f:
with open(Prefix + ".NodeStat", 'w') as f: #output file
    for tree in TreeList: #Main loop. Loop over all trees
        print("Tree:" + str(tree.index) + "/" + str(nbtree)) #print the output
        nodes=list(tree.nodes(order='timeasc')) #Get the node list in time ascending order (start with sample, end with root)
        #print(nodes)
        nodes.sort()#sort the nodes to git the nodes in sorting order
        #print(nodes)
        parent=getParentArray(nodes)#Get the array of the parent of each node (allow reconstucting the tree)
        uniq, uniqInd = np.unique(parent, return_inverse=True)#Get unique ID and sorted indices of parents. This is the main trick to save time. In the exemples from the tskit github page, kelleher et al use a function like "parent = np.zeros(ts.num_nodes, dtype=np.int32); for u in range(ts.num_nodes): parent[u] = tree.parent(u)". The issue with this is that it create and manipulate a huge array that it is mostly empty (each tree involves only a small subset of nodes). To manipulate only the nodes form the focal tree, I first obtain the parent of each nodes. Then I extract only uniq value, so I obtain unique node. Then, to store the relation between the nodes (parent-child) in an "UniqInd" indice array. For each node, the the uniqINd arrayu indicate what is its parent in the uniqArray. So, the UniqInd array contains a kind of new ID for each parent nodes (so excluding sample node). These new ID are comprise between 0 (the parent node with the smaller age) and "z", the parent node with the oldest age, which is the root. We will use these new ID to traverse trees. Importantly, we need to modify the new ID by adding to them the number of samples (the ID "5" become "15" if there are 10 different samples). This allow traversing function to traverse iteratively the parents and not loop over samples again and again (if there are 10 samples [0-9], the parent new ID must start at "1O" and not at 0 (IDs 0 to 9 are the sample nodes)
        time=getParentTime(uniq) #Get the time of each parent nodes
        sampleMean, sampleHarmoMean=GetSampleSum(uniqInd+samples) #traverse the tree and get the means
        np.savetxt(f, np.rot90([uniq,sampleMean,sampleHarmoMean, np.repeat(tree.interval[0], len(uniq)), np.repeat(tree.interval[1], len(uniq)), time]), delimiter=" ") #save output
        
print("Done ! Well done !")
