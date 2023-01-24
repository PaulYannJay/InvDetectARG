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
if NSample=="All":
    ts_subset = ts.simplify(filter_sites=False)
else:
    subset = range(0, int(NSample))
    ts_subset = ts.simplify(subset, filter_sites=False)

with open(Positions) as file:
        lines = [line.strip() for line in file]

for position in lines:
    print(position)
    tree=ts_subset.at(int(position))
    tree.draw_svg(
            path=Prefix +  "_N=" + NSample + "Pos=" + str(position) + ".svg",
            size=(1200,800),
            y_axis=True, y_label=" ",  # optional: show a time scale on the left
            y_ticks=(1000,2000, 3000,4000,5000,6000,7000,8000,9000,10000,11000,12000,13000,14000,15000)
            )

