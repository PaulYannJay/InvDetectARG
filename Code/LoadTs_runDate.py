import sys
import os
import cyvcf2
import tsinfer
import getopt
import tskit
import tsdate
from tqdm import tqdm

def main(argv):
    global TreeFile
    global Prefix 
    try:
        #opts, args = getopt.getopt(argv,"hg:o:w:k:t:k:",["ifile=","ofile="])
        opts, args = getopt.getopt(argv,"hi:o:")
    except getopt.GetoptError:
        print('Bug1 ! Usage: LoadTs_OutEdgeStat.py -i <file.trees> -o <OutputPrefix>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('Usage: LoadTs_OutEdgeStat.py -i <file.trees> -o <OutputPrefix>')
            sys.exit()
        elif opt in ("-i"):
            TreeFile = arg
        elif opt in ("-o"):
            Prefix = arg
    print('Input file is ', TreeFile)
    print('Output prefix is ', Prefix)

main(sys.argv[1:])

ts = tskit.load(TreeFile)
tsSimple=ts.simplify(keep_unary=False)
dated_ts = tsdate.date(tsSimple, Ne=12500, mutation_rate=2.5e-8)
dated_ts.dump("../Output/" + Prefix + ".dated.trees")
