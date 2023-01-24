import sys
import os
import cyvcf2
import tsinfer
import getopt
import json
from tqdm import tqdm
from tskit import MISSING_DATA

def main(argv):
    global vcfFile
    global Prefix 
    try:
        #opts, args = getopt.getopt(argv,"hg:o:w:k:t:k:",["ifile=","ofile="])
        opts, args = getopt.getopt(argv,"hi:o:p:")
    except getopt.GetoptError:
        print('Bug1 ! Usage: LoadVCF_InferARG.py -i <file.vcf> -o <OutputPrefix> ')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('Usage: LoadVCF_InferARG.py -i <file.vcf> -o <OutputPrefix>')
            sys.exit()
        elif opt in ("-i"):
            vcfFile = arg
        elif opt in ("-p"):
            InversionGeno = arg
        elif opt in ("-o"):
            Prefix = arg
    print('Input file is ', vcfFile)
    print('Output prefix is ', Prefix)


def add_diploid_sites(vcf, samples):
    """
    Read the sites in the vcf and add them to the samples object.
    """
    # You may want to change the following line, e.g. here we allow
    # "*" (a spanning deletion) to be a valid allele state
    #allele_chars = set("ATGCatgc*")
    allele_chars = set("ATGCatgc*")
    pos = 0
    for variant in vcf:  # Loop over variants, each assumed at a unique site
        if pos == variant.POS:
            print(f"Duplicate entries at position {pos}, ignoring all but the first")
            continue
        else:
            pos = variant.POS
        if any([not phased for _, _, phased in variant.genotypes]):
            raise ValueError("Unphased genotypes for variant at position", pos)
        alleles = [variant.REF.upper()] + [v.upper() for v in variant.ALT]
        #ancestral = variant.INFO.get("AA", ".")  # "." means unknown
        ## some VCFs (e.g. from 1000G) have many values in the AA field: take the 1st
        #ancestral = ancestral.split("|")[0].upper()
        #if ancestral == "." or ancestral == "":
        #    ancestral_allele = MISSING_DATA
        #    # alternatively, you could specify `ancestral = variant.REF.upper()`
        #else:
        #    ancestral_allele = alleles.index(ancestral)
        ## Check we have ATCG alleles
        #for a in alleles:
        #    if len(set(a) - allele_chars) > 0:
        #        print(f"Ignoring site at pos {pos}: allele {a} not in {allele_chars}")
        #        continue
        # Map original allele indexes to their indexes in the new alleles list.
        genotypes = [g for row in variant.genotypes for g in row[0:2]]
        #samples.add_site(pos, genotypes, alleles, ancestral_allele = 0) ###NOTE: this work only with the last version of tsinfer, the 0.3.0 SO I removed this section
        samples.add_site(pos, genotypes, alleles)


#def add_haploid_individuals(vcf, samples): #For every individuals, create two haploid individual with the metadata "Name" + "_Haplotype number". Also check if they have the focal inversion
#    for name in vcf.samples:
#        if name + "_1" in SampleWInv:
#            samples.add_individual(ploidy=1, metadata={"name": name + "_1", "Inverted" : "T"})
#        else:
#            samples.add_individual(ploidy=1, metadata={"name": name + "_1", "Inverted" : "F"})
#        if name + "_2" in SampleWInv:
#            samples.add_individual(ploidy=1, metadata={"name": name + "_2", "Inverted" : "T"})
#        else:
#            samples.add_individual(ploidy=1, metadata={"name": name + "_2", "Inverted" : "F"})


def chromosome_length(vcf):
    assert len(vcf.seqlens) == 1
    return vcf.seqlens[0]

if __name__ == "__main__":
        main(sys.argv[1:])

#with open(InversionGeno) as file:
#    SampleWInv = [line.strip() for line in file] #Store the list of haploid genome with inversion in a list

vcf = cyvcf2.VCF(vcfFile)

with tsinfer.SampleData(
    path=Prefix+".samples", 
    #sequence_length=chromosome_length(vcf)
) as samples:
    #add_haploid_individuals(vcf, samples)
    add_diploid_sites(vcf, samples)

print(
    "Sample file created for {} samples ".format(samples.num_samples)
    + "({} individuals) ".format(samples.num_individuals)
    + "with {} variable sites.".format(samples.num_sites),
    flush=True,
)

# Do the inference
ts = tsinfer.infer(samples)
print(
    "Inferred tree sequence: {} trees over {} Mb ({} edges)".format(
        ts.num_trees, ts.sequence_length / 1e6, ts.num_edges
    )
)

ts.dump("../Output/" + Prefix + ".trees")
