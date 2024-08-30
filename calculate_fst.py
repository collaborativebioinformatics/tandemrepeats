import sys
import pandas as pd
from itertools import combinations
import argparse


def parse_args(args):
    parser = argparse.ArgumentParser(prog="bench", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input", type=str, help="Allele count tsv created by population_ac_by_length.py")
    parser.add_argument("-o", "--output", type=str, default='/dev/stdout',
                        help="Ouput TSV file containing fst and pairwise fst")
    return parser.parse_args()


# I am follwoing equations in the following paper to calcuate Fst.
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3141729/#APP1title
# we calculate the frequency of an allele give all other alleles per locus in pk
# weights are also calculated from the allele counts
# the function takes the allele count table and populations you are trying to find markers 
def calculate_fst(allele_ac, pops):
    all_FST = []
    for locus, df in allele_ac.groupby(['chrom', 'start', 'end']):
        
        Hs = 0
        Ht_1 = 0
        Ht_2 = 0
        
        number_of_haps = df[[f"AC_{p}" for p in pops]].sum()
        all_haps = sum(number_of_haps)
        #for pop in ['EAS', 'AMR', 'AFR', 'SAS']:
        for pop in pops:
            pk = df[f"AC_{pop}"] / number_of_haps[f"AC_{pop}"]
            qk = 1 - pk
            w = number_of_haps[f"AC_{pop}"] /all_haps
            Hs += 2 * w * pk * qk
            Ht_1 += w * pk 
            Ht_2 += w * qk 
        
        Ht = 2 * Ht_1 * Ht_2
        
        Fst = 1 - (Hs/Ht)
        all_FST.append(Fst)
      
    return pd.concat(all_FST)




if __name__ == '__main__':
    args = parse_args(sys.argv[1:])
    inputTable = pd.read_csv(args.input, sep='\t')
    populations = ['SAS', 'EAS', "AMR", 'AFR']
    inputTable['fst'] = calculate_fst(inputTable, populations)


    for pop1, pop2 in combinations(populations, 2):
        inputTable[f'fst_{pop1}_{pop2}'] = calculate_fst(inputTable,[pop1,pop2])

    inputTable.to_csv(args.output, sep="\t") 