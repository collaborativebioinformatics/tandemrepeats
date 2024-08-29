"""
Merge putatively redundant alleles per-locus
"""
import sys
import json
import argparse

import tdb
import truvari

import numpy as np
import pandas as pd


def build_matrix(alleles, threshold=0.98):
    """
    Compare all-against-all to build a boolean adjacency matrix
    """
    n_entries = len(alleles)
    match_matrix = np.zeros((n_entries, n_entries), dtype=bool)
    for i in range(n_entries - 1):
        for j in range(i + 1, n_entries):
            szsim, _ = truvari.sizesim(len(alleles[i]), len(alleles[j]))
            state = True
            if szsim < threshold:
                state = False
            else:
                sqsim = truvari.seqsim(alleles[i], alleles[j])
                if sqsim < threshold:
                    state = False
            match_matrix[i, j] = state
            match_matrix[j, i] = state
    return match_matrix


def find_matching_sets(matrix, locus_id):
    """
    Creates a lookup of which alleles match
    returns the new allele numbers lookup as dict and a list of original
    allele numbers to keep
    """
    n = len(matrix)
    visited = [False] * n
    matched_sets = []

    def dfs(item, current_set):
        """
        Depth first search to find chain of matches
        """
        visited[item] = True
        current_set.append(item)  # alt alleles start at number 1
        for other in range(n):
            if matrix[item][other] and not visited[other]:
                dfs(other, current_set)

    for i in range(n):
        if not visited[i]:
            current_set = []
            dfs(i, current_set)
            matched_sets.append(current_set)

    # This just keeps the first allele
    # A better strategy would be to keep the most frequently observed allele
    to_keep = [(locus_id, idx[0]) for idx in matched_sets]
    # Create a lookup of old allele number to the new allele number
    to_rename = {(locus_id, old_num): new_num
                 for new_num, entry_set in enumerate(matched_sets)
                 for old_num in entry_set}

    return to_rename, to_keep

def table_updater(table, all_to_rename, all_to_keep):
    """
    Given an allele or a sample table, subset to locus/allele that need to be kept and
    rename the remaining allele numbers.
    All updating happens in-place
    """
    table.set_index(['LocusID', 'allele_number'], inplace=True)
    keep = table.index.isin(all_to_keep)
    table.drop(table.index[~keep], inplace=True)
    table.index = table.index.map(lambda idx: (idx[0], all_to_rename.get(idx, idx[1])))
    table.reset_index(inplace=True)


def merge_main(args):
    """
    Create a new tdb from multiple input calls
    """
    parser = argparse.ArgumentParser(prog="tdb create", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-o", "--output", metavar="OUT", required=True,
                        help="Output tdb directory")
    parser.add_argument("--threshold", type=float, default=0.98,
                        help="Similarity threshold")
    parser.add_argument("--debug", action="store_true",
                        help="Verbose logging")
    parser.add_argument("input", metavar="IN",
                        help="Input tdb")
    args = parser.parse_args(args)

    d = tdb.load_tdb(args.input, lfilters=[('LocusID', '=', 9)])
    print(d)

    stats = {'n_loci': 0,
             'n_collap_loci': 0,
             'n_alleles': 0,
             'n_collap_alleles': 0}
    all_to_rename = {}
    all_to_keep = []
    for (grp,), alleles in d['allele'].groupby(['LocusID']):
        matrix = build_matrix(list(alleles['sequence']), args.threshold)
        to_rename, to_keep = find_matching_sets(matrix, grp)
        a1 = len(alleles)
        a2 = len(to_keep)
        stats['n_alleles'] += a1
        stats['n_loci'] += a2
        if a1 != a2:
            stats['n_collap_loci'] += 1
            stats['n_collap_alleles'] += a1 - a2
        all_to_rename.update(to_rename)
        all_to_keep.extend(to_keep)

    all_to_rename = pd.Series(all_to_rename)
    print('what')
    print(all_to_rename)
    print()
    table_updater(d['allele'], all_to_rename, all_to_keep)
    
    for samp in d['sample'].values():
        table_updater(samp, all_to_rename, all_to_keep)

    tdb.save_tdb(d, "out.tdb")
    print(d)
    print(json.dumps(stats, indent=4))

if __name__ == '__main__':
    merge_main(sys.argv[1:])

