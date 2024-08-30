"""
Annotates tdb loci with status of intersection with a gene / exon
"""
import os
import sys
import argparse

import tdb
import pandas as pd
import pyarrow.parquet as pq
from pybedtools import BedTool

def parse_args(args):
    parser = argparse.ArgumentParser(prog="bench", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("tdb", type=str, help="TDB file")
    parser.add_argument("gtf", type=str, help="GTF file")
    parser.add_argument("-o", "--output", type=str, default='/dev/stdout',
                        help="Ouput TSV file")
    return parser.parse_args()

# abstract tdb.gene_annotate(tdb[lids], gtf)
if __name__ == '__main__':
    args = parse_args(sys.argv[1:])
    names = tdb.get_tdb_filenames(args.tdb)
    annotations = BedTool(args.gtf)

    # Hard coded filtering- turn off for production
    loci = pq.read_table(names['locus'], filters=[('chrom', '=', 'chr20')]).to_pandas()
    regions_df = loci[['chrom', 'start', 'end']]

    query_regions = regions = BedTool.from_dataframe(regions_df)

    # Extract gene features
    genes = annotations.filter(lambda x: x[2] == 'gene')

    # Extract exon features
    exons = annotations.filter(lambda x: x[2] == 'exon')

    # Intersect with genes
    intersects_with_genes = regions.intersect(genes, wa=True, wb=True)
    genes_intervals = set((line.chrom, int(line.start), int(line.end)) for line in intersects_with_genes)

    # Intersect with exons
    intersects_with_exons = regions.intersect(exons, wa=True, wb=True)
    exons_intervals = set((line.chrom, int(line.start), int(line.end)) for line in intersects_with_exons)

    test = regions_df.set_index(['chrom', 'start', 'end'])
    test['hits_gene'] = test.index.isin(genes_intervals)
    test['hits_exon'] = test.index.isin(exons_intervals)

    loci.set_index(['chrom', 'start', 'end'], inplace=True)
    output = loci.join(test)
    output.to_csv(args.output, sep='\t', index=False)

