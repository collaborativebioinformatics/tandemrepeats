"""
Annotates tdb loci with status of intersection with a gene / exon
"""
import re
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
    parser.add_argument("-c", "--chrom", type=str, help="Chromsome name to parse")
    parser.add_argument("-l", "--locusids", type=str, help="LocusIDs to parse, multiple ids via e.g. (123,528,529)")
    parser.add_argument("-o", "--output", type=str, default='/dev/stdout',
                        help="Ouput TSV file")
    return parser.parse_args()

def get_gene_names(regions, annotations):
    # Convert the result to a DataFrame
    regions_bed = BedTool.from_dataframe(regions)
    loci_intersects_with_genes = regions_bed.intersect(annotations, wa=True, wb=True)

    columns = columns = ['chrom', 'start', 'end', 'gene_chrom', 'unk', 'type', 'gene_start', 'gene_end', 'unk2', 'strand', 'unk3', 'gtf_anno']
    intersected_df = loci_intersects_with_genes.to_dataframe(names=columns)

    # Reset index to work with indiv. columns
    intersected_df = intersected_df.reset_index()
    def fetch_gene_name(anno):
        match = re.search(r'gene_name "([^"]+)"', anno)
        if match:
            return match.group(1)
        return None

    intersected_df['gene_name'] = intersected_df['gtf_anno'].apply(fetch_gene_name)
    view = intersected_df[['chrom','start', 'end', 'gene_name']].drop_duplicates()
    ret = view.groupby(['chrom', 'start', 'end'])['gene_name'].agg(set).apply(lambda x: ','.join(x))
    return ret

def gtf_annotate(loci, bedtool_gtf):
    """
    Returns a dataframe of LocusID, hits_gene, hits_exon
    """
    regions_df = loci[['chrom', 'start', 'end']]

    regions = BedTool.from_dataframe(regions_df)

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
    
    gene_names = get_gene_names(regions_df[['chrom', 'start', 'end']], annotations)
    test = regions_df.set_index(['chrom', 'start', 'end'])
    test['hits_gene'] = test.index.isin(genes_intervals)
    test['hits_exon'] = test.index.isin(exons_intervals)
    test['gene_name'] = gene_names

    loci.set_index(['chrom', 'start', 'end'], inplace=True)
    output = loci.join(test)
    return output

if __name__ == '__main__':
    args = parse_args(sys.argv[1:])
    names = tdb.get_tdb_filenames(args.tdb)
    annotations = BedTool(args.gtf)
    filters = []
    if args.chrom:
        filters.append(('chrom', 'in', args.chrom.split(',')))
    if args.locusids:
        filters.append(('LocusID', 'in', args.locusids.split(',')))

    loci = pq.read_table(names['locus'], filters=filters).to_pandas()
    output = gtf_annotate(loci, annotations)
    output.to_csv(args.output, sep='\t', index=False)
