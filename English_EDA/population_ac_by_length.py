"""
Calculates allele counts by length overall and for each subpopulation
Expects a tsv metadata file which has columns of:
    Sample name - identical sample name as found in the tdb
    Superpopulation code - identifier for sample's population
Samples without an entry in the Sample_name will be binned into an unknown ('UNK') population
"""
import sys
import argparse
import tdb
import pandas as pd


def parse_args(args):
    parser = argparse.ArgumentParser(prog="pop_ac", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("tdb", type=str, help="TDB file")
    parser.add_argument("metadata", type=str, help="TSV file")
    parser.add_argument("-c", "--chrom", type=str, help="Chromsome name to parse")
    parser.add_argument("-l", "--locusids", type=str, help="LocusIDs to parse, multiple ids via e.g. (123,528,529)")
    parser.add_argument("-s", "--samples", type=str, help="Samples to process, multiple via commas, default all")
    parser.add_argument("-o", "--output", type=str, default='/dev/stdout',
                        help="Ouput TSV file")
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args(sys.argv[1:])

    db_fn = args.tdb
    meta_fn = args.metadata

    filters = []
    if args.chrom:
        filters.append(('chrom', 'in', args.chrom.split(',')))
    if args.locusids:
        filters.append(('LocusID', 'in', args.locusids.split(',')))
    samples = None
    if args.samples:
        samples = args.samples.split(',')

    print("Loading metadata")
    available_samples = set(tdb.get_tdb_samplenames(db_fn))
    if samples:
        available_samples = samples.intersection(available_samples)

    metadata = (pd.read_csv(meta_fn, sep='\t')
                    .where(lambda x: x["Sample name"].isin(available_samples))
                    .dropna())

    view = metadata[['Sample name', 'Superpopulation code']].copy()

    pop_lookup = view.set_index(['Sample name'])['Superpopulation code'].to_dict()
    n_unk = 0
    for i in available_samples:
        if i not in pop_lookup:
            n_unk += 1
            pop_lookup[i] = 'UNK'

    print("Analyzing %d of %d samples" % (len(pop_lookup), len(available_samples)))

    print("Ancestry Distribution")
    print(view['Superpopulation code'].value_counts())
    print(f"UNK\t{n_unk}")

    print("Loading TDB")
    data = tdb.load_tdb(db_fn, samples=samples,
                        lfilters=filters)

    #LocusID allele_length   chrom   start   end     is_ref  AC      AF
    print("Getting AC by length")
    ac_table = tdb.allele_count_length(data)
    ac_table.set_index(['LocusID', 'allele_length'], inplace=True)

    # Columns for each population's observed allele count
    for i in view['Superpopulation code'].unique():
        ac_table[f'AC_{i}'] = 0
    if n_unk:
        ac_table[f'AC_UNK'] = 0

    indexed_allele = data['allele'].set_index(['LocusID', 'allele_number'])

    for sample_name, sample_table in data['sample'].items():
        m_popname = pop_lookup[sample_name]
        print(f"Analyzing {m_popname} sample {sample_name}")
        m_table = sample_table.set_index(['LocusID', 'allele_number']).join(indexed_allele, how='left')
        m_table.reset_index(inplace=True)
        m_table.set_index(['LocusID', 'allele_length'], inplace=True)

        ac_table.loc[m_table.index, f'AC_{m_popname}'] += 1
        
    ac_table.to_csv(args.output, sep='\t', index=False)
