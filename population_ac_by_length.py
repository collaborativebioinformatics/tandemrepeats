import tdb
import pandas as pd

db_fn = "/Users/english/code/references/hprc_105.tdb"
meta_fn = "/Users/english/code/references/hprc_105.tdb/igsr_samples.tsv"

print("Loading metadata")
available_samples = tdb.get_tdb_samplenames(db_fn)
metadata = (pd.read_csv(meta_fn, sep='\t')
                .where(lambda x: x["Sample name"].isin(available_samples))
                .dropna())

view = metadata[['Sample name', 'Superpopulation code']].copy()
pop_lookup = view.set_index(['Sample name'])['Superpopulation code'].to_dict()

print("Analyzing %d of %d samples" % (len(pop_lookup), len(available_samples)))

print("Ancestry Distribution")
print(view['Superpopulation code'].value_counts())

print("Loading TDB")
data = tdb.load_tdb(db_fn,
                    samples=list(view['Sample name']),
                    lfilters=[("chrom", "=", "chr20")])

#LocusID allele_length   chrom   start   end     is_ref  AC      AF
print("Getting AC by length")
ac_table = tdb.allele_count_length(data)
ac_table.set_index(['LocusID', 'allele_length'], inplace=True)

# Columns for each population's observed allele count
for i in view['Superpopulation code'].unique():
    ac_table[f'AC_{i}'] = 0

indexed_allele = data['allele'].set_index(['LocusID', 'allele_number'])

for sample_name, sample_table in data['sample'].items():
    m_popname = pop_lookup[sample_name]
    print(f"Analyzing {m_popname} sample {sample_name}")
    m_table = sample_table.set_index(['LocusID', 'allele_number']).join(indexed_allele, how='left')
    m_table.reset_index(inplace=True)
    m_table.set_index(['LocusID', 'allele_length'], inplace=True)

    ac_table.loc[m_table.index, f'AC_{m_popname}'] += 1
    

ac_table.to_csv("result.txt", sep='\t', index=False)
    
    

