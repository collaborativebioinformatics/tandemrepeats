import tdb
import pandas as pd
from pybedtools import BedTool

# The tdb to analyze
tdb_fn = "../hprc_105.tdb/"

# The results of tdb query len_poly_score
len_poly_fn = "../hprc_105_len_poly_score.txt"

# The results of Fst query by allele length
fst_fn = "../result_fst.tsv"

# The metadata with "Sample name" and "Superpopulation code"
meta_fn = "../igsr_samples.tsv"

lp = pd.read_csv("../hprc_105_len_poly_score.txt", sep='\t')
fst = pd.read_csv("../result_fst.tsv", sep='\t')

lp.set_index(['chrom', 'start', 'end'], inplace=True)
fst.set_index(['chrom', 'start', 'end'], inplace=True)

fst['len_poly_score'] = lp['len_poly_score']

# Add annotations from GENCODE GTF (Release 38)
# gtf file on https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz 
annotations = BedTool('/gencode.v38.annotation.gtf')

# Reset index
fst_reset = fst.reset_index()
regions_df = fst_reset[['chrom', 'start', 'end']]

# Create a BedTool object for regions
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

# Define a function to classify the regions
def classify_region(row):
    region = (row['chrom'], row['start'], row['end'])
    if region in exons_intervals:
        return 'exon'
    elif region in genes_intervals:
        return 'gene'
    else:
        return 'intergenic'

# Apply the function to the dataframe
regions_df['region_type'] = regions_df.apply(classify_region, axis=1)
regions_df.set_index(['chrom', 'start', 'end'], inplace=True)

# Add region_type info back to the fst
fst = fst.join(regions_df['region_type'], how='left')

# Candidate loci
keep_fs = fst['fst'] >= 0.20
print(keep_fs.sum(), 'fst alleles')

keep_lp = fst['len_poly_score'] >= 20
print(keep_lp.sum(), 'lpoly loci')

keep_gene = fst['region_type'] == "gene"
print(keep_gene.sum(), 'gene')

keep_all = keep_fs & keep_lp & keep_gene

candidate_alleles = fst[keep_all]
candidate_loci = candidate_alleles.reset_index()[['chrom', 'start', 'end']].drop_duplicates()

print(len(candidate_alleles), 'total alleles')
print(len(candidate_loci), 'total loci')

# PCA with 7 TRs of interest within genes
candidate_loci_ids = loci[loci.index.isin(idx_candidate_loci)]['LocusID']

data = tdb.load_tdb("/Users/english/code/references/hprc_105.tdb",
                    lfilters=[('LocusID', 'in', candidate_loci_ids)])
                    
alleles = data['allele'].set_index(['LocusID', 'allele_number'])[[]]

all_lengths = []
for i, sample in enumerate(data['sample']):
    m_samp = data['sample'][sample]
    m_samp = m_samp[m_samp['spanning_reads'] >= 10].set_index(['LocusID', 'allele_number'])
    alleles[sample] = alleles.index.isin(m_samp.index).astype(int)
    # de-fragment the frame every once in a while
    if i % 50 == 0:
        alleles = alleles.copy()

alleles = alleles.fillna(0)

keep = alleles.mean(axis=1) < 1000
len(keep), keep.sum()

pca = PCA(n_components=2)
values = alleles[keep].values.T
X_r = pca.fit(values).transform(values)

m_X = pd.DataFrame(X_r, columns=["PC1", "PC2"])
m_X['sample'] = alleles.columns

# Add Population information for the plot
meta = pd.read_csv("/Users/english/code/references/hprc_105.tdb/igsr_samples.tsv", sep='\t')
slookup = dict(zip(meta["Sample name"], meta["Superpopulation code"]))

m_X['Super Pop'] = m_X['sample'].map(slookup).fillna('UNK')

pc1, pc2 = pca.explained_variance_ratio_
pc1 = round(pc1*100, 1)
pc2 = round(pc2*100, 1)

order = sorted(m_X['Super Pop'].unique())
p = sb.scatterplot(data=m_X, x="PC1", y="PC2", hue="Super Pop",
                   hue_order=order)
plt.legend(loc='lower center', bbox_to_anchor=(0.5, 1.15), ncol=3)
_ = p.set(title=f"PCA of TR from {len(data['sample'])} Samples ({len(candidate_loci)} loci within a gene)",
         xlabel=f"PC1 ({pc1}%)", ylabel=f"PC2  ({pc2}%)")
plt.savefig("PCA_Loci_within_genes.png")
