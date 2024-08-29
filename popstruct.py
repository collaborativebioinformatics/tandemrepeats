import tdb
import numpy as np
import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt
import scipy as scipy
from matplotlib import cm, colors

from collections import defaultdict

data = tdb.load_tdb("repo_utils/test_files/databases/hprc_105.tdb/",
                    lfilters=[("chrom", "=", "chr4")]) # subsetting to run faster

metadata = (pd.read_csv("repo_utils/test_files/databases/igsr_samples.tsv", sep='\t')
                .where(lambda x: x["Sample name"].isin(data["sample"].keys()))
                .dropna()
                .set_index(["Sample name"]))

allele_lengths = data["allele"].set_index(["LocusID", "allele_number"])["allele_length"]

rows = []
pops = []
for samp, table in data["sample"].items():
    if samp not in metadata.index:
        continue
    pops.extend([metadata.loc[samp]["Superpopulation code"]])# * 2)
    view = table.set_index(["LocusID", "allele_number"])
    view["allele_length"] = allele_lengths
    # 1::2 - get column separation - because of hets...? or non-ref I should say.. but lose ::2's almost perfect 
    # row separation
    # I think I should be taking the more common allele.
    rows.append(view.iloc[::2]["allele_length"].reset_index(level=1, drop=True))
    #just take the first allele - or every other, whatever
    #rows.append(view.iloc[1::2]["allele_length"].reset_index(level=1, drop=True))

# Subset to Loci with >= 10 alleles
loci_ac = data["allele"]["LocusID"].value_counts().where(lambda x: x >= 20).dropna()

heatmap = np.empty((len(rows), len(loci_ac)))
heatmap[:] = np.nan
heatmap.shape

allen_idx_lookup = dict(zip(loci_ac.index, range(len(loci_ac))))

for samp_idx, r in enumerate(rows):
    for locus, observed_len in r.items():
        if locus not in allen_idx_lookup: continue
        heatmap[samp_idx, allen_idx_lookup[locus]] = observed_len

#keep = heatmap.ptp(0) >= 20

sub_heat = heatmap#[:, keep] # only keep the most variable by length sites

center = np.nanmean(sub_heat, axis=0)
shifted_hm = (sub_heat - center) / sub_heat.ptp(0)
shifted_hm = shifted_hm[:, sub_heat.mean(axis=0).argsort()] # sort by loci's average length
# remove unobserved loci, I guess?
# m_mins = m_mins[~np.isnan(m_mins)]

set(pops)
lut = dict(zip(set(pops), ['darksalmon', 'palegreen', 'deepskyblue', 'violet']))
row_colors = [lut[_] for _ in pops]

row_colors_allele = 'rb' * (len(pops) // 2) #['darksalmon', 'palegreen', 'deepskyblue', 'violet']))
#row_colors = [lut[_] for _ in pops]

clustermap = sb.clustermap(np.nan_to_num(shifted_hm), 
              col_cluster=False, 
              row_colors=row_colors, 
              method="complete",
              metric="correlation", # cityblock, braycurtis - i liked how this looked
              cmap=cm.RdBu_r)

clustermap.savefig("clustermap.png")

allele_freq = tdb.allele_count_length(data)

allele_freq = allele_freq[allele_freq["AF"] >= 0.02]

al = (allele_lengths.reset_index()
          .sort_values(["LocusID", "allele_number", "allele_length"])
          .drop_duplicates(["LocusID", "allele_length"])
          .set_index(["LocusID", "allele_length"]))
af = allele_freq.reset_index().set_index(["LocusID", "allele_length"])

af['allele_number'] = al
af = af.reset_index().set_index(["LocusID", "allele_number"])["AF"]

rows = []
pops = []
for samp, table in data["sample"].items():
    if samp not in metadata.index:
        continue
    pops.extend([metadata.loc[samp]["Superpopulation code"]])# * 2)
    view = table.set_index(["LocusID", "allele_number"])
    view["allele_freq"] = af
    view = view[~view["allele_freq"].isna()]
    # I think I should be taking the more common allele.
    rows.append(view["allele_freq"].reset_index().groupby(["LocusID"])["allele_freq"].max())
    #just take the first allele - or every other, whatever
    #rows.append(view.iloc[1::2]["allele_length"].reset_index(level=1, drop=True))

loci_ac = af.reset_index()["LocusID"].value_counts().where(lambda x: (x >= 5)).dropna()


heatmap2 = np.empty((len(rows), len(loci_ac)))
heatmap2[:] = np.nan

allen_idx_lookup = dict(zip(loci_ac.index, range(len(loci_ac))))
for samp_idx, r in enumerate(rows):
    for locus, observed_len_af in r.items():
        if locus not in allen_idx_lookup: continue
        heatmap2[samp_idx, allen_idx_lookup[locus]] = observed_len_af

g = sb.clustermap(np.nan_to_num(heatmap2), 
                  col_cluster=False, 
                  row_colors=row_colors, 
                  method="complete",
                  metric="correlation", # cityblock, braycurtis - i liked how this looked
                  cmap=cm.RdBu_r)
g.ax_col_dendrogram.set_visible(False)

g.savefig("allelefq_clustermap.png")

df = np.nan_to_num(heatmap2)

#FloatingPointError: NaN dissimilarity value with metric=correlation
g = sb.clustermap(df, method='complete', metric="canberra")

g.savefig("allelefq_clustermap2.png")

#den = scipy.cluster.hierarchy.dendrogram(g.dendrogram_row.linkage,
den = scipy.cluster.hierarchy.dendrogram(g.dendrogram_col.linkage,
                                         labels = loci_ac.index,
                                         color_threshold=0.60)  

class Clusters(dict):
    def _repr_html_(self):
        html = '<table style="border: 0;">'
        for c in self:
            hx = rgb2hex(colorConverter.to_rgb(c))
            html += '<tr style="border: 0;">' \
            '<td style="background-color: {0}; ' \
                       'border: 0;">' \
            '<code style="background-color: {0};">'.format(hx)
            html += c + '</code></td>'
            html += '<td style="border: 0"><code>' 
            html += repr(self[c]) + '</code>'
            html += '</td></tr>'
        
        html += '</table>'
        
        return html

def get_cluster_classes(den, label='ivl'):
    cluster_idxs = defaultdict(list)
    for c, pi in zip(den['color_list'], den['icoord']):
        for leg in pi[1:3]:
            i = (leg - 5.0) / 10.0
            if abs(i - int(i)) < 1e-5:
                cluster_idxs[c].append(int(i))

    cluster_classes = Clusters()
    for c, l in cluster_idxs.items():
        i_l = [den[label][i] for i in l]
        cluster_classes[c] = i_l

    return cluster_classes

clusters = get_cluster_classes(den)

