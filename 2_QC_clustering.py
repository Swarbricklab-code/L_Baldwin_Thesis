# Script for the QC and clustering of single cell RNA seq data
# By Louise Baldwin, based on scripts from Sergio Erdal Irac
# takes merged h5ad as input, merged (but now clustered) h5ad as output

###################
# Set up
###################

# import packages
import numpy as np
import pandas as pd
import scanpy as sc
import os
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rcParams


# directories
in_file = ("data/processed/merged_scrublet.h5ad")
results_file = ("data/processed/merged_withscrub.h5ad")
figdir = ("outs/2_QC_clustering_withscrub/figures/")
tabdir = ("outs/2_QC_clustering_withscrub/tables/")
os.makedirs(figdir, exist_ok=True)
os.makedirs(tabdir, exist_ok=True)
# set parameters for scanpy
# verbosity: errors (0), warnings (1), info (2), hints (3), detailed traceback (4)
# change default figdir to desired figdir
sc.settings.verbosity = 3             
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')
#sc.settings.figdir='/share/ScratchGeneral/loubal/projects/MSC/mouse-single-cell/outs/QC/figures/'
sc.settings.figdir=figdir

# read in file
adata=sc.read_h5ad(in_file)

#####################
# QC and clustering
#####################

sc.pl.highest_expr_genes(adata, n_top=20, save='.png')
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.var['mt'] = adata.var_names.str.startswith('mt-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True, save='.png')

sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', save='_counts_vs_mt.png')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts',save='_counts_vs_genes.png')

#use quantile to filter number of genes and counts
upper_lim = np.quantile(adata.obs.n_genes_by_counts.values, .98)
lower_lim = np.quantile(adata.obs.n_genes_by_counts.values, .02)
print(f'{lower_lim} to {upper_lim}')
#adata = adata[adata.obs.n_genes_by_counts < 7000, :] #example if you wanted to pick a number yourself
adata = adata[(adata.obs.n_genes_by_counts < upper_lim) & (adata.obs.n_genes_by_counts > lower_lim)]
adata = adata[adata.obs.pct_counts_mt < 15]

#normalise and find highly variable genes
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata, save=".png")

# save as filtered but unannotated
adata.raw = adata
#regress out mt, scale and cluster
adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
adata.write(results_file)
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.tl.leiden(adata)

#for troubleshooting
#adata=sc.read_h5ad(results_file)

#calculate and plot counts per cluster
# figure out how to change x and y axis labels to percentage and leiden cluster
counts = adata.obs.groupby(['Treatment', 'leiden']).count().reset_index()
def map_per(x):
    samp, y = x
    tot = counts[counts.Treatment == samp].total_counts.sum()
    return (y/tot)*100
counts['per'] = counts[['Treatment', 'total_counts']].apply(map_per, axis = 1)

counts
plt.figure(figsize = (20,3))

ax = sns.barplot(x = 'leiden', y = 'per', hue = 'Treatment', data = counts)
plt.savefig(("percentage_cells_per_leiden01_cluster.png"))

#############################
# Finding cluster markers
#############################
sc.settings.verbosity = 2
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, save='_DEG.png')

# calculate statistics for marker genes
## NOTE: This is automatically t-test overstim var. Should this be wilcoxon or logrank?
sc.tl.rank_genes_groups(adata, groupby='leiden', key_added='rank_genes_leiden')


pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(60)
results = adata.uns['rank_genes_groups']
results['names'].dtype.names

#create table with DEG statistics
out = np.array([[0,0,0,0,0]])
for group in results['names'].dtype.names:
    out = np.vstack((out, np.vstack((results['names'][group],
                                     results['scores'][group],
                                     results['pvals_adj'][group],
                                     results['logfoldchanges'][group],
                                     np.array([group] * len(results['names'][group])).astype('object'))).T))

# # print table to csv
# #### NOTE: check were these csv's end up
# #### hopefully this works - remove if not
pd.DataFrame(out).to_csv((tabdir+'_rank_genes_groups.csv'), index=False)

# #check cols and rows of output
# out.shape

# # filter DEG statistics
markers = pd.DataFrame(out[1:], columns = ['Gene', 'scores', 'pval_adj', 'logfc', 'cluster'])
markers = sc.get.rank_genes_groups_df(adata, None)
markers = markers[(markers.pvals_adj < 0.05) & (markers.logfoldchanges > .5)]
markers

markers.to_csv((tabdir+'_markers_rank_genes_groups.csv'), index=False)


#save results to file.
adata.write_h5ad(results_file)

#troubleshooting: why wont it make a umap
# adata=sc.read_h5ad(results_file)
# # plot some umaps
# rcParams['figure.figsize']=(7,7)

# # sc.pl.umap(adata, color=['leiden','ReactionID','Treatment','Tissue','Model'], legend_loc='on data', size=1,
# #  show=False, add_outline=True,  frameon=False, return_fig=True, ncols=2, wspace=0.5 title='leiden_0.1', save='umap.pdf')

# sc.pl.umap(adata, color=['leiden','ReactionID','Treatment','Tissue','Model'], size=1,
#  show=False, add_outline=True,  frameon=False, return_fig=True, ncols=2, wspace=0.5, title='leiden_0.1', save=('umap.png'))

# adata.write_h5ad(results_file)