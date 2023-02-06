# Script for the QC, merging and clustering of single cell RNA seq data
# Written by Louise Baldwin, based on scripts from Sergio Erdal Irac

##############################
# Set up #
##############################

## data location
# https://github.com/Swarbricklab-datasets/mouse_immunotherapy.git

# import packages
import numpy as np
import pandas as pd
import scanpy as sc
import os

# set param for scanpy
# verbosity: errors (0), warnings (1), info (2), hints (3), detailed traceback (4)
sc.settings.verbosity = 3             
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

# make directory structure
cellranger = ("data/raw/single-cell/")
results_file = ("data/processed/updated_merged.h5ad")
tabdir = ("outs/1_merge_reactions/tables/")
os.makedirs(tabdir, exist_ok=True)


############################
#  Data import and merge #
############################

# Read in data

# samples = ["Ms_i_4T1_LN01", "Ms_i_4T1_LN02", "Ms_i_4T1_LN03", "Ms_i_4T1_LN04", "Ms_i_4T1_LN05",
#             "Ms_i_4T1_T01", "Ms_i_4T1_T02", "Ms_i_4T1_T03", "Ms_i_4T1_T04", "Ms_i_4T1_T05",
#             "Ms_i_EMT6_LN01", "Ms_i_EMT6_LN02", "Ms_i_EMT6_LN03", "Ms_i_EMT6_LN04", "Ms_i_EMT6_LN05",
#             "Ms_i_EMT6_T01", "Ms_i_EMT6_T02A", "Ms_i_EMT6_T02B", "Ms_i_EMT6_T03_1", "Ms_i_EMT6_T03_2", "Ms_i_EMT6_T04_1", "Ms_i_EMT6_T04_2", "Ms_i_EMT6_T05"]

# # for sample in samples:
# #     print(sample)
# for sample in samples:
#     sample=sc.read(cellranger+sample,
#     var_names="gene_symbols",
#     cache=True)

adata1 = sc.read_10x_mtx(cellranger + "Ms_i_4T1_LN01", # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)                              # write a cache file for faster subsequent reading

adata2 = sc.read_10x_mtx(cellranger + "Ms_i_4T1_LN02", # the directory with the `.mtx` file
        var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)   

adata3 = sc.read_10x_mtx(cellranger + "Ms_i_4T1_LN03", # the directory with the `.mtx` file
        var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)   

adata4 = sc.read_10x_mtx(cellranger +"Ms_i_4T1_LN04",  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)   

adata5 = sc.read_10x_mtx(cellranger + "Ms_i_4T1_LN05",  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)   

adata6 = sc.read_10x_mtx(cellranger + "Ms_i_4T1_T01",  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)   

adata7 = sc.read_10x_mtx(cellranger + "Ms_i_4T1_T02",  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)   

adata8 = sc.read_10x_mtx(cellranger + "Ms_i_4T1_T03",  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)   

adata9 = sc.read_10x_mtx(cellranger + "Ms_i_4T1_T04",  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)   

adata10 = sc.read_10x_mtx(cellranger + "Ms_i_4T1_T05",  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)   

adata11 = sc.read_10x_mtx(cellranger + "Ms_i_EMT6_LN01",  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)   

adata12 = sc.read_10x_mtx(cellranger + "Ms_i_EMT6_LN02",  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)   

adata13 = sc.read_10x_mtx(cellranger + "Ms_i_EMT6_LN03",  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)   

adata14 = sc.read_10x_mtx(cellranger + "Ms_i_EMT6_LN04",  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)   

adata15 = sc.read_10x_mtx(cellranger + "Ms_i_EMT6_LN05",  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True) 

adata16 = sc.read_10x_mtx(cellranger + "Ms_i_EMT6_T01",  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)   

adata17 = sc.read_10x_mtx(cellranger +"Ms_i_EMT6_T02A",  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)

adata18 = sc.read_10x_mtx(cellranger + "Ms_i_EMT6_T02B",  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)

adata19 = sc.read_10x_mtx(cellranger + "Ms_i_EMT6_T03_1",  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)   

adata20 = sc.read_10x_mtx(cellranger + "Ms_i_EMT6_T03_2",  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)  

adata21 = sc.read_10x_mtx(cellranger + "Ms_i_EMT6_T04_1",  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)   
    
adata22 = sc.read_10x_mtx(cellranger + "Ms_i_EMT6_T04_2",  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)   

adata23 = sc.read_10x_mtx(cellranger + "Ms_i_EMT6_T05",  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)   

# create merged object with metadata
# add ReactionID
adata1.obs['ReactionID'] = 'Ms_i_4T1_LN01'
adata2.obs['ReactionID'] = 'Ms_i_4T1_LN02'
adata3.obs['ReactionID'] = 'Ms_i_4T1_LN03'
adata4.obs['ReactionID'] = 'Ms_i_4T1_LN04'
adata5.obs['ReactionID'] = 'Ms_i_4T1_LN05'
adata6.obs['ReactionID'] = 'Ms_i_4T1_T01'
adata7.obs['ReactionID'] = 'Ms_i_4T1_T02'
adata8.obs['ReactionID'] = 'Ms_i_4T1_T03'
adata9.obs['ReactionID'] = 'Ms_i_4T1_T04'
adata10.obs['ReactionID'] = 'Ms_i_4T1_T05'
adata11.obs['ReactionID'] = 'Ms_i_EMT6_LN01'
adata12.obs['ReactionID'] = 'Ms_i_EMT6_LN02'
adata13.obs['ReactionID'] = 'Ms_i_EMT6_LN03'
adata14.obs['ReactionID'] = 'Ms_i_EMT6_LN04'
adata15.obs['ReactionID'] = 'Ms_i_EMT6_LN05'
adata16.obs['ReactionID'] = 'Ms_i_EMT6_T01'
adata17.obs['ReactionID'] = 'Ms_i_EMT6_T02A'
adata18.obs['ReactionID'] = 'Ms_i_EMT6_T02B'
adata19.obs['ReactionID'] = 'Ms_i_EMT6_T03_1'
adata20.obs['ReactionID'] = 'Ms_i_EMT6_T03_2'
adata21.obs['ReactionID'] = 'Ms_i_EMT6_T04_1'
adata22.obs['ReactionID'] = 'Ms_i_EMT6_T04_2'
adata23.obs['ReactionID'] = 'Ms_i_EMT6_T05'

#add treatment arm
adata1.obs['Treatment'] = 'Baseline'
adata2.obs['Treatment'] = 'PD1+CTLA4'
adata3.obs['Treatment'] = 'CTLA4'
adata4.obs['Treatment'] = 'PD1'
adata5.obs['Treatment'] = 'Control IgG'
adata6.obs['Treatment'] = 'Baseline'
adata7.obs['Treatment'] = 'PD1+CTLA4'
adata8.obs['Treatment'] = 'PD1'
adata9.obs['Treatment'] = 'CTLA4'
adata10.obs['Treatment'] = 'Control IgG'
adata11.obs['Treatment'] = 'Baseline'
adata12.obs['Treatment'] = 'PD1+CTLA4'
adata13.obs['Treatment'] = 'PD1'
adata14.obs['Treatment'] = 'CTLA4'
adata15.obs['Treatment'] = 'Control IgG'
adata16.obs['Treatment'] = 'Baseline'
adata17.obs['Treatment'] = 'PD1+CTLA4'
adata18.obs['Treatment'] = 'PD1+CTLA4'
adata19.obs['Treatment'] = 'CTLA4'
adata20.obs['Treatment'] = 'CTLA4'
adata21.obs['Treatment'] = 'PD1'
adata22.obs['Treatment'] = 'PD1'
adata23.obs['Treatment'] = 'Control IgG'

#add tissue
adata1.obs['Tissue'] = 'Lymph node'
adata2.obs['Tissue'] = 'Lymph node'
adata3.obs['Tissue'] = 'Lymph node'
adata4.obs['Tissue'] = 'Lymph node'
adata5.obs['Tissue'] = 'Lymph node'
adata6.obs['Tissue'] = 'Primary tumour'
adata7.obs['Tissue'] = 'Primary tumour'
adata8.obs['Tissue'] = 'Primary tumour'
adata9.obs['Tissue'] = 'Primary tumour'
adata10.obs['Tissue'] = 'Primary tumour'
adata11.obs['Tissue'] = 'Lymph node'
adata12.obs['Tissue'] = 'Lymph node'
adata13.obs['Tissue'] = 'Lymph node'
adata14.obs['Tissue'] = 'Lymph node'
adata15.obs['Tissue'] = 'Lymph node'
adata16.obs['Tissue'] = 'Primary tumour'
adata17.obs['Tissue'] = 'Primary tumour'
adata18.obs['Tissue'] = 'Primary tumour'
adata19.obs['Tissue'] = 'Primary tumour'
adata20.obs['Tissue'] = 'Primary tumour'
adata21.obs['Tissue'] = 'Primary tumour'
adata22.obs['Tissue'] = 'Primary tumour'
adata23.obs['Tissue'] = 'Primary tumour'

#add model
adata1.obs['Model'] = '4T1'
adata2.obs['Model'] = '4T1'
adata3.obs['Model'] = '4T1'
adata4.obs['Model'] = '4T1'
adata5.obs['Model'] = '4T1'
adata6.obs['Model'] = '4T1'
adata7.obs['Model'] = '4T1'
adata8.obs['Model'] = '4T1'
adata9.obs['Model'] = '4T1'
adata10.obs['Model'] = '4T1'
adata11.obs['Model'] = 'EMT6'
adata12.obs['Model'] = 'EMT6'
adata13.obs['Model'] = 'EMT6'
adata14.obs['Model'] = 'EMT6'
adata15.obs['Model'] = 'EMT6'
adata16.obs['Model'] = 'EMT6'
adata17.obs['Model'] = 'EMT6'
adata18.obs['Model'] = 'EMT6'
adata19.obs['Model'] = 'EMT6'
adata20.obs['Model'] = 'EMT6'
adata21.obs['Model'] = 'EMT6'
adata22.obs['Model'] = 'EMT6'
adata23.obs['Model'] = 'EMT6'

#add StudyID
adata1.obs['StudyID'] = 'LA24'
adata2.obs['StudyID'] = 'LA24'
adata3.obs['StudyID'] = 'LA24'
adata4.obs['StudyID'] = 'LA24'
adata5.obs['StudyID'] = 'LA24'
adata6.obs['StudyID'] = 'LA24'
adata7.obs['StudyID'] = 'LA24'
adata8.obs['StudyID'] = 'LA24'
adata9.obs['StudyID'] = 'LA24'
adata10.obs['StudyID'] = 'LA24'
adata11.obs['StudyID'] = 'LA25'
adata12.obs['StudyID'] = 'LA25'
adata13.obs['StudyID'] = 'LA25'
adata14.obs['StudyID'] = 'LA25'
adata15.obs['StudyID'] = 'LA25'
adata16.obs['StudyID'] = 'LA25'
adata17.obs['StudyID'] = 'LA25'
adata18.obs['StudyID'] = 'LA25'
adata19.obs['StudyID'] = 'LA25'
adata20.obs['StudyID'] = 'LA25'
adata21.obs['StudyID'] = 'LA25'
adata22.obs['StudyID'] = 'LA25'
adata23.obs['StudyID'] = 'LA25'

#make variable names unique
#skipping this step will over-normalise the reactions
adata1.var_names_make_unique()
adata2.var_names_make_unique()
adata3.var_names_make_unique()
adata4.var_names_make_unique()
adata5.var_names_make_unique()

adata6.var_names_make_unique()
adata7.var_names_make_unique()
adata8.var_names_make_unique()
adata9.var_names_make_unique()
adata10.var_names_make_unique()
adata11.var_names_make_unique()
adata12.var_names_make_unique()

adata13.var_names_make_unique()
adata14.var_names_make_unique()
adata15.var_names_make_unique()
adata16.var_names_make_unique()
adata17.var_names_make_unique()
adata18.var_names_make_unique()
adata19.var_names_make_unique()
adata20.var_names_make_unique()
adata21.var_names_make_unique()
adata22.var_names_make_unique()
adata23.var_names_make_unique()

# merge into one adata using concatenate
adata = adata1.concatenate(adata2, adata3, adata4, adata5, adata6, adata7, adata8, adata9, adata10, adata11,
    adata12, adata13, adata14, adata15, adata16, adata17, adata18, adata19, adata20, adata21,
     adata22, adata23)

# save merged adata
adata.write(results_file)

# print some metrics
print(adata.obs['Treatment'].value_counts())
print(adata.obs['Model'].value_counts())
print(adata.obs['Tissue'].value_counts())
print(adata.obs['ReactionID'].value_counts())

# generate cells per variable and write to df

#test_stats = df['A'].value_counts().to_frame()
treatment_stats = adata.obs['Treatment'].value_counts().to_frame()
model_stats = adata.obs['Model'].value_counts().to_frame()
tissue_stats = adata.obs['Tissue'].value_counts().to_frame()
reaction_stats = adata.obs['ReactionID'].value_counts().to_frame()

# save df as csv
#test_stats.to_csv('test.csv')
##Figure out a way to make these save somewhere other than logs##
treatment_stats.to_csv((tabdir+'cells_per_treatment.csv'))
model_stats.to_csv(tabdir+'cells_per_model.csv')
tissue_stats.to_csv(tabdir+'cells_per_tissue.csv')
reaction_stats.to_csv(tabdir+'cells_per_reaction.csv')

# make variable names unique in adata
adata.var_names_make_unique()


#####################################
# QC and Filtering #
#####################################


