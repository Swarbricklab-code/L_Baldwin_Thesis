# The Great (Immune) Escape: Unpacking breast cancer immune evasion cell by cell
Repo for scripts associated with the above thesis by Louise Baldwin

## Project overview
The immune checkpoint inhibitors, or immunotherapies, have been a paradigm shift in cancer treatment. Blocking the immune-inhibitory interactions of CTLA4 with it's ligand B7, or PD1 with it's ligands PDL1/PDL2, can yield remarkable and durable anti-tumour immune reponses in some cancers. However, these therapies are yet to achieve the same outcomes in triple negative breast cancer (TNBC).

As part of my doctoral studies, I utilised single cell sequencing technologies to investigate mechanisms of TNBC immune evasion in two murine models. To analyse both the single cell RNA sequencing and matched single cell VDJ T cell receptor (TCR) sequencing, I generated a number of scripts. The scripts used for the analysis of this dataset can be found in this repo. 

## Computational overview
The analysis of this project is described below. Each script requires an input file and generates an output file.

| Script                               | In_file                             | Results_file                         
| :---                                 |    :----:                           |                                         ---:| 
| 1_merge_data.py                      | Cellranger outputs                  | updated_merged.h5ad                         | 
| 1a_Scrublet.ipynb                    | updated_merged.h5ad                 | merged_scrublet.h5ad                        | 
| 2_QC_clustering.py                   | merged_scrublet.h5ad                | merged_withscrub.h5ad                       | 
| 3_umap.ipynb                         | merged_withscrub.h5ad               | clustered_merged_umap.h5ad                  | 
| 4_major_level_annotation             | clustered_merged_umap.h5ad          | annotated.h5ad                              | 
| BBKNN_just_Tcells.ipynb              | annotated.h5ad                      | Subset_Tcells_BBKNN.h5ad                    | 
| BBKNN_just_Tcells_annotation.ipynb   | Subset_Tcells_BBKNN.h5ad            | Subset_Tcells_BBKNN_annotated.h5ad          | 
| Tcells_mergeVDJ.ipynb                | Subset_Tcells_BBKNN_annotated.h5ad  | Tcells_withTCR.h5ad                         | 
| TCR_anaysis.ipynb                    | Tcells_withTCR.h5ad                 | Plots, figures                              |

### Script descriptions
**1_merge_data.py:** concatenate cellranger outputs into one .adata object

**1a_Scrublet.ipynb:** predict and remove doublets from the data

**2_QC_clustering.py:** Run QC, filter and cluster data

**3_umap.ipynb:** Calculate UMAP

**4_major_level_annotation:** Annotate clusters at the major (lineage) level

**BBKNN_just_Tcells.ipynb:** Subset T cells and batch correct using BBKNN

**BBKNN_just_Tcells_annotation.ipynb:** Annotate T cells at the minor level

**Tcells_mergeVDJ.ipynb:** Append TCR data to T cell object

**TCR_anaysis.ipynb:** Carry out TCR QC and analysis

## Link to thesis
The use of these scripts is described further in the methods (section 2.11), and the outputs are found throughout chapter 5.
