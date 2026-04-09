# Repository for: [Manuscript Title]

This repository contains the code, data pointers, and workflow used to reproduce the analyses and figures for:

Author A, Author B, Author C. Year. *Full manuscript title*. Journal. DOI

## Abstract

Metabolic niche partitioning links genome-encoded functions to community assembly and ecological roles, yet delineating these niches directly from genome-centric metagenomes remains limited by the lack of standardized workflows. Here, we present a semi-automated, genome-centric framework for defining metabolic niches in microbial consortia and demonstrate its application in a model sulfate-reducing, benzene-degrading enrichment culture.
We sequenced the solid and liquid phases of the culture and reconstructed metagenome-assembled genomes (MAGs) using MuDoGeR. To capture functional diversity, we integrated annotations from multiple biological knowledge domains, including genes involved in anaerobic aromatic degradation, KEGG modules and reactions, genome-scale metabolic reconstructions, and PROKKA, resulting in multiple matrix-style datasets. MAGs were grouped by shared functional traits using K-means and hierarchical clustering, followed by ANOVA and chi-square statistical filtering, depending on dataset type. Spearman correlation further refined trait associations and revealed functional co-occurrence patterns, which were visualized as undirected, weighted networks in Cytoscape. Network topology was then used to identify functional hubs and organize them into metabolic niches guided by biochemical pathway knowledge.
We recovered 127 MAGs (120 Bacteria, 7 Archaea) and dereplicated them into 52 species-level operational taxonomic units. We identified XX distinct functional network hubs, which were resolved into XX metabolic niches, including aromatic compound degradation and carbon fixation niches. While most MAGs were linked to multiple niches, suggesting functional redundancy, a few performed bottleneck roles such as acetogenesis.
This workflow provides a systematic, scalable approach to defining metabolic niches from metagenomic data and to illuminating microbial cooperation, niche overlap, and key functional contributors within complex microbial ecosystems.

## Overview

This repository provides an end-to-end workflow for defining metabolic niches from genome-centric metagenomic data. The workflow begins with raw reads or genome FASTA files and proceeds through genome recovery, functional annotation, dataset generation, clustering, statistical testing, correlation analysis, and network construction.

The framework is modular. Users may generate one or more functional potential datasets depending on their study design and available data. At least one dataset matrix is required for downstream clustering, statistics, and network analysis.

## Workflow summary

The workflow consists of three phases:

 # 1. Genome recovery and dataset generation

Genome FASTA files are obtained either from using the Multi-domain Genome Recovery tool- MuDoGeR (DOI: 10.1111/1755-0998.13904) or from any alternative MAG recovery tool. These genomes are then functionally annotated using one or more complementary approaches, including targeted pathway/gene detection, KEGG module profiling, genome-scale metabolic reconstruction, and whole-genome annotation. These steps produce one or more matrix-style datasets describing the functional potential of each genome.

2. Clustering, statistics, and correlation analysis

Each selected dataset is cleaned and evaluated for cluster number optimization. Genomes are then grouped using K-means + MDS and hierarchical clustering. Statistical testing is applied according to dataset type, with chi-square used for binary datasets and ANOVA used for proportional datasets. Feature-feature, feature-cluster, and cluster-cluster correlations are then calculated.

3. Network construction

Filtered correlation outputs are merged into edge lists and node attribute tables for network construction in Cytoscape. These networks are used to identify functional hubs and support metabolic niche interpretation.

Datasets generated in this workflow

Depending on the selected annotation strategy, the workflow can generate one or more of the following datasets:

Targeted gene/pathway presence-absence matrix
Targeted pathway completion/proportion matrix
KEGG module presence-absence matrix
KEGG module completion/proportion matrix
Genome reconstruction versus KEGG module matrix
GEM quality assessment table

Users may analyze one dataset or several datasets independently.

Inputs

The workflow accepts either:

raw metagenomic reads for users who wish to reconstruct MAGs, or
genome FASTA files for users starting directly from genomes

Additional optional inputs include:

pathway/gene definition files
KEGG module or reaction mappings
sample or genome metadata
user-defined configuration files

The minimum requirement for downstream analysis is at least one dataset matrix generated from the genomes.

Repository structure
# [Manuscript Title]

![Python](https://img.shields.io/badge/Python-3.11-blue)
![R](https://img.shields.io/badge/R-4.3.2-blue)
![Snakemake](https://img.shields.io/badge/Snakemake-8.0-brightgreen)
![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)

Repository for the code, workflow, and data pointers used to reproduce the analyses and figures for:

**Author A, Author B, Author C. Year.** *Full manuscript title*. Journal. DOI

---

## Table of contents

- [Abstract](#abstract)
- [Overview](#overview)
- [Workflow at a glance](#workflow-at-a-glance)
- [Datasets generated](#datasets-generated)
- [Repository structure](#repository-structure)
- [Inputs, intermediate files, and outputs](#inputs-intermediate-files-and-outputs)
- [How to run](#how-to-run)
- [Expected outputs](#expected-outputs)
- [Reproducibility](#reproducibility)
- [Data availability](#data-availability)
- [Software requirements](#software-requirements)
- [Setup](#setup)
- [License](#license)
- [Contact](#contact)

---

## Abstract

Metabolic niche partitioning links genome-encoded functions to community assembly and ecological roles, yet delineating these niches directly from genome-centric metagenomes remains limited by the lack of standardized workflows. Here, we present a semi-automated, genome-centric framework for defining metabolic niches in microbial consortia and demonstrate its application in a model sulfate-reducing, benzene-degrading enrichment culture.

We sequenced the solid and liquid phases of the culture and reconstructed metagenome-assembled genomes (MAGs) using MuDoGeR. To capture functional diversity, we integrated annotations from multiple biological knowledge domains, including genes involved in anaerobic aromatic degradation, KEGG modules and reactions, genome-scale metabolic reconstructions, and PROKKA, resulting in multiple matrix-style datasets. MAGs were grouped by shared functional traits using K-means and hierarchical clustering, followed by ANOVA and chi-square statistical filtering, depending on dataset type. Spearman correlation further refined trait associations and revealed functional co-occurrence patterns, which were visualized as undirected, weighted networks in Cytoscape. Network topology was then used to identify functional hubs and organize them into metabolic niches guided by biochemical pathway knowledge.

We recovered 127 MAGs (120 Bacteria, 7 Archaea) and dereplicated them into 52 species-level operational taxonomic units. We identified XX distinct functional network hubs, which were resolved into XX metabolic niches, including aromatic compound degradation and carbon fixation niches. While most MAGs were linked to multiple niches, suggesting functional redundancy, a few performed bottleneck roles such as acetogenesis.

This workflow provides a systematic, scalable approach to defining metabolic niches from metagenomic data and to illuminating microbial cooperation, niche overlap, and key functional contributors within complex microbial ecosystems.

---

## Overview

This repository contains an end-to-end workflow for defining metabolic niches from genome-centric metagenomic data.

The workflow begins with either:
- raw metagenomic reads, or
- genome FASTA files

It then proceeds through:
- genome recovery
- functional annotation
- dataset generation
- clustering
- statistical testing
- correlation analysis
- network construction

The framework is modular. Users may generate one or more dataset matrices depending on their goals. At least one matrix is required for downstream analysis.

---

## Workflow at a glance

The workflow has three phases.

### 1. Genome recovery and dataset generation
Genome FASTA files are obtained from MuDoGeR or another MAG recovery workflow. These genomes are then annotated using one or more approaches, including targeted pathway/gene detection, KEGG module profiling, genome-scale metabolic reconstruction, and whole-genome annotation.

### 2. Clustering, statistics, and correlation analysis
Each selected dataset matrix is cleaned, evaluated for cluster number optimization, clustered using K-means + MDS and hierarchical clustering, and analyzed using ANOVA or chi-square depending on data type. Feature-feature, feature-cluster, and cluster-cluster correlations are then computed.

### 3. Network construction
Filtered outputs are merged into edge lists and node attribute tables for network construction in Cytoscape.

---

## Datasets generated

Depending on the selected annotation strategy, the workflow can generate one or more of the following:

1. **Targeted gene/pathway presence-absence matrix**  
2. **Targeted pathway completion/proportion matrix**  
3. **KEGG module presence-absence matrix**  
4. **KEGG module completion/proportion matrix**  
5. **Genome reconstruction versus KEGG module matrix**  
6. **GEM quality assessment table**

These datasets can be analyzed independently or in parallel.

---

## Repository structure
.
├── README.md
├── config/
├── env/
├── inputs/
│   ├── raw_reads/
│   ├── genomes/
│   ├── metadata/
│   ├── pathway_definitions/
│   ├── kegg_definitions/
│   └── external_data_pointers/
├── scripts/
│   ├── 00_mag_recovery/
│   ├── 01_dataset_generation/
│   ├── 02_downstream_analysis/
│   └── 03_network_visualization/
├── workflow/
│   ├── manifests/
│   ├── notebooks/
│   └── pipeline_notes/
├── intermediate/
│   ├── 00_mag_recovery/
│   ├── 01_dataset_generation/
│   └── 02_downstream_analysis/
├── outputs/
│   ├── datasets/
│   ├── clustering_results/
│   ├── statistical_results/
│   ├── correlation_results/
│   ├── network_files/
│   └── summary_tables/
├── figures/
├── manuscript/
└── LICENSE

Folder organization
inputs/

Contains all starting materials required to run the workflow, including raw reads, genome FASTA files, metadata, curated pathway definitions, KEGG mappings, and pointers to large external datasets.

intermediate/

Contains stepwise outputs produced throughout the workflow, including MAG recovery results, annotation outputs, cleaned datasets, clustering outputs, statistical results, and correlation tables. These files are retained to support reproducibility and troubleshooting.

outputs/

Contains the main final products of the workflow, including processed dataset matrices, clustering summaries, statistical outputs, correlation results, edge lists, node tables, and network-ready files.

Downstream analysis scripts

The downstream analysis begins with at least one selected dataset matrix and proceeds through the following scripts:

01_removing_zero_columns_K_Optimization.py
Cleans the dataset and identifies candidate K values.
02_kmeans_mds_clustering.py
Performs K-means clustering and multidimensional scaling.
03_hierarchical_clustering.py
Performs hierarchical clustering.
04_ANOVA_testing.py
Used for proportion/completion datasets.
04_chi-square_testing.py
Used for binary datasets.
05_spearman_feature_vs_feature_correlation.py
Computes correlations among features.
06_spearman_feature_vs_cluster_correlation.py
Computes correlations between features and cluster groups.
07_spearman_cluster_vs_cluster_correlation.py
Computes correlations among clustering outputs.
08_merge_filtered_files.py
Merges filtered results into edge lists and node attribute tables for network construction.
Running the workflow

A typical execution order is:

Dataset generation
MAG recovery
targeted annotation
pathway completion calculation
KEGG module detection
KEGG module completion calculation
PROKKA annotation
filtering of hypothetical proteins
CarveMe reconstruction
memote quality assessment
Downstream analysis
cleaning and K optimization
K-means + MDS clustering
hierarchical clustering
ANOVA or chi-square testing, depending on dataset type
feature-feature correlation
feature-cluster correlation
cluster-cluster correlation
merge filtered files for network construction
Network visualization
import edge and node tables into Cytoscape for network construction and export
Expected outputs

Typical outputs include:

one or more functional/metabolic dataset matrices
cleaned matrices for clustering
K optimization summaries
K-means and hierarchical clustering outputs
ANOVA or chi-square result tables
feature-feature, feature-cluster, and cluster-cluster correlation tables
merged edge lists and node attribute tables
Cytoscape-ready network files
figures and summary tables used in the manuscript
Reproducibility

This repository reproduces:

genome-resolved functional/metabolic datasets derived from MAGs
clustering of genomes by shared functional potential
statistical filtering of traits across cluster groups
correlation-based network inputs for metabolic niche interpretation
manuscript figures and summary outputs
Data availability

Raw data are available from:

[repository / accession / DOI]

Because some inputs are too large or restricted, they are not stored directly in this repository. Instructions for obtaining these files are provided in inputs/.

Processed data required to reproduce the published analyses and figures are included in outputs/ where appropriate.

Software requirements

This workflow was tested with:

R 4.3.2
Python 3.11
Snakemake 8.0
BLAST 2.14

Additional dependencies and package versions are listed in:

env/environment.yml
renv.lock or equivalent lock files

External tools used in the full workflow may include:

MuDoGeR
standena
METABOLIC
PROKKA
CarveMe
memote
Cytoscape

Equivalent tools may be substituted if they generate compatible outputs.

Setup
Conda
conda env create -f env/environment.yml
conda activate manuscript-workflow

License

Code is released under the MIT License.
Data are released under CC-BY 4.0 unless otherwise stated.

Contact

For questions about the workflow, open a GitHub issue or contact:

Name – email@domain



