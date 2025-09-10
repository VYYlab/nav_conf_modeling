# nav_conf_modeling

## Introduction

This repository accompanies the manuscript *Exploring voltage-gated sodium channel conformations and protein-protein interactions using AlphaFold2*.  
The goal of the study is to evaluate the capability of AlphaFold2 to sample multiple conformational states of voltage-gated sodium (NaV) channels, and to assess the accuracy of AlphaFold Multimer in modeling interactions between the NaV α-subunit and its key protein partners, including auxiliary β-subunits and calmodulin (CaM).  

The work addresses five main research questions:  
1. Does AlphaFold2 generate a single state or multiple distinct conformations when modeling full NaV channels?  
2. What is the range and extent of conformational states that can be sampled?  
3. Which modeling settings (e.g., number of recycles, use of custom templates) influence the states obtained?  
4. Can AlphaFold Multimer reproduce experimentally observed NaV α–partner complexes?  
5. Does the presence of protein partners alter the conformational landscape of the NaV α-subunit?  

This repository provides a **step-by-step guide** to reproduce the main results of the paper. It contains all required input files, scripts, and instructions for running the modeling, calculating conformational metrics, and analyzing the resulting ensembles.

The workflow is divided into four stages:  
1. **Prepare inputs**  
2. **Define distances and regional subsets**  
3. **Generate models with ColabFold**  
4. **Analyze models (distances and pLDDT subsets)**  

---

## 1. Preparing the inputs

All twenty input files used in the study are provided in this repository under `inputs/`. These correspond to the twenty study cases described in the manuscript, which include NaV1.1, NaV1.2, NaV1.4, NaV1.5, NaV1.7, NaV1.9, and NaX α-subunits, as well as complexes with β-subunits (β1–β4) and calmodulin.  

![Regions of interest in NaV channels](figures/Figure1.png)

All sequences were retrieved from UniProt (see manuscript for accession codes). For multimeric complexes, the sequences of the partners are included in the input FASTA file on the same line, separated by a colon (`:`). This is the syntax required by AlphaFold/ColabFold to interpret the proteins as interacting partners.

---

## 2. Defining distances and subsets

In order to assign NaV channel models to specific conformational states, we analyzed structural features known to change during gating and inactivation. These include distances between specific residues and atoms located in flexible regions of the channel that are sensitive to conformational transitions.

**(Placeholder for figure)**  

The reference table of atom pairs is provided in `distances/atom_pairs.tsv`. This file defines which residues and atoms should be used for distance calculations. Each row specifies the pair and includes a `tag` field, which is used as the column name in the resulting dataset.

In addition, subsets of residues are defined for pLDDT analysis, corresponding to regions of interest (e.g., individual VSDs, the pore domain, the IFM motif). These subsets allow us to quantify the confidence of AlphaFold2 predictions in specific structural elements.

---

## 3. Generating the models

Models are generated using [ColabFold](https://github.com/sokrypton/ColabFold), specifically with the **subsampled MSA approach** described in the paper. Example command:

```bash
colabfold_batch --num-models 5 --model-type auto --msa-mode mmseqs2_uniref_env \
--num-seeds 20 \
--templates --max-seq 256 --max-extra-seq 512 \
--num-recycle 6 --save-recycles inputs/nav1.7_alphaonly.fasta outputs/nav1.7_alphaonly/ ```
