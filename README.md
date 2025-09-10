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


All sequences were retrieved from UniProt (see manuscript for accession codes). For multimeric complexes, the sequences of the partners are included in the input FASTA file on the same line, separated by a colon (`:`). This is the syntax required by AlphaFold/ColabFold to interpret the proteins as interacting partners.

---

## 2. Defining distances and subsets

In order to assign NaV channel models to specific conformational states, we analyzed structural regions known to change during Nav conformational cycle and that have already been described in experimental strucutres. These include distances between specific residues and atoms located in  regions of the channel that are sensitive to conformational transitions:

- **Voltage-sensing domains (VSDs I–IV):**  
  Distance between the Cα of the first gating charge in the S4 helix (GC1-S4) and the Cα of the hydrophobic constriction site residue in the S2 helix (HC-S2).  
  - Large values → activated (“up”) VSD state.  
  - Small values → deactivated (“down”) VSD state.

- **Activation gate (AG):**  
  Distances between opposing Cα atoms in the S6 helices that form the hydrophobic intracellular gate.  
  - AG1 = S6I–S6III pair  
  - AG2 = S6II–S6IV pair  
  These values indicate whether the pore is open or closed.

- **Fast inactivation (IFM motif):**  
  Distance between the Cα of the phenylalanine in the IFM motif and the Cα of the aspartate in its binding site in the pore domain.  
  - Low values → IFM motif bound, channel in fast-inactivated state.  

- **Selectivity filter (SF):**  
  Distance between the Cα of the lysine and aspartate residues of the DEKA motif.  
  - Variations reflect dilation or constriction of the selectivity filter, associated with slow inactivation.

To benchmark the models, the same distances were also measured in available experimental structures of human NaV channels and selected non-human NaV channels representing distinct conformational states.


![Regions of interest in NaV channels](figures/Figure1.png)

---

### Structure of the reference table

All distances and residue subsets used in this study are summarized in `references/distances-and-subsets.csv`. Each row corresponds to a modeling case (e.g., NaV1.1, NaV1.7) and provides both the atom pairs for distance calculations and the residue ranges for pLDDT analysis.

**File structure (columns):**

- **Modeling case**: Channel subtype (e.g., NaV1.1, NaV1.2, NaV1.7).  
- **VSDI, VSDII, VSDIII, VSDIV**: Atom pairs defining the GC1–HC distances for each VSD.  
- **AG1, AG2**: Atom pairs defining the activation gate distances (S6I–S6III and S6II–S6IV).  
- **IFM**: Atom pair between IFM phenylalanine and the aspartate in the pore binding site.  
- **SF**: Atom pair within the DEKA motif to monitor selectivity filter dilation.  
- **plddt full**: Residue ranges covering the entire α-subunit, used for overall pLDDT averaging.  
- **plddt VSDI–VSDIV**: Residue ranges defining each VSD for regional pLDDT analysis.  
- **plddt AG, plddt IFM, plddt SF**: Residue ranges for the activation gate, IFM motif, and selectivity filter.  
- **plddt partner**: Residue ranges for the partner protein (β-subunit or CaM) when present; `NA` if not applicable.  
- **List of reference pdbs**: PDB identifiers of experimental structures used as references for the corresponding case; `NA` if not available.

This table provides the central reference for both structural metrics (distances) and regional confidence analysis (pLDDT subsets) used throughout the study.
---

## 3. Generating the models

Models are generated using [local ColabFold](https://github.com/YoshitakaMo/localcolabfold), specifically with the **subsampled MSA approach** described in the paper [High-throughput prediction of protein conformational distributions with subsampled AlphaFold2
](https://www.nature.com/articles/s41467-024-46715-9). Example command:

```bash
colabfold_batch --num-models 5 --model-type auto --msa-mode mmseqs2_uniref_env \
--num-seeds 20 \
--templates --max-seq 256 --max-extra-seq 512 \
--num-recycle 6 --save-recycles inputs/nav1.7_alphaonly.fasta outputs/nav1.7_alphaonly/ ```
