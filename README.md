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
--num-recycle 6 --save-recycles inputs/hNaV1.7-alpha-only.fa outputs/hNaV1.7-alpha-only/ ```
```
**Explanation of arguments:**

- `--num-models 5`  
  Number of AlphaFold2 models to generate per random seed.

- `--model-type auto`  
  Automatically selects the best model configuration depending on the input (monomer or multimer).

- `--msa-mode mmseqs2_uniref_env`  
  Specifies the database and method used to generate multiple sequence alignments (MSAs).

- `--num-seeds 20`  
  Number of random seeds for initialization. Each seed explores different stochastic trajectories of the network.

- `--templates`  
  Enables the use of template structures if available in the database or provided by the user.

- `--max-seq 256`  
  Maximum number of sequences randomly subsampled from the full MSA.  
  **Key parameter for the subsampled MSA approach.**

- `--max-extra-seq 512`  
  Maximum number of additional sequences selected around each cluster center.  
  **Together with `--max-seq`, this parameter controls the degree of subsampling.**

- `--num-recycle 6`  
  Number of recycling iterations (re-feeding the structure back into the network to refine predictions).

- `--save-recycles`  
  Saves the intermediate model at each recycle step (from recycle 0 to recycle 6).

- `inputs/nav1.7_alphaonly.fasta`  
  Input FASTA file containing the sequence(s) to model.

- `outputs/nav1.7_alphaonly/`  
  Output directory where all generated models and metadata will be stored.

---

>**Important note on model counts:**  
>This run uses **20 seeds × 5 models = 100 models per case**.  
>Since intermediate models from recycle 0 through 6 are saved (`--save-recycles`), each model produces **7 structures**.  
>Thus, the total output is **100 × 7 = 700 models** stored in the output directory defined by the user.
>This is the number of models we generated in all cases for our study. 

### 3.1 Running with custom templates

In addition to the default template mode (where AlphaFold automatically selects templates from a curated version of the PDB), we also tested the effect of **state-specific custom templates**. As described in the manuscript, this analysis was performed for **NaV1.7**, focusing on VSDII, which has several experimental structures captured in the deactivated conformation. In particular, two chimeric structures were used as templates:

- **PDB 6N4R** (Xu et al., 2019)  
- **PDB 7K48** (Wisedchaisri et al., 2020)

Both are NaVAb-based chimeras containing the human NaV1.7 VSDII sequence, stabilized in the deactivated state by peptide toxins. 

To run ColabFold in **custom template mode**, we used exactly the same options as in Section 3, but added the flag:

- `--custom-template-path $PATH`  
  where `$PATH` points to the folder containing the desired template PDBs (in this repository, we provide them under `custom_templates/`).

**Example execution:**

```bash
colabfold_batch --num-models 5 --model-type auto --msa-mode mmseqs2_uniref_env \
--num-seeds 20 \
--templates --max-seq 256 --max-extra-seq 512 \
--num-recycle 6 --save-recycles \
--custom-template-path ./custom_templates/ \
inputs/hNaV1.7-alpha-only.fa outputs/hNaV1.7-alpha-only_customtemplates/
```

## 4. Analyze models

This section illustrates how to compute (i) **distance coordinates** and (ii) **regional pLDDT subsets** from a folder containing all the PDB models for a given case. We keep the running example from above and assume the output directory is: outputs/hNaV1.7-alpha-only/

A typical layout (abbreviated) after a ColabFold run with `--save-recycles` is:

```bash
ll outputs/nav1.7_alphaonly/
# total ...
# -rw-r--r--  model_000_seed_000.r0.pdb
# -rw-r--r--  model_000_seed_000.r1.pdb
# ...
# -rw-r--r--  model_004_seed_019.r6.pdb
```

We provide two scripts under `scripts/`:
- `calculate_distances.py` → computes atom–atom distances defined in a **metrics** file.  
- `subset_plddt.py` → computes average pLDDT on residue **subsets** defined in a **regions** file.  

Both scripts take:
- `--source` → the directory with all the PDBs to analyze.  
- A second argument file:  
  - `--metrics <file>` specifying atom pairs for distance calculations of residue ranges for subset pLDDT calculation

> **Note on tags.** Each line in the input TXT defines a **tag** that becomes the column name in the resulting CSV. Use stable, descriptive tags (e.g., `VSDI`, `VSDII`, `AG1`, `IFM`, `SF`, `plddt_VSDI`, …).

### 4.1 Distance coordinates

Create a plain-text **metrics** file that lists atom pairs (chain, residue index, atom name) per line with a trailing **tag**. Below is an example for **NaV1.7**, derived directly from the reference table (see `references/distances-and-subsets.csv`, row `NaV1.7`):

**`references/nav17-alphaonly_distances.txt`**

```
chain1 res1 atom1 chain2 res2 atom2 tag  
A 163  CA    A 214  CA    VSDI  
A 776  CA    A 824  CA    VSDII  
A 1226 CA    A 1276 CA    VSDIII  
A 1547 CA    A 1599 CA    VSDIV  
A 398  CA    A 1442 CA    AG1  
A 953  CA    A 1745 CA    AG2  
A 1462 CA    A 1635 CA    IFM  
A 361  CA    A 1395 CA    SF  
```

Run the distance calculation:

```bash
python scripts/calculate_distances.py \
  --source outputs/hNaV1.7-alpha-only/ \
  --metrics references/nav17-alphaonly_distances.txt
  ```

**Output**: a CSV file outputs/hNaV1.7-alpha-only/distances.csv with one row per PDB and one column per tag. Example columns: filename,VSDI,VSDII,VSDIII,VSDIV,AG1,AG2,IFM,SF.

### 4.2 pLDDT subsets

Create a plain-text **regions** file listing residue ranges (comma-separated if multiple segments) and a **tag** for each subset. Again, for **NaV1.7** we take ranges from the reference table:

- `plddt full` → `A17-430,A701-991,A1163-1897`  
- `plddt VSDI` → `A114-228`  
- `plddt VSDII` → `A717-838`  
- `plddt VSDIII` → `A1169-1293`  
- `plddt VSDIV` → `A1492-1616`  
- `plddt AG` → `A397-399,A952-954,A1441-1443,A1744-1746`  
- `plddt IFM` → `A1461-1463`  
- `plddt SF` → `A359-363,A913-918,A1392-1397,A1685-1689`  

**`references/nav17-alphaonly_ranges.txt`**

```
range tag
A17-430,A701-991,A1163-1897 plddt_full
A114-228 plddt_VSDI
A717-838 plddt_VSDII
A1169-1293 plddt_VSDIII
A1492-1616 plddt_VSDIV
A397-399,A952-954,A1441-1443,A1744-1746 plddt_AG
A1461-1463 plddt_IFM
A359-363,A913-918,A1392-1397,A1685-1689 plddt_SF
```

Run the subset pLDDT calculation:

```bash
python scripts/subset_plddt.py \
  --source outputs/hNaV1.7-alpha-only/ \
  --metrics references/nav17-alphaonly_ranges.txt
```

**Output**: a CSV file outputs/hNaV1.7-alpha-only/plddt_subsets.csv with columns like filename,plddt_full,plddt_VSDI,...,plddt_SF.

### 4.3 Merging results and example output

Join distance and pLDDT tables on filename to obtain a single analysis table suitable for plotting or further statistics (e.g., in Python/pandas, R, or your tool of choice). A typical merged CSV resembles:
```
filename,plddt_full,plddt_VSDI,plddt_VSDII,plddt_VSDIII,plddt_VSDIV,plddt_AG,plddt_IFM,plddt_SF,VSDI,VSDII,VSDIII,VSDIV,AG1,AG2,IFM,SF
...
```

This mirrors the structure of example_outputs/NaV1.1-alphaOnly.csv provided in the repository (columns may vary with your chosen tags).

>**Tip**. Keep your tag names consistent across cases (e.g., always VSDI, VSDII, AG1, IFM, SF, plddt_VSDI, …) so merged datasets are directly comparable.
