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
- **plddt full**: Residue ranges covering the entire α-subunit excluding unstructured flexible loops, used for overall pLDDT averaging.  
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
> **Note on input sequences.**  
> All sequences used in this study are provided in  
> `references/input-sequences.csv`, including protein names, UniProt accessions,  
> and the exact sequences used for each modeling case.

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

This section illustrates how to compute:
1. **Distance coordinates** from atom pairs  
2. **Regional pLDDT subsets**  
3. **RMSD to reference structures**  
4. **ipTM scores** (from `.json` metadata or a `.log` file)

We keep the running example and assume all PDB models for a case are stored under:
```
outputs/hNaV1.7-alpha-only/
```

A typical layout (abbreviated) after a ColabFold run with `--save-recycles` is:

```bash
ll outputs/hNaV1.7-alpha-only/
# total ...
# -rw-r--r--  model_000_seed_000.r0.pdb
# -rw-r--r--  model_000_seed_000.r1.pdb
# ...
# -rw-r--r--  model_004_seed_019.r6.pdb
```

### Available analysis scripts
All scripts live in `scripts/`:

| Script | What it computes | Inputs | Output |
|--------|-----------------|--------|--------|
| `calculate_distances.py` | Atom–atom distances using a metrics file | `--source`, `--metrics` | `distances.csv` |
| `subset_plddt.py` | Average pLDDT on residue subsets | `--source`, `--regions` | `subset_plddt.csv` |
| `rmsd.py` | RMSD vs. one or more reference PDBs | `<models_dir> <REF1> [REF2…]` | `rmsd.csv` |
| `iptm.py` | Extracts ipTM scores from ColabFold JSON dump | `--source <folder_with_json>` | `iptm.csv` |
| `iptm_log.py` | Extracts ipTM and recycle index from ColabFold `.log` | `--source <logfile>` | `iptm_by_log.csv` |

#### Python packages
- numpy
- pandas
- tqdm
- tap
- matplotlib
- pyrosetta (for calculate_distances and subset_plddt)
- chimerax (for rmsd.py)


### Example usage

#### Distances
```bash
python scripts/calculate_distances.py \
    --source outputs/hNaV1.7-alpha-only \
    --metrics metrics.txt
```
*See structure of `metrics.txt` in section **4.1 Distance coordinates** below.*


#### pLDDT subsets
```bash
python scripts/subset_plddt.py \
    --source outputs/hNaV1.7-alpha-only \
    --regions regions.txt
```

*See structure of `regions.txt` in section **4.2 pLDDT subsets** below.*

#### RMSD to reference structure(s)
```bash
# One reference
python scripts/rmsd.py outputs/hNaV1.7-alpha-only 7DTD

# Multiple references
python scripts/rmsd.py outputs/hNaV1.7-alpha-only 7DTD 6J8E
```

#### Extract ipTM values
From JSON metadata:
```bash
python scripts/iptm.py --source outputs/hNaV1.7-alpha-only
```

From a run log:
```bash
python scripts/iptm_log.py --source outputs/hNaV1.7-alpha-only/log.txt
```

### About tags
Each line in the input TXT defines a **tag** that becomes the column name in the resulting CSV. Use stable, descriptive tags (e.g., `VSDI`, `VSDII`, `AG1`, `IFM`, `SF`, `plddt_VSDI`, …).

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

This mirrors the structure of example_outputs/NaV1.1-alphaOnly.csv provided in the repository (columns may vary with your chosen tag names).

>**Tip**. Keep your tag names consistent across cases (e.g., always VSDI, VSDII, AG1, IFM, SF, plddt_VSDI, …) so merged datasets from different study cases are directly comparable.



## 5. Correlation and clustering analysis

This repository includes ready-to-run Jupyter notebooks and the **final analysis tables** used in the paper. The folder layout is:

```
correlation_analysis/
├─ AlphaOnly/
│  ├─ Correlations-final.ipynb
│  └─ data/   # final per-model tables (CSV) for α-subunits only
├─ Betas/
│  ├─ Correlations-final.ipynb
│  └─ data/   # final per-model tables (CSV) for α+β complexes and corresponding α-subunits only as reference
└─ CaM/
   ├─ Correlations-final.ipynb
   └─ data/   # final per-model tables (CSV) for α+CaM complexes and corresponding α-subunits only as reference

clustering_analysis/
├─ clustering-final.ipynb
└─ data/   # combined tables used for joint clustering/UMAP
```

- The **`data/`** subfolders contain the **final CSV tables** produced in this study (e.g., `NaV1.7-alphaOnly.csv`, `NaV1.1-beta1.csv`). Each row is a single model; columns include distance-based state coordinates (e.g., VSDI–VSDIV, AG1/AG2, IFM, SF), confidence metrics (global and subset **pLDDT**), **RMSDs to experimental reference structures**, and metadata such as `filename` (from which the recycle index can be extracted).
- The notebooks (`*final.ipynb`) are configured to **run out-of-the-box** using the corresponding `data/` CSVs to reproduce the **figures and statistics** reported in the paper (correlations, mutual information heatmaps, recycle trends, etc.). Minimal editing should be required beyond adjusting paths if you move the repository.


## 6. Data availability & model retrieval

- The complete **3D model ensembles** for all cases will be available on **Dryad**: *[DOI link placeholder]*. Each Dryad package mirrors this repository’s case naming (e.g., `hNaV1.7-alpha-only`).
- Using the CSVs in the folders above, the community can **search for specific models** (desired state coordinates and confidence) and then retrieve the corresponding PDB from the Dryad bundle.

**Example workflow (NaV1.7 α-subunit only):**
1. Download the Dryad bundle for *NaV1.7 α-only* models.
2. Open `correlation_analysis/AlphaOnly/data/NaV1.7-alphaOnly.csv`.
3. Filter by column values to locate models with the features you need, e.g.:
   - `VSDII < 13` (more deactivated VSDII)
   - `IFM < 10` (IFM bound)
   - `AG_area` between `160` and `190` Å² (narrow/closed–like gate)
   - `plddt_VSDII > 80` (high confidence in VSDII)
4. Use the `filename` column to find the exact PDB file in the Dryad package (e.g., `model_000_seed_012.r0.pdb`).

> Tip: keep tags/column names consistent across cases (e.g., `VSDI–VSDIV`, `AG1/AG2`, `IFM`, `SF`, `AG_area`, `plddt_*`, `recycle`) so that filtering and cross-case comparisons are straightforward.


