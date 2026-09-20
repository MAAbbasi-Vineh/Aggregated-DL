# Development of a Novel Aggregated Deep Learning Framework for Small Biological Datasets Using Overlapping Subsequences

A hybrid CNN–LSTM–Attention–Residual deep learning framework for the classification of full-length genomic regulatory sequences from small biological datasets. This repository contains all scripts, datasets, pre-computed embeddings, saved models, and supplementary data associated with the study.

---

## Table of Contents

- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [How to Use This Repository](#how-to-use-this-repository)
  - [Step 1 — Overlapping Subsequence Augmentation](#step-1--overlapping-subsequence-augmentation)
  - [Step 2 — Train and Evaluate the Aggregated Model](#step-2--train-and-evaluate-the-aggregated-model)
  - [Step 3 — Further Analyses](#step-3--further-analyses)
- [Datasets](#datasets)
- [Embeddings and Saved Models](#embeddings-and-saved-models)
- [Requirements](#requirements)
- [Installation](#installation)
- [Model Architecture](#model-architecture)
- [Citation](#citation)
- [License](#license)
- [Contact](#contact)

---

## Overview

This framework addresses a core challenge in computational genomics: classifying full-length regulatory sequences (200 bp) when only a small number of labeled examples are available (50–100 sequences per class). It extends the overlapping subsequence augmentation strategy introduced by Abbasi-Vineh et al. (2025) ([DOI: 10.1038/s41598-025-12796-9](https://doi.org/10.1038/s41598-025-12796-9)) by introducing a novel aggregation mechanism that transfers features learned at the subsequence level back to the original full-length sequences.

**Key components:**

- **Overlapping subsequence augmentation** — each 200 bp sequence is padded with 40-nt poly-N at both ends and decomposed into overlapping 40-nt windows using variable overlaps (5–20 nt) and a 15-consecutive-nucleotide commonality criterion, producing ~240 unique subsequences per original sequence.
- **Hybrid CNN–LSTM–Attention–Residual model** — hierarchically extracts local motif features (CNN), models positional dependencies (bidirectional LSTM), and dynamically weights informative positions (CNN-Attention and LSTM-Attention).
- **Coverage-aware feature aggregation** — subsequence-level features are aggregated to the original full-length sequence using binary masking (training) and attention-weighted averaging (evaluation).
- **Dual hybrid loss function** — jointly optimises subsequence-level and sequence-level cross-entropy losses (weighting factor α = 0.7).

The framework was evaluated on three independent genomic datasets: divergent organisms (6 classes, 100 sequences/class), chloroplast genomes (6 classes, 50 sequences/class), and regulatory sequences with and without Shine–Dalgarno motifs (2 classes, 100 sequences/class).

---

## Repository Structure

```
Aggregated-DL/
│
├── Script_for_Augmentation_Overlapping_Subsequences.ipynb   # Step 1: generate augmented subsequences
├── Script_for_Aggregated_DL_Main_Model.ipynb                # Step 2: train and evaluate the model
├── Primary_Scripts_Aggregated_DL.ipynb                      # Step 3: further analyses and visualisation
│
├── Datasets.rar                                             # All three datasets (FASTA + CSV)
├── Final saved models.rar                                   # Pre-trained model weights and metadata
├── Figures.zip                                              # All figures generated in the study
├── Supplementary Data.zip                                   # Supplementary tables and results
│
├── README.md
└── LICENSE
```

---

## How to Use This Repository

The pipeline consists of three sequential steps, each implemented in a dedicated notebook.

---

### Step 1 — Overlapping Subsequence Augmentation

**Notebook:** [`Script_for_Augmentation_Overlapping_Subsequences.ipynb`](https://github.com/parkingvarsson/Aggregated-DL/blob/main/Script_for_Augmentation_Overlapping_Subsequences.ipynb)

Run this notebook **first**, before the main model. It reads FASTA files of original 200-nt sequences and generates the augmented 40-nt overlapping subsequences required as input to Step 2.

**What it does:**
1. Adds 40-nt poly-N padding to both ends of each original sequence.
2. Generates all unique 40-nt subsequences using variable overlaps and a 15-nt commonality criterion.
3. Writes the subsequences to a CSV file (`sequence, label` format, one row per subsequence).

**Adjustable parameters:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `k` | `40` | Length of each augmented subsequence (nt) |
| `min_overlap` | `5` | Minimum overlap between consecutive windows (nt) |
| `max_overlap` | `20` | Maximum overlap between consecutive windows (nt) |
| Commonality threshold | `15` | Minimum consecutive shared nucleotides for inclusion |

**Usage:**
- Update `input_fasta` and `output_file` to your actual paths in the final cell.
- Run **once per class** (one FASTA file → one CSV file).
- Repeat for all classes before proceeding to Step 2.

---

### Step 2 — Train and Evaluate the Aggregated Model

**Notebook:** [`Script_for_Aggregated_DL_Main_Model.ipynb`](https://github.com/parkingvarsson/Aggregated-DL/blob/main/Script_for_Aggregated_DL_Main_Model.ipynb)

Run this notebook after Step 1. It loads the augmented CSV files (from Step 1) and the original 200-nt CSV files, trains the hybrid model under stratified k-fold cross-validation, and evaluates performance at both the subsequence and original sequence levels.

**The notebook is structured sequentially — run all cells top to bottom. The `main()` function at the end of the notebook executes the full pipeline automatically.**

| Section | Content |
|---------|---------|
| 1. Dependencies and Hyperparameters | Imports, device, all training parameters |
| 2. Data Loading and Encoding | Sequence mapping, validation, one-hot encoding |
| 3. PyTorch Dataset | `SequenceDataset` wrapper |
| 4. Model Architecture | CNN blocks, attention modules, full model |
| 5. Training | Dual hybrid loss, validation step |
| 6. Evaluation | Subsequence-level and original sequence-level (attention-weighted aggregation) |
| 7. Model Serialisation | Saves weights, metadata, mappings, embeddings |
| 8. Main Pipeline | `main()` — called automatically at the end |

**Before running, update the data paths in Section 1:**

```python
DATA_FOLDER_40nt  = '/content/drive/MyDrive/.../40nt'   # CSV files from Step 1
DATA_FOLDER_200nt = '/content/drive/MyDrive/.../200nt'  # original 200-nt CSV files
SAVE_DIR          = '/content/drive/MyDrive/.../output'
```

**Key training hyperparameters:**

| Parameter | Value |
|-----------|-------|
| Subsequence length | 40 nt |
| Batch size | 256 |
| Epochs | 100 |
| Early stopping patience | 10 |
| Optimizer | AdamW (weight decay = 0.01) |
| Scheduler | OneCycleLR |
| Loss weighting (α) | 0.7 |
| Cross-validation | Stratified 3-fold |

---

### Step 3 — Further Analyses

**Notebook:** [`Primary_Scripts_Aggregated_DL.ipynb`](https://github.com/parkingvarsson/Aggregated-DL/blob/main/Primary_Scripts_Aggregated_DL.ipynb)

Use this notebook after completing Steps 1 and 2 for all downstream and complementary analyses, including:

- ROC and precision-recall curves
- Confusion matrices (subsequence and sequence level)
- Calibration analysis
- Gradient-based saliency mapping and positional importance scores
- Gene identity matrix and pairwise cosine similarity
- Positional SD-motif frequency analysis and sequence logo generation
- Aggregation validity experiments (Leave-One-Gene-Out, cross-class mapping disruption, nucleotide corruption, similarity-based splitting, Leave-One-Species-Out)

---

## Datasets

All datasets are provided in `Datasets.rar`. Each dataset contains original 200-nt sequences (FASTA and CSV) and the corresponding augmented 40-nt subsequences (CSV) produced by Step 1. All sequences are 200 bp upstream regions of coding sequences (CDS) retrieved from NCBI (accessed February 2025).

### Dataset 1 — Divergent Organisms (6 classes, 100 sequences/class)

| Class | Organism | NCBI Accession |
|-------|----------|----------------|
| Archaeon | — | NZ_CP145900 |
| Bacterium | — | CP178715 |
| Fungus | — | NC_018946 |
| Plant | *Arabidopsis thaliana* (chr. 4) | CP002687 |
| Protist | — | JAPFFF010000003 |
| Virus | — | PV416404 |

### Dataset 2 — Chloroplast Genomes (6 classes, 50 sequences/class)

| Class | Organism | NCBI Accession |
|-------|----------|----------------|
| Arabidopsis | *Arabidopsis thaliana* | NC_000932.1 |
| Chlamydomonas | *Chlamydomonas reinhardtii* | NC_005353.1 |
| Chlorella | *Chlorella vulgaris* | NC_001865.1 |
| Nicotiana | *Nicotiana tabacum* | MZ707522.1 |
| Porphyridium | *Porphyridium purpureum* | NC_023133 |
| Triticum | *Triticum aestivum* | NC_002762.1 |

### Dataset 3 — Shine–Dalgarno Motifs (2 classes, 100 sequences/class)

Chloroplast regulatory sequences from the same species, with and without canonical Shine–Dalgarno (SD) motifs in the 20 bp immediately upstream of the CDS start codon. Both classes derive from the same genomic background to minimise confounding compositional differences between classes.

---

## Embeddings and Saved Models

Pre-trained model weights, training metadata, sequence mappings, and pre-computed CNN embeddings (FC3 feature vectors, 512-dimensional) for all three datasets are provided in `Final saved models.rar`.
> **Note:** One of the saved model files exceeds GitHub's file size limit and therefore could not be uploaded to this repository. Please contact us by email to request the file.

Embeddings can be loaded directly for downstream analyses — clustering, transfer learning, or similarity analysis — without re-running the full training pipeline:

```python
import numpy as np

embeddings = np.load('dataset1_embeddings.npy', allow_pickle=True).item()
# Returns a dict: {sequence_id: numpy array (512-dimensional FC3 feature vector)}


```

---

## Requirements

Python ≥ 3.10. A GPU is strongly recommended (the study used a Tesla T4 GPU on Google Colab).

```bash
pip install torch numpy pandas scikit-learn matplotlib seaborn scipy biopython
```

| Package | Purpose |
|---------|---------|
| `torch >= 2.0` | Model training and inference |
| `numpy` | Array operations |
| `pandas` | Data loading |
| `scikit-learn` | Cross-validation and metrics |
| `matplotlib`, `seaborn` | Visualisation |
| `scipy` | Saliency smoothing |
| `biopython` | FASTA parsing (Step 1 only) |

---

## Installation

**Option 1 — Clone the repository:**

```bash
git clone https://github.com/parkingvarsson/Aggregated-DL.git
cd Aggregated-DL
```

**Option 2 — Open directly in Google Colab:**

Click the badge to open the main model notebook directly in Colab:

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/parkingvarsson/Aggregated-DL/blob/main/Script_for_Aggregated_DL_Main_Model.ipynb)

When running on Colab, first mount Google Drive and enable a **T4 GPU** runtime (Runtime → Change runtime type → T4 GPU):

```python
from google.colab import drive
drive.mount('/content/drive')
```

---

## Model Architecture

```
Input (40-nt one-hot encoded, 5 channels: A, T, C, G, N)
    │
    ▼
CNN Block 1 — Conv1D (128 filters, k=3) + BN + MaxPool + Dropout(0.3) + Residual + CNN-Attention
    │
    ▼
CNN Block 2 — Multi-kernel Conv1D (k=3,5,7; 256 filters each → 768 ch) + BN + MaxPool + Dropout(0.5) + Residual + CNN-Attention
    │
    ▼
CNN Block 3 — Multi-kernel Conv1D (k=3,5,7; 512 filters each → 1536 ch) + BN + MaxPool + Dropout(0.3) + Residual + CNN-Attention
    │
    ▼
Positional Encoding
    │
    ▼
Bidirectional LSTM (2 layers, hidden=256/direction, total=512) + Dropout(0.3)
    │
    ▼
LSTM-Attention → context vector + attention weights (stored for aggregation)
    │
    ▼
FC1 (256 units) + BN + Dropout(0.3) + ReLU
FC2 (512 units) + Dropout(0.5) + ReLU
FC3 (512 units) + Dropout(0.3) + ReLU   ← embeddings extracted here
    │
    ▼
FC4 → Softmax → class probabilities
```

For original sequence-level classification, LSTM-Attention weights from all 240 subsequences are normalised and used to compute a weighted average of FC3 feature vectors (Equation 4, manuscript), which is then passed through FC4 for the final prediction.

---

## Citation

If you use this code, data, or models in your research, please cite:

```
[Citation will be added upon publication]
```

The data augmentation method this framework extends is described in:

> Abbasi-Vineh et al. (2026). *Scientific Reports*. DOI: 10.1038/s41598-026-69140-y

---

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

---

## Contact

For technical questions, requests for additional data or analysis details, or information about the latest version of the framework, please contact the authors:

* **Dr. Mohammad Ali Abbasi-Vineh:** [maa.vineh@gmail.com](mailto:maa.vineh@gmail.com); [Ali.abbasi@mpimp-golm.mpg.de](mailto:Ali.abbasi@mpimp-golm.mpg.de)
* **Dr. Naser Farrokhi:** [n_farrokhi@sbu.ac.ir](mailto:n_farrokhi@sbu.ac.ir)
* **Dr. Pär K. Ingvarsson:** [par.ingvarsson@slu.se](mailto:par.ingvarsson@slu.se)

> **Note:** Optimised versions of this framework, informed by ongoing functional and experimental studies, are under active development. If you plan to apply the framework to your own datasets, we encourage you to contact the authors for information about the latest updates and recommended implementation practices.

