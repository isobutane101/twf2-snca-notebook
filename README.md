# TWF2 knockdown and reduced *SNCA* expression

Analysis code for **"Association Between TWF2 Knockdown and Reduced SNCA Expression: A Screen of 300 CRISPRi Perturbations in Human Embryonic Stem Cells."**

Aditya S. Desai¹, Daniel Plymire²
¹ Lakeside School, Seattle, WA, USA · ² University of Texas at Arlington, Arlington, TX, USA

This repository contains the single notebook that produces every number, table and figure in the manuscript: [`VCC_SNCA_FullScreen_v9.ipynb`](VCC_SNCA_FullScreen_v9.ipynb).

## Summary

Elevated α-synuclein, encoded by *SNCA*, is central to Parkinson's disease pathology, and lowering *SNCA* expression is an actively pursued therapeutic strategy. A gene whose inhibition reduces *SNCA* **without broadly disturbing the transcriptome** is a more tractable starting point than one that reorganizes it.

We screened all 300 CRISPRi perturbations in the Arc Institute Virtual Cell Challenge 2025 dataset on two axes at once — the effect on *SNCA*, and the number of other genes that change.

Twenty-eight perturbations significantly altered *SNCA*; twelve reduced it, four by at least 50%, and only two by at least 75%. Under criteria fixed before the data were loaded, **TWF2 was the only perturbation of the 300 satisfying both.**

| Perturbation | Cells | *SNCA* reduction | Adjusted p | Genes meeting criteria |
|---|---:|---:|---:|---:|
| SNCA | 386 | 95.2% | 1.4 × 10⁻¹⁴⁰ | 97 |
| **TWF2** | **1,008** | **92.4%** | **< 1 × 10⁻³⁰⁰** | **14** |
| EP300 | 2,217 | 69.4% | < 1 × 10⁻³⁰⁰ | 2,855 |
| DOT1L | 1,484 | 57.1% | 1.2 × 10⁻¹⁶⁸ | 2,472 |
| SOX2 | 806 | 46.2% | 4.3 × 10⁻⁵⁹ | 2,529 |

Median across all 300 perturbations: 136 genes altered.

EP300 and DOT1L are chromatin regulators; both would appear as hits in a screen measuring *SNCA* alone. Measuring the transcriptome-wide footprint alongside the target effect is what separates them from TWF2.

> These are **transcript-level associations** in a stem cell line. They identify TWF2 as a candidate for experimental investigation. They do not establish mechanism, biological replication, or therapeutic suitability. See [Limitations](#limitations).

## Dataset

Arc Institute Virtual Cell Challenge 2025 — CRISPRi in H1 human embryonic stem cells, 10x Genomics Flex chemistry.

```
gs://arc-institute-virtual-cell-atlas/virtual-cell-challenge/2025/
```

| | |
|---|---|
| Total cells | 491,046 |
| Non-targeting controls | 114,528 |
| CRISPRi target genes | 300 |
| Processing batches | 48 (`Flex_1_01`–`Flex_3_16`) |
| Genes | 18,080 → 12,760 after filtering (detected in > 5% of control cells) |

The three released files (training: 221,273 cells / 150 targets; validation: 98,927 / 50; test: 170,846 / 100) are combined. All 48 batches appear in all three files, confirming they are partitions of one experiment rather than separate experiments; file membership is treated as a label, not a batch variable.

## Method

**Prespecified criteria** (fixed before the data were loaded):

- A perturbation qualifies if it reduces *SNCA* by ≥ 50% at adjusted p < 0.10 **and** alters no more than 20 genes.
- Differential expression: adjusted p < 0.10 and |log2FC| > 0.40, following Replogle et al.

**Three analyses with different independence assumptions**, which gave concordant estimates of the TWF2 effect on *SNCA*:

| Analysis | log2FC | Adjusted p | Genes altered |
|---|---:|---:|---:|
| Cell-level Wilcoxon rank-sum, batch-matched | −3.72 | < 1 × 10⁻³⁰⁰ | 14 |
| Paired within-batch pseudobulk, 48 batches | −3.36 | 1.2 × 10⁻²⁶ | 52 |
| DESeq2, design `~ batch + condition` | −3.49 | 5.5 × 10⁻²³¹ | 21 |

A reduction in *SNCA* was observed in **30 of 30** batches containing ≥ 20 TWF2 cells (sign test p = 1.9 × 10⁻⁹).

Two implementation notes that matter: the rank-sum test is **tie-corrected** — single-cell counts are heavily tied at zero and omitting the correction inflates significance — and effect sizes follow the scanpy definition (base-2 log of the ratio of back-transformed group means). The engine is verified against `scipy.stats.mannwhitneyu` and against scanpy with tie correction enabled (§8.1).

**Robustness** (§11): null calibration by splitting controls in half, control subsampling to 1,000 cells across five seeds, threshold sensitivity grids, matched-cell-number comparison, and leave-one-batch-out. At matched cell number the difference in gene counts persists, indicating it is not explained by statistical power.

## Notebook contents

| § | |
|---|---|
| 0–3 | Connect (Drive, GCS), prespecified config, environment, data acquisition |
| 4–5 | Dataset inventory — what replicate structure actually exists — and QC |
| 6–7 | Pseudobulk construction (from raw counts, before normalization); preprocessing |
| 8 | Differential expression engine + correctness check against scipy and scanpy |
| 9 | Full screen: every perturbation, every gene; null calibration |
| 10 | Replicate-level (pseudobulk) analysis; agreement with cell level; DESeq2 shortlist |
| 11 | Robustness: subsampling, threshold sensitivity, power diagnostic |
| 12–13 | Knockdown efficiency across all targets; the screen, applied |
| 14–16 | Guide-level concordance; gnomAD constraint annotation; TWF2 vs *SNCA* responses |
| 17–19 | Figures; single source of truth for the manuscript; open decisions |

## Outputs

Running the notebook writes two directories:

- **`results/`** — 30+ CSV/JSON tables, numbered by stage: `04_cell_level_screen.csv`, `05_combined_screen.csv`, `07_threshold_sensitivity.csv`, `09_screen_with_criteria.csv`, `12_per_batch_replication.csv`, and `manuscript_numbers.json`, which holds every value quoted in the paper.
- **`figures/`** — `Fig1_landscape`, `Fig2_expression_context`, `Fig3_volcanoes`, `Fig4_criteria_plot`, `Fig5_constraint`, `Fig6_threshold_sensitivity`, `Fig7_power_diagnostic`, `Fig8_pseudobulk_concordance`, `Fig9_per_batch_replication`.

## Reproducing

The notebook is written for Google Colab with GPU. Open it, connect a GPU runtime, and run cells in order — §0 handles Drive and Google Cloud authentication, §3 pulls the dataset from the public bucket.

The complete screen runs in **≈ 15 minutes on a single NVIDIA A100**. All analyses use a fixed random seed (**42**); every subsampling seed is recorded in the notebook and in `results/`.

```
Python 3.13.15   numpy 2.1.3    pandas 3.0.5    scipy 1.16.3
scanpy 1.12.3    anndata 0.13.2  PyDESeq2        PyTorch 2.11.0
```

## Limitations

- The 48 batches are **processing units within a single experiment, not independent biological replicates.** Blocking on batch removes cell-level pseudoreplication and controls for processing variation, but no analysis here demonstrates independent biological replication.
- **Guide-level validation is not possible in this dataset.** Guide identity is recorded as a jointly delivered pair, and every cell for a given target received both guides. Only 17 of 300 targets carry more than one guide value, and neither TWF2 nor *SNCA* is among them.
- **The gene count is not a specificity measure.** It depends on thresholds, control set size and perturbation toxicity, which is why threshold-sensitivity and matched-sample-size analyses are reported alongside it.
- Nineteen of 300 perturbations had too few cells per batch for replicate-level analysis and are reported at cell level only.
- No mechanism is proposed. TWF2 encodes twinfilin-2, an actin-monomer-binding protein; transcriptional, post-transcriptional and indirect cytoskeletal-signaling explanations are all consistent with these data.

## Acknowledgements

We thank the Arc Institute for making the Virtual Cell Challenge dataset publicly available, and Alice Yu for guidance throughout this project. Computational resources were provided by Google Colab Pro+.

The authors received no funding for this work and declare no conflicts of interest.
