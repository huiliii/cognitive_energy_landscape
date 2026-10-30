# Code Availability

This repository contains the code used in our paper. The analysis pipeline is organized into three components: control energy computation, GAMLSS modeling, and gene expression analysis.

## 1. Control Energy

### 1.1 Extract Cognitive Activations

- **Extract keywords:** `keywords_extract.py`
- **Keywords → activation:** `neurosyth_activation.py`, `activation_state.py`
- **Activation → control energy:** `control_energy.py`

## 2. GAMLSS Model

- **Fit model:** `GAMLSS_plot.R`
- **Distribution selection:** `GAMLSS_model_selection.R`

## 3. Gene Expression Analysis

- **Toy example and AHBA profile construction:** [`CODE/Gene_Expression/gene_expression_example.py`](CODE/Gene_Expression/gene_expression_example.py)
- **Workflow, parameters, gene lists, and usage:** [Gene expression documentation](CODE/Gene_Expression/README.md)
- **Process-associated gene sets:** [`kang2011_gene_sets.csv`](CODE/Gene_Expression/kang2011_gene_sets.csv), based on Kang et al. (2011), Supplementary Table 13.

The default example uses synthetic data. An optional real-data workflow extracts regional AHBA expression with abagen and documents gene matching, aggregation, and common normalization.
