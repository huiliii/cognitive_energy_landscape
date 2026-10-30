# Code Availability

This repository contains the code used in our paper. The analysis pipeline is organized into three components: control energy computation, GAMLSS modeling, and transcriptional analysis.

## 1. Control Energy

### 1.1 Extract Cognitive Activations

- **Extract keywords:** `keywords_extract.py`
- **Keywords → activation:** `neurosyth_activation.py`, `activation_state.py`
- **Activation → control energy:** `control_energy.py`

## 2. GAMLSS Model

- **Fit model:** `GAMLSS_plot.R`
- **Distribution selection:** `GAMLSS_model_selection.R`

## 3. Transcriptional Analysis

- **Extract regional AHBA expression:** [`extract_ahba.py`](CODE/Transcriptional_Analysis/extract_ahba.py)
- **Construct six process-specific profiles:** [`build_process_profiles.m`](CODE/Transcriptional_Analysis/build_process_profiles.m)
- **Gene lists and workflow documentation:** [Transcriptional Analysis](CODE/Transcriptional_Analysis/README.md)

The study workflow covers AHBA extraction, process-associated gene matching, regional averaging, and common normalization. File locations are configurable; machine-specific directories are removed.
