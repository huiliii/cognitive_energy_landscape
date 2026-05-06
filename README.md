# Code Availability

This repository contains the code used in our paper. The analysis pipeline is organized into two main components: control energy computation and GAMLSS modeling.

## 1. Control Energy

### 1.1 Extract Cognitive Activations

- **Extract keywords:** `keywords_extract.py`
- **Keywords → activation:** `neurosyth_activation.py`, `activation_state.py`
- **Activation → control energy:** `control_energy.py`

## 2. GAMLSS Model

- **Fit model:** `GAMLSS_plot.R`
- **Distribution selection:** `GAMLSS_model_selection.R`
