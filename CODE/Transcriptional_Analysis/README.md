# 3. Transcriptional Analysis

This section provides the study workflow for extracting adult Allen Human Brain Atlas (AHBA) expression and constructing regional profiles for six neurodevelopmental processes. File locations are supplied as inputs; no machine-specific directories are required.

## 3.1 Extract regional AHBA expression

The original extraction used `abagen.get_expression_data` with a volumetric AAL2 atlas and default parameters. The script adds a processing report for reproducibility.

```bash
pip install abagen==0.1.3
python CODE/Transcriptional_Analysis/extract_ahba.py \
  --atlas /path/to/AAL2_mni.nii.gz \
  --output-dir transcriptional_output
```

Outputs are `abagen_expression.csv` (regions × genes, with a label column) and `abagen_report.txt`. Use a Python environment compatible with abagen v0.1.3.

The defaults include probe reannotation, intensity-based filtering at 0.5, differential-stability probe selection, corrected MNI sample coordinates, 2-mm sample-matching tolerance, scaled robust sigmoid sample and gene normalization, and averaging within donors followed by averaging across donors. No atlas-info table, hemisphere mirroring, or missing-region interpolation was supplied. See the [versioned abagen documentation](https://abagen.readthedocs.io/en/stable/generated/abagen.get_expression_data.html) and the generated report for details.

## 3.2 Construct process-specific regional profiles

The accompanying `kang2011_gene_sets.csv` contains the six process-associated lists from the study's transcription of Supplementary Table 13 of [Kang et al. (2011)](https://doi.org/10.1038/nature10523). These are genes correlated with markers of myelination, axonal development, dendrite development, neuron differentiation, neuron migration, and synapse development. Original gene symbols are retained without automatic alias conversion.

Run the MATLAB function with your processed expression matrix:

```matlab
addpath('CODE/Transcriptional_Analysis');
build_process_profiles( ...
    'transcriptional_output/abagen_expression.csv', ...
    'CODE/Transcriptional_Analysis/kang2011_gene_sets.csv', ...
    'transcriptional_output/process_profiles.mat');
```

The workflow:

1. Retains the original AAL2 rows `[1:28 33:94]` (90 regions; MATLAB indexing), after checking the input labels.
2. Deduplicates each process-associated gene list and matches symbols exactly to the expression columns.
3. Averages matched gene expression within each region for each process.
4. Constructs a reference profile from the unique union of matched genes across all six processes.
5. Divides all profiles by the same reference mean across retained regions.

The MAT output records raw gene-expression submatrices, process means, normalized profiles, retained region labels, matched and unmatched genes, and the common scaling factor. Profile columns are: union of six sets, synapse development, neuron migration, neuron differentiation, dendrite development, axonal development, and myelination. Columns 2–7 provide the six regional control-input weight vectors.

Adult-derived expression profiles remain fixed across developmental periods; they are not age-specific expression measurements.

## 3.3 Historical intermediate processing and missing regions

The archived extraction CSV contains 15,633 genes; the downstream matrix used in the study contains 15,184, with identical expression values for retained columns. The intermediate exclusion step is not recorded in the available extraction code. For exact reproduction, use the archived downstream matrix and report its gene membership; a fresh extraction should not be assumed identical.

The archived intermediate matrix contains missing expression at AAL2 label 92, while its saved final process-profile values match contralateral label 91. The available extraction script does not record this edit. The supplied function therefore stops on nonfinite expression and requires missing-region handling to be explicitly documented and applied, rather than silently imputing values.

## 3.4 Spatial null comparison

The process profiles are inputs to the study's control-energy calculation. The saved analysis compares observed CE with CE obtained from 1,000 spatial spin-null gene-expression maps. A standardized null-relative effect can be summarized as:

```text
SES = (mean(CE_null) - CE_observed) / SD(CE_null)
```

Positive values indicate lower energy for the observed spatial profile. Empirical spin-test significance is evaluated separately. This section provides profile construction; it does not generate spin permutations or rerun the CE solver. Independent random shuffling of regions is not equivalent to a spatial spin model.

## References

- Kang et al. (2011). [Spatio-temporal transcriptome of the human brain](https://doi.org/10.1038/nature10523), Supplementary Table 13.
- Markello et al. (2021). [Standardizing workflows in imaging transcriptomics with the abagen toolbox](https://doi.org/10.7554/eLife.72129).
