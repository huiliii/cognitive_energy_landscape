# 3. Gene Expression Analysis

This example shows how regional AHBA expression is combined into six neurodevelopmental process profiles and scaled into regional control-input weights. The default run uses synthetic data and requires no AHBA download. It is an educational example, not a complete reproduction of the manuscript's control-energy or spin-test analysis.

## Run the toy example

Requires Python 3, NumPy, and pandas:

```bash
pip install numpy pandas
python CODE/Gene_Expression/gene_expression_example.py --output-dir gene_expression_output
```

The example creates 12 synthetic regions and two synthetic genes per process. The synthetic genes are unrelated to the accompanying real gene lists.

## Use a preprocessed AHBA matrix

Provide a CSV whose first column contains atlas region labels and whose remaining columns contain numeric gene-expression values, with gene symbols as headers:

```bash
python CODE/Gene_Expression/gene_expression_example.py \
  --expression-csv /path/to/regional_expression.csv \
  --region-labels /path/to/connectome_region_labels.csv \
  --output-dir ahba_process_profiles
```

The optional region-label file has a `label` column in the exact order used by the structural connectomes. Select regions by their atlas labels rather than assuming row positions. The manuscript MATLAB script retains rows `[1:28, 33:94]` from its original 120-row AAL2 matrix (MATLAB indexing); verify the atlas labels before applying that selection to a different input.

## Optional AHBA extraction

Install the study's reported abagen version in a compatible Python environment:

```bash
pip install abagen==0.1.3
python CODE/Gene_Expression/gene_expression_example.py \
  --atlas /path/to/AAL2_mni.nii.gz \
  --region-labels /path/to/connectome_region_labels.csv \
  --output-dir ahba_process_profiles
```

This downloads public adult postmortem AHBA measurements. Extraction explicitly specifies probe reannotation, intensity filtering at 0.5, differential-stability probe selection, corrected MNI coordinates, 2-mm matching tolerance, scaled robust sigmoid normalization across genes and samples, and averaging within donors followed by averaging across donors. No hemisphere mirroring, missing-region interpolation, or atlas-info table is requested. A processing report is saved. Consult the [abagen parameter documentation](https://abagen.readthedocs.io/en/stable/generated/abagen.get_expression_data.html) for exact behavior.

## Gene sets and aggregation

`kang2011_gene_sets.csv` contains the six ranked process-associated lists transcribed from Supplementary Table 13 of [Kang et al. (2011)](https://doi.org/10.1038/nature10523), using the source table retained with the study analysis. The table lists genes correlated with process-associated markers. These lists should not be described as exhaustive pathway annotations. Original symbols are retained; no automatic alias conversion is performed.

For each process, the script:

1. Removes duplicate gene symbols within that process and matches symbols exactly against the expression matrix.
2. Records matched and unmatched symbols and gene counts.
3. Averages matched gene expression within each retained region.
4. Divides all six profiles by one common factor: the regional mean of the expression profile averaged over the **unique union** of matched genes across the six sets.

Thus, for expression matrix `X`, matched gene set `G_p`, and their union `U`:

```python
profile_p = X[sorted(G_p)].mean(axis=1)
common_scale = X[sorted(U)].mean(axis=1).mean()
weights_p = profile_p / common_scale
```

Gene sets can overlap; a gene may contribute to multiple process profiles but is counted once in the union. Adult-derived spatial profiles are fixed across developmental periods. They are not age-specific expression measurements.

## Outputs and reproducibility boundaries

- `process_profiles.csv`: regional means before common scaling.
- `control_weights.csv`: the six commonly scaled regional weight vectors.
- `gene_matching.csv` and `gene_counts.csv`: retained and excluded symbols and counts.
- `union_profile.csv` and `run_info.txt`: normalization reference and run details.
- `illustrative_diagonal_B.npz`: diagonal input matrices illustrating one conventional NCT configuration. These are not used to calculate manuscript CE in this example.

Missing or nonfinite expression values trigger an error. Any interpolation, exclusion, or contralateral substitution must be explicitly implemented and documented before running the aggregation. This example does not silently reproduce undocumented edits to the historical matrix or impute missing regions.

Spatial spin permutations, their empirical p-values, and CE calculations are not generated here. Randomly shuffling regional rows is not a substitute for a spatially constrained spin null model.

## References

- Kang et al. (2011). [Spatio-temporal transcriptome of the human brain](https://doi.org/10.1038/nature10523). Source gene lists: Supplementary Table 13.
- Markello et al. (2021). [Standardizing workflows in imaging transcriptomics with the abagen toolbox](https://doi.org/10.7554/eLife.72129).
