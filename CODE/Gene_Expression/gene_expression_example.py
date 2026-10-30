"""Build process-specific regional profiles. Default run uses synthetic data only."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

PROCESSES = ['Myelination', 'Axon development', 'Dendrite development',
             'Neuron differentiation', 'Neuron migration', 'Synapse development']


def extract_ahba(atlas, output):
    """Optional real-data extraction; downloads AHBA data via abagen."""
    import abagen
    expression, report = abagen.get_expression_data(
        str(atlas), ibf_threshold=0.5, probe_selection='diff_stability',
        donor_probes='aggregate', sim_threshold=None, lr_mirror=None,
        missing=None, tolerance=2, sample_norm='srs', gene_norm='srs',
        norm_matched=True, norm_structures=False, region_agg='donors',
        agg_metric='mean', corrected_mni=True, reannotated=True,
        donors='all', return_report=True)
    expression.index.name = 'label'
    expression.to_csv(output / 'regional_expression.csv')
    (output / 'abagen_report.txt').write_text(report)
    return expression


def synthetic_inputs():
    """Small deterministic illustration; no biological measurements or spin tests."""
    rng = np.random.default_rng(42)
    expression = pd.DataFrame(rng.uniform(.1, .9, (12, 12)),
                              index=np.arange(1, 13),
                              columns=[f'TOY_GENE_{i}' for i in range(12)])
    genes = pd.DataFrame([(p, f'TOY_GENE_{2*i+j}')
                          for i, p in enumerate(PROCESSES) for j in range(2)],
                         columns=['process', 'gene'])
    return expression, genes


def build_profiles(expression, genes):
    if expression.index.has_duplicates or expression.columns.has_duplicates:
        raise ValueError('Region labels and gene names must be unique.')
    genes = genes[['process', 'gene']].dropna().drop_duplicates()
    genes['matched'] = genes.gene.isin(expression.columns)
    profiles = pd.DataFrame(index=expression.index)
    selected = set()
    for process in PROCESSES:
        matched = sorted(set(genes.loc[(genes.process == process) & genes.matched, 'gene']))
        if not matched:
            raise ValueError(f'No matching genes for {process}')
        values = expression[matched].to_numpy(dtype=float)
        if not np.isfinite(values).all():
            raise ValueError(f'Missing/nonfinite expression for {process}. '
                             'Resolve and document missing-region handling explicitly.')
        profiles[process] = values.mean(axis=1)
        selected.update(matched)
    # Union counts a gene once, even if it belongs to several process sets.
    union_profile = expression[sorted(selected)].mean(axis=1)
    scale = float(union_profile.mean())
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError('Common scaling factor must be finite and positive.')
    weights = profiles / scale
    return profiles, weights, genes, union_profile, scale


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group()
    source.add_argument('--expression-csv', type=Path,
                        help='Preprocessed region-by-gene CSV; first column is region ID.')
    source.add_argument('--atlas', type=Path, help='MNI volumetric atlas for optional AHBA extraction.')
    parser.add_argument('--gene-sets', type=Path,
                        default=Path(__file__).with_name('kang2011_gene_sets.csv'))
    parser.add_argument('--region-labels', type=Path,
                        help='CSV with a label column in desired connectome order.')
    parser.add_argument('--output-dir', type=Path, default=Path('gene_expression_output'))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    if args.expression_csv or args.atlas:
        genes = pd.read_csv(args.gene_sets)
        expression = (pd.read_csv(args.expression_csv, index_col=0) if args.expression_csv
                      else extract_ahba(args.atlas, args.output_dir))
        mode = 'AHBA input'
    else:
        expression, genes = synthetic_inputs()
        mode = 'SYNTHETIC TOY DATA'
    if args.region_labels:
        labels = pd.read_csv(args.region_labels)['label'].tolist()
        expression = expression.loc[labels]
    profiles, weights, matches, union, scale = build_profiles(expression, genes)
    profiles.to_csv(args.output_dir / 'process_profiles.csv', index_label='label')
    weights.to_csv(args.output_dir / 'control_weights.csv', index_label='label')
    matches.to_csv(args.output_dir / 'gene_matching.csv', index=False)
    union.to_csv(args.output_dir / 'union_profile.csv', index_label='label')
    counts = matches.groupby('process', sort=False).agg(
        listed_genes=('gene', 'size'), matched_genes=('matched', 'sum'))
    counts.to_csv(args.output_dir / 'gene_counts.csv')
    # Each column contains regional weights; diagonal B is illustrated only.
    matrices = {f'B_{i+1}': np.diag(weights[p]) for i, p in enumerate(PROCESSES)}
    np.savez(args.output_dir / 'illustrative_diagonal_B.npz', **matrices)
    (args.output_dir / 'run_info.txt').write_text(
        f'Mode: {mode}\nRegions: {len(expression)}\nCommon scale: {scale:.12g}\n'
        'B_1..B_6 follow the process order in gene_counts.csv.\n'
        'The diagonal B matrices illustrate a conventional NCT input configuration; '
        'they do not reproduce the study control-energy computation.\n')
    print(mode)
    print(counts.to_string())
    print(f'Common scale: {scale:.6g}; outputs: {args.output_dir}')


if __name__ == '__main__':
    main()
