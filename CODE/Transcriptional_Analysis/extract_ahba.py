"""Extract regional AHBA expression with the study's default abagen workflow."""
import argparse
from pathlib import Path
import abagen


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--atlas', required=True, type=Path,
                        help='Path to the volumetric AAL2 atlas in MNI space.')
    parser.add_argument('--output-dir', type=Path, default=Path('transcriptional_output'))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    # Original analysis: abagen.get_expression_data(AAL2_atlas).
    # return_report adds provenance without changing the extraction parameters.
    expression, report = abagen.get_expression_data(str(args.atlas), return_report=True)
    expression.index.name = 'label'
    expression.to_csv(args.output_dir / 'abagen_expression.csv')
    (args.output_dir / 'abagen_report.txt').write_text(report)
    print(f'Saved {expression.shape[0]} regions and {expression.shape[1]} genes.')


if __name__ == '__main__':
    main()
