#!/usr/bin/env python3
"""
End-to-End Validation Pipeline for CircleSeeker

This script:
1. Generates simulated eccDNA data (or uses existing)
2. Runs CircleSeeker pipeline on the simulated data
3. Compares results with ground truth
4. Reports recall, precision, F1, FP, FN rates

Usage:
    python run_validation.py                    # Run full pipeline
    python run_validation.py --skip-simulation  # Use existing simulation data
    python run_validation.py --skip-pipeline    # Use existing pipeline results
"""

import os
import sys
import json
import subprocess
from pathlib import Path
from datetime import datetime

# Add project root to path
PROJECT_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from eccdna_simulator import SimulationConfig, run_simulation
from validation_metrics import run_validation


def run_circleseeker_pipeline(
    reads_fasta: Path,
    reference_fasta: Path,
    output_dir: Path,
    threads: int = 4,
    minimap2_preset: str = "sr",
) -> bool:
    """
    Run CircleSeeker pipeline on simulated data.

    Note: This is a simplified version focused on the U/M/C classification path.
    The main CircleSeeker pipeline uses minimap2 twice:
      - `minimap2_align` (PAF -> TSV) for candidate alignment (default preset: `sr`)
      - `minimap2` (SAM/BAM) for inference/Cyrcular (default preset: `map-hifi`)
    This simulation validator only needs the former to compute recall/precision.
    """
    print("\n" + "=" * 60)
    print("Running CircleSeeker Pipeline")
    print("=" * 60)

    output_dir.mkdir(parents=True, exist_ok=True)

    # For simulation validation, we run the alignment and classification steps
    # The full pipeline requires TideHunter which needs real HiFi reads

    try:
        # Step 1: Run minimap2 alignment
        print("\n[1/3] Running minimap2 alignment...")
        alignment_tsv = output_dir / "alignment.tsv"

        # Build minimap2 command
        cmd = [
            "minimap2",
            "-x", minimap2_preset,
            "-N", "200",
            "--secondary=yes",
            "-c",  # Output CIGAR
            "--cs",  # Output cs tag
            "-t", str(threads),
            str(reference_fasta),
            str(reads_fasta),
        ]

        # Run alignment and convert PAF to TSV format
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=600
        )

        if result.returncode != 0:
            print(f"  Warning: minimap2 failed: {result.stderr}")
            # Create empty results for validation
            _create_empty_results(output_dir)
            return False

        # Parse PAF output and convert to alignment TSV
        paf_output = result.stdout
        _convert_paf_to_alignment_tsv(paf_output, alignment_tsv)
        print(f"  Alignment written to: {alignment_tsv}")

        # Step 2: Run UM classification
        print("\n[2/3] Running U/M classification...")
        from circleseeker.modules.um_classify import UMeccClassifier

        classifier = UMeccClassifier()
        uecc_df, mecc_df, unclassified_df = classifier.run(alignment_tsv)

        # Save results
        um_output_dir = output_dir / "um_classify"
        um_output_dir.mkdir(parents=True, exist_ok=True)

        if not uecc_df.empty:
            uecc_path = output_dir / "uecc.csv"
            uecc_df.to_csv(uecc_path, index=False)
            print(f"  UeccDNA: {uecc_df['query_id'].nunique()} (saved to {uecc_path})")
        else:
            print("  UeccDNA: 0")

        if not mecc_df.empty:
            mecc_path = output_dir / "mecc.csv"
            mecc_df.to_csv(mecc_path, index=False)
            print(f"  MeccDNA: {mecc_df['query_id'].nunique()} (saved to {mecc_path})")
        else:
            print("  MeccDNA: 0")

        if not unclassified_df.empty:
            unclass_path = output_dir / "unclassified.csv"
            unclassified_df.to_csv(unclass_path, index=False)
            print(f"  Unclassified: {unclassified_df['query_id'].nunique()}")

        # Step 3: Run Cecc detection on unclassified
        print("\n[3/3] Running CeccDNA detection...")
        if not unclassified_df.empty:
            from circleseeker.modules.cecc_build import CeccBuild

            unclass_csv = output_dir / "unclassified.csv"
            if unclass_csv.exists():
                cecc_path = output_dir / "cecc.csv"

                try:
                    runner = CeccBuild()
                    cecc_df = runner.run_pipeline(
                        input_csv=unclass_csv,
                        output_csv=cecc_path,
                    )
                    if cecc_df is not None and not cecc_df.empty:
                        print(f"  CeccDNA: {cecc_df['query_id'].nunique()} (saved to {cecc_path})")
                    else:
                        print("  CeccDNA: 0")
                except Exception as e:
                    print(f"  CeccDNA detection error: {e}")
                    import traceback
                    traceback.print_exc()
        else:
            print("  CeccDNA: 0 (no unclassified reads)")

        print("\nPipeline completed!")
        return True

    except subprocess.TimeoutExpired:
        print("  Error: Pipeline timeout")
        return False
    except Exception as e:
        print(f"  Error running pipeline: {e}")
        import traceback
        traceback.print_exc()
        return False


def _convert_paf_to_alignment_tsv(paf_output: str, output_path: Path):
    """Convert PAF format to alignment TSV (BLAST outfmt 6 style)."""
    import pandas as pd

    records = []
    for line in paf_output.strip().split('\n'):
        if not line:
            continue

        fields = line.split('\t')
        if len(fields) < 12:
            continue

        # PAF format:
        # 0: query_id, 1: query_len, 2: q_start, 3: q_end,
        # 4: strand, 5: target, 6: target_len, 7: t_start, 8: t_end,
        # 9: matches, 10: block_len, 11: mapq

        query_id = fields[0]
        query_len = int(fields[1])
        q_start = int(fields[2]) + 1  # Convert to 1-based
        q_end = int(fields[3])
        strand = fields[4]
        subject_id = fields[5]
        t_start = int(fields[7]) + 1  # Convert to 1-based
        t_end = int(fields[8])
        matches = int(fields[9])
        alignment_length = int(fields[10])
        mapq = int(fields[11]) if len(fields) > 11 else 0

        # Calculate identity
        identity = (matches / alignment_length * 100) if alignment_length > 0 else 0
        mismatches = alignment_length - matches

        # For BLAST format, negative strand has s_start > s_end
        # PAF always has t_start < t_end regardless of strand
        if strand == "+":
            s_start = t_start
            s_end = t_end
            sstrand = "plus"
        else:
            # Swap for negative strand (BLAST convention)
            s_start = t_end
            s_end = t_start
            sstrand = "minus"

        records.append([
            query_id, subject_id, identity, alignment_length,
            mismatches, 0,  # gap_opens
            q_start, q_end, s_start, s_end,
            0, 0,  # evalue, bit_score
            sstrand,
            mapq,
        ])

    # Write to TSV
    columns = [
        "query_id", "subject_id", "identity", "alignment_length",
        "mismatches", "gap_opens", "q_start", "q_end",
        "s_start", "s_end", "evalue", "bit_score", "sstrand", "mapq"
    ]

    df = pd.DataFrame(records, columns=columns)
    df.to_csv(output_path, sep='\t', index=False, header=False)


def _create_empty_results(output_dir: Path):
    """Create empty result files for validation when pipeline fails."""
    import pandas as pd

    columns = ["eccDNA_id", "eccdna_type", "reads", "chr", "start0", "end0",
               "strand", "length", "copy_number"]

    for name in ["uecc.csv", "mecc.csv", "cecc.csv"]:
        pd.DataFrame(columns=columns).to_csv(output_dir / name, index=False)


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Run end-to-end validation for CircleSeeker"
    )
    parser.add_argument(
        "--skip-simulation",
        action="store_true",
        help="Skip simulation, use existing data"
    )
    parser.add_argument(
        "--skip-pipeline",
        action="store_true",
        help="Skip pipeline, use existing results"
    )
    parser.add_argument(
        "--simulation-dir",
        type=Path,
        default=Path("tests/simulation/simulation_data"),
        help="Simulation data directory"
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=Path("tests/simulation/pipeline_results"),
        help="Pipeline results directory"
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Number of threads"
    )
    parser.add_argument(
        "--minimap2-preset",
        type=str,
        default="sr",
        help="minimap2 preset to use for simulation validation (default: sr)",
    )
    parser.add_argument(
        "--num-uecc",
        type=int,
        default=100,
        help="Number of UeccDNA to simulate"
    )
    parser.add_argument(
        "--num-mecc",
        type=int,
        default=20,
        help="Number of MeccDNA to simulate"
    )
    parser.add_argument(
        "--num-cecc",
        type=int,
        default=20,
        help="Number of CeccDNA to simulate"
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for simulation reproducibility (default: 42)",
    )

    args = parser.parse_args()

    # Change to project root
    os.chdir(PROJECT_ROOT)

    print("=" * 60)
    print("CircleSeeker Validation Pipeline")
    print(f"Started: {datetime.now().isoformat()}")
    print("=" * 60)

    # Step 1: Run simulation
    if not args.skip_simulation:
        config = SimulationConfig(
            output_dir=str(args.simulation_dir),
            num_uecc=args.num_uecc,
            num_mecc=args.num_mecc,
            num_cecc=args.num_cecc,
            seed=args.seed,
        )
        run_simulation(config)
    else:
        print("\n[Skipping simulation - using existing data]")
        if not args.simulation_dir.exists():
            print(f"Error: Simulation directory not found: {args.simulation_dir}")
            sys.exit(1)

    # Step 2: Run pipeline
    if not args.skip_pipeline:
        reads_fasta = args.simulation_dir / "simulated_reads.fa"
        reference_fasta = args.simulation_dir / "reference.fa"

        success = run_circleseeker_pipeline(
            reads_fasta=reads_fasta,
            reference_fasta=reference_fasta,
            output_dir=args.results_dir,
            threads=args.threads,
            minimap2_preset=args.minimap2_preset,
        )

        if not success:
            print("\nWarning: Pipeline had issues. Validation may be incomplete.")
    else:
        print("\n[Skipping pipeline - using existing results]")
        if not args.results_dir.exists():
            print(f"Error: Results directory not found: {args.results_dir}")
            sys.exit(1)

    # Step 3: Run validation
    report_path = args.results_dir / "validation_report.json"
    run_validation(args.simulation_dir, args.results_dir, report_path)

    print("\n" + "=" * 60)
    print(f"Validation Complete: {datetime.now().isoformat()}")
    print("=" * 60)


if __name__ == "__main__":
    main()
