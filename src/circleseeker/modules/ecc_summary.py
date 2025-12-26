#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ecc_summary.py - eccDNA Analysis Summary Report Generator
Generates HTML and text reports from eccDNA analysis results
"""

import json
import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime
from typing import Dict, Optional, Tuple, Union
import argparse
import logging
from circleseeker.utils.logging import get_logger
import sys


class EccSummary:
    """Generate comprehensive eccDNA analysis reports."""

    def __init__(self, sample_name: str, output_dir: Path, logger: Optional[logging.Logger] = None):
        """Initialize the summary generator.

        Args:
            sample_name: Sample name for the analysis
            output_dir: Directory to save output files
            logger: Optional logger instance
        """
        self.sample_name = sample_name
        self.output_dir = Path(output_dir)
        self.logger = logger or get_logger(self.__class__.__name__)

        # Data containers
        self.read_stats = {}
        self.ctcr_stats = {}
        self.eccdna_stats = {}
        self.overlap_stats = {}
        self.version = "v2.1.0"

        # Internal bookkeeping for compatibility methods
        self._original_fasta: Optional[Path] = None
        self._processed_csv: Optional[Path] = None

    def collect_read_statistics(self, fasta_path: Path, processed_csv: Path) -> Dict:
        """Collect read classification statistics.

        Args:
            fasta_path: Path to original reads FASTA file
            processed_csv: Path to processed reads CSV file

        Returns:
            Dictionary containing read statistics
        """
        self.logger.info(f"Collecting read statistics from {fasta_path} and {processed_csv}")

        # Count total reads from FASTA
        total_reads = 0
        total_length = 0
        current_len = 0
        try:
            with open(fasta_path, "r") as f:
                for line in f:
                    if line.startswith(">"):
                        if current_len:
                            total_length += current_len
                            current_len = 0
                        total_reads += 1
                    else:
                        current_len += len(line.strip())
            total_length += current_len
            self.logger.info(f"Total reads counted: {total_reads}")
        except FileNotFoundError:
            self.logger.error(f"FASTA file not found: {fasta_path}")
            total_reads = 0
            total_length = 0
        except Exception as e:
            self.logger.error(f"Error reading FASTA file: {e}")
            total_reads = 0
            total_length = 0
        # Read classification from processed CSV
        ctcr_subtypes = {"CtcR-perfect": 0, "CtcR-hybrid": 0, "CtcR-inversion": 0}
        ctcr_total = 0

        try:
            df = pd.read_csv(processed_csv)

            # Count CtcR subtypes
            for subtype in ctcr_subtypes.keys():
                count = len(df[df["readClass"] == subtype])
                ctcr_subtypes[subtype] = count
                ctcr_total += count

            self.logger.info(f"CtcR reads counted: {ctcr_total}")

        except FileNotFoundError:
            self.logger.warning(f"Processed CSV not found: {processed_csv}")
        except Exception as e:
            self.logger.error(f"Error reading processed CSV: {e}")

        # Calculate other reads
        other_reads = total_reads - ctcr_total

        # Store statistics
        self.read_stats = {
            "total_reads": total_reads,
            "total_sequences": total_reads,
            "total_length": total_length,
            "ctcr_reads": ctcr_total,
            "ctcr_percentage": round(ctcr_total / total_reads * 100, 1) if total_reads > 0 else 0,
            "other_reads": other_reads,
            "other_percentage": round(other_reads / total_reads * 100, 1) if total_reads > 0 else 0,
        }

        # Store CtcR subtype statistics
        self.ctcr_stats = {}
        for subtype, count in ctcr_subtypes.items():
            key = subtype.replace("CtcR-", "ctcr_") + "_count"
            self.ctcr_stats[key] = count

            # Percentage of CtcR
            pct_key = subtype.replace("CtcR-", "ctcr_") + "_pct_ctcr"
            self.ctcr_stats[pct_key] = round(count / ctcr_total * 100, 1) if ctcr_total > 0 else 0

            # Percentage of total
            total_pct_key = subtype.replace("CtcR-", "ctcr_") + "_pct_total"
            self.ctcr_stats[total_pct_key] = (
                round(count / total_reads * 100, 1) if total_reads > 0 else 0
            )

        # Additional compatibility keys
        self.ctcr_stats["total_reads"] = total_reads
        self.ctcr_stats["ctcr_reads"] = ctcr_total
        self.ctcr_stats["circular_reads"] = ctcr_subtypes.get("CtcR-perfect", 0)

        return {**self.read_stats, **self.ctcr_stats}

    # ------------------------------------------------------------------
    # Compatibility wrappers for legacy pipeline interface
    # ------------------------------------------------------------------
    def process_fasta(self, fasta_path: Union[Path, str]) -> Dict:
        """Legacy wrapper to record FASTA statistics."""
        path = Path(fasta_path)
        self._original_fasta = path

        total_sequences = 0
        total_length = 0
        current_len = 0

        try:
            with open(path, "r") as handle:
                for line in handle:
                    line = line.strip()
                    if not line:
                        continue
                    if line.startswith(">"):
                        if current_len:
                            total_length += current_len
                            current_len = 0
                        total_sequences += 1
                    else:
                        current_len += len(line)
                total_length += current_len
        except FileNotFoundError:
            self.logger.warning(f"FASTA file not found: {path}")
        except Exception as exc:
            self.logger.error(f"Failed to read FASTA {path}: {exc}")

        average_length = round(total_length / total_sequences, 2) if total_sequences else 0

        self.read_stats.update(
            {
                "total_sequences": total_sequences,
                "total_length": total_length,
                "average_length": average_length,
            }
        )

        return self.read_stats

    def process_processed_csv(self, processed_csv: Union[Path, str]) -> Dict:
        """Legacy wrapper to collect read classification statistics."""
        path = Path(processed_csv)
        self._processed_csv = path

        if self._original_fasta and self._original_fasta.exists():
            return self.collect_read_statistics(self._original_fasta, path)

        self.logger.debug(
            "Processed CSV provided before FASTA; delaying read statistics collection"
        )
        return self.ctcr_stats

    def process_merged_csv(self, merged_csv: Union[Path, str]) -> Dict:
        """Legacy wrapper to collect eccDNA statistics."""
        return self.collect_eccdna_statistics(Path(merged_csv))

    def process_overlap_stats(self, overlap_json: Union[Path, str]) -> Dict:
        """Legacy wrapper to collect overlap statistics."""
        return self.collect_overlap_statistics(Path(overlap_json))

    def collect_eccdna_statistics(self, merged_csv: Path) -> Dict:
        """Collect eccDNA statistics from merged output.

        Args:
            merged_csv: Path to merged eccDNA CSV file

        Returns:
            Dictionary containing eccDNA statistics
        """
        self.logger.info(f"Collecting eccDNA statistics from {merged_csv}")

        try:
            df = pd.read_csv(merged_csv)

            # Function to calculate statistics for a subset
            def calc_stats(data):
                if len(data) == 0:
                    return {
                        "count": 0,
                        "min_length": 0,
                        "max_length": 0,
                        "mean_length": 0,
                        "median_length": 0,
                        "mode_length": 0,
                        "std_length": 0,
                    }

                lengths = data["Length"].values

                # Calculate mode (most frequent length)
                mode_result = pd.Series(lengths).mode()
                mode_val = mode_result.iloc[0] if len(mode_result) > 0 else lengths[0]

                return {
                    "count": len(data),
                    "min_length": int(lengths.min()),
                    "max_length": int(lengths.max()),
                    "mean_length": round(lengths.mean(), 1),
                    "median_length": round(np.median(lengths), 1),
                    "mode_length": int(mode_val),
                    "std_length": round(lengths.std(), 1) if len(lengths) > 1 else 0,
                }

            # Calculate for all eccDNA
            all_stats = calc_stats(df)

            # Calculate for each type
            # Note: UeccDNA* and CeccDNA* include both confirmed and inferred
            uecc_df = df[df["eccDNA_type"] == "UeccDNA"]
            mecc_df = df[df["eccDNA_type"] == "MeccDNA"]
            cecc_df = df[df["eccDNA_type"] == "CeccDNA"]

            uecc_stats = calc_stats(uecc_df)
            mecc_stats = calc_stats(mecc_df)
            cecc_stats = calc_stats(cecc_df)

            # Store with appropriate keys
            self.eccdna_stats = {
                # All eccDNA
                "all_count": all_stats["count"],
                "all_min_length": f"{all_stats['min_length']:,}",
                "all_max_length": f"{all_stats['max_length']:,}",
                "all_mean_length": f"{all_stats['mean_length']:,.1f}",
                "all_median_length": f"{all_stats['median_length']:,.1f}",
                "all_mode_length": f"{all_stats['mode_length']:,}",
                "all_std_length": f"{all_stats['std_length']:,.1f}",
                # UeccDNA*
                "uecc_count": uecc_stats["count"],
                "uecc_min_length": f"{uecc_stats['min_length']:,}",
                "uecc_max_length": f"{uecc_stats['max_length']:,}",
                "uecc_mean_length": f"{uecc_stats['mean_length']:,.1f}",
                "uecc_median_length": f"{uecc_stats['median_length']:,.1f}",
                "uecc_mode_length": f"{uecc_stats['mode_length']:,}",
                "uecc_std_length": f"{uecc_stats['std_length']:,.1f}",
                # MeccDNA
                "mecc_count": mecc_stats["count"],
                "mecc_min_length": f"{mecc_stats['min_length']:,}",
                "mecc_max_length": f"{mecc_stats['max_length']:,}",
                "mecc_mean_length": f"{mecc_stats['mean_length']:,.1f}",
                "mecc_median_length": f"{mecc_stats['median_length']:,.1f}",
                "mecc_mode_length": f"{mecc_stats['mode_length']:,}",
                "mecc_std_length": f"{mecc_stats['std_length']:,.1f}",
                # CeccDNA*
                "cecc_count": cecc_stats["count"],
                "cecc_min_length": f"{cecc_stats['min_length']:,}",
                "cecc_max_length": f"{cecc_stats['max_length']:,}",
                "cecc_mean_length": f"{cecc_stats['mean_length']:,.1f}",
                "cecc_median_length": f"{cecc_stats['median_length']:,.1f}",
                "cecc_mode_length": f"{cecc_stats['mode_length']:,}",
                "cecc_std_length": f"{cecc_stats['std_length']:,.1f}",
            }

            self.logger.info(f"Processed {all_stats['count']} total eccDNA sequences")

        except FileNotFoundError:
            self.logger.error(f"Merged CSV file not found: {merged_csv}")
            self.eccdna_stats = self._get_empty_eccdna_stats()
        except Exception as e:
            self.logger.error(f"Error processing merged CSV: {e}")
            self.eccdna_stats = self._get_empty_eccdna_stats()

        return self.eccdna_stats

    def _get_empty_eccdna_stats(self) -> Dict:
        """Return empty eccDNA statistics structure."""
        stats = {}
        for prefix in ["all", "uecc", "mecc", "cecc"]:
            stats[f"{prefix}_count"] = 0
            for suffix in [
                "min_length",
                "max_length",
                "mean_length",
                "median_length",
                "mode_length",
                "std_length",
            ]:
                stats[f"{prefix}_{suffix}"] = "0"
        return stats

    def collect_overlap_statistics(self, overlap_json: Path) -> Dict:
        """Collect overlap statistics from JSON file.

        Args:
            overlap_json: Path to overlap statistics JSON file

        Returns:
            Dictionary containing overlap statistics
        """
        self.logger.info(f"Collecting overlap statistics from {overlap_json}")

        try:
            with open(overlap_json, "r") as f:
                data = json.load(f)

            # Extract inference statistics
            simple_stats = data.get("inferred_simple", {})
            chimeric_stats = data.get("inferred_chimeric", {})

            self.overlap_stats = {
                "inferred_uecc_count": simple_stats.get("total", 0),
                "uecc_overlap_count": simple_stats.get("overlapping", 0),
                "uecc_overlap_pct": f"{simple_stats.get('overlapping_pct', 0):.2f}",
                "inferred_cecc_count": chimeric_stats.get("total", 0),
                "cecc_overlap_count": chimeric_stats.get("overlapping", 0),
                "cecc_overlap_pct": f"{chimeric_stats.get('overlapping_pct', 0):.2f}",
            }

            self.logger.info("Loaded overlap statistics successfully")

        except FileNotFoundError:
            self.logger.warning(f"Overlap JSON file not found: {overlap_json}")
            self.overlap_stats = {
                "inferred_uecc_count": 0,
                "uecc_overlap_count": 0,
                "uecc_overlap_pct": "0.00",
                "inferred_cecc_count": 0,
                "cecc_overlap_count": 0,
                "cecc_overlap_pct": "0.00",
            }
        except Exception as e:
            self.logger.error(f"Error reading overlap JSON: {e}")
            self.overlap_stats = {
                "inferred_uecc_count": 0,
                "uecc_overlap_count": 0,
                "uecc_overlap_pct": "0.00",
                "inferred_cecc_count": 0,
                "cecc_overlap_count": 0,
                "cecc_overlap_pct": "0.00",
            }

        return self.overlap_stats

    def generate_html_report(self, output_path: Optional[Path] = None) -> str:
        """Generate HTML report.

        Args:
            output_path: Optional path to save HTML file

        Returns:
            HTML content as string
        """
        if output_path is None:
            output_path = self.output_dir / f"{self.sample_name}_report.html"

        # Get current date and time
        analysis_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        # Prepare template data
        template_data = {
            "sample_name": self.sample_name,
            "analysis_date": analysis_date,
            "version": self.version,
            **self.read_stats,
            **self.ctcr_stats,
            **self.eccdna_stats,
            **self.overlap_stats,
        }

        # HTML template
        html_content = self._get_html_template()

        # Replace placeholders
        for key, value in template_data.items():
            # Format numbers with commas for display
            if isinstance(value, int) and key not in ["ctcr_percentage", "other_percentage"]:
                value = f"{value:,}"
            html_content = html_content.replace(f"{{{key}}}", str(value))

        # Save HTML file
        with open(output_path, "w", encoding="utf-8") as f:
            f.write(html_content)

        self.logger.info(f"HTML report saved to {output_path}")
        return html_content

    def generate_text_summary(self, output_path: Optional[Path] = None) -> str:
        """Generate text summary report.

        Args:
            output_path: Optional path to save text file

        Returns:
            Text summary as string
        """
        if output_path is None:
            output_path = self.output_dir / f"{self.sample_name}_summary.txt"

        # Get current date and time
        analysis_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        # Build text report
        lines = []
        lines.append("=" * 80)
        lines.append(f"eccDNA Analysis Report - {self.sample_name}")
        lines.append("=" * 80)
        lines.append("")

        # Read Classification
        lines.append("1. READ CLASSIFICATION")
        lines.append("-" * 22)
        lines.append(f"Total Reads: {self.read_stats['total_reads']:,}")
        lines.append(
            f"CtcR Reads: {self.read_stats['ctcr_reads']:,} ({self.read_stats['ctcr_percentage']}%)"
        )
        perfect = self.ctcr_stats
        lines.append(
            f"  - CtcR-perfect: {perfect['ctcr_perfect_count']:,} "
            f"({perfect['ctcr_perfect_pct_ctcr']}% of CtcR)"
        )
        lines.append(
            f"  - CtcR-hybrid: {perfect['ctcr_hybrid_count']:,} "
            f"({perfect['ctcr_hybrid_pct_ctcr']}% of CtcR)"
        )
        lines.append(
            f"  - CtcR-inversion: {perfect['ctcr_inversion_count']:,} "
            f"({perfect['ctcr_inversion_pct_ctcr']}% of CtcR)"
        )
        lines.append(
            f"Other Reads: {self.read_stats['other_reads']:,} "
            f"({self.read_stats['other_percentage']}%)"
        )
        lines.append("")

        # eccDNA Statistics
        lines.append("2. eccDNA STATISTICS")
        lines.append("-" * 20)
        lines.append("Type        Count    MinLen    MaxLen    MeanLen    MedianLen")
        lines.append("-" * 65)

        # Format eccDNA statistics table
        for prefix, name in [
            ("all", "All"),
            ("uecc", "UeccDNA*"),
            ("mecc", "MeccDNA"),
            ("cecc", "CeccDNA*"),
        ]:
            count = self.eccdna_stats[f"{prefix}_count"]
            min_len = self.eccdna_stats[f"{prefix}_min_length"]
            max_len = self.eccdna_stats[f"{prefix}_max_length"]
            mean_len = self.eccdna_stats[f"{prefix}_mean_length"]
            median_len = self.eccdna_stats[f"{prefix}_median_length"]

            # Format for alignment
            line = f"{name:<11} {count:<8} {min_len:<9} {max_len:<9} {mean_len:<10} {median_len}"
            lines.append(line)

        lines.append("")

        # Inference Analysis
        lines.append("3. INFERENCE ANALYSIS")
        lines.append("-" * 21)
        ov = self.overlap_stats
        lines.append(
            f"Inferred UeccDNA: {ov['inferred_uecc_count']} sequences, "
            f"{ov['uecc_overlap_count']} ({ov['uecc_overlap_pct']}%) overlap"
        )
        lines.append(
            f"Inferred CeccDNA: {ov['inferred_cecc_count']} sequences, "
            f"{ov['cecc_overlap_count']} ({ov['cecc_overlap_pct']}%) overlap"
        )
        lines.append("")

        # Footer
        lines.append(f"Generated by CircleSeeker {self.version}")
        lines.append(f"Date: {analysis_date}")
        lines.append("=" * 80)

        # Join lines
        text_content = "\n".join(lines)

        # Save text file
        with open(output_path, "w", encoding="utf-8") as f:
            f.write(text_content)

        self.logger.info(f"Text summary saved to {output_path}")
        return text_content

    def _get_html_template(self) -> str:
        """Return the HTML template."""
        return """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>eccDNA Analysis Report - {sample_name}</title>
    <style>
        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }

        body {
            font-family: Arial, 'Helvetica Neue', sans-serif;
            background: #fafafa;
            min-height: 100vh;
            padding: 30px 20px;
            line-height: 1.6;
            color: #333;
            font-size: 14px;
        }

        .container {
            max-width: 1200px;
            margin: 0 auto;
            background: white;
            border: 1px solid #e0e0e0;
        }

        h1 {
            background: #f8f9fa;
            color: #2c3e50;
            padding: 30px;
            text-align: center;
            font-size: 24px;
            font-weight: bold;
            border-bottom: 2px solid #dee2e6;
            margin-bottom: 30px;
        }

        h2 {
            color: #34495e;
            padding: 15px 20px;
            font-size: 18px;
            border: 1px solid #dee2e6;
            margin: 30px 40px 25px 40px;
            font-weight: bold;
            background: #f8f9fa;
            border-radius: 5px;
        }

        h2:first-of-type {
            margin-top: 0;
        }

        h3 {
            color: #495057;
            font-size: 14px;
            margin-bottom: 10px;
            font-weight: bold;
            text-transform: uppercase;
        }

        h4 {
            color: #34495e;
            font-size: 16px;
            margin-bottom: 15px;
            font-weight: bold;
        }

        /* Statistics Grid */
        .stats-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
            gap: 20px;
            margin: 0 40px 30px 40px;
        }

        .stat-card {
            background: #fafbfc;
            padding: 20px;
            border: 1px solid #dee2e6;
            text-align: center;
        }

        .stat-card h3 {
            color: #6c757d;
            font-size: 12px;
            margin-bottom: 10px;
            letter-spacing: 0.5px;
        }

        .stat-value {
            font-size: 28px;
            font-weight: bold;
            margin: 8px 0;
        }

        .stat-percentage {
            color: #6c757d;
            font-size: 12px;
        }

        /* Tables */
        table {
            width: calc(100% - 80px);
            margin: 20px 40px;
            border-collapse: collapse;
            border: 1px solid #dee2e6;
            font-size: 14px;
        }

        thead {
            background: #f8f9fa;
        }

        th {
            padding: 12px 15px;
            text-align: center;
            color: #495057;
            font-weight: bold;
            font-size: 14px;
            border-bottom: 2px solid #dee2e6;
        }

        tbody tr {
            background: white;
        }

        tbody tr:nth-child(even) {
            background-color: #fafbfc;
        }

        td {
            padding: 10px 15px;
            border-bottom: 1px solid #e9ecef;
            font-size: 14px;
            color: #495057;
            text-align: center;
        }

        tbody tr:last-child td {
            border-bottom: none;
        }

        /* Summary and Info Boxes */
        .summary-box {
            background: #f8f9fa;
            border-left: 3px solid #adb5bd;
            padding: 15px 20px;
            margin: 20px 40px;
            font-size: 14px;
        }

        .summary-box strong {
            color: #495057;
            font-weight: bold;
        }

        .note-box {
            background: #fffefb;
            border-left: 3px solid #dec181;
            padding: 15px 20px;
            margin: 20px 40px;
            font-style: italic;
            color: #5a5a5a;
            font-size: 14px;
        }

        .note-box strong {
            color: #856404;
            font-style: normal;
            font-weight: bold;
        }

        /* Inference Summary Section - Compact Version */
        .inference-summary {
            background: #fafbfc;
            margin: 30px 40px;
            padding: 20px;
            border: 1px solid #dee2e6;
        }

        .inference-details {
            margin-top: 15px;
            padding: 15px;
            background: white;
            border: 1px solid #e9ecef;
        }

        .inference-item {
            margin-bottom: 8px;
            font-size: 14px;
            color: #495057;
        }

        .inference-item:last-child {
            margin-bottom: 0;
        }

        .inference-item strong {
            color: #495057;
            font-weight: bold;
        }

        .inference-count {
            color: #2c3e50;
            font-weight: bold;
        }

        .inference-note {
            color: #6c757d;
            font-size: 12px;
        }

        /* Footer */
        .footer {
            background: #f8f9fa;
            color: #495057;
            text-align: center;
            padding: 20px;
            border-top: 2px solid #dee2e6;
            font-size: 13px;
        }

        .footer p {
            margin: 3px 0;
        }

        /* Special Markers */
        sup {
            color: #dc3545;
            font-weight: bold;
            font-size: 12px;
        }

        /* Color coding for different values */
        .color-primary { color: #5a6c7d; }
        .color-danger { color: #a94442; }
        .color-secondary { color: #777777; }

        /* Responsive Design */
        @media (max-width: 768px) {
            body { font-size: 13px; }
            h1 { font-size: 20px; padding: 20px; }
            h2 { font-size: 16px; margin: 20px; padding: 15px 20px 10px 20px; }
            .stats-grid { margin: 0 20px 20px 20px; }
            table { margin: 20px; width: calc(100% - 40px); font-size: 13px; }
            .summary-box, .note-box, .inference-summary { margin: 20px; }
        }

        /* Print Styles */
        @media print {
            body { background: white; font-size: 12px; }
            .container { border: none; }
        }
    </style>
</head>
<body>
    <div class="container">
        <h1>🧬 eccDNA Analysis Report - {sample_name}</h1>

        <h2>📊 1. Read Classification Statistics</h2>
        <div class="stats-grid">
            <div class="stat-card">
                <h3>Total Reads</h3>
                <div class="stat-value color-primary">{total_reads}</div>
                <div class="stat-percentage">100% of dataset</div>
            </div>
            <div class="stat-card">
                <h3>CtcR Reads</h3>
                <div class="stat-value color-danger">{ctcr_reads}</div>
                <div class="stat-percentage">{ctcr_percentage}% of total</div>
            </div>
            <div class="stat-card">
                <h3>Other Reads</h3>
                <div class="stat-value color-secondary">{other_reads}</div>
                <div class="stat-percentage">{other_percentage}% of total</div>
            </div>
        </div>

        <h2>🔍 2. CtcR Subtype Distribution</h2>
        <table>
            <thead>
                <tr>
                    <th>CtcR Subtype</th>
                    <th>Count</th>
                    <th>Percentage of CtcR</th>
                    <th>Percentage of Total Reads</th>
                </tr>
            </thead>
            <tbody>
                <tr>
                    <td><strong>CtcR-perfect</strong></td>
                    <td>{ctcr_perfect_count}</td>
                    <td>{ctcr_perfect_pct_ctcr}%</td>
                    <td>{ctcr_perfect_pct_total}%</td>
                </tr>
                <tr>
                    <td><strong>CtcR-hybrid</strong></td>
                    <td>{ctcr_hybrid_count}</td>
                    <td>{ctcr_hybrid_pct_ctcr}%</td>
                    <td>{ctcr_hybrid_pct_total}%</td>
                </tr>
                <tr>
                    <td><strong>CtcR-inversion</strong></td>
                    <td>{ctcr_inversion_count}</td>
                    <td>{ctcr_inversion_pct_ctcr}%</td>
                    <td>{ctcr_inversion_pct_total}%</td>
                </tr>
            </tbody>
        </table>

        <h2>📈 3. eccDNA Summary Statistics</h2>

        <div class="summary-box">
            <strong>🔬 eccDNA Type Definitions:</strong><br>
            <strong>Unique-locus eccDNAs (UeccDNAs)</strong> derive from a single contiguous genomic region and have one unambiguous origin locus.<br>
            <strong>Multi-locus eccDNAs (MeccDNAs)</strong> are composed predominantly of repetitive sequences, resulting in alignments to multiple genomic locations.<br>
            <strong>Chimeric eccDNAs (CeccDNAs)</strong> contain two or more segments from non-contiguous genomic loci, indicating a fusion of distant sequences.
        </div>

        <table>
            <thead>
                <tr>
                    <th>eccDNA Type</th>
                    <th>Count</th>
                    <th>Min Length (bp)</th>
                    <th>Max Length (bp)</th>
                    <th>Mean Length (bp)</th>
                    <th>Median Length (bp)</th>
                    <th>Mode Length (bp)</th>
                    <th>Std Dev</th>
                </tr>
            </thead>
            <tbody>
                <tr>
                    <td><strong>All eccDNA</strong></td>
                    <td>{all_count}</td>
                    <td>{all_min_length}</td>
                    <td>{all_max_length}</td>
                    <td>{all_mean_length}</td>
                    <td>{all_median_length}</td>
                    <td>{all_mode_length}</td>
                    <td>{all_std_length}</td>
                </tr>
                <tr>
                    <td><strong>UeccDNA<sup>*</sup></strong></td>
                    <td>{uecc_count}</td>
                    <td>{uecc_min_length}</td>
                    <td>{uecc_max_length}</td>
                    <td>{uecc_mean_length}</td>
                    <td>{uecc_median_length}</td>
                    <td>{uecc_mode_length}</td>
                    <td>{uecc_std_length}</td>
                </tr>
                <tr>
                    <td><strong>MeccDNA</strong></td>
                    <td>{mecc_count}</td>
                    <td>{mecc_min_length}</td>
                    <td>{mecc_max_length}</td>
                    <td>{mecc_mean_length}</td>
                    <td>{mecc_median_length}</td>
                    <td>{mecc_mode_length}</td>
                    <td>{mecc_std_length}</td>
                </tr>
                <tr>
                    <td><strong>CeccDNA<sup>*</sup></strong></td>
                    <td>{cecc_count}</td>
                    <td>{cecc_min_length}</td>
                    <td>{cecc_max_length}</td>
                    <td>{cecc_mean_length}</td>
                    <td>{cecc_median_length}</td>
                    <td>{cecc_mode_length}</td>
                    <td>{cecc_std_length}</td>
                </tr>
            </tbody>
        </table>

        <div class="note-box">
            <strong>📌 Important Note:</strong> UeccDNA<sup>*</sup> and CeccDNA<sup>*</sup> represent combined counts integrating both confirmed and inferred eccDNA sequences for comprehensive analysis.
        </div>

        <div class="inference-summary">
            <h4>🔮 Inference Analysis Details</h4>
            <div class="inference-details">
                <div class="inference-item">
                    <strong>Inferred UeccDNA:</strong> <span class="inference-count">{inferred_uecc_count}</span> sequences, <span class="inference-note">{uecc_overlap_count} ({uecc_overlap_pct}%) overlap with confirmed UeccDNA</span>
                </div>
                <div class="inference-item">
                    <strong>Inferred CeccDNA:</strong> <span class="inference-count">{inferred_cecc_count}</span> sequences, <span class="inference-note">{cecc_overlap_count} ({cecc_overlap_pct}%) overlap with confirmed CeccDNA</span>
                </div>
            </div>
        </div>

        <div class="footer">
            <p><strong>Report generated by CircleSeeker {version}</strong></p>
            <p>Analysis Date: {analysis_date}</p>
            <p>For questions regarding this analysis, please refer to the CircleSeeker documentation</p>
        </div>
    </div>
</body>
</html>"""

    def process_all(
        self, fasta_path: Path, processed_csv: Path, merged_csv: Path, overlap_json: Path
    ) -> Tuple[str, str]:
        """Process all data and generate reports.

        Args:
            fasta_path: Path to original reads FASTA
            processed_csv: Path to processed reads CSV
            merged_csv: Path to merged eccDNA CSV
            overlap_json: Path to overlap statistics JSON

        Returns:
            Tuple of (HTML content, text content)
        """
        self.logger.info("Starting eccDNA summary report generation")

        # Collect all statistics
        self.collect_read_statistics(fasta_path, processed_csv)
        self.collect_eccdna_statistics(merged_csv)
        self.collect_overlap_statistics(overlap_json)

        # Generate reports
        html_content = self.generate_html_report()
        text_content = self.generate_text_summary()

        self.logger.info("Report generation completed successfully")

        return html_content, text_content


def main():
    """Command-line interface for ecc_summary."""
    parser = argparse.ArgumentParser(
        description="Generate eccDNA analysis summary reports",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Required arguments
    parser.add_argument(
        "-s",
        "--sample-name",
        required=True,
        dest="sample_name",
        help="Sample name for the analysis",
    )

    parser.add_argument("-f", "--fasta", required=True, help="Path to original reads FASTA file")

    parser.add_argument("-p", "--processed", required=True, help="Path to processed reads CSV file")

    parser.add_argument("-m", "--merged", required=True, help="Path to merged eccDNA CSV file")

    parser.add_argument("-j", "--json", required=True, help="Path to overlap statistics JSON file")

    # Optional arguments
    parser.add_argument(
        "-o",
        "--output-dir",
        default=".",
        help="Output directory for reports (default: current directory)",
    )

    parser.add_argument("--html-only", action="store_true", help="Generate only HTML report")

    parser.add_argument("--text-only", action="store_true", help="Generate only text summary")

    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")

    args = parser.parse_args()

    # Setup logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(level=log_level, format="%(asctime)s - %(levelname)s - %(message)s")

    # Convert paths
    fasta_path = Path(args.fasta)
    processed_csv = Path(args.processed)
    merged_csv = Path(args.merged)
    overlap_json = Path(args.json)
    output_dir = Path(args.output_dir)

    # Validate input files
    for file_path, name in [
        (fasta_path, "FASTA"),
        (processed_csv, "Processed CSV"),
        (merged_csv, "Merged CSV"),
        (overlap_json, "Overlap JSON"),
    ]:
        if not file_path.exists():
            logging.error(f"{name} file not found: {file_path}")
            sys.exit(1)

    # Create output directory if needed
    output_dir.mkdir(parents=True, exist_ok=True)

    # Initialize summary generator
    summary = EccSummary(args.sample_name, output_dir)

    try:
        # Process all data
        summary.process_all(fasta_path, processed_csv, merged_csv, overlap_json)

        # Generate requested reports
        if not args.text_only:
            html_path = output_dir / f"{args.sample_name}_report.html"
            logging.info(f"HTML report saved: {html_path}")

        if not args.html_only:
            text_path = output_dir / f"{args.sample_name}_summary.txt"
            logging.info(f"Text summary saved: {text_path}")

        logging.info("✅ Summary report generation completed successfully!")

    except Exception as e:
        logging.error(f"Error during report generation: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
