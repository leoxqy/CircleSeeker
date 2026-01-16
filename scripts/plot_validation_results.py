#!/usr/bin/env python3
"""Plot batch validation results for CircleSeeker."""

import json
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# Set style
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 11
plt.rcParams['axes.titlesize'] = 13
plt.rcParams['axes.labelsize'] = 12

def load_results(json_path: str) -> dict:
    """Load validation results from JSON file."""
    with open(json_path, 'r') as f:
        return json.load(f)

def plot_metrics_by_scale(data: dict, output_dir: Path):
    """Plot Recall, Precision, F1 by scale."""
    by_scale = data['aggregate']['by_scale']
    scales = sorted([int(k) for k in by_scale.keys()])

    recalls = [by_scale[str(s)]['recall'] * 100 for s in scales]
    precisions = [by_scale[str(s)]['precision'] * 100 for s in scales]
    f1s = [by_scale[str(s)]['f1'] * 100 for s in scales]

    fig, ax = plt.subplots(figsize=(10, 6))

    x = np.arange(len(scales))
    width = 0.25

    bars1 = ax.bar(x - width, recalls, width, label='Recall', color='#2ecc71', edgecolor='black', linewidth=0.5)
    bars2 = ax.bar(x, precisions, width, label='Precision', color='#3498db', edgecolor='black', linewidth=0.5)
    bars3 = ax.bar(x + width, f1s, width, label='F1 Score', color='#e74c3c', edgecolor='black', linewidth=0.5)

    ax.set_xlabel('Number of eccDNAs')
    ax.set_ylabel('Percentage (%)')
    ax.set_title('CircleSeeker v0.10.3 Validation: Performance Metrics by Scale')
    ax.set_xticks(x)
    ax.set_xticklabels([f'{s:,}' for s in scales], rotation=45, ha='right')
    ax.legend(loc='lower left')
    ax.set_ylim(99.9, 100.02)

    # Add value labels on bars
    for bars in [bars1, bars2, bars3]:
        for bar in bars:
            height = bar.get_height()
            if height < 100:
                ax.annotate(f'{height:.3f}',
                           xy=(bar.get_x() + bar.get_width() / 2, height),
                           xytext=(0, 3),
                           textcoords="offset points",
                           ha='center', va='bottom', fontsize=8, rotation=90)

    plt.tight_layout()
    plt.savefig(output_dir / 'metrics_by_scale.png', dpi=150, bbox_inches='tight')
    plt.savefig(output_dir / 'metrics_by_scale.pdf', bbox_inches='tight')
    plt.close()
    print(f"Saved: metrics_by_scale.png/pdf")

def plot_fp_fn_by_scale(data: dict, output_dir: Path):
    """Plot FP and FN counts by scale."""
    by_scale = data['aggregate']['by_scale']
    scales = sorted([int(k) for k in by_scale.keys()])

    fps = [by_scale[str(s)]['fp'] for s in scales]
    fns = [by_scale[str(s)]['fn'] for s in scales]
    avg_fps = [by_scale[str(s)]['avg_fp'] for s in scales]
    avg_fns = [by_scale[str(s)]['avg_fn'] for s in scales]

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Left: Total FP/FN
    ax1 = axes[0]
    x = np.arange(len(scales))
    width = 0.35

    bars1 = ax1.bar(x - width/2, fps, width, label='False Positives', color='#e74c3c', edgecolor='black', linewidth=0.5)
    bars2 = ax1.bar(x + width/2, fns, width, label='False Negatives', color='#3498db', edgecolor='black', linewidth=0.5)

    ax1.set_xlabel('Number of eccDNAs')
    ax1.set_ylabel('Total Count (10 runs)')
    ax1.set_title('Total FP/FN Across 10 Runs per Scale')
    ax1.set_xticks(x)
    ax1.set_xticklabels([f'{s:,}' for s in scales], rotation=45, ha='right')
    ax1.legend()

    # Add value labels
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            if height > 0:
                ax1.annotate(f'{int(height)}',
                           xy=(bar.get_x() + bar.get_width() / 2, height),
                           xytext=(0, 3),
                           textcoords="offset points",
                           ha='center', va='bottom', fontsize=9)

    # Right: Average FP/FN per run
    ax2 = axes[1]

    bars3 = ax2.bar(x - width/2, avg_fps, width, label='Avg FP per Run', color='#e74c3c', alpha=0.7, edgecolor='black', linewidth=0.5)
    bars4 = ax2.bar(x + width/2, avg_fns, width, label='Avg FN per Run', color='#3498db', alpha=0.7, edgecolor='black', linewidth=0.5)

    ax2.set_xlabel('Number of eccDNAs')
    ax2.set_ylabel('Average Count per Run')
    ax2.set_title('Average FP/FN per Run')
    ax2.set_xticks(x)
    ax2.set_xticklabels([f'{s:,}' for s in scales], rotation=45, ha='right')
    ax2.legend()

    plt.tight_layout()
    plt.savefig(output_dir / 'fp_fn_by_scale.png', dpi=150, bbox_inches='tight')
    plt.savefig(output_dir / 'fp_fn_by_scale.pdf', bbox_inches='tight')
    plt.close()
    print(f"Saved: fp_fn_by_scale.png/pdf")

def plot_accuracy_heatmap(data: dict, output_dir: Path):
    """Plot accuracy as a summary figure."""
    by_scale = data['aggregate']['by_scale']
    scales = sorted([int(k) for k in by_scale.keys()])
    overall = data['aggregate']['overall']

    fig, ax = plt.subplots(figsize=(12, 6))

    # Create data for plotting
    accuracies = [by_scale[str(s)]['accuracy'] * 100 for s in scales]
    tps = [by_scale[str(s)]['tp'] for s in scales]
    gts = [by_scale[str(s)]['gt'] for s in scales]

    # Bar plot
    colors = ['#27ae60' if acc == 100 else '#f39c12' if acc >= 99.99 else '#e74c3c' for acc in accuracies]
    bars = ax.bar(range(len(scales)), accuracies, color=colors, edgecolor='black', linewidth=0.5)

    ax.set_xlabel('Number of eccDNAs')
    ax.set_ylabel('Accuracy (%)')
    ax.set_title(f'CircleSeeker v0.10.3 Validation Summary\n'
                 f'Overall: {overall["accuracy"]*100:.4f}% Accuracy | '
                 f'Recall: {overall["recall"]*100:.4f}% | '
                 f'Precision: {overall["precision"]*100:.4f}%')
    ax.set_xticks(range(len(scales)))
    ax.set_xticklabels([f'{s:,}' for s in scales], rotation=45, ha='right')
    ax.set_ylim(99.9, 100.05)

    # Add annotations
    for i, (bar, acc, tp, gt) in enumerate(zip(bars, accuracies, tps, gts)):
        ax.annotate(f'{acc:.4f}%\n({tp:,}/{gt:,})',
                   xy=(bar.get_x() + bar.get_width() / 2, bar.get_height()),
                   xytext=(0, 5),
                   textcoords="offset points",
                   ha='center', va='bottom', fontsize=9)

    # Add legend for colors
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#27ae60', edgecolor='black', label='100% Accuracy'),
        Patch(facecolor='#f39c12', edgecolor='black', label='≥99.99% Accuracy'),
    ]
    ax.legend(handles=legend_elements, loc='lower right')

    plt.tight_layout()
    plt.savefig(output_dir / 'accuracy_summary.png', dpi=150, bbox_inches='tight')
    plt.savefig(output_dir / 'accuracy_summary.pdf', bbox_inches='tight')
    plt.close()
    print(f"Saved: accuracy_summary.png/pdf")

def plot_combined_figure(data: dict, output_dir: Path):
    """Create a combined publication-ready figure."""
    by_scale = data['aggregate']['by_scale']
    scales = sorted([int(k) for k in by_scale.keys()])
    overall = data['aggregate']['overall']

    fig = plt.figure(figsize=(14, 10))

    # Create grid
    gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.25)

    # Panel A: Accuracy by scale
    ax1 = fig.add_subplot(gs[0, 0])
    accuracies = [by_scale[str(s)]['accuracy'] * 100 for s in scales]
    colors = ['#27ae60' if acc == 100 else '#f39c12' for acc in accuracies]
    bars = ax1.bar(range(len(scales)), accuracies, color=colors, edgecolor='black', linewidth=0.5)
    ax1.set_xlabel('Number of eccDNAs')
    ax1.set_ylabel('Accuracy (%)')
    ax1.set_title('A. Detection Accuracy by Scale', fontweight='bold', loc='left')
    ax1.set_xticks(range(len(scales)))
    ax1.set_xticklabels([f'{s//1000}k' if s >= 1000 else str(s) for s in scales], rotation=45, ha='right')
    ax1.set_ylim(99.9, 100.02)
    ax1.axhline(y=100, color='gray', linestyle='--', alpha=0.5)

    # Panel B: Recall, Precision, F1
    ax2 = fig.add_subplot(gs[0, 1])
    recalls = [by_scale[str(s)]['recall'] * 100 for s in scales]
    precisions = [by_scale[str(s)]['precision'] * 100 for s in scales]
    f1s = [by_scale[str(s)]['f1'] * 100 for s in scales]

    x = np.arange(len(scales))
    width = 0.25
    ax2.bar(x - width, recalls, width, label='Recall', color='#2ecc71', edgecolor='black', linewidth=0.5)
    ax2.bar(x, precisions, width, label='Precision', color='#3498db', edgecolor='black', linewidth=0.5)
    ax2.bar(x + width, f1s, width, label='F1', color='#9b59b6', edgecolor='black', linewidth=0.5)
    ax2.set_xlabel('Number of eccDNAs')
    ax2.set_ylabel('Percentage (%)')
    ax2.set_title('B. Performance Metrics', fontweight='bold', loc='left')
    ax2.set_xticks(x)
    ax2.set_xticklabels([f'{s//1000}k' if s >= 1000 else str(s) for s in scales], rotation=45, ha='right')
    ax2.legend(loc='lower left', fontsize=9)
    ax2.set_ylim(99.9, 100.02)

    # Panel C: FP/FN counts
    ax3 = fig.add_subplot(gs[1, 0])
    fps = [by_scale[str(s)]['fp'] for s in scales]
    fns = [by_scale[str(s)]['fn'] for s in scales]

    width = 0.35
    ax3.bar(x - width/2, fps, width, label='False Positives', color='#e74c3c', edgecolor='black', linewidth=0.5)
    ax3.bar(x + width/2, fns, width, label='False Negatives', color='#3498db', edgecolor='black', linewidth=0.5)
    ax3.set_xlabel('Number of eccDNAs')
    ax3.set_ylabel('Count (Total of 10 Runs)')
    ax3.set_title('C. Error Distribution', fontweight='bold', loc='left')
    ax3.set_xticks(x)
    ax3.set_xticklabels([f'{s//1000}k' if s >= 1000 else str(s) for s in scales], rotation=45, ha='right')
    ax3.legend(loc='upper left', fontsize=9)

    # Panel D: Summary statistics as text
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.axis('off')

    summary_text = f"""
    Validation Summary
    ══════════════════════════════════════

    Configuration:
    • Scales tested: 100 to 100,000 eccDNAs
    • Runs per scale: 10 (seeds 42-51)
    • eccDNA composition: 70% Unique, 15% Multimap, 15% Cluster
    • Total ground truth: {overall['gt']:,} eccDNAs

    Overall Performance:
    • Successful runs: 80/80 (100%)
    • True Positives: {overall['tp']:,}
    • False Positives: {overall['fp']} (avg {overall['avg_fp']:.2f}/run)
    • False Negatives: {overall['fn']} (avg {overall['avg_fn']:.2f}/run)

    • Accuracy:  {overall['accuracy']*100:.4f}%
    • Recall:    {overall['recall']*100:.4f}%
    • Precision: {overall['precision']*100:.4f}%
    • F1 Score:  {overall['f1']*100:.4f}%

    ══════════════════════════════════════
    CircleSeeker v0.10.3 | minimap2 v2.30
    """

    ax4.text(0.05, 0.95, summary_text, transform=ax4.transAxes,
             fontsize=10, verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='#f8f9fa', edgecolor='#dee2e6'))
    ax4.set_title('D. Summary Statistics', fontweight='bold', loc='left')

    plt.savefig(output_dir / 'validation_combined.png', dpi=200, bbox_inches='tight', facecolor='white')
    plt.savefig(output_dir / 'validation_combined.pdf', bbox_inches='tight')
    plt.close()
    print(f"Saved: validation_combined.png/pdf")

def main():
    script_dir = Path(__file__).parent
    results_file = script_dir / 'batch_validation_results.json'
    output_dir = script_dir

    if not results_file.exists():
        print(f"Error: Results file not found: {results_file}")
        return

    print(f"Loading results from: {results_file}")
    data = load_results(results_file)

    print("Generating plots...")
    plot_metrics_by_scale(data, output_dir)
    plot_fp_fn_by_scale(data, output_dir)
    plot_accuracy_summary(data, output_dir)
    plot_combined_figure(data, output_dir)

    print("\nAll plots generated successfully!")

def plot_accuracy_summary(data: dict, output_dir: Path):
    """Alias for plot_accuracy_heatmap."""
    plot_accuracy_heatmap(data, output_dir)

if __name__ == '__main__':
    main()
