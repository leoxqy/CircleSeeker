# CircleSeeker: Advanced Circular DNA Analysis for PacBio HiFi Data

[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
[![Project Status: Active](https://img.shields.io/badge/status-active-success.svg)](https://github.com/leoxqy/CircleSeeker/)

**CircleSeeker** is a specialized bioinformatics pipeline designed for the identification, classification, and characterization of extrachromosomal circular DNA (eccDNA) from PacBio HiFi long-read sequencing data.

---

## Overview

Traditional methods for eccDNA identification are often platform-specific. Approaches for Next-Generation Sequencing (NGS) data primarily rely on **Split-Read** reconstruction algorithms, while earlier methods for Nanopore long reads focused on directly identifying tandem repeat structures (**CTC reads**).

**CircleSeeker** innovatively combines the strengths of these two strategies, specifically optimized for the unique advantages of PacBio HiFi data—long read lengths and high accuracy. It first identifies potential circular structures by detecting tandem repeats within reads and then employs a Split-Read based strategy to precisely reconstruct large eccDNAs (>20 kb). This hybrid approach ensures the accurate and efficient recovery of complete eccDNA sequences.

Following identification, CircleSeeker automatically classifies each eccDNA based on its alignment pattern and chimerism, providing a powerful tool for investigating genome instability.

## Key Features

-   **Hybrid Detection Strategy**: Merges tandem repeat detection with Split-Read assembly, tailored to leverage the strengths of PacBio HiFi data for maximum recall and precision.
-   **Comprehensive eccDNA Classification**: Sorts eccDNA into biologically meaningful categories:
    -   **UeccDNA (Unique-locus eccDNA):** Originates from a single, contiguous genomic locus.
    -   **MeccDNA (Multi-locus eccDNA):** Composed of repetitive sequences (e.g., tandem repeats) that map to multiple genomic locations, but lacks a chimeric structure.
    -   **CeccDNA (Chimeric eccDNA):** A chimeric circle formed by the fusion of two or more distant, non-contiguous genomic segments.
    -   **MCeccDNA (Multi-locus Chimeric eccDNA):** A hybrid circle exhibiting both multi-locus and chimeric features.
    -   **XeccDNA:** Sequences that cannot be classified into the above categories, potentially representing contamination, other unknown circular molecules, or complex unclassified structures.
-   **Automated & Efficient Pipeline**: Integrates multiple third-party tools (e.g., TideHunter, BLAST, minimap2, samtools) into a single, cohesive workflow with multi-threading support.
-   **Checkpoint & Resumption**: Creates a checkpoint file after each major step, allowing the pipeline to be resumed from the last completed point in case of interruption.
-   **Detailed Analysis Reports**: Automatically generates an HTML summary report, comprehensive CSV tables, and final FASTA files for easy interpretation and downstream analysis.

## Installation

We currently recommend installing CircleSeeker by downloading a release package and installing it with Conda. More installation options will be available in the future.

### Option 1: From Release File (Recommended)

1.  **Navigate to the Releases Page**: Visit the project's [Releases page](https://github.com/leoxqy/CircleSeeker/releases).
2.  **Download the Package**: Download the latest `.tar.bz2` package (e.g., `circle_seeker-1.0.1-pyh2b07374_0.tar.bz2`).
3.  **Install with Conda**: Open your terminal, navigate to the download directory, and run:
    ```bash
    # Replace the filename with the version you downloaded
    conda install circle_seeker-1.0.1-pyh2b07374_0.tar.bz2
    ```

### Option 2: From Source (For Developers)

1.  **Clone the Repository**:
    ```bash
    git clone [https://github.com/leoxqy/CircleSeeker.git](https://github.com/leoxqy/CircleSeeker.git)
    cd CircleSeeker
    ```
2.  **Install with pip**:
    ```bash
    pip install .
    ```
    **Important**: If installing from source, ensure the following dependencies are in your system's `$PATH`:
    `minimap2`, `samtools`, `blast+`, `seqkit`, `mosdepth`, `bcftools`, `tidehunter`.

## Usage

### Basic Command
```bash
CircleSeeker \
    -i <input.fasta> \
    -p <output_prefix> \
    -r <reference.fasta> \
    -t <threads>
```

### Arguments

#### Required Arguments
| Argument            | Description                              |
| :------------------ | :--------------------------------------- |
| `-i`, `--input`     | Input PacBio HiFi reads in FASTA format. |
| `-p`, `--prefix`    | Prefix for all output files.             |
| `-r`, `--reference` | Reference genome in FASTA format.        |

#### Optional Arguments
| Argument              | Description                                  | Default           |
| :-------------------- | :------------------------------------------- | :---------------- |
| `-o`, `--output`      | Output directory path.                       | Current directory |
| `-t`, `--num_threads` | Number of threads to use.                    | `8`               |
| `--enable_X`          | Enable the detection of XeccDNA.             | `False`           |
| `--force`             | Force re-run and overwrite existing outputs. | `False`           |
| `--start_from`        | Skip to a specific step number to start.     | `None`            |

## Example Workflow

1.  **Prepare Input Files**
    -   PacBio HiFi Reads: `HeLa_rep1.fasta`
    -   Reference Genome: `hg38.p13.fa`

2.  **Run CircleSeeker**
    ```bash
    CircleSeeker \
        -i HeLa_rep1.fasta \
        -p HeLa_rep1 \
        -r hg38.p13.fa \
        -t 16 \
        --enable_X
    ```

3.  **Check the Results**
    -   Review the output files in the specified directory.
    -   Open `HeLa_rep1.report.html` for a visual summary of the analysis.
    -   Confirm successful completion by checking the `HeLa_rep1.checkpoint` file.

## Output Files

The pipeline generates several output files. The most important ones include:

-   `<prefix>.Final.Confirmed.{Uecc,Mecc,Cecc,MCecc}.{csv,fasta}`: Tables and sequences for eccDNA types confirmed by coverage analysis.
-   `<prefix>.Final.Confirmed.Xecc.fasta`: Confirmed XeccDNA sequences (if `--enable_X` is used).
-   `<prefix>.Final.Inferred.Uecc.{csv,fasta}`: **UeccDNA inferred using the Split-Reads method.**
-   `<prefix>.report.{html,txt}`: A formatted summary report of the analysis.
-   `<prefix>.summary.csv`: A master table containing statistics on read classifications.
-   `<prefix>.checkpoint`: The checkpoint file used for pipeline resumption.

## Citation

If you use CircleSeeker in your research, please cite our paper:

> [Placeholder for citation information - to be updated upon publication]

## Contributing

We welcome contributions and new feature ideas. To contribute:
1.  Fork this repository and create a new branch.
2.  Make your changes and test them thoroughly.
3.  Submit a Pull Request for review.

Please see `CONTRIBUTING.md` for more details.

## License

CircleSeeker is licensed under the [GNU General Public License v2 (GPLv2)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html).

## Authors & Support

-   **Yaoxin Zhang** (yxzhang@ncgr.com)
-   **Leo Xinqi Yu** (leoxqy@hotmail.com)

For any questions, bug reports, or feature requests, please [open an issue](https://github.com/leoxqy/CircleSeeker/issues) on GitHub.