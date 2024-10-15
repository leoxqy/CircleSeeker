# CircleSeeker: Advanced Circular DNA Analysis Tool

CircleSeeker is a specialized bioinformatics tool designed for advanced circular DNA analysis using PacBio HiFi sequencing data. This comprehensive software processes high-throughput long-read sequences to efficiently identify and characterize various types of extrachromosomal circular DNA (eccDNA) structures, leveraging the high accuracy and long read lengths of PacBio HiFi technology.

## Features

- Efficient processing of PacBio HiFi sequencing data
- Identification and classification of various eccDNA types (UeccDNA, MeccDNA, CeccDNA, XeccDNA)
- Integration with popular bioinformatics tools (TideHunter, BLAST, minimap2, samtools)
- Multi-threaded processing for improved performance
- Comprehensive analysis report generation

## Installation

CircleSeeker can be installed using conda:

```bash
conda install -c bioconda circle_seeker
```

Alternatively, you can install it from source:

```bash
git clone https://github.com/leoxqy/CircleSeeker.git
cd CircleSeeker
pip install .
```

## Dependencies

CircleSeeker requires the following tools and libraries:

- Python 3.6+
- numpy >=1.22
- pandas
- pysam
- biopython
- minimap2
- samtools
- BLAST
- bedtools
- seqkit
- TideHunter
- pbtk

These dependencies are automatically managed when installing via conda.

## Usage

Basic usage of CircleSeeker:

```bash
CircleSeeker -i <input_fasta> -p <output_prefix> -r <reference_genome> [-o <output_folder>] [-t <num_threads>]
```

Arguments:
- `-i, --input`: Input FASTA file (required)
- `-p, --prefix`: Prefix for all generated files (required)
- `-r, --reference`: Reference genome FASTA file (required)
- `-o, --output`: Output folder (optional)
- `-t, --num_threads`: Number of threads to use (default: 8)

For a detailed workflow, please refer to the [CircleSeeker Tutorial](CircleSeeker完整工作流程教程.md).

## Output

CircleSeeker generates several output files, including:

- Confirmed eccDNA sequences (UeccDNA, MeccDNA, CeccDNA, XeccDNA)
- CSV files with detailed information about identified eccDNA
- A comprehensive analysis report

## License

CircleSeeker is licensed under the GNU General Public License v2 (GPLv2).

## Authors

- Yaoxin Zhang (yxzhang@ncgr.com)
- Leo Xinqi Yu (leoxqy@hotmail.com)

## Citation

If you use CircleSeeker in your research, please cite:

[Placeholder for citation information]

## Support

For issues, bug reports, or feature requests, please use the [GitHub issue tracker](https://github.com/leoxqy/CircleSeeker/issues).

## Contributing

We welcome contributions to CircleSeeker. Please refer to our [contribution guidelines](CONTRIBUTING.md) for more information.

## Acknowledgments

We would like to thank the developers of the tools and libraries that CircleSeeker depends on, as well as the research community for their valuable feedback and support.
