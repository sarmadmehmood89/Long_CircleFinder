
# Long-CircleFinder

**Long-CircleFinder** is a computational tool for detecting and characterizing extrachromosomal circular DNA (eccDNA) from long-read sequencing data. The tool identifies overamplified regions in the genome, extracts supplementary reads, constructs a graph to identify eccDNA cycles, and refines structural information using Partial Order Alignment (POA). This makes it ideal for studying eccDNA’s role in oncogene amplification, genomic instability, and drug resistance.

## Table of Contents

- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Scripts](#scripts)
  - [long_circlefinder.py](#long_circlefinderpy)
  - [eccDNA_generator.py](#eccdna_generatorpy)
  - [simulate_reads.py](#simulate_readspy)
- [Outputs](#outputs)
- [Example Commands](#example-commands)
- [License](#license)
- [Citation](#citation)

---

## Features

- Detects **overamplified regions** in long-read sequencing data based on sequencing depth.
- Extracts **supplementary reads** that may indicate structural rearrangements.
- Constructs a **graph** and identifies cycles representing eccDNA.
- Utilizes **POA** to generate accurate consensus sequences and refine breakpoints.

## Installation

### Requirements

To run **Long-CircleFinder**, you will need:

- **Python >= 3.7**
- **Pysam**: for handling SAM/BAM files.
  ```bash
  pip install pysam
  ```
- **Minimap2**: a long-read aligner, used for consensus sequence alignment.
  - Install via package manager or [from source](https://github.com/lh3/minimap2).
- **SAMtools**: for BAM file processing.
  - Install via package manager or [from source](https://github.com/samtools/samtools).
- **Wtdbg2 and WTPOA-CNS**: for POA and consensus generation.
  - Install [from source](https://github.com/ruanjue/wtdbg2).

### Installation Steps

1. Clone this repository:
   ```bash
   git clone https://github.com/your-username/Long-CircleFinder.git
   cd Long-CircleFinder
   ```

2. Install Python dependencies:
   ```bash
   pip install -r requirements.txt
   ```

3. Ensure that **Minimap2**, **SAMtools**, **Wtdbg2**, and **WTPOA-CNS** are accessible from your system PATH.

## Usage

```bash
python long_circlefinder.py --bam <path_to_bam_file> --output_dir <output_directory> --reference_genome <path_to_reference_genome> --min_support <min_support_value> --depth <sequencing_depth>
```

### Arguments
- `--bam`: Path to the BAM file containing aligned long-read sequencing data.
- `--output_dir`: Directory to store output files.
- `--reference_genome`: Path to the reference genome in FASTA format.
- `--min_support`: Minimum support value for supplementary reads.
- `--depth`: Sequencing depth (optional; calculated automatically if not provided).

## Scripts

### `long_circlefinder.py`

The **main script** for detecting eccDNA using long-read sequencing data. It performs the following tasks:
1. Identifies **overamplified regions** based on sequencing depth.
2. Extracts **supplementary reads** with split alignments.
3. Constructs a **graph** of connected regions.
4. Detects cycles in the graph, suggesting eccDNA structures.
5. Uses **POA** to generate consensus sequences and refines breakpoints.

### `eccDNA_generator.py`

Generates random DNA fragments from a genome FASTA file to simulate eccDNA. It reads specified chromosomes, avoids regions from a BED file, and generates FASTA files for each simulated eccDNA sequence.

**Usage:**
```bash
python eccDNA_generator.py --genome_file <path_to_genome_fasta> --bed_file <path_to_bed_file> --num_fragments_per_chromosome <number> --total_eccDNA <total_molecules> --output_directory <output_dir>
```

### `simulate_reads.py`

Simulates sequencing reads for eccDNA FASTA files using **NanoSim**. It generates synthetic reads based on coverage levels and the specified read length.

**Usage:**
```bash
python simulate_reads.py --input_dir <path_to_eccDNA_fasta_files> --output_dir <output_dir> --model_prefix <path_to_NanoSim_model> --read_length <read_length>
```

## Outputs

Long-CircleFinder generates the following files in the output directory:

1. **Overamplified Regions**: A list of genomic regions with high sequencing depth.
   - File: `<output_dir>/overamplified_regions.txt`
2. **Supplementary Reads**: Details of extracted supplementary reads.
   - File: `<output_dir>/supplementary_reads.txt`
3. **Detected Cycles (eccDNA Candidates)**: Graph cycles indicating potential eccDNA structures.
   - Files:
     - Cycles: `<output_dir>/cycles.txt`
     - Graph: `<output_dir>/graph.txt`
4. **Consensus Sequences**: Consensus sequences for each detected eccDNA.
   - File: `<output_dir>/consensus_sequences.fa`
5. **Refined Breakpoints**: Information on the exact breakpoints for each eccDNA.
   - File: `<output_dir>/refined_breakpoints.txt`

## Example Commands

### Running the Main Detection Script

```bash
python long_circlefinder.py --bam data/sample_data.bam --output_dir results --reference_genome reference.fa --min_support 3 --depth 30
```

### Generating Random Fragments for eccDNA Simulation

```bash
python eccDNA_generator.py --genome_file reference_genome.fa --bed_file regions_to_avoid.bed --num_fragments_per_chromosome 5 --total_eccDNA 10 --output_directory eccDNA_sequences
```

### Simulating Reads from eccDNA FASTA Files

```bash
python simulate_reads.py --input_dir eccDNA_sequences --output_dir simulated_reads --model_prefix /path/to/NanoSim/model --read_length 20000
```

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for more details.

## Citation

If you use **Long-CircleFinder** in your research, please cite:


