
# Long-CircleFinder

**Long-CircleFinder** is a tool for detecting and characterizing extrachromosomal circular DNA (eccDNA) from long-read sequencing data. This tool combines sequencing depth analysis, supplementary read extraction, graph-based cycle detection, and Partial Order Alignment (POA) to identify and refine eccDNA structures.

## Features
- Detects **overamplified regions** based on sequencing depth.
- Extracts **supplementary reads** to establish structural connections.
- Constructs a **graph** and identifies cycles indicative of eccDNA.
- Uses **POA** to generate accurate consensus sequences and refine breakpoints.

## Installation

### Requirements
- **Python >= 3 
- **Pysam**: `pip install pysam`
- **Minimap2** and **SAMtools**: Install via package manager or from source.
- **Wtdbg2** (includes **WTPOA-CNS**): [Install from source](https://github.com/ruanjue/wtdbg2).

### Clone the Repository
```bash
git clone https://github.com/your-username/Long-CircleFinder.git
cd Long-CircleFinder
pip install -r requirements.txt
```

## Usage
```bash
python long_circlefinder.py --bam <path_to_bam_file> --output_dir <output_directory> --reference_genome <path_to_reference_genome> --min_support <min_support_value> --depth <sequencing_depth>
```

### Example
```bash
python long_circlefinder.py --bam data/sample_data.bam --output_dir results/ --reference_genome reference.fa --min_support 3 --depth 30
```

## Outputs
- **Overamplified Regions**: `<output_dir>/overamplified_regions.txt`
- **Supplementary Reads**: `<output_dir>/supplementary_reads.txt`
- **Detected Cycles**: `<output_dir>/cycles.txt`
- **Consensus Sequences**: `<output_dir>/consensus_sequences.fa`
- **Refined Breakpoints**: `<output_dir>/refined_breakpoints.txt`

## eccDNA Generator

The `eccDNA_generator.py` script generates eccDNA sequences by selecting random DNA fragments from specified chromosomes in a genome, while avoiding regions specified in a BED file. This tool can be used to simulate eccDNA structures by combining random fragments into contiguous sequences.

### Usage
```bash
python eccDNA_generator.py --genome_file <path_to_genome_fasta> --bed_file <path_to_bed_file> --num_fragments_per_chromosome 5 --total_eccDNA 10 --output_directory eccDNA_sequences

## Simulate Reads

The `simulate_reads.py` script simulates sequencing reads from eccDNA FASTA files using **NanoSim**. It generates synthetic reads for each eccDNA sequence based on specified coverage levels, useful for testing eccDNA analysis workflows.

### Usage
```bash
python simulate_reads.py --input_dir <path_to_eccDNA_FASTA_files> --output_dir <output_directory> --model_prefix <path_to_NanoSim_model_directory> --read_length <read_length>



## License
This project is licensed under the MIT License.

## Citation
If you use **Long-CircleFinder** in your research, please cite:

