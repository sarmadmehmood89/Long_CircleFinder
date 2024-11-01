
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

## License
This project is licensed under the MIT License.

## Citation
If you use **Long-CircleFinder** in your research, please cite:

