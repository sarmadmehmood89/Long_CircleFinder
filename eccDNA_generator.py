import random
import os
import argparse
import sys

def ensure_directory(directory):
    """Ensure that a directory exists; if it doesn't, create it."""
    if not os.path.exists(directory):
        os.makedirs(directory)

def read_genome(fasta_file, chromosome):
    """Read a specific chromosome from a FASTA file and return its sequence."""
    sequence = ''
    current_chromosome = None
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if line[1:] == chromosome:
                    current_chromosome = chromosome
                else:
                    current_chromosome = None
            elif current_chromosome:
                sequence += line
    return sequence

def read_bed_file(bed_file, chromosome):
    """Read a BED file and return a list of tuples with start and end positions to avoid for a specific chromosome."""
    excluded_regions = []
    with open(bed_file, 'r') as file:
        for line in file:
            parts = line.strip().split()
            if parts[0] == chromosome:
                start = int(parts[1])
                end = int(parts[2])
                excluded_regions.append((start, end))
    return excluded_regions

def is_valid_fragment(start, length, excluded_regions, used_regions, genome_length):
    """Check if the fragment overlaps with any excluded or previously used region."""
    end = start + length
    if end > genome_length or length < 1000:
        return False
    for region in excluded_regions + used_regions:
        ex_start, ex_end = region
        if start < ex_end and end > ex_start:
            return False
    return True

def pick_random_fragments(chromosome, genome_sequence, num_fragments, excluded_regions):
    """Pick multiple random fragments from a chromosome while avoiding excluded and overlapping regions."""
    fragments = []
    used_regions = []  # List to track regions already used for fragments
    genome_length = len(genome_sequence)
    attempts = 0
    max_attempts = 1000
    while len(fragments) < num_fragments and attempts < max_attempts:
        length = int(random.gauss(150000, 50000))  # Gaussian distribution with mean 150,000 and std 50,000
        length = max(1000, length)  # Ensure length is at least 1000 base pairs
        start = random.randint(0, genome_length - length)
        if is_valid_fragment(start, length, excluded_regions, used_regions, genome_length):
            fragment_sequence = genome_sequence[start:start + length]
            fragments.append((chromosome, start, start + length - 1, length, fragment_sequence))
            used_regions.append((start, start + length - 1))
        attempts += 1
    return fragments

def save_to_fasta(sequence, header, directory, filename):
    """Save the given DNA sequence to a FASTA file in the specified directory."""
    with open(os.path.join(directory, filename), 'w') as f:
        f.write(f'>{header}\n')
        f.write(sequence + '\n')

def main():
    parser = argparse.ArgumentParser(description="Generate eccDNA sequences from randomly chosen fragments across chromosomes.")
    parser.add_argument("--genome_file", required=True, help="Path to the genome FASTA file")
    parser.add_argument("--bed_file", required=True, help="Path to the BED file")
    parser.add_argument("--num_fragments_per_chromosome", type=int, default=5, help="Number of fragments per chromosome to generate")
    parser.add_argument("--total_eccDNA", type=int, default=10, help="Total number of eccDNA molecules to generate")
    parser.add_argument("--output_directory", default="eccDNA_sequences", help="Output directory to save eccDNA sequences and coordinates")

    args = parser.parse_args()

    chromosomes = ['chr' + str(i) for i in range(1, 22)] + ['chrX']
    ensure_directory(args.output_directory)

    chromosome_fragments = {}
    for chromosome in chromosomes:
        print(f"Processing {chromosome}...")
        sequence = read_genome(args.genome_file, chromosome)
        excluded_regions = read_bed_file(args.bed_file, chromosome)
        chromosome_fragments[chromosome] = pick_random_fragments(chromosome, sequence, args.num_fragments_per_chromosome, excluded_regions)

    with open(os.path.join(args.output_directory, "eccDNA_coordinates.txt"), 'w') as coord_file:
        for i in range(args.total_eccDNA):
            eccdna_fragments = [random.choice(fragments) for fragments in chromosome_fragments.values()]
            eccdna_fragments = random.sample(eccdna_fragments, k=random.randint(2, 22))  # Pick 2-22 fragments
            header = f'eccDNA_B{len(eccdna_fragments)}_{i+1} Total Size: {sum(frag[3] for frag in eccdna_fragments)}bp'
            coord_file.write(f"# {header} Coordinates:\n")
            for frag in eccdna_fragments:
                coord_file.write(f"{frag[0]}\t{frag[1]}\t{frag[2]}\tFragment Size: {frag[3]}bp\n")
            coord_file.write("\n")  # Space between eccDNA entries

            # Save individual FASTA file for each eccDNA
            fasta_filename = f"eccDNA_B{len(eccdna_fragments)}_{i+1}.fasta"
            eccdna_sequence = ''.join(frag[4] for frag in eccdna_fragments)
            save_to_fasta(eccdna_sequence, header, args.output_directory, fasta_filename)

    print("All eccDNA generation complete.")

if __name__ == "__main__":
    main()

