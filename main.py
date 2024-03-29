import pysam
from statistics import mean, stdev
from collections import defaultdict
import logging
import sys

# Constants
MIN_REGION_LENGTH = 20000
STD_DEV_MULTIPLIER = 2  # Multiplier for the standard deviation to define overamplification
KNOWN_CHROMOSOMES = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(message)s')

def get_chromosomes_from_bam(bam_file):
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        chromosomes = [chrom for chrom in bam.references if chrom in KNOWN_CHROMOSOMES]
    return chromosomes


def calculate_coverage(bam_file, chromosome):
    coverage = defaultdict(int)
    with pysam.AlignmentFile(bam_file, "rb") as file:
        for read in file.fetch(chromosome):
            # Exclude reads with mapping quality <= 0
            if read.mapping_quality <= 0:
                continue
            for position in range(read.reference_start, read.reference_end):
                coverage[position] += 1
    return coverage

def identify_overamplified_regions(chromosome, coverage):
    coverage_values = list(coverage.values())
    if not coverage_values:
        return []

    coverage_mean = mean(coverage_values)
    coverage_std = stdev(coverage_values)
    overamplified_threshold = coverage_mean + (STD_DEV_MULTIPLIER * coverage_std)

    regions = []
    start = None
    for position in sorted(coverage.keys()):
        if coverage[position] > overamplified_threshold and coverage[position] >= 1:
            if start is None:
                start = position
        else:
            if start is not None:
                end_position = position - 1
                if all(coverage.get(p, 0) >= 1 for p in range(start, end_position + 1)):
                    if end_position - start >= MIN_REGION_LENGTH:
                        average_coverage = mean([coverage[p] for p in range(start, end_position + 1)])
                        regions.append((chromosome, start, end_position, average_coverage))
                start = None

    if start is not None and all(coverage.get(p, 0) >= 1 for p in range(start, max(coverage.keys()) + 1)):
        end_position = max(coverage.keys())
        if end_position - start >= MIN_REGION_LENGTH:
            average_coverage = mean([coverage[p] for p in range(start, end_position + 1)])
            regions.append((chromosome, start, end_position, average_coverage))

    return regions

def get_reads_with_supplementary_alignments(bam_file):
    supplementary_reads = []
    with pysam.AlignmentFile(bam_file, "rb") as file:
        for read in file.fetch(until_eof=True):
            # Exclude secondary alignments (flag 256)
            if read.flag & 256:
                continue
            # Include reads that have an SA tag, this will include primary and supplementary alignments
            if read.has_tag('SA'):
                supplementary_reads.append(read)
    return supplementary_reads

def map_regions_to_indices(all_overamplified_regions):
    region_to_index = {}
    index_to_region = {}
    idx = 0
    for region in all_overamplified_regions:
        key = (region[0], region[1], region[2])  # chromosome, start, end as key
        region_to_index[key] = idx
        index_to_region[idx] = key
        idx += 1
    return region_to_index, index_to_region

def generate_graph_from_supplementary_alignments(all_overamplified_regions, supplementary_reads, region_to_index):
    G = defaultdict(set)
    for read in supplementary_reads:
        for region in all_overamplified_regions:
            if read.reference_start in range(region[1], region[2] + 1) and read.reference_name == region[0]:
                primary_idx = region_to_index[(region[0], region[1], region[2])]
                sa_tags = read.get_tag('SA').strip(';').split(';')
                for sa_tag in sa_tags:
                    sa_parts = sa_tag.split(',')
                    if len(sa_parts) >= 2:
                        sa_chrom, sa_pos = sa_parts[0], int(sa_parts[1])
                        sa_region = next((reg for reg in all_overamplified_regions if sa_chrom == reg[0] and sa_pos in range(reg[1], reg[2] + 1)), None)
                        if sa_region:
                            sa_idx = region_to_index[(sa_region[0], sa_region[1], sa_region[2])]
                            if primary_idx != sa_idx:
                                G[primary_idx].add(sa_idx)
                                G[sa_idx].add(primary_idx)
    return G

def dfs_visit(node, graph, path, visited, all_cycles):
    if node in path:
        cycle_start_index = path.index(node)
        all_cycles.append(path[cycle_start_index:] + [node])
        return
    if node in visited:
        return
    visited.add(node)
    path.append(node)
    for neighbor in graph[node]:
        dfs_visit(neighbor, graph, path.copy(), visited, all_cycles)

def find_cycles(graph):
    visited = set()
    all_cycles = []
    for node in graph:
        if node not in visited:
            dfs_visit(node, graph, [], visited, all_cycles)
    return all_cycles

def main(bam_file):
    chromosomes = get_chromosomes_from_bam(bam_file)
    all_overamplified_regions = []
    for chromosome in chromosomes:
        coverage = calculate_coverage(bam_file, chromosome)
        overamplified_regions = identify_overamplified_regions(chromosome, coverage)
        all_overamplified_regions.extend(overamplified_regions)

    # Logging overamplified regions
    logging.info("Identified overamplified regions:")
    for region in all_overamplified_regions:
        logging.info(f"Chromosome {region[0]}: Region {region[1]}-{region[2]}, Average Coverage: {region[3]}")

    supplementary_reads = get_reads_with_supplementary_alignments(bam_file)
    region_to_index, index_to_region = map_regions_to_indices(all_overamplified_regions)
    G = generate_graph_from_supplementary_alignments(all_overamplified_regions, supplementary_reads, region_to_index)

    # Assuming print_graph_with_coordinates is defined to log the graph structure
    #print_graph_with_coordinates(G, index_to_region)  # Uncomment if implemented

    cycles = find_cycles(G)
    logging.info(f"Total number of cycles detected: {len(cycles)}")
    for i, cycle in enumerate(cycles, 1):
        cycle_nodes = " -> ".join(str(index_to_region[node]) for node in cycle)
        logging.info(f"Cycle {i}: {cycle_nodes}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        logging.error("Usage: python script.py <BAM_FILE>")
        sys.exit(1)
    bam_file_path = sys.argv[1]
    main(bam_file_path)
