import argparse
import os
import pysam
import subprocess
import time
from collections import defaultdict

# Logging helpers
def log_time_start(message):
    start_time = time.time()
    start_time_str = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(start_time))
    print(f"{message} started at {start_time_str}")
    return start_time

def log_time_end(start_time, message):
    end_time = time.time()
    end_time_str = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(end_time))
    elapsed_time = end_time - start_time
    print(f"{message} ended at {end_time_str} - Elapsed Time: {elapsed_time:.2f} seconds")

# Function to estimate sequencing depth and calculate min_supp
def estimate_depth_and_minsupp(bam_file_path, min_support, depth, output_dir):
    start_time = log_time_start("Estimating sequencing depth and calculating min_supp")

    if not min_support:
        if depth:
            minsupp = round(depth / 10.0) + 2
        else:
            mapped_length = 0
            reflength = sum(pysam.AlignmentFile(bam_file_path, 'rb').lengths)
            with pysam.AlignmentFile(bam_file_path, 'rb') as bamfile:
                for read in bamfile.fetch(until_eof=True):
                    if not read.is_unmapped:
                        mapped_length += read.query_length
            depth = mapped_length / reflength
            minsupp = round(depth / 10.0) + 2
    else:
        minsupp = min_support
        if not depth:
            mapped_length = 0
            reflength = sum(pysam.AlignmentFile(bam_file_path, 'rb').lengths)
            with pysam.AlignmentFile(bam_file_path, 'rb') as bamfile:
                for read in bamfile.fetch(until_eof=True):
                    if not read.is_unmapped:
                        mapped_length += read.query_length
            depth = mapped_length / reflength

    log_time_end(start_time, "Estimating sequencing depth and calculating min_supp")

    # Write depth and minsupp to a file
    depth_file = os.path.join(output_dir, f"{os.path.basename(bam_file_path)}_depth_minsupp.txt")
    with open(depth_file, 'w') as f:
        f.write(f"Estimated Depth: {depth}\n")
        f.write(f"Min Supp: {minsupp}\n")

    return depth, minsupp

# Function to identify overamplified regions based on 2x the average depth
def identify_overamplified_regions(bam_file, avg_coverage, threshold_multiplier=2, min_region_length=1000):
    start_time = log_time_start("Identifying overamplified regions")
    
    overamplified_regions = []
    threshold = threshold_multiplier * avg_coverage

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        current_region_start = None

        for pileupcolumn in bam.pileup():
            chrom = pileupcolumn.reference_name
            pos = pileupcolumn.pos
            coverage = pileupcolumn.n

            if coverage >= threshold:
                if current_region_start is None:
                    current_region_start = pos
            else:
                if current_region_start is not None:
                    region_length = pos - current_region_start
                    if region_length >= min_region_length:
                        overamplified_regions.append((chrom, current_region_start, pos - 1))
                    current_region_start = None

        if current_region_start is not None:
            region_length = pileupcolumn.pos - current_region_start
            if region_length >= min_region_length:
                overamplified_regions.append((pileupcolumn.reference_name, current_region_start, pileupcolumn.pos))

    log_time_end(start_time, "Identifying overamplified regions")
    return overamplified_regions

# Function to extract supplementary reads from BAM file
def get_supplementary_reads(bam_file_path):
    start_time = log_time_start("Extracting supplementary reads")

    supplementary_reads = []
    bam = pysam.AlignmentFile(bam_file_path, "rb")
    header = bam.header

    for read in bam.fetch():
        if read.has_tag('SA'):
            supplementary_reads.append(read.to_dict())  # Convert read to dict for serialization

    bam.close()

    log_time_end(start_time, "Extracting supplementary reads")
    return supplementary_reads, header

# Function to generate a graph from supplementary alignments with at least 3 supporting supplementary reads
def generate_graph_from_supplementary_alignments(supplementary_reads, overamplified_regions, region_to_index, bam_header):
    start_time = log_time_start("Generating graph from supplementary alignments")

    G = defaultdict(lambda: defaultdict(lambda: {'count': 0, 'orientations': [], 'connections': []}))
    for read_dict in supplementary_reads:
        read = pysam.AlignedSegment.from_dict(read_dict, bam_header)
        primary_region = next((reg for reg in overamplified_regions if read.reference_name == reg[0] and read.reference_start in range(reg[1], reg[2] + 1)), None)
        if primary_region:
            primary_idx = region_to_index[primary_region]
            sa_tags = read.get_tag('SA').strip(';').split(';')
            for sa_tag in sa_tags:
                sa_parts = sa_tag.split(',')
                sa_chrom, sa_pos = sa_parts[0], int(sa_parts[1])
                sa_region = next((reg for reg in overamplified_regions if sa_chrom == reg[0] and int(sa_pos) in range(reg[1], reg[2] + 1)), None)
                if sa_region:
                    sa_idx = region_to_index[sa_region]
                    connection_detail = f"{sa_chrom}:{sa_pos}"
                    if primary_idx != sa_idx:
                        G[primary_idx][sa_idx]['count'] += 1
                        G[primary_idx][sa_idx]['orientations'].append((primary_region[0], sa_parts[2]))
                        G[primary_idx][sa_idx]['connections'].append(connection_detail)

    # Filter out connections with fewer than 3 supplementary reads
    filtered_graph = defaultdict(lambda: defaultdict(lambda: {'count': 0, 'orientations': [], 'connections': []}))
    for node1 in G:
        for node2 in G[node1]:
            if G[node1][node2]['count'] >= 3:  # Only retain connections with 3 or more supplementary reads
                filtered_graph[node1][node2] = G[node1][node2]

    log_time_end(start_time, "Generating graph from supplementary alignments")
    return filtered_graph

# Function to find cycles in a graph (DFS-based)
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
    start_time = log_time_start("Cycle detection")

    visited = set()
    all_cycles = []
    nodes = list(graph.keys())
    for node in nodes:
        if node not in visited:
            dfs_visit(node, graph, [], visited, all_cycles)

    log_time_end(start_time, "Cycle detection")
    return all_cycles

# Function to write cycles to file, separating longer cycles (>3 nodes) and shorter ones
def write_cycles_to_file(cycles, overamplified_regions, cycle_file, intermediate_file):
    start_time = log_time_start("Writing cycles to file")

    with open(cycle_file, 'w') as cycles_file, open(intermediate_file, 'w') as intermediate_file:
        for cycle in cycles:
            # Get distinct regions (chrom:start-end)
            distinct_nodes = set()
            for node in cycle:
                chrom, start, end = overamplified_regions[node]
                distinct_nodes.add(f"{chrom}:{start}-{end}")

            if len(distinct_nodes) > 3:  # Only write cycles with more than 3 distinct nodes to cycles_file
                cycles_file.write("Cycle:\n")
                for node in cycle:
                    chrom, start, end = overamplified_regions[node]
                    cycles_file.write(f"{chrom}:{start}-{end}\n")
                cycles_file.write("\n")
            else:  # Write shorter cycles to intermediate connections
                intermediate_file.write("Intermediate Connection:\n")
                for node in cycle:
                    chrom, start, end = overamplified_regions[node]
                    intermediate_file.write(f"{chrom}:{start}-{end}\n")
                intermediate_file.write("\n")

    print(f"Cycles with more than 3 nodes written to {cycle_file}")
    print(f"Shorter cycles written to {intermediate_file}")

    log_time_end(start_time, "Writing cycles to file")

# Function to extract reads strictly from a specific region (either start or end)
def extract_reads_for_poa(bam_file, chrom, pos, window_size=500):
    reads = []
    bam = pysam.AlignmentFile(bam_file, "rb")

    # Fetch reads near the position (either start or end of the region)
    for read in bam.fetch(chrom, max(0, pos - window_size), pos + window_size):
        if read.mapping_quality > 20:
            reads.append(read)

    bam.close()
    return reads

# Function to perform POA and map consensus sequence back to the genome
def poa_and_mapping(consensus_prefix, reference_genome, reads):
    start_time = log_time_start("Performing POA and mapping consensus")

    # Write reads to a FASTA file for POA
    fasta_file = f"{consensus_prefix}_reads.fasta"
    with open(fasta_file, 'w') as fasta_out:
        for i, read in enumerate(reads):
            fasta_out.write(f">{read.query_name}_{i}\n{read.query_sequence}\n")

    consensus_fasta = f"{consensus_prefix}_consensus.fa".replace(" ", "_")
    mapped_bams = []

    # POA step (strictly using the reads from the two ends)
    subprocess.call(['wtdbg2', '-i', fasta_file, '-fo', consensus_prefix, '-t', '4'])
    subprocess.call(['wtpoa-cns', '-i', f"{consensus_prefix}.ctg.lay.gz", '-fo', consensus_fasta])

    # Process the consensus and map to the reference genome
    with open(consensus_fasta, 'r') as fasta_file:
        seq_id = None
        seq_data = []
        for line in fasta_file:
            if line.startswith(">"):
                if seq_id:  # Save previous sequence if exists
                    bam_file = map_consensus_sequence(seq_id, seq_data, reference_genome, consensus_prefix)
                    mapped_bams.append(bam_file)
                seq_id = line.strip()[1:]  # Extract the sequence ID
                seq_data = []
            else:
                seq_data.append(line.strip())
        if seq_id:  # Map the last sequence
            bam_file = map_consensus_sequence(seq_id, seq_data, reference_genome, consensus_prefix)
            mapped_bams.append(bam_file)

    log_time_end(start_time, "Performing POA and mapping consensus")
    return mapped_bams

# Helper function to map a single consensus sequence to the genome
def map_consensus_sequence(seq_id, seq_data, reference_genome, consensus_prefix):
    sequence_fasta = f"{consensus_prefix}_{seq_id}.fa".replace(" ", "_")
    with open(sequence_fasta, 'w') as out_fasta:
        out_fasta.write(f">{seq_id}\n")
        out_fasta.write("\n".join(seq_data) + "\n")

    sorted_bam = f"{consensus_prefix}_{seq_id}_sorted.bam".replace(" ", "_")
    minimap2_cmd = ['minimap2', '-ax', 'map-ont', reference_genome, sequence_fasta]
    samtools_sort_cmd = ['samtools', 'sort', '-o', sorted_bam]
    samtools_index_cmd = ['samtools', 'index', sorted_bam]

    print(f"Mapping consensus sequence {seq_id} to reference genome with minimap2")
    
    # Run minimap2 and pipe output to samtools sort
    with subprocess.Popen(minimap2_cmd, stdout=subprocess.PIPE) as minimap2_proc:
        with subprocess.Popen(samtools_sort_cmd, stdin=minimap2_proc.stdout) as samtools_proc:
            minimap2_proc.stdout.close()
            samtools_proc.communicate()

    print(f"Sorting and indexing BAM file {sorted_bam}")
    subprocess.call(samtools_index_cmd)

    return sorted_bam

# Function to process all potential connections strictly between the two ends of regions
def process_all_connections(bam_file, reference_genome, output_dir, regions):
    start_time = log_time_start("Processing all connections")

    poa_workspace = os.path.join(output_dir, "poa_workspace")
    os.makedirs(poa_workspace, exist_ok=True)

    refined_connections_file = os.path.join(output_dir, "refined_connections.txt").replace(" ", "_")
    with open(refined_connections_file, 'w') as refined_file:
        for i in range(len(regions) - 1):
            chrom1, start1, end1 = regions[i]
            chrom2, start2, end2 = regions[i + 1]
            
            # Start of chrom1 to Start of chrom2
            reads_start1_start2 = extract_reads_for_poa(bam_file, chrom1, start1) + extract_reads_for_poa(bam_file, chrom2, start2)
            if reads_start1_start2:
                process_poa_and_log(refined_file, reads_start1_start2, poa_workspace, reference_genome)

    print(f"Refined connections file written to {refined_connections_file}")

    log_time_end(start_time, "Processing all connections")

# Helper function to perform POA, map the reads, and log the result
def process_poa_and_log(refined_file, reads, poa_workspace, reference_genome):
    if not reads:
        print("No reads found for the connection.")
        return

    poa_prefix = os.path.join(poa_workspace, "poa_connection").replace(" ", "_")
    fasta_file = f"{poa_prefix}_reads.fasta"
    
    with open(fasta_file, 'w') as fasta_out:
        for i, read in enumerate(reads):
            fasta_out.write(f">{read.query_name}_{i}\n{read.query_sequence}\n")

    refined_bams = poa_and_mapping(poa_prefix, reference_genome, reads)

    for bam_file in refined_bams:
        refined_coordinates, supporting_reads = extract_refined_coordinates_and_supporting_reads(bam_file, poa_prefix)
        if len(refined_coordinates) >= 2:
            refined_start = refined_coordinates[0]
            refined_end = refined_coordinates[-1]

            start_orientation = refined_start[3]
            end_orientation = refined_end[3]

            refined_file.write(
                f"{refined_start[0]}:{refined_start[1]}-{refined_start[2]} ({start_orientation}) -> "
                f"{refined_end[0]}:{refined_end[1]}-{refined_end[2]} ({end_orientation}), "
                f"Supporting Reads: {supporting_reads}\n"
            )

# Function to extract refined connection points and supporting reads from the consensus
def extract_refined_coordinates_and_supporting_reads(bam_file, poa_prefix):
    refined_coordinates = []
    bam = pysam.AlignmentFile(bam_file, "rb")
    
    for read in bam.fetch():
        if not read.is_unmapped:
            chrom = read.reference_name
            start = read.reference_start
            end = read.reference_end
            strand = '-' if read.is_reverse else '+'
            refined_coordinates.append((chrom, start, end, strand))
    
    bam.close()
    reads_fasta_file = f"{poa_prefix}_reads.fasta".replace(" ", "_")
    supporting_reads = count_reads_in_fasta(reads_fasta_file)

    return refined_coordinates, supporting_reads

# Function to count the number of reads in the FASTA file
def count_reads_in_fasta(fasta_file):
    count = 0
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith(">"):
                count += 1
    return count

# Main execution function
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Identify overamplified regions, detect cycles, and perform POA.")
    parser.add_argument("--bam", required=True, help="Path to the BAM file")
    parser.add_argument("--output_dir", required=True, help="Output directory for logs and files")
    parser.add_argument("--reference_genome", required=True, help="Path to the reference genome for mapping consensus")
    parser.add_argument("--min_support", type=int, help="Minimum support value")
    parser.add_argument("--depth", type=float, help="Sequencing depth")

    args = parser.parse_args()

    bam_file = args.bam
    output_dir = args.output_dir
    reference_genome = args.reference_genome
    min_support = args.min_support
    depth = args.depth

    os.makedirs(output_dir, exist_ok=True)

    overall_start_time = log_time_start(f"Processing {bam_file}")

    # Estimate average depth and calculate min_supp
    avg_coverage, minsupp = estimate_depth_and_minsupp(bam_file, min_support, depth, output_dir)
    
    # Identify overamplified regions
    overamplified_regions = identify_overamplified_regions(bam_file, avg_coverage)
    bam_basename = os.path.splitext(os.path.basename(bam_file))[0]

    # Write overamplified regions to a file
    output_regions_file = os.path.join(output_dir, f"{bam_basename}_overamplified_regions.txt")
    with open(output_regions_file, 'w') as out_file:
        for region in overamplified_regions:
            out_file.write(f"{region[0]}\t{region[1]}\t{region[2]}\n")
    print(f"Overamplified regions written to {output_regions_file}")

    # Extract supplementary reads
    supplementary_reads, bam_header = get_supplementary_reads(bam_file)

    # Create index mapping for overamplified regions
    region_to_index = {region: idx for idx, region in enumerate(overamplified_regions)}

    # Generate graph from supplementary alignments
    graph = generate_graph_from_supplementary_alignments(supplementary_reads, overamplified_regions, region_to_index, bam_header)

    # Write graph to file
    output_graph_file = os.path.join(output_dir, f"{bam_basename}_graph.txt")
    with open(output_graph_file, 'w') as graph_file:
        for node1 in graph:
            for node2 in graph[node1]:
                count = graph[node1][node2]['count']
                orientations = graph[node1][node2]['orientations']
                connections = graph[node1][node2]['connections']
                graph_file.write(f"Region {node1} -> Region {node2}: Count: {count}, Orientations: {orientations}, Connections: {connections}\n")
    print(f"Graph written to {output_graph_file}")

    # Detect cycles in the graph
    cycles = find_cycles(graph)

    # Write cycles to file and handle intermediate connections
    output_cycles_file = os.path.join(output_dir, f"{bam_basename}_cycles.txt")
    output_intermediate_file = os.path.join(output_dir, f"{bam_basename}_intermediate_connections.txt")
    write_cycles_to_file(cycles, overamplified_regions, output_cycles_file, output_intermediate_file)

    # Extract reads around breakpoints in cycles for POA (only for cycles in cycles.txt)
    cycle_regions = [overamplified_regions[node] for cycle in cycles if len(set(cycle)) > 3 for node in cycle]

    # Perform POA and refined connections
    process_all_connections(bam_file, reference_genome, output_dir, cycle_regions)

    log_time_end(overall_start_time, f"Processing {bam_file}")

