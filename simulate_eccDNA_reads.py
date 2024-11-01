import os
import subprocess
import multiprocessing
import argparse

def ensure_directory(directory):
    """Ensure that a directory exists; if it doesn't, create it."""
    if not os.path.exists(directory):
        os.makedirs(directory)
        print(f"Directory created: {directory}")

def calculate_num_reads(eccDNA_length, coverage, read_length):
    """Calculate the number of reads required for the specified coverage."""
    return (coverage * eccDNA_length) // read_length

def simulate_reads(args):
    """Wrapper function for simulating reads that handles exceptions and logs progress."""
    eccDNA_file, eccDNA_name, eccDNA_length, model_prefix, output_base_dir, coverage, read_length, idx = args
    try:
        print(f"Starting simulation for {eccDNA_name} at coverage {coverage}X, process index {idx}")
        num_reads = calculate_num_reads(eccDNA_length, coverage, read_length)
        output_dir = os.path.join(output_base_dir, f"{eccDNA_name}_coverage_{coverage}X")
        ensure_directory(output_dir)
        output_prefix = os.path.join(output_dir, "simulated_reads")
        command = [
            "python", "./NanoSim/src/simulator.py", "genome",
            "-dna_type", "circular",
            "-rg", eccDNA_file,
            "-c", model_prefix,
            "-b", "guppy",
            "-o", output_prefix,
            "-n", str(num_reads),
            "-t", "16"  # Number of threads per subprocess
        ]
        subprocess.run(command, check=True)
        print(f"Simulation completed for {eccDNA_name} at coverage {coverage}X")
    except Exception as e:
        print(f"Error during simulation for {eccDNA_name}: {e}")

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Simulate sequencing reads from eccDNA files using NanoSim.")
    parser.add_argument("--input_dir", required=True, help="Directory containing eccDNA FASTA files.")
    parser.add_argument("--output_dir", required=True, help="Output directory for simulated reads.")
    parser.add_argument("--model_prefix", default="/path/to/default/model", help="Prefix for the NanoSim model directory.")
    parser.add_argument("--read_length", type=int, default=20000, help="Length of reads to simulate.")
    return parser.parse_args()

def main():
    args = parse_arguments()

    ensure_directory(args.output_dir)
    eccDNA_files = [f for f in os.listdir(args.input_dir) if f.endswith('.fasta')]
    eccDNA_files.sort()  # Sort to ensure consistent processing order

    # Coverage levels can be adjusted or passed as an argument if variability is needed
    coverage_levels = [100, 200, 500]

    # Prepare arguments for multiprocessing
    tasks = []
    for idx, eccDNA_file in enumerate(eccDNA_files):
        eccDNA_path = os.path.join(args.input_dir, eccDNA_file)
        eccDNA_name = eccDNA_file.replace('.fasta', '')
        with open(eccDNA_path, 'r') as f:
            f.readline()  # Skip the header
            eccDNA_sequence = f.read().replace('\n', '')
            eccDNA_length = len(eccDNA_sequence)

        coverage = coverage_levels[idx % len(coverage_levels)]
        tasks.append((eccDNA_path, eccDNA_name, eccDNA_length, args.model_prefix, args.output_dir, coverage, args.read_length, idx))

    # Create a process pool and execute simulations
    with multiprocessing.Pool(processes=min(len(tasks), multiprocessing.cpu_count())) as pool:
        pool.map(simulate_reads, tasks)

    print("Read simulation complete. Check the output directories for details.")

if __name__ == "__main__":
    main()

