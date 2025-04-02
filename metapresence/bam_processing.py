"""
Functions for processing BAM files and calculating metrics with memory efficiency optimizations.
"""

import logging
import pysam
import numpy as np
from itertools import islice
import gc

# Get logger
logger = logging.getLogger('metapresence')


def process_genome(genome_name, contigs, sorted_bam, contigs_lengths_d, genome_contigs, av_read_len, args):
    """
    Process a single genome and calculate its metrics.

    Args:
        genome_name: Name of the genome to process.
        contigs: List of contigs for this genome.
        sorted_bam: Path to the sorted BAM file.
        contigs_lengths_d: Dictionary with contig lengths.
        genome_contigs: Dictionary mapping genomes to their contigs.
        av_read_len: Average read length.
        args: Command line arguments.

    Returns:
        Dictionary with metrics for the genome, or None if no reads mapped.
    """
    logger.debug(f'Processing genome: {genome_name}')

    # Initialize data structures
    genome_data = collect_genome_data(genome_name, contigs, sorted_bam, contigs_lengths_d, args)

    # If no reads mapped, skip this genome
    if genome_data["read_count"] == 0:
        logger.debug(f'Genome {genome_name} has zero mapped reads')
        return None

    # Calculate metrics for the current genome
    metrics = calculate_metrics(genome_name, genome_data, av_read_len, args.unpaired)

    # Clean up large data structures to free memory
    del genome_data
    gc.collect()

    return {genome_name: metrics}


def collect_genome_data(genome_name, contigs, sorted_bam, contigs_lengths_d, args):
    """
    Collect data from BAM file for a single genome with streaming approach.

    Args:
        genome_name: Name of the genome.
        contigs: List of contigs belonging to this genome.
        sorted_bam: Path to the sorted BAM file.
        contigs_lengths_d: Dictionary with contig lengths.
        args: Command line arguments.

    Returns:
        Dictionary containing genome data including read positions, coverage, etc.
    """
    # Open BAM file for streaming
    sam = pysam.AlignmentFile(sorted_bam.split()[0], "rb")

    genome_covbases_lengths = {}
    noread_contig_length = 0
    genome_readcount = 0
    genome_pointer = 0

    # Use lists for positions (will be converted to arrays only when needed)
    allpositions = [0]  # Start with 0 position
    allpositions2 = [0]  # For second mate in paired reads

    # Process each contig for this genome
    for contig in contigs:
        # Check if contig exists in the length dictionary
        if contig not in contigs_lengths_d:
            if not args.quiet:
                logger.warning(f'Contig {contig} of sequence {genome_name} was not present in the fasta file(s). Skipping sequence...')
            genome_readcount = 0
            break

        # Process this contig with streaming
        contig_data = process_contig_streaming(contig, genome_pointer, sam, contigs_lengths_d,
                                              genome_covbases_lengths, genome_name, noread_contig_length)

        # Update genome data with contig results
        genome_pointer = contig_data["genome_pointer"]
        noread_contig_length = contig_data["noread_contig_length"]
        genome_readcount += contig_data["read_count"]

        # Add positions to master lists (only store if needed for metrics calculation)
        if "positions1" in contig_data and contig_data["positions1"]:
            allpositions.extend(contig_data["positions1"])
        if "positions2" in contig_data and contig_data["positions2"]:
            allpositions2.extend(contig_data["positions2"])

    sam.close()

    # Return collected genome data
    return {
        "covbases_lengths": genome_covbases_lengths,
        "read_count": genome_readcount,
        "positions1": allpositions,
        "positions2": allpositions2
    }


def process_contig_streaming(contig, genome_pointer, sam, contigs_lengths_d, genome_covbases_lengths,
                            genome_name, noread_contig_length, batch_size=10000):
    """
    Process a single contig using a streaming approach to reduce memory usage.

    Args:
        contig: Name of the contig to process.
        genome_pointer: Current pointer position in the genome.
        sam: Open BAM file handle.
        contigs_lengths_d: Dictionary with contig lengths.
        genome_covbases_lengths: Dictionary to update with coverage data.
        genome_name: Name of the genome.
        noread_contig_length: Accumulated length of contigs with no reads.
        batch_size: Size of read batches to process at once.

    Returns:
        Dictionary with contig processing results.
    """
    # Initialize data structures
    contig_length = contigs_lengths_d[contig]

    # Use a sparse representation for coverage changes
    # Instead of an array of zeros, use a dictionary to store only non-zero changes
    coverage_changes = {}
    read_count = 0
    pair_taken = set()
    positions1 = []
    positions2 = []

    try:
        # Stream reads in batches to reduce memory usage
        read_batch = []
        fetch_iterator = sam.fetch(contig=contig)

        while True:
            # Get the next batch of reads
            read_batch = list(islice(fetch_iterator, batch_size))
            if not read_batch:
                break

            # Process this batch
            for read in read_batch:
                read_count += 1
                start_pos = read.reference_start
                end_pos = start_pos + read.query_length

                # Track positions for calculating read distribution
                if read.query_name not in pair_taken:
                    positions1.append(start_pos + genome_pointer)
                    pair_taken.add(read.query_name)
                else:
                    positions2.append(start_pos + genome_pointer)

                # Update coverage change points using sparse representation
                if start_pos < contig_length:
                    coverage_changes[start_pos] = coverage_changes.get(start_pos, 0) + 1

                if end_pos < contig_length:
                    coverage_changes[end_pos] = coverage_changes.get(end_pos, 0) - 1

            # Clear the batch to free memory
            read_batch = []

    except ValueError:
        # This can happen if the contig is not in the BAM file
        pass

    # Update genome pointer
    new_genome_pointer = genome_pointer + contig_length

    # Handle case where no reads mapped to this contig
    new_noread_contig_length = noread_contig_length

    if read_count == 0:
        new_noread_contig_length += contig_length
        if genome_name in genome_covbases_lengths:
            genome_covbases_lengths[genome_name][0] += contig_length
            genome_covbases_lengths[genome_name][1] += contig_length
            genome_covbases_lengths[genome_name][2] += 0
        else:
            genome_covbases_lengths[genome_name] = [contig_length, contig_length, 0]
    else:
        # Calculate coverage statistics for this contig using sparse representation
        calculate_contig_coverage_sparse(coverage_changes, contig_length, genome_covbases_lengths, genome_name)

    # Clear large data structures
    pair_taken.clear()

    return {
        "read_count": read_count,
        "genome_pointer": new_genome_pointer,
        "noread_contig_length": new_noread_contig_length,
        "positions1": positions1,
        "positions2": positions2
    }


def calculate_contig_coverage_sparse(coverage_changes, contig_length, genome_covbases_lengths, genome_name):
    """
    Calculate coverage statistics for a contig using a sparse representation of coverage changes.

    Args:
        coverage_changes: Dictionary of position -> change value for coverage changes.
        contig_length: Length of the contig.
        genome_covbases_lengths: Dictionary to update with coverage data.
        genome_name: Name of the genome.
    """
    # Convert sparse representation to a sorted list of positions and changes
    change_positions = sorted(coverage_changes.keys())

    # Initialize counters
    current_coverage = 0
    total_depth = 0
    zerodepth = 0
    last_pos = 0

    # Process each change point
    for pos in change_positions:
        # Calculate coverage for the region from last_pos to pos
        if current_coverage == 0:
            zerodepth += (pos - last_pos)
        else:
            total_depth += current_coverage * (pos - last_pos)

        # Update coverage at this position
        current_coverage += coverage_changes[pos]
        last_pos = pos

    # Process the final region from last_pos to contig_length
    if current_coverage == 0:
        zerodepth += (contig_length - last_pos)
    else:
        total_depth += current_coverage * (contig_length - last_pos)

    # Update the coverage data dictionary
    if genome_name in genome_covbases_lengths:
        genome_covbases_lengths[genome_name][0] += zerodepth
        genome_covbases_lengths[genome_name][1] += contig_length
        genome_covbases_lengths[genome_name][2] += total_depth
    else:
        genome_covbases_lengths[genome_name] = [zerodepth, contig_length, total_depth]


def calculate_metrics(genome_name, genome_data, av_read_len, unpaired):
    """
    Calculate metrics for a genome based on collected data.

    Args:
        genome_name: Name of the genome.
        genome_data: Dictionary containing collected genome data.
        av_read_len: Average read length.
        unpaired: Whether reads are unpaired.

    Returns:
        String containing tab-separated metrics.
    """
    covbases_lengths = genome_data["covbases_lengths"][genome_name]
    allpositions = genome_data["positions1"]
    allpositions2 = genome_data["positions2"]
    read_count = genome_data["read_count"]

    # Add fake read at the end of the genome
    allpositions.append(covbases_lengths[1] - av_read_len)
    allpositions2.append(covbases_lengths[1] - av_read_len)

    # Calculate FUG1 (read distribution metric for first mate)
    fug1 = calculate_fug_memory_efficient(allpositions, covbases_lengths[1])

    # Calculate FUG2 (read distribution metric for second mate)
    if unpaired:
        fug2 = 'unpaired'
    else:
        fug2 = calculate_fug_memory_efficient(allpositions2, covbases_lengths[1])

    # Calculate coverage and breadth
    coverage = covbases_lengths[2] / covbases_lengths[1]
    breadth = (covbases_lengths[1] - covbases_lengths[0]) / covbases_lengths[1]

    # Calculate BER (Breadth-Expected breadth Ratio)
    expected_breadth = 1 - np.exp(-0.883 * coverage)
    ber = breadth / expected_breadth if expected_breadth > 0 else 0

    # Format the metrics as a tab-separated string
    metrics = f"{covbases_lengths[1]}\t{coverage}\t{breadth}\t{ber}\t{fug1}\t{fug2}\t{read_count}"
    logger.debug(f'Calculated metrics for {genome_name}: coverage={coverage:.4f}, breadth={breadth:.4f}, BER={ber:.4f}')

    return metrics


def calculate_fug_memory_efficient(positions, genome_length, chunk_size=100000):
    """
    Calculate Fraction of Unexpected Gaps (FUG) metric with memory efficiency.

    Args:
        positions: List of read positions.
        genome_length: Length of the genome.
        chunk_size: Size of chunks for processing large datasets.

    Returns:
        FUG value as float, or 'nan' if calculation is not possible.
    """
    # Check if we have enough positions to calculate
    if len(positions) <= 1:
        return 'nan'

    window_size = round(genome_length / len(positions))

    if window_size <= 0:
        return 'nan'

    # For large position lists, process in chunks to save memory
    if len(positions) > chunk_size:
        logger.debug(f"Processing {len(positions)} positions in chunks of {chunk_size}")

        # Sort positions (needed for calculating distances)
        positions.sort()

        # Calculate distances in chunks
        distances = []
        for i in range(0, len(positions) - 1, chunk_size):
            end_idx = min(i + chunk_size, len(positions) - 1)
            # Process current chunk
            if i == 0:
                # First chunk - use adjacent positions
                chunk_positions = positions[i:end_idx + 1]
                chunk_distances = np.diff(chunk_positions)
            else:
                # Include the last position from previous chunk
                chunk_positions = [positions[i-1]] + positions[i:end_idx + 1]
                chunk_distances = np.diff(chunk_positions)

            distances.extend(chunk_distances)

        # Convert to numpy array for histogram
        distances = np.array(distances)
    else:
        # For smaller datasets, use vectorized operations
        positions = np.array(positions)
        positions.sort()  # Ensure positions are sorted
        distances = np.diff(positions)

    # Calculate histogram with memory efficiency
    max_distance = np.max(distances)

    # Create histogram bins starting from window_size
    if max_distance < window_size:
        # No distances above window_size
        return 1.0  # Perfect distribution

    hist_range = range(window_size, int(max_distance) + 1)

    # Calculate histogram
    hist, _ = np.histogram(distances, bins=hist_range)
    frequency = hist / (len(positions) - 1)

    # Calculate FUG
    indices = np.arange(len(hist))
    expected_distance = np.sum(frequency * indices)
    fug = (window_size - expected_distance) / window_size

    # Clean up large arrays
    del distances
    gc.collect()

    return fug


def get_average_read_length(sorted_bam, sample_size=1000):
    """
    Calculate the average read length from a BAM file.

    Args:
        sorted_bam: Path to the sorted BAM file.
        sample_size: Number of reads to sample.

    Returns:
        Average read length as an integer.
    """
    logger.debug(f'Calculating average read length from {sorted_bam}')

    try:
        save = pysam.set_verbosity(0)  # to avoid missing index error coming up
        sam = pysam.AlignmentFile(sorted_bam.split()[0], "rb")
        pysam.set_verbosity(save)

        totlen = 0
        c = 0

        # Stream reads for calculating average length
        for read in sam.fetch():
            totlen += read.query_length
            c += 1
            if c == sample_size:
                break

        sam.close()

        if c == 0:
            logger.warning('No reads found in BAM file')
            return 100  # Default value

        av_read_len = round(totlen / c)
        logger.debug(f'Average read length: {av_read_len}')
        return av_read_len

    except Exception as e:
        logger.error(f'Error calculating average read length: {str(e)}')
        return 100  # Default value