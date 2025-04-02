"""
Core functionality for metapresence with memory optimization.
"""

import os
import sys
import multiprocessing
import logging
import gc
from .metrics import parse_metrics_criteria
from .utils.parsing import parse_single_fasta, parse_multi_fasta, chunks
from .utils.plot import create_metrics_plot
from .bam_processing import process_genome, get_average_read_length

# Get logger
logger = logging.getLogger('metapresence')

try:
    import numpy as np
except ImportError:
    logger.error('Error: numpy is not installed.')
    sys.exit(1)

try:
    import pysam
except ImportError:
    logger.error('Error: pysam is not installed.')
    sys.exit(1)


# Helper function for multiprocessing - must be at module level
def _process_genome_wrapper(args_tuple):
    """Wrapper function for process_genome to make it compatible with multiprocessing."""
    genome_list, sorted_bam, contigs_lengths_d, genome_contigs, av_read_len, params = args_tuple

    results = {}
    for genome in genome_list:
        genome_result = process_genome(genome, genome_contigs[genome], sorted_bam,
                                      contigs_lengths_d, genome_contigs, av_read_len, params)
        if genome_result:
            results.update(genome_result)

    # Force garbage collection to free memory
    gc.collect()

    return results


def process_alignment(args):
    """
    Process the alignment file and generate metrics and abundance outputs.

    Args:
        args: Command line arguments.
    """
    # Extract args for easier reference
    scaffold_input = args.input_fasta
    sorted_bam = args.indexed_sorted_bam
    processes = args.p

    # Parse fasta files
    contigs_lengths_d, genome_contigs = load_reference_data(scaffold_input, args, processes)

    # Read BAM file and calculate average read length
    logger.info('Reading BAM file')
    try:
        av_read_len = get_average_read_length(sorted_bam)
    except ValueError:
        logger.error(f'Cannot find index for {sorted_bam}. Index and BAM file should be placed in the same directory.')
        sys.exit(1)

    # Process genomes and calculate metrics with memory efficiency
    metrics_dict = process_genomes_with_multiprocessing(genome_contigs, sorted_bam, contigs_lengths_d,
                                                      av_read_len, args, processes)

    # Generate output files
    generate_output_files(metrics_dict, args)

    # Final cleanup
    gc.collect()

    logger.info('All processing complete')


def load_reference_data(scaffold_input, args, processes):
    """
    Load reference genome data from FASTA files with memory optimization.

    Args:
        scaffold_input: Path to input FASTA file or directory.
        args: Command line arguments.
        processes: Number of processes to use for parsing.

    Returns:
        Tuple of (contigs_lengths_d, genome_contigs).
    """
    if os.path.isdir(scaffold_input):
        logger.info('Reading fasta files')
        fastas = [os.path.join(scaffold_input, x) for x in os.listdir(scaffold_input)]
        logger.debug(f'Found {len(fastas)} fasta files')

        # Make sure we create enough chunks to utilize available processes
        actual_processes = min(processes, len(fastas))
        chunked_fastas = chunks(fastas, max(1, len(fastas) // actual_processes))
        parseFasta_args = [(x, args.all_contigs) for x in chunked_fastas]

        contigs_lengths_d = {}
        genome_contigs = {}

        # Process files using available processes
        logger.debug(f'Processing fasta files using {len(chunked_fastas)} chunks with {actual_processes} processes')

        with multiprocessing.Pool(processes=actual_processes) as pool:
            for result_chunk in pool.map(parse_multi_fasta, parseFasta_args):
                # Update dictionaries with chunk results
                for contig, length in result_chunk[0].items():
                    contigs_lengths_d[contig] = length

                for genome, contigs_list in result_chunk[1].items():
                    genome_contigs[genome] = contigs_list

                # Force garbage collection
                gc.collect()

        if args.all_contigs:
            genome_contigs = {y: [y] for x in genome_contigs for y in genome_contigs[x]}

        logger.debug(f'Found {len(contigs_lengths_d)} contigs in {len(genome_contigs)} genomes')

    elif os.path.isfile(scaffold_input):
        logger.info('Reading fasta file')
        contigs_lengths_d, genome_contigs = parse_single_fasta(scaffold_input, args.input_bins, args.all_contigs)
        logger.debug(f'Found {len(contigs_lengths_d)} contigs in {len(genome_contigs)} genomes')

    else:
        logger.error(f'Cannot open {scaffold_input}')
        sys.exit(1)

    logger.info('Fasta processing complete')
    return contigs_lengths_d, genome_contigs


def process_genomes_with_multiprocessing(genome_contigs, sorted_bam, contigs_lengths_d, av_read_len, args, processes):
    """
    Process genomes using multiprocessing with memory optimization.

    Args:
        genome_contigs: Dictionary mapping genomes to their contigs.
        sorted_bam: Path to the sorted BAM file.
        contigs_lengths_d: Dictionary with contig lengths.
        av_read_len: Average read length.
        args: Command line arguments.
        processes: Number of processes to use.

    Returns:
        Dictionary with metrics for each genome.
    """
    # Get list of genomes
    genomes = list(genome_contigs.keys())
    logger.debug(f'Processing {len(genomes)} genomes')

    # Calculate optimal number of chunks based on available processes
    # Make sure we create at least as many chunks as processes
    num_chunks = processes
    if args.low_memory and len(genomes) > processes * 10:
        # For low memory mode with many genomes, create more chunks
        num_chunks = min(len(genomes), processes * 2)

    # Create chunks with approximately equal size
    chunk_size = max(1, len(genomes) // num_chunks)
    chunked_genomes = [genomes[i:i+chunk_size] for i in range(0, len(genomes), chunk_size)]

    # Make sure we're not creating more chunks than processes
    while len(chunked_genomes) > processes * 2:
        # Merge the smallest chunks
        chunked_genomes.sort(key=len)
        chunked_genomes[1].extend(chunked_genomes[0])
        chunked_genomes.pop(0)

    # Ensure we're not creating more chunks than genomes
    actual_chunks = min(len(chunked_genomes), len(genomes))
    actual_processes = min(processes, actual_chunks)

    logger.debug(f'Split {len(genomes)} genomes into {actual_chunks} chunks for {actual_processes} processes')

    # Prepare arguments for the multiprocessing pool
    pool_args = [(genome_chunk, sorted_bam, contigs_lengths_d, genome_contigs, av_read_len, args)
                for genome_chunk in chunked_genomes]

    # Process using available processes
    metrics_dict = {}
    logger.info(f'Calculating metrics using {actual_processes} processes')

    # Use a pool with the calculated number of processes
    with multiprocessing.Pool(processes=actual_processes) as pool:
        for batch_idx, result_batch in enumerate(pool.imap(_process_genome_wrapper, pool_args)):
            if result_batch:
                metrics_dict.update(result_batch)

            # Log progress
            logger.debug(f'Processed batch {batch_idx+1}/{len(pool_args)} containing {len(result_batch) if result_batch else 0} genomes')

            # Force garbage collection after each batch
            gc.collect()

    logger.info('BAM processing complete')

    # Sort by genome name
    metrics_dict = {x: metrics_dict[x] for x in sorted(metrics_dict.keys())}

    return metrics_dict


def generate_output_files(metrics_dict, args):
    """
    Generate output files based on calculated metrics.

    Args:
        metrics_dict: Dictionary with metrics for each genome.
        args: Command line arguments.
    """
    logger.info('Preparing output files')

    # Write metrics file
    metrics_file_path = args.o + "_metrics.tsv"
    logger.info(f'Writing metrics to {metrics_file_path}')
    present_coverage = write_metrics_file(metrics_file_path, metrics_dict, args)

    # Write abundances file
    abundances_file_path = args.o + "_abundances.tsv"
    logger.info(f'Writing abundances to {abundances_file_path}')
    write_abundances_file(abundances_file_path, present_coverage)

    # Generate plot if requested
    if args.plot_metrics:
        # Determine output format for static plots
        plot_format = getattr(args, 'plot_format', 'png')
        plot_path = f"{args.o}_scatterplot.{plot_format}"

        logger.info(f'Generating metrics plot at {plot_path}')

        # Create the appropriate type of plot
        create_metrics_plot(metrics_file_path, plot_path, args.min_reads, args.unpaired,
                           interactive=getattr(args, 'interactive', False))


def write_metrics_file(metrics_file_path, metrics_dict, args):
    """
    Write metrics to a TSV file.

    Args:
        metrics_file_path: Path to the output metrics file.
        metrics_dict: Dictionary with metrics for each genome.
        args: Command line arguments.

    Returns:
        Dictionary mapping genome names to coverage values for present genomes.
    """
    present_coverage = {}

    with open(metrics_file_path, "w") as metrics_file:
        metrics_file.write('genome\tlength\tcoverage\tbreadth\tBER\tFUG1\tFUG2\tread_count\n')

        # Process genomes in batches to limit memory usage for large datasets
        batch_size = 1000
        genomes = list(metrics_dict.keys())

        for i in range(0, len(genomes), batch_size):
            batch_genomes = genomes[i:i+batch_size]

            for genome in batch_genomes:
                values = metrics_dict[genome].strip().split('\t')
                cov, ber, fug1 = float(values[1]), float(values[3]), float(values[4])

                if values[5] == 'unpaired' or args.unpaired:
                    fug2 = None
                else:
                    try:
                        fug2 = float(values[5])
                    except ValueError:
                        fug2 = None

                nreads = int(values[6])

                if parse_metrics_criteria(cov, ber, fug1, fug2, nreads, args.fug_criterion,
                                         args.min_reads, args.max_for_fug, args.ber_threshold,
                                         args.fug_threshold, args.unpaired):
                    present_coverage[genome] = cov
                    logger.debug(f'Genome {genome} is present with coverage {cov}')
                else:
                    logger.debug(f'Genome {genome} is not considered present')

                metrics_file.write(f'{genome}\t{metrics_dict[genome]}\n')

            # Force garbage collection after each batch
            gc.collect()

    return present_coverage


def write_abundances_file(abundances_file_path, present_coverage):
    """
    Write abundances to a TSV file.

    Args:
        abundances_file_path: Path to the output abundances file.
        present_coverage: Dictionary mapping genome names to coverage values for present genomes.
    """
    with open(abundances_file_path, "w") as abundances_file:
        abundances_file.write('genome\trelative_abundance_%\n')

        if present_coverage:
            totcov = sum(present_coverage.values())

            # Sort by abundance (highest first) for better readability
            sorted_genomes = sorted(present_coverage.keys(),
                                   key=lambda g: present_coverage[g],
                                   reverse=True)

            for genome in sorted_genomes:
                abundance = present_coverage[genome] / totcov * 100
                abundances_file.write(f'{genome}\t{abundance}\n')
                logger.debug(f'Genome {genome} abundance: {abundance:.2f}%')
        else:
            logger.warning('No genomes were determined to be present in the sample')