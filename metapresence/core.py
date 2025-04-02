"""
Core functionality for metapresence.
"""

import os
import sys
import multiprocessing
import logging
from .metrics import parse_metrics_criteria
from .utils.parsing import parse_single_fasta, parse_multi_fasta, chunks
from .utils.plot import create_metrics_plot

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
def _parse_bam_wrapper(args_tuple):
    """Wrapper function for parse_bam to make it compatible with multiprocessing."""
    genome_list, sorted_bam, contigs_lengths_d, genome_contigs, av_read_len, params = args_tuple
    return parse_bam(genome_list, sorted_bam, contigs_lengths_d, genome_contigs, av_read_len, params)


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
    if os.path.isdir(scaffold_input):
        logger.info('Reading fasta files')
        fastas = [os.path.join(scaffold_input, x) for x in os.listdir(scaffold_input)]
        logger.debug(f'Found {len(fastas)} fasta files')

        chunked_fastas = chunks(fastas, processes)
        parseFasta_args = [(x, args.all_contigs) for x in chunked_fastas]

        with multiprocessing.Pool(processes=processes) as pool:
            logger.debug(f'Processing fasta files with {processes} processes')
            all_parsing = pool.map(parse_multi_fasta, parseFasta_args)

        contigs_lengths_d = {x: result[0][x] for result in all_parsing for x in result[0]}
        genome_contigs = {x: result[1][x] for result in all_parsing for x in result[1]}

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

    # Read BAM file and calculate average read length
    logger.info('Reading BAM file')
    try:
        save = pysam.set_verbosity(0)  # to avoid missing index error coming up
        sam = pysam.AlignmentFile(sorted_bam.split()[0], "rb")
        pysam.set_verbosity(save)

        totlen = 0
        c = 0
        for read in sam.fetch():
            totlen += read.query_length
            c += 1
            if c == 1000:
                break

        av_read_len = round(totlen / c)
        logger.debug(f'Average read length: {av_read_len}')
        sam.close()

    except ValueError:
        logger.error(f'Cannot find index for {sorted_bam}. Index and BAM file should be placed in the same directory.')
        sys.exit(1)

    # Process BAM file and calculate metrics
    genomes = [x for x in genome_contigs]
    logger.debug(f'Processing {len(genomes)} genomes')

    chunked_genomes = chunks(genomes, processes)

    # Prepare arguments for the multiprocessing pool
    pool_args = [(genome_chunk, sorted_bam, contigs_lengths_d, genome_contigs, av_read_len, args)
                 for genome_chunk in chunked_genomes]

    logger.info(f'Calculating metrics using {processes} processes')
    with multiprocessing.Pool(processes=processes) as pool:
        all_metrics = pool.map(_parse_bam_wrapper, pool_args)

    logger.info('BAM processing complete')
    logger.info('Preparing output files')

    # Generate output files
    metrics_dict = {}
    for result in all_metrics:
        if result:  # Check if result is not None
            for x in result:
                metrics_dict[x] = result[x]

    metrics_dict = {x: metrics_dict[x] for x in sorted(metrics_dict.keys())}

    # Write metrics file
    metrics_file_path = args.o + "_metrics.tsv"
    logger.info(f'Writing metrics to {metrics_file_path}')

    with open(metrics_file_path, "w") as metrics_file:
        metrics_file.write('genome\tlength\tcoverage\tbreadth\tBER\tFUG1\tFUG2\tread_count\n')

        present_coverage = {}
        for genome in metrics_dict:
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

    # Write abundances file
    abundances_file_path = args.o + "_abundances.tsv"
    logger.info(f'Writing abundances to {abundances_file_path}')

    with open(abundances_file_path, "w") as abundances_file:
        abundances_file.write('genome\trelative_abundance_%\n')

        if present_coverage:
            totcov = sum(present_coverage.values())
            for genome in present_coverage:
                abundance = present_coverage[genome] / totcov * 100
                abundances_file.write(f'{genome}\t{abundance}\n')
                logger.debug(f'Genome {genome} abundance: {abundance:.2f}%')
        else:
            logger.warning('No genomes were determined to be present in the sample')

    # Generate plot if requested
    if args.plot_metrics:
        plot_path = args.o + "_scatterplot.png"
        logger.info(f'Generating metrics plot at {plot_path}')
        create_metrics_plot(metrics_file_path, plot_path, args.min_reads, args.unpaired)

    logger.info('All processing complete')


def parse_bam(list_of_genomes, sorted_bam, contigs_lengths_d, genome_contigs, av_read_len, args):
    """
    Parse BAM file and calculate metrics for each genome.

    Args:
        list_of_genomes: List of genome names to process.
        sorted_bam: Path to the sorted BAM file.
        contigs_lengths_d: Dictionary with contig lengths.
        genome_contigs: Dictionary mapping genomes to their contigs.
        av_read_len: Average read length.
        args: Command line arguments.

    Returns:
        Dictionary with metrics for each genome.
    """
    local_logger = logging.getLogger('metapresence')
    printing_dict = {}
    sam = pysam.AlignmentFile(sorted_bam.split()[0], "rb")

    for seq in list_of_genomes:
        local_logger.debug(f'Processing genome: {seq}')
        genome_covbases_lengths = {}
        noread_contig_length = 0
        genome_readcount = 0
        genome_pointer = 0
        allpositions = [0]
        allpositions2 = [0]

        # Retrieve mapping positions from BAM file
        for contig in genome_contigs[seq]:
            if contig not in contigs_lengths_d:
                if not args.quiet:
                    local_logger.warning(f'Contig {contig} of sequence {seq} was not present in the fasta file(s). Skipping sequence...')
                genome_readcount = 0
                break

            genome_change = [0] * contigs_lengths_d[contig]
            rn = 0
            pair_taken = set()

            try:
                for read in sam.fetch(contig=contig):
                    rn += 1
                    start_pos = read.reference_start
                    end_pos = start_pos + read.query_length

                    if read.query_name not in pair_taken:
                        allpositions.append(start_pos + genome_pointer)
                        pair_taken.add(read.query_name)
                    else:
                        allpositions2.append(start_pos + genome_pointer)

                    if start_pos < len(genome_change):
                        genome_change[start_pos] += 1
                    else:
                        continue

                    if end_pos < len(genome_change):
                        genome_change[end_pos] -= 1
                    else:
                        continue

            except ValueError:
                pass

            genome_pointer += contigs_lengths_d[contig]

            if rn == 0:
                noread_contig_length += contigs_lengths_d[contig]
                if seq in genome_covbases_lengths:
                    genome_covbases_lengths[seq][0] += contigs_lengths_d[contig]
                    genome_covbases_lengths[seq][1] += contigs_lengths_d[contig]
                    genome_covbases_lengths[seq][2] += 0
                    continue
                else:
                    genome_covbases_lengths[seq] = [contigs_lengths_d[contig], contigs_lengths_d[contig], 0]
                    continue

            # Calculate contig depth
            current_coverage = 0
            total_depth = 0
            zerodepth = 0

            for position in range(contigs_lengths_d[contig]):
                current_coverage = current_coverage + genome_change[position]
                if current_coverage == 0:
                    zerodepth += 1
                    continue
                total_depth += current_coverage

            if seq in genome_covbases_lengths:
                genome_covbases_lengths[seq][0] += zerodepth
                genome_covbases_lengths[seq][1] += contigs_lengths_d[contig]
                genome_covbases_lengths[seq][2] += total_depth
            else:
                genome_covbases_lengths[seq] = [zerodepth, contigs_lengths_d[contig], total_depth]

            genome_readcount += rn

        # Skip genomes with zero mapped reads
        if genome_readcount == 0:
            local_logger.debug(f'Genome {seq} has zero mapped reads')
            continue

        # Add fake read at the end of the genome
        allpositions.append(genome_covbases_lengths[seq][1] - av_read_len)
        allpositions2.append(genome_covbases_lengths[seq][1] - av_read_len)

        # Calculate metrics for the current genome
        window_size = round(genome_covbases_lengths[seq][1] / (len(allpositions)))

        if window_size > 0:
            distances = np.array(allpositions[1:]) - np.array(allpositions[:-1])
            hist = np.histogram(distances, bins=range(window_size, np.max(distances) + 1))[0]
            frequency = hist / (len(allpositions) - 1)
            gennext_read_ratio = (window_size - (np.sum(frequency * np.arange(len(frequency))))) / window_size
        else:
            gennext_read_ratio = 'nan'

        if not args.unpaired:
            window_size = round(genome_covbases_lengths[seq][1] / (len(allpositions2)))
            if window_size > 0:
                distances = np.array(allpositions2[1:]) - np.array(allpositions2[:-1])
                hist, edges = np.histogram(distances, bins=range(window_size, np.max(distances) + 1))
                frequency = hist / (len(allpositions2) - 1)
                gennext_read_ratio2 = (window_size - (np.sum(frequency * np.arange(len(frequency))))) / window_size
            else:
                gennext_read_ratio2 = 'nan'
        else:
            gennext_read_ratio2 = 'unpaired'

        coverage = genome_covbases_lengths[seq][2] / genome_covbases_lengths[seq][1]
        breadth = (genome_covbases_lengths[seq][1] - genome_covbases_lengths[seq][0]) / genome_covbases_lengths[seq][1]
        beb_ratio = breadth / (1 - np.exp(-0.883 * coverage))

        printing_dict[seq] = f"{genome_covbases_lengths[seq][1]}\t{coverage}\t{breadth}\t{beb_ratio}\t{gennext_read_ratio}\t{gennext_read_ratio2}\t{genome_readcount}"
        local_logger.debug(f'Calculated metrics for {seq}: coverage={coverage:.4f}, breadth={breadth:.4f}, BER={beb_ratio:.4f}')

    sam.close()
    return printing_dict