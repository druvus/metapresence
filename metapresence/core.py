"""
Core functionality for metapresence.
"""

import os
import sys
import multiprocessing
from .metrics import parse_metrics_criteria
from .utils.parsing import parse_single_fasta, parse_multi_fasta, chunks
from .utils.plot import create_metrics_plot

try:
    import numpy as np
except ImportError:
    print('Error: numpy is not installed.')
    sys.exit(1)

try:
    import pysam
except ImportError:
    print('Error: pysam is not installed.')
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
        print('Reading fasta files')
        fastas = [os.path.join(scaffold_input, x) for x in os.listdir(scaffold_input)]
        chunked_fastas = chunks(fastas, processes)
        parseFasta_args = [(x, args.all_contigs) for x in chunked_fastas]

        with multiprocessing.Pool(processes=processes) as pool:
            all_parsing = pool.map(parse_multi_fasta, parseFasta_args)

        contigs_lengths_d = {x: result[0][x] for result in all_parsing for x in result[0]}
        genome_contigs = {x: result[1][x] for result in all_parsing for x in result[1]}

        if args.all_contigs:
            genome_contigs = {y: [y] for x in genome_contigs for y in genome_contigs[x]}

    elif os.path.isfile(scaffold_input):
        print('Reading fasta file')
        contigs_lengths_d, genome_contigs = parse_single_fasta(scaffold_input, args.input_bins, args.all_contigs)

    else:
        print(f'Error: cannot open {scaffold_input}\naborted')
        sys.exit(1)

    print('\tdone\n')

    # Read BAM file and calculate average read length
    print('Reading bam file')
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
        sam.close()

    except ValueError:
        print(f'Error: cannot find index for {sorted_bam}. Index and bam file should be placed in the same directory.\nAborted')
        sys.exit(1)

    # Process BAM file and calculate metrics
    genomes = [x for x in genome_contigs]
    chunked_genomes = chunks(genomes, processes)

    # Prepare arguments for the multiprocessing pool
    pool_args = [(genome_chunk, sorted_bam, contigs_lengths_d, genome_contigs, av_read_len, args)
                 for genome_chunk in chunked_genomes]

    with multiprocessing.Pool(processes=processes) as pool:
        all_metrics = pool.map(_parse_bam_wrapper, pool_args)

    print('\tdone\n')
    print('preparing output files\n')

    # Generate output files
    metrics_dict = {}
    for result in all_metrics:
        if result:  # Check if result is not None
            for x in result:
                metrics_dict[x] = result[x]

    metrics_dict = {x: metrics_dict[x] for x in sorted(metrics_dict.keys())}

    # Write metrics file
    with open(args.o + "_metrics.tsv", "w") as metrics_file:
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

            metrics_file.write(f'{genome}\t{metrics_dict[genome]}\n')

    # Write abundances file
    with open(args.o + "_abundances.tsv", "w") as abundances_file:
        abundances_file.write('genome\trelative_abundance_%\n')

        if present_coverage:
            totcov = sum(present_coverage.values())
            for genome in present_coverage:
                abundances_file.write(f'{genome}\t{present_coverage[genome] / totcov * 100}\n')

    # Generate plot if requested
    if args.plot_metrics:
        create_metrics_plot(args.o + "_metrics.tsv", args.o + "_scatterplot.png",
                          args.min_reads, args.unpaired)

    print('All done!\n')


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
    printing_dict = {}
    sam = pysam.AlignmentFile(sorted_bam.split()[0], "rb")

    for seq in list_of_genomes:
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
                    print(f'Warning: contig {contig} of sequence {seq} was not present in the fasta file(s). Skipping sequence...')
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

    sam.close()
    return printing_dict