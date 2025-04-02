"""
Parsing utilities for metapresence.
"""

import os
import sys
import logging
from ..metrics import parse_metrics_criteria

# Get logger
logger = logging.getLogger('metapresence')

try:
    import numpy as np
except ImportError:
    logger.error('Error: numpy is not installed.')
    sys.exit(1)


def chunks(l, nproc):
    """
    Divide a list into n sublists.

    Args:
        l: The list to divide.
        nproc: Number of chunks.

    Returns:
        List of sublists.
    """
    if not l or nproc <= 0:
        return []

    if nproc > len(l):
        logger.debug(f"Reducing number of processes from {nproc} to {len(l)} to match list length")
        nproc = len(l)

    right_div = len(l) // nproc
    nmore = len(l) % nproc

    result = [l[i * (right_div + 1):i * (right_div + 1) + right_div + 1] for i in range(nmore)]
    result += [l[nmore * (right_div + 1) + i * right_div:nmore * (right_div + 1) + i * right_div + right_div]
              for i in range(nproc - nmore) if nmore < len(l)]

    logger.debug(f"Split list of {len(l)} items into {len(result)} chunks")
    return result


def parse_single_fasta(fastafile, input_bins, all_contigs):
    """
    Parse a single fasta file.

    Args:
        fastafile: Path to the fasta file.
        input_bins: Path to the bins file.
        all_contigs: Whether to treat all contigs as independent.

    Returns:
        Tuple of dictionaries: (contigs_lengths, genome_contigs)
    """
    logger.debug(f"Parsing single fasta file: {fastafile}")
    genome_contigs = {}

    if not all_contigs and input_bins:
        logger.debug(f"Reading bins from: {input_bins}")
        with open(input_bins) as f:
            for line_num, line in enumerate(f, 1):
                a = line.strip().split('\t')
                if len(a) < 2:
                    logger.warning(f"Line {line_num} in {input_bins} doesn't have tab-separated values, skipping")
                    continue
                if a[1] in genome_contigs:
                    genome_contigs[a[1]].append(a[0])
                else:
                    genome_contigs[a[1]] = [a[0]]
        logger.debug(f"Read {len(genome_contigs)} bin entries")

    contigs_lengths_d = {}

    with open(fastafile) as f:
        contig_name = ''
        contig_count = 0
        total_length = 0

        for line in f:
            if line[0] == '>':
                contig_count += 1
                header = line[1:].strip().split()
                contig_name = header[0]
                if all_contigs:
                    genome_contigs[contig_name] = [contig_name]
                contigs_lengths_d[contig_name] = 0
            else:
                contigs_lengths_d[contig_name] += len(line.strip())
                total_length += len(line.strip())

    logger.debug(f"Read {contig_count} contigs with total length {total_length:,} bp")
    return contigs_lengths_d, genome_contigs


def parse_multi_fasta(args):
    """
    Parse multiple fasta files.

    Args:
        args: Tuple of (file_list, all_contigs).

    Returns:
        Tuple of dictionaries: (contigs_lengths, genome_contigs)
    """
    ff, all_contigs = args
    contigs_lengths_d, genome_contigs = {}, {}
    local_logger = logging.getLogger('metapresence')

    total_contigs = 0
    total_length = 0

    for fa in ff:
        local_logger.debug(f"Parsing file: {fa}")
        try:
            with open(fa) as f:
                fa_name = fa.split('/')[-1]
                contig_name = None
                file_contigs = 0
                file_length = 0

                for line in f:
                    if line[0] == '>':
                        file_contigs += 1
                        total_contigs += 1
                        contig_name = line.strip().split()[0][1:]
                        contigs_lengths_d[contig_name] = 0

                        if not all_contigs:
                            if fa_name in genome_contigs:
                                genome_contigs[fa_name].append(contig_name)
                            else:
                                genome_contigs[fa_name] = [contig_name]
                        else:
                            genome_contigs[contig_name] = [contig_name]
                    else:
                        seq_len = len(line.strip())
                        contigs_lengths_d[contig_name] += seq_len
                        file_length += seq_len
                        total_length += seq_len

                local_logger.debug(f"File {fa_name}: {file_contigs} contigs, {file_length:,} bp")
        except Exception as e:
            local_logger.error(f"Error parsing file {fa}: {str(e)}")

    local_logger.debug(f"Processed {len(ff)} files with {total_contigs} contigs and {total_length:,} bp total")
    return contigs_lengths_d, genome_contigs


def parse_metrics(args):
    """
    Parse metrics file to generate abundance file.

    Args:
        args: Command line arguments containing metrics file path, output path, and threshold parameters.
    """
    logger.info(f"Parsing metrics file: {args.metrics}")
    present_coverage = {}

    try:
        with open(args.metrics) as f:
            # Skip header
            header = f.readline().strip()
            logger.debug(f"Header: {header}")

            for line_num, line in enumerate(f, 2):
                parts = line.strip().split('\t')
                if len(parts) < 8:  # Make sure we have all required fields
                    logger.warning(f"Line {line_num} doesn't have enough columns, skipping")
                    continue

                genome = parts[0]
                cov = float(parts[2])
                ber = float(parts[4])

                if parts[5] == 'nan' or parts[5] == 'unpaired':
                    logger.debug(f"Genome {genome} has FUG1 value of {parts[5]}")
                    fug1 = 0
                else:
                    fug1 = float(parts[5])

                if args.unpaired or parts[6] == 'unpaired':
                    fug2 = None
                else:
                    if parts[6] == 'nan':
                        logger.debug(f"Genome {genome} has FUG2 value of {parts[6]}")
                        fug2 = 0
                    else:
                        fug2 = float(parts[6])

                nreads = int(parts[7])

                logger.debug(f"Genome {genome}: cov={cov:.4f}, BER={ber:.4f}, reads={nreads}")

                if parse_metrics_criteria(cov, ber, fug1, fug2, nreads, args.fug_criterion,
                                         args.min_reads, args.max_for_fug, args.ber_threshold,
                                         args.fug_threshold, args.unpaired):
                    present_coverage[genome] = cov
                    logger.debug(f"Genome {genome} is considered present")
                else:
                    logger.debug(f"Genome {genome} is not considered present")

    except Exception as e:
        logger.error(f"Error reading metrics file: {str(e)}")
        sys.exit(1)

    logger.info(f"Found {len(present_coverage)} present genomes")
    logger.info(f"Writing abundance results to {args.output_abundances}")

    try:
        with open(args.output_abundances, "w") as out_file:
            out_file.write('genome\trelative_abundance_%\n')

            if present_coverage:
                totcov = sum(present_coverage.values())
                for genome in present_coverage:
                    abundance = present_coverage[genome] / totcov * 100
                    out_file.write(f'{genome}\t{abundance}\n')
                    logger.debug(f"Genome {genome}: {abundance:.4f}%")
            else:
                logger.warning("No genomes were determined to be present based on the criteria")

    except Exception as e:
        logger.error(f"Error writing to output file: {str(e)}")
        sys.exit(1)

    logger.info("Abundance calculation complete")