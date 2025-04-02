"""
Parsing utilities for metapresence.
"""

import os
import sys
from ..metrics import parse_metrics_criteria

try:
    import numpy as np
except ImportError:
    print('Error: numpy is not installed.')
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
        nproc = len(l)

    right_div = len(l) // nproc
    nmore = len(l) % nproc

    result = [l[i * (right_div + 1):i * (right_div + 1) + right_div + 1] for i in range(nmore)]
    result += [l[nmore * (right_div + 1) + i * right_div:nmore * (right_div + 1) + i * right_div + right_div]
              for i in range(nproc - nmore) if nmore < len(l)]

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
    genome_contigs = {}

    if not all_contigs and input_bins:
        with open(input_bins) as f:
            for line in f:
                a = line.strip().split('\t')
                if len(a) < 2:
                    continue
                if a[1] in genome_contigs:
                    genome_contigs[a[1]].append(a[0])
                else:
                    genome_contigs[a[1]] = [a[0]]

    contigs_lengths_d = {}

    with open(fastafile) as f:
        contig_name = ''
        for line in f:
            if line[0] == '>':
                header = line[1:].strip().split()
                contig_name = header[0]
                if all_contigs:
                    genome_contigs[contig_name] = [contig_name]
                contigs_lengths_d[contig_name] = 0
            else:
                contigs_lengths_d[contig_name] += len(line.strip())

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

    for fa in ff:
        with open(fa) as f:
            fa_name = fa.split('/')[-1]
            contig_name = None

            for line in f:
                if line[0] == '>':
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
                    contigs_lengths_d[contig_name] += len(line.strip())

    return contigs_lengths_d, genome_contigs


def parse_metrics(args):
    """
    Parse metrics file to generate abundance file.

    Args:
        args: Command line arguments containing metrics file path, output path, and threshold parameters.
    """
    present_coverage = {}

    with open(args.metrics) as f:
        # Skip header
        f.readline()

        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 8:  # Make sure we have all required fields
                continue

            genome = parts[0]
            cov = float(parts[2])
            ber = float(parts[4])
            fug1 = float(parts[5]) if parts[5] != 'nan' and parts[5] != 'unpaired' else 0

            if args.unpaired or parts[6] == 'unpaired':
                fug2 = None
            else:
                fug2 = float(parts[6]) if parts[6] != 'nan' else 0

            nreads = int(parts[7])

            if parse_metrics_criteria(cov, ber, fug1, fug2, nreads, args.fug_criterion,
                                     args.min_reads, args.max_for_fug, args.ber_threshold,
                                     args.fug_threshold, args.unpaired):
                present_coverage[genome] = cov

    with open(args.output_abundances, "w") as new:
        new.write('genome\trelative_abundance_%\n')

        if present_coverage:
            totcov = sum(present_coverage.values())
            for genome in present_coverage:
                new.write(f'{genome}\t{present_coverage[genome] / totcov * 100}\n')