"""
Command-line interface for metapresence.
"""

import argparse
import sys
import os
from .core import process_alignment
from .utils.merge import merge_abundances
from .utils.parsing import parse_metrics
from .utils.logging_config import setup_logging


def main():
    """Main entry point for the metapresence application."""
    parser = argparse.ArgumentParser(description="Calculation of metrics for evaluating the distribution of aligned reads onto nucleotidic sequences.")
    parser.add_argument("input_fasta", help="Either a folder of fasta files containing the sequences to analyze, or a (multi)fasta file containing the sequences to analyze. For large reference sizes, an input folder and more than one process are recommended for speeding-up.")
    parser.add_argument("indexed_sorted_bam", help="Indexed sorted bam file, with bam and index (.bai) in the same directory")
    parser.add_argument("-input_bins", metavar='input_bins.txt', help="A text file with each line listing a contig and a bin name, tab-seperated. This file is needed only if <input_fasta> is a fasta file and not a folder of fasta files.")
    parser.add_argument("-p", type=int, metavar='processes', default=1, help="number of processes to be executed in parallel. Default: 1")
    parser.add_argument("-o", metavar='output_prefix', default='metapout', help="the prefix of the output files. Default: metapout")
    parser.add_argument("-ber_threshold", metavar='[float]', type=float, default=0.8, help="Breadth-Expected breadth Ratio threshold. All genomes with BER value below the threshold are considered absent. Default: 0.8")
    parser.add_argument("-fug_threshold", metavar='[float]', type=float, default=0.5, help="Fraction of Unexpected Gaps threshold. All genomes with coverage lower than <max_for_fug> and with FUG value - for both mates - below the threshold are considered absent. Default: 0.5")
    parser.add_argument("-min_reads", type=int, metavar='[int]', default=80, help="Number of mapped reads on a given genome below which it is considered absent. Default: 80")
    parser.add_argument("-max_for_fug", type=float, metavar='[float]', default=0.1, help="Coverage value above which only BER metric is used. Default=0.1")
    parser.add_argument("-fug_criterion", metavar='["all","any","mean"]', default="all", help="Write < all > if a present species must have the FUG values for both the group of mates above the threshold, < any > if only one FUG value, < mean > if the mean FUG value. Irrelevant if --unpaired is set. Default=all")

    # File format options
    parser.add_argument("--stream", help="Use streaming mode for large BAM files to reduce memory usage. May be slower but uses less memory. Default: FALSE", action="store_true")
    parser.add_argument("--batch_size", type=int, default=10000, help="Number of reads to process in each batch when using streaming mode. Default: 10000")

    # Basic options
    parser.add_argument("--unpaired", help='--unpaired: set this flag if the aligned reads are not paired. Default: FALSE', action="store_true")
    parser.add_argument("--all_contigs", help='--all_contigs: set this flag if each contiguous sequence in < input_fasta > should be evaluated independently and be reported in the output files. If this flag is set < -input_bins > is irrelevant. Default: FALSE', action="store_true")

    # Visualization options
    parser.add_argument("--plot_metrics", help='--plot_metrics: set this flag if a scatterplot of the metric values has to be generated. Default: FALSE', action="store_true")
    parser.add_argument("--interactive", help='--interactive: create interactive HTML plots instead of static images. Requires Plotly. Default: FALSE', action="store_true")
    parser.add_argument("--plot_format", choices=["png", "svg", "pdf"], default="png", help="Format for static plots. Default: png")

    # Memory management options
    parser.add_argument("--low_memory", help="Optimize for low memory usage. May be slower but uses less memory. Default: FALSE", action="store_true")

    # Logging options
    parser.add_argument("--quiet", help='--quiet: suppress informational messages. Default: FALSE', action="store_true")
    parser.add_argument("--verbose", help='--verbose: enable verbose logging. Default: FALSE', action="store_true")
    parser.add_argument('--version', '-v', '-version', action='version', version=f'metapresence 1.0.0')

    args = parser.parse_args()

    # Set up logging
    logger = setup_logging(verbose=args.verbose, quiet=args.quiet)

    # Check input validity
    if not os.path.isdir(args.input_fasta) and not args.input_bins and not args.all_contigs:
        logger.error('If <input_fasta> is not a folder of fasta files, <-input_bins> must be specified or <--all_contigs> must be set.')
        sys.exit(1)

    if args.fug_criterion not in ['all', 'any', 'mean']:
        logger.error(f'"{args.fug_criterion}" is not a valid argument for -fug_criterion. Write either "all", "any" or "mean"')
        sys.exit(1)

    # Check for required packages based on options
    if args.interactive:
        try:
            import plotly.graph_objects
            logger.info("Using interactive plots with Plotly")
        except ImportError:
            logger.warning("Plotly is not installed. Interactive plots will not be generated.")
            logger.info("Install plotly with: pip install plotly")
            args.interactive = False

    if args.low_memory:
        logger.info("Running in low memory mode - this may be slower but will use less RAM")

    if args.stream:
        logger.info(f"Using streaming mode for BAM processing with batch size {args.batch_size}")

    process_alignment(args)


def merge():
    """Entry point for merging abundance outputs."""
    parser = argparse.ArgumentParser(description='Merge together multiple abundance outputs from metapresence.py into a single table')
    parser.add_argument("-abundance_folder", help="Folder containing abundance outputs from metapresence.py", required=True, metavar='')
    parser.add_argument("-output_table", help="Output file to store merged abundances", required=True, metavar='')

    # Visualization options for merged data
    parser.add_argument("--visualize", help="Generate visualization of merged abundances", action="store_true")
    parser.add_argument("--interactive", help="Create interactive HTML visualization. Requires Plotly", action="store_true")

    # Standard options
    parser.add_argument("--quiet", help='--quiet: suppress informational messages. default: FALSE', action="store_true")
    parser.add_argument("--verbose", help='--verbose: enable verbose logging. default: FALSE', action="store_true")
    parser.add_argument("--low_memory", help="Optimize for low memory usage. May be slower but uses less memory. Default: FALSE", action="store_true")

    args = parser.parse_args()

    # Set up logging
    logger = setup_logging(verbose=args.verbose, quiet=args.quiet)

    # Check for required packages based on options
    if args.interactive:
        try:
            import plotly.graph_objects
        except ImportError:
            logger.warning("Plotly is not installed. Interactive plots will not be generated.")
            logger.info("Install plotly with: pip install plotly")
            args.interactive = False

    merge_abundances(args.abundance_folder, args.output_table, low_memory=args.low_memory)

    # TODO: Add visualization of merged abundance data if requested
    if args.visualize:
        try:
            from .utils.visualize_merged import create_merged_visualization
            create_merged_visualization(args.output_table, interactive=args.interactive)
        except ImportError as e:
            logger.error(f"Could not generate visualization: {str(e)}")


def parse():
    """Entry point for parsing metrics output."""
    parser = argparse.ArgumentParser(description='Parse metrics file to generate abundance file. Useful to test different metrics thresholds without parsing the bam file.')
    parser.add_argument('-metrics', help='Output of metapresence.py containing the metric values (*_metrics.tsv).', required=True)
    parser.add_argument('-output_abundances', help='Output file to store the abundance values.', required=True)
    parser.add_argument("-ber_threshold", metavar='[float]', type=float, default=0.8, help="Breadth-Expected breadth Ratio threshold. All genomes with BER value below the threshold are considered absent. Default: 0.8")
    parser.add_argument("-fug_threshold", metavar='[float]', type=float, default=0.5, help="Fraction of Unexpected Gaps threshold. All genomes with coverage lower than <max_for_fug> and with FUG value - for both mates - below the threshold are considered absent. Default: 0.5")
    parser.add_argument("-min_reads", type=int, metavar='[int]', default=80, help="Number of mapped reads on a given genome below which it is considered absent. Default: 80")
    parser.add_argument("-max_for_fug", type=float, metavar='[float]', default=0.1, help="Coverage value above which only BER metric is used. Default=0.1")
    parser.add_argument("-fug_criterion", metavar='["all","any","mean"]', default="all", help="Write < all > if a present species must have the FUG values for both the group of mates above the threshold, < any > if only one FUG value, < mean > if the mean FUG value. Irrelevant if --unpaired is set. Default=all")

    # Standard options
    parser.add_argument("--unpaired", help='--unpaired: set this flag if the aligned reads are not paired. Default: FALSE', action="store_true")
    parser.add_argument("--plot_metrics", help='--plot_metrics: set this flag if a scatterplot of the metric values has to be generated. Requires Matplotlib. Default: FALSE', action="store_true")
    parser.add_argument("--interactive", help='--interactive: create interactive HTML plots instead of static images. Requires Plotly. Default: FALSE', action="store_true")
    parser.add_argument("--quiet", help='--quiet: suppress informational messages. default: FALSE', action="store_true")
    parser.add_argument("--verbose", help='--verbose: enable verbose logging. default: FALSE', action="store_true")
    parser.add_argument("--low_memory", help="Optimize for low memory usage. May be slower but uses less memory. Default: FALSE", action="store_true")

    args = parser.parse_args()

    # Set up logging
    logger = setup_logging(verbose=args.verbose, quiet=args.quiet)

    if args.fug_criterion not in ['all', 'any', 'mean']:
        logger.error(f'"{args.fug_criterion}" is not a valid argument for -fug_criterion. Write either "all", "any" or "mean"')
        sys.exit(1)

    # Check for required packages based on options
    if args.interactive:
        try:
            import plotly.graph_objects
        except ImportError:
            logger.warning("Plotly is not installed. Interactive plots will not be generated.")
            logger.info("Install plotly with: pip install plotly")
            args.interactive = False

    parse_metrics(args)


if __name__ == "__main__":
    main()