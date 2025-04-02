"""
Plotting utilities for metapresence.
"""

import sys
import logging
import os
from .interactive_plot import create_interactive_metrics_plot

# Get logger
logger = logging.getLogger('metapresence')


def create_metrics_plot(metrics_file, output_plot, min_reads, unpaired, interactive=False):
    """
    Create a scatterplot of BER vs FUG metrics.

    Args:
        metrics_file: Path to the metrics file.
        output_plot: Path to the output plot file.
        min_reads: Minimum read threshold.
        unpaired: Whether reads are unpaired.
        interactive: Whether to create an interactive plot.
    """
    if interactive:
        # Change output file extension to .html for interactive plots
        base_name = os.path.splitext(output_plot)[0]
        output_html = base_name + '.html'

        # Create interactive plot
        create_interactive_metrics_plot(metrics_file, output_html, min_reads, unpaired)
        return

    # Traditional static plot with matplotlib
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        logger.error('Cannot generate plot: matplotlib is not installed')
        logger.info('Install matplotlib with: pip install matplotlib')
        return

    logger.info(f"Creating metrics plot from {metrics_file}")
    logger.debug(f"Plot settings: min_reads={min_reads}, unpaired={unpaired}")

    fug = []
    ber = []
    genome_labels = []
    coverage_values = []

    try:
        logger.debug(f"Reading metrics from {metrics_file}")
        with open(metrics_file) as f:
            for line_num, line in enumerate(f, 1):
                parts = line.strip().split('\t')

                if line_num == 1 or len(parts) < 8:
                    continue

                genome = parts[0]

                # Skip entries with too few reads
                read_count = float(parts[7])
                if read_count < min_reads:
                    logger.debug(f"Skipping {genome} with only {read_count} reads (threshold: {min_reads})")
                    continue

                # Get BER value
                try:
                    ber_value = float(parts[4])
                    ber.append(ber_value)
                except ValueError:
                    logger.warning(f"Invalid BER value for {genome}: {parts[4]}")
                    continue

                # Get FUG value(s)
                try:
                    if not unpaired and parts[6] != 'unpaired' and parts[6] != 'nan' and parts[5] != 'nan':
                        fug1 = float(parts[5])
                        fug2 = float(parts[6])
                        fug_value = (fug1 + fug2) / 2
                        logger.debug(f"Using mean FUG for {genome}: {fug_value:.4f} (FUG1={fug1:.4f}, FUG2={fug2:.4f})")
                    elif parts[5] != 'nan' and parts[5] != 'unpaired':
                        fug_value = float(parts[5])
                        logger.debug(f"Using FUG1 for {genome}: {fug_value:.4f}")
                    else:
                        logger.debug(f"Skipping {genome} with invalid FUG values: FUG1={parts[5]}, FUG2={parts[6]}")
                        continue

                    fug.append(fug_value)
                    genome_labels.append(genome)

                    # Get coverage value for coloring
                    coverage_values.append(float(parts[2]))

                except ValueError:
                    logger.warning(f"Invalid FUG values for {genome}: FUG1={parts[5]}, FUG2={parts[6]}")
                    continue
    except Exception as e:
        logger.error(f"Error reading metrics file: {str(e)}")
        return

    if len(ber) == 0 or len(fug) == 0:
        logger.warning("No valid data points found for plotting")
        return

    logger.debug(f"Plotting {len(ber)} data points")

    # Create the plot
    plt.figure(figsize=(10, 8))

    # Create scatter plot with coverage-based coloring
    scatter = plt.scatter(ber, fug, alpha=0.7, c=coverage_values, cmap='viridis',
                         s=[max(20, min(300, c*20)) for c in coverage_values])

    # Add colorbar
    cbar = plt.colorbar(scatter)
    cbar.set_label('Coverage')

    plt.xlabel('BER (Breadth-Expected Breadth Ratio)')
    plt.ylabel('FUG (Fraction of Unexpected Gaps)' if unpaired else 'Mean FUG')

    plt.title('BER vs FUG Metrics')
    plt.grid(True, alpha=0.3)
    plt.xlim(0, max(1.1, max(ber) * 1.1))
    plt.ylim(0, max(1.1, max(fug) * 1.1))

    # Add reference lines for thresholds
    plt.axhline(y=0.5, color='red', linestyle='--', alpha=0.5, label='Default FUG threshold (0.5)')
    plt.axvline(x=0.8, color='red', linestyle='--', alpha=0.5, label='Default BER threshold (0.8)')

    # Add legend
    plt.legend(loc='lower right')

    # Add annotations for the highest coverage genomes (top 5)
    if len(genome_labels) > 0:
        top_indices = sorted(range(len(coverage_values)), key=lambda i: coverage_values[i], reverse=True)[:5]
        for idx in top_indices:
            plt.annotate(genome_labels[idx].split('.')[0],
                        (ber[idx], fug[idx]),
                        xytext=(5, 5),
                        textcoords='offset points',
                        fontsize=8)

    plt.tight_layout()

    try:
        logger.info(f"Saving plot to {output_plot}")
        plt.savefig(output_plot, dpi=300)
        plt.close()
        logger.info(f"Plot saved to {output_plot}")
    except Exception as e:
        logger.error(f"Error saving plot: {str(e)}")

    # Clean up
    plt.close('all')