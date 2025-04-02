"""
Plotting utilities for metapresence.
"""

import sys


def create_metrics_plot(metrics_file, output_plot, min_reads, unpaired):
    """
    Create a scatterplot of BER vs FUG metrics.

    Args:
        metrics_file: Path to the metrics file.
        output_plot: Path to the output plot file.
        min_reads: Minimum read threshold.
        unpaired: Whether reads are unpaired.
    """
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        print('Warning: not able to generate the plot. Matplotlib is not installed.')
        print('Install matplotlib with: pip install matplotlib')
        return

    fug = []
    ber = []

    try:
        with open(metrics_file) as f:
            for line in f:
                parts = line.strip().split('\t')

                if parts[0] == 'genome' or len(parts) < 8:
                    continue

                # Skip entries with too few reads
                if float(parts[7]) < min_reads:
                    continue

                # Get BER value
                ber_value = float(parts[4])
                ber.append(ber_value)

                # Get FUG value(s)
                if not unpaired and parts[6] != 'unpaired' and parts[6] != 'nan' and parts[5] != 'nan':
                    fug_value = (float(parts[5]) + float(parts[6])) / 2
                elif parts[5] != 'nan' and parts[5] != 'unpaired':
                    fug_value = float(parts[5])
                else:
                    # Skip entries with invalid FUG
                    continue

                fug.append(fug_value)
    except Exception as e:
        print(f"Error reading metrics file: {str(e)}")
        return

    if len(ber) == 0 or len(fug) == 0:
        print("No valid data points found for plotting.")
        return

    # Create the plot
    plt.figure(figsize=(10, 8))
    plt.scatter(ber, fug, alpha=0.7)

    plt.xlabel('BER (Breadth-Expected Breadth Ratio)')
    plt.ylabel('FUG (Fraction of Unexpected Gaps)' if unpaired else 'Mean FUG')

    plt.title('BER vs FUG Metrics')
    plt.grid(True, alpha=0.3)
    plt.xlim(0, max(1.1, max(ber) * 1.1))
    plt.ylim(0, max(1.1, max(fug) * 1.1))

    # Add reference lines
    plt.axhline(y=0.5, color='red', linestyle='--', alpha=0.5)  # Default FUG threshold
    plt.axvline(x=0.8, color='red', linestyle='--', alpha=0.5)  # Default BER threshold

    plt.tight_layout()

    try:
        plt.savefig(output_plot)
        plt.close()
        print(f"Plot saved to {output_plot}")
    except Exception as e:
        print(f"Error saving plot: {str(e)}")