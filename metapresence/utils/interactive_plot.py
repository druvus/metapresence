"""
Interactive plotting utilities for metapresence using Plotly.
"""

import os
import logging
import json
import sys
import traceback

# Get logger
logger = logging.getLogger('metapresence')


def create_interactive_metrics_plot(metrics_file, output_html, min_reads, unpaired):
    """
    Create an interactive scatterplot of BER vs FUG metrics using Plotly.

    Args:
        metrics_file: Path to the metrics file.
        output_html: Path to the output HTML file.
        min_reads: Minimum read threshold.
        unpaired: Whether reads are unpaired.
    """
    try:
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
        import numpy as np
    except ImportError:
        logger.error('Cannot generate interactive plot: plotly is not installed')
        logger.info('Install plotly with: pip install plotly')
        return

    logger.info(f"Creating interactive metrics plot from {metrics_file}")
    logger.debug(f"Plot settings: min_reads={min_reads}, unpaired={unpaired}")

    # Data structures to store values for plotting
    metrics_data = []

    try:
        if not os.path.exists(metrics_file):
            logger.error(f"Metrics file not found: {metrics_file}")
            return

        logger.debug(f"Reading metrics from {metrics_file}")
        with open(metrics_file) as f:
            # Get header
            header = f.readline().strip().split('\t')
            logger.debug(f"Header: {header}")

            # Process each line (each genome)
            line_count = 0
            valid_data_count = 0

            for line in f:
                line_count += 1
                parts = line.strip().split('\t')

                if len(parts) < 8:
                    logger.warning(f"Line has insufficient columns (expected 8, got {len(parts)}): {line.strip()}")
                    continue

                # Extract values
                try:
                    genome = parts[0]
                    length = float(parts[1])
                    coverage = float(parts[2])
                    breadth = float(parts[3])
                    ber = float(parts[4])

                    # Parse FUG values
                    if parts[5] != 'nan' and parts[5] != 'unpaired':
                        fug1 = float(parts[5])
                    else:
                        fug1 = None
                        logger.debug(f"Invalid FUG1 value for {genome}: {parts[5]}")

                    if not unpaired and parts[6] != 'nan' and parts[6] != 'unpaired':
                        fug2 = float(parts[6])
                    else:
                        fug2 = None
                        if not unpaired:
                            logger.debug(f"Invalid FUG2 value for {genome}: {parts[6]}")

                    read_count = int(parts[7])

                    # Skip entries with too few reads
                    if read_count < min_reads:
                        logger.debug(f"Skipping {genome} with only {read_count} reads (threshold: {min_reads})")
                        continue

                    # Calculate mean FUG for paired data
                    if fug1 is not None and fug2 is not None and not unpaired:
                        fug_mean = (fug1 + fug2) / 2
                    elif fug1 is not None:
                        fug_mean = fug1
                    else:
                        logger.debug(f"Skipping {genome} with no valid FUG values")
                        continue  # Skip if no valid FUG values

                    # Store all genome data for the plot
                    metrics_data.append({
                        'genome': genome,
                        'length': length,
                        'coverage': coverage,
                        'breadth': breadth,
                        'ber': ber,
                        'fug1': fug1,
                        'fug2': fug2,
                        'fug_mean': fug_mean,
                        'read_count': read_count
                    })
                    valid_data_count += 1

                except (ValueError, IndexError) as e:
                    logger.warning(f"Error parsing line: {line.strip()} - {str(e)}")

            logger.info(f"Read {line_count} lines from metrics file, found {valid_data_count} valid data points")

    except Exception as e:
        logger.error(f"Error reading metrics file: {str(e)}")
        logger.debug(traceback.format_exc())
        return

    if not metrics_data:
        logger.warning("No valid data points found for plotting. Check your min_reads threshold and metrics file.")
        return

    logger.debug(f"Plotting {len(metrics_data)} data points")

    # Create a standard non-interactive plot first as a fallback
    try:
        import matplotlib.pyplot as plt

        plt.figure(figsize=(10, 8))

        # Extract data
        ber_values = [d['ber'] for d in metrics_data]
        fug_values = [d['fug_mean'] for d in metrics_data]
        coverage_values = [d['coverage'] for d in metrics_data]

        # Create scatter plot
        plt.scatter(ber_values, fug_values, c=coverage_values, cmap='viridis')
        plt.colorbar(label='Coverage')
        plt.xlabel('BER (Breadth-Expected Breadth Ratio)')
        plt.ylabel('FUG (Fraction of Unexpected Gaps)')
        plt.title('BER vs FUG Metrics')

        # Add reference lines
        plt.axhline(y=0.5, color='red', linestyle='--')
        plt.axvline(x=0.8, color='red', linestyle='--')

        # Save the plot
        fallback_plot = os.path.splitext(output_html)[0] + "_fallback.png"
        plt.savefig(fallback_plot)
        plt.close()

        logger.info(f"Created fallback static plot at {fallback_plot}")
    except Exception as e:
        logger.debug(f"Could not create fallback plot: {str(e)}")

    # Create interactive plot with Plotly
    try:
        # Save data for debugging
        debug_data_path = os.path.splitext(output_html)[0] + "_debug_data.json"
        with open(debug_data_path, 'w') as f:
            json.dump(metrics_data, f, indent=2)
        logger.debug(f"Saved debug data to {debug_data_path}")

        # Create figure with 2 subplots
        fig = make_subplots(
            rows=1, cols=2,
            subplot_titles=("BER vs FUG Metrics", "Coverage vs Breadth"),
            specs=[[{"type": "scatter"}, {"type": "scatter"}]]
        )

        # Extract data for plotting
        genome_names = [d['genome'] for d in metrics_data]
        ber_values = [d['ber'] for d in metrics_data]
        fug_values = [d['fug_mean'] for d in metrics_data]
        coverage_values = [d['coverage'] for d in metrics_data]
        breadth_values = [d['breadth'] for d in metrics_data]
        read_counts = [d['read_count'] for d in metrics_data]

        # Log data ranges to help with debugging
        logger.debug(f"BER range: {min(ber_values)} to {max(ber_values)}")
        logger.debug(f"FUG range: {min(fug_values)} to {max(fug_values)}")
        logger.debug(f"Coverage range: {min(coverage_values)} to {max(coverage_values)}")

        # Calculate marker sizes based on read count
        marker_sizes = []
        for count in read_counts:
            # Log scale for marker size to avoid extremely large bubbles
            size = 5 + 5 * np.log10(max(1, count / min_reads))
            marker_sizes.append(size)

        # Create colorscale for the points based on coverage
        colorscale = 'Viridis'

        # Add BER vs FUG scatter plot
        fig.add_trace(
            go.Scatter(
                x=ber_values,
                y=fug_values,
                mode='markers',
                marker=dict(
                    size=marker_sizes,
                    color=coverage_values,
                    colorscale=colorscale,
                    colorbar=dict(title="Coverage"),
                    showscale=True
                ),
                text=[f"Genome: {g}<br>BER: {b:.4f}<br>FUG: {f:.4f}<br>Coverage: {c:.4f}<br>Reads: {r}"
                      for g, b, f, c, r in zip(genome_names, ber_values, fug_values, coverage_values, read_counts)],
                hoverinfo='text',
                name='Genomes'
            ),
            row=1, col=1
        )

        # Add Coverage vs Breadth scatter plot
        fig.add_trace(
            go.Scatter(
                x=coverage_values,
                y=breadth_values,
                mode='markers',
                marker=dict(
                    size=marker_sizes,
                    color=ber_values,
                    colorscale='Plasma',
                    colorbar=dict(title="BER"),
                    showscale=True
                ),
                text=[f"Genome: {g}<br>Coverage: {c:.4f}<br>Breadth: {b:.4f}<br>BER: {ber:.4f}<br>Reads: {r}"
                      for g, c, b, ber, r in zip(genome_names, coverage_values, breadth_values, ber_values, read_counts)],
                hoverinfo='text',
                name='Genomes'
            ),
            row=1, col=2
        )

        # Add reference lines for BER vs FUG plot
        fig.add_shape(
            type="line",
            x0=0.8, y0=0, x1=0.8, y1=1,
            line=dict(color="red", width=2, dash="dash"),
            row=1, col=1
        )
        fig.add_shape(
            type="line",
            x0=0, y0=0.5, x1=1, y1=0.5,
            line=dict(color="red", width=2, dash="dash"),
            row=1, col=1
        )

        # Add annotations for the reference lines
        fig.add_annotation(
            x=0.8, y=0.05,
            text="BER threshold (0.8)",
            showarrow=False,
            font=dict(color="red"),
            row=1, col=1
        )
        fig.add_annotation(
            x=0.2, y=0.5,
            text="FUG threshold (0.5)",
            showarrow=False,
            font=dict(color="red"),
            row=1, col=1
        )

        # Add reference line for expected breadth curve in Coverage vs Breadth plot
        x_curve = np.linspace(0, max(coverage_values) * 1.1, 100)
        y_curve = 1 - np.exp(-0.883 * x_curve)

        fig.add_trace(
            go.Scatter(
                x=x_curve,
                y=y_curve,
                mode='lines',
                line=dict(color='red', dash='dash'),
                name='Expected Breadth'
            ),
            row=1, col=2
        )

        # Update layout
        fig.update_layout(
            title_text="Metapresence Metrics Visualization",
            height=700,
            width=1200,
            legend=dict(
                orientation="h",
                yanchor="bottom",
                y=-0.2,
                xanchor="center",
                x=0.5
            )
        )

        # Update x and y axis labels
        fig.update_xaxes(title_text="BER (Breadth-Expected Breadth Ratio)", row=1, col=1)
        fig.update_yaxes(title_text="FUG (Fraction of Unexpected Gaps)", row=1, col=1)
        fig.update_xaxes(title_text="Coverage", row=1, col=2)
        fig.update_yaxes(title_text="Breadth", row=1, col=2)

        # Add buttons for enabling/disabling traces
        fig.update_layout(
            updatemenus=[
                dict(
                    type="buttons",
                    buttons=[
                        dict(
                            label="Reset View",
                            method="relayout",
                            args=[{"xaxis.autorange": True, "yaxis.autorange": True}]
                        )
                    ],
                    direction="left",
                    pad={"r": 10, "t": 10},
                    showactive=False,
                    x=0.1,
                    xanchor="left",
                    y=1.1,
                    yanchor="top"
                )
            ]
        )

        # Save the plot to HTML file
        fig.write_html(output_html)
        logger.info(f"Interactive plot saved to {output_html}")

        # Create a simplified self-contained dashboard
        create_simplified_dashboard(metrics_data, output_html)

    except Exception as e:
        logger.error(f"Error creating interactive plot: {str(e)}")
        logger.debug(traceback.format_exc())


def create_simplified_dashboard(metrics_data, output_base):
    """
    Create a simplified, self-contained dashboard HTML file.

    Args:
        metrics_data: List of dictionaries with metrics for each genome.
        output_base: Base path for output HTML file.
    """
    if not metrics_data:
        logger.warning("No data to create dashboard")
        return

    output_dir = os.path.dirname(output_base)
    base_name = os.path.splitext(os.path.basename(output_base))[0]

    # Create simplified dashboard HTML file
    dashboard_html = os.path.join(output_dir, f"{base_name}_simple_dashboard.html")

    try:
        # Serialize the data directly into the HTML
        data_json = json.dumps(metrics_data)

        html_content = f"""<!DOCTYPE html>
<html>
<head>
    <title>Metapresence Simple Dashboard</title>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 0; padding: 20px; }}
        #plot {{ width: 90%; height: 600px; margin: 20px auto; }}
        .controls {{ padding: 20px; background: #f5f5f5; margin-bottom: 20px; }}
        h1 {{ color: #333; text-align: center; }}
    </style>
</head>
<body>
    <h1>Metapresence Metrics Dashboard</h1>

    <div class="controls">
        <h3>Debug Information</h3>
        <p>Data points: <span id="data-count">0</span></p>
        <p>BER range: <span id="ber-range">N/A</span></p>
        <p>FUG range: <span id="fug-range">N/A</span></p>
    </div>

    <div id="plot"></div>

    <script>
        // Embedded data
        const data = {data_json};

        // Update debug info
        document.getElementById('data-count').textContent = data.length;

        if (data.length > 0) {{
            const berValues = data.map(d => d.ber);
            const fugValues = data.map(d => d.fug_mean);

            document.getElementById('ber-range').textContent =
                `${{Math.min(...berValues).toFixed(4)}} to ${{Math.max(...berValues).toFixed(4)}}`;
            document.getElementById('fug-range').textContent =
                `${{Math.min(...fugValues).toFixed(4)}} to ${{Math.max(...fugValues).toFixed(4)}}`;

            // Create plot
            const plotData = {{
                x: berValues,
                y: fugValues,
                mode: 'markers',
                type: 'scatter',
                text: data.map(d => `Genome: ${{d.genome}}<br>BER: ${{d.ber.toFixed(4)}}<br>FUG: ${{d.fug_mean.toFixed(4)}}<br>Coverage: ${{d.coverage.toFixed(4)}}`),
                marker: {{
                    size: data.map(d => 10 + 10 * Math.log10(Math.max(1, d.read_count / 80))),
                    color: data.map(d => d.coverage),
                    colorscale: 'Viridis',
                    colorbar: {{title: 'Coverage'}}
                }}
            }};

            const layout = {{
                title: 'BER vs FUG Metrics',
                xaxis: {{title: 'BER', range: [0, Math.max(1.1, Math.max(...berValues) * 1.1)]}},
                yaxis: {{title: 'FUG', range: [0, Math.max(1.1, Math.max(...fugValues) * 1.1)]}},
                shapes: [
                    {{
                        type: 'line',
                        x0: 0.8, y0: 0, x1: 0.8, y1: 1,
                        line: {{color: 'red', width: 2, dash: 'dash'}}
                    }},
                    {{
                        type: 'line',
                        x0: 0, y0: 0.5, x1: 1, y1: 0.5,
                        line: {{color: 'red', width: 2, dash: 'dash'}}
                    }}
                ],
                annotations: [
                    {{
                        x: 0.8, y: 0.05,
                        text: "BER threshold (0.8)",
                        showarrow: false,
                        font: {{color: "red"}}
                    }},
                    {{
                        x: 0.2, y: 0.5,
                        text: "FUG threshold (0.5)",
                        showarrow: false,
                        font: {{color: "red"}}
                    }}
                ]
            }};

            Plotly.newPlot('plot', [plotData], layout);

            console.log('Plot created with ' + data.length + ' data points');
        }} else {{
            document.getElementById('plot').innerHTML = '<p style="text-align: center; color: red;">No data available for plotting</p>';
            console.error('No data available for plotting');
        }}
    </script>
</body>
</html>
"""

        with open(dashboard_html, 'w') as f:
            f.write(html_content)

        logger.info(f"Simple dashboard saved to {dashboard_html}")

    except Exception as e:
        logger.error(f"Error creating simple dashboard: {str(e)}")
        logger.debug(traceback.format_exc())


def create_additional_interactive_plots(metrics_data, output_base):
    """
    Create additional interactive visualizations based on the metrics data.

    Args:
        metrics_data: List of dictionaries with metrics for each genome.
        output_base: Base path for output HTML files.
    """
    if not metrics_data:
        logger.warning("No data for additional plots")
        return

    try:
        import plotly.graph_objects as go
        import plotly.express as px
        import numpy as np
    except ImportError:
        logger.warning("Plotly not installed, skipping additional plots")
        return

    # Get base name without extension
    output_dir = os.path.dirname(output_base)
    base_name = os.path.splitext(os.path.basename(output_base))[0]

    # Create abundance summary visualization
    if metrics_data:
        try:
            # Sort by coverage (highest first)
            sorted_data = sorted(metrics_data, key=lambda x: x['coverage'], reverse=True)

            # Take top 20 genomes (or all if less than 20)
            top_limit = min(20, len(sorted_data))
            top_genomes = sorted_data[:top_limit]

            if not top_genomes:
                logger.warning("No genomes to plot for abundance chart")
                return

            # Calculate total coverage
            total_coverage = sum(d['coverage'] for d in sorted_data)

            if total_coverage == 0:
                logger.warning("Total coverage is zero, cannot create abundance chart")
                return

            # Calculate relative abundance
            for genome in top_genomes:
                genome['abundance'] = (genome['coverage'] / total_coverage) * 100

            # Create pie chart of relative abundances
            fig = go.Figure(data=[go.Pie(
                labels=[d['genome'] for d in top_genomes],
                values=[d['abundance'] for d in top_genomes],
                hole=0.4,
                textinfo='label+percent',
                insidetextorientation='radial',
                marker=dict(
                    colors=px.colors.qualitative.Plotly
                )
            )])

            fig.update_layout(
                title_text=f"Top {top_limit} Genomes by Relative Abundance",
                annotations=[dict(text='Relative<br>Abundance', x=0.5, y=0.5, font_size=20, showarrow=False)]
            )

            # Save abundance plot
            abundance_html = os.path.join(output_dir, f"{base_name}_abundance.html")
            fig.write_html(abundance_html)
            logger.info(f"Abundance visualization saved to {abundance_html}")

        except Exception as e:
            logger.error(f"Error creating abundance plot: {str(e)}")
            logger.debug(traceback.format_exc())