"""
Utilities for merging abundance outputs.
"""

import os
import logging
import sys

# Get logger
logger = logging.getLogger('metapresence')


def merge_abundances(abundance_folder, output_table):
    """
    Merge multiple abundance outputs into a single table.

    Args:
        abundance_folder: Folder containing abundance outputs.
        output_table: Path to the output table file.
    """
    if not os.path.isdir(abundance_folder):
        logger.error(f"Error: {abundance_folder} is not a valid directory")
        sys.exit(1)

    logger.info(f"Scanning directory: {abundance_folder}")
    files = [x for x in sorted(os.listdir(abundance_folder))
             if os.path.isfile(os.path.join(abundance_folder, x))]

    if not files:
        logger.error(f"Error: No files found in {abundance_folder}")
        sys.exit(1)

    logger.info(f"Found {len(files)} files to merge")
    logger.debug(f"Files to merge: {', '.join(files)}")

    printing_dict = {}
    skipped_files = 0

    for n, file in enumerate(files):
        file_path = os.path.join(abundance_folder, file)
        try:
            logger.debug(f"Processing file: {file}")
            with open(file_path) as f:
                # Check for and skip header
                header = f.readline().strip()
                if not header.startswith('genome'):
                    logger.warning(f"File {file} may not be an abundance file (unexpected header: {header})")

                genome_count = 0
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) < 2:
                        logger.warning(f"Line in {file} doesn't have enough columns, skipping")
                        continue

                    genome = parts[0]
                    abundance = parts[1]
                    genome_count += 1

                    if genome in printing_dict:
                        printing_dict[genome].append(abundance)
                    else:
                        printing_dict[genome] = ['0.0'] * n
                        printing_dict[genome].append(abundance)

            # Ensure all genomes have the same number of values
            for genome in printing_dict:
                if len(printing_dict[genome]) < n + 1:
                    printing_dict[genome].append('0.0')

            logger.debug(f"Added {genome_count} genomes from {file}")

        except Exception as e:
            logger.error(f"Error processing file {file}: {str(e)}")
            skipped_files += 1

    if skipped_files > 0:
        logger.warning(f"Skipped {skipped_files} files due to errors")

    if not printing_dict:
        logger.error("No data found in files. Cannot create merged table.")
        sys.exit(1)

    # Get total unique genomes
    total_genomes = len(printing_dict)
    logger.info(f"Found {total_genomes} unique genomes across all files")

    # Sort genomes alphabetically
    printing_dict = {genome: printing_dict[genome] for genome in sorted(printing_dict.keys())}

    logger.info(f"Writing merged table to {output_table}")
    try:
        with open(output_table, 'w') as new:
            new.write(f'fasta\t{"\t".join(files)}\n')

            for genome in printing_dict:
                new.write(f'{genome}\t{"\t".join(printing_dict[genome])}\n')
    except Exception as e:
        logger.error(f"Error writing output file: {str(e)}")
        sys.exit(1)

    logger.info(f"Successfully merged {len(files)} files into {output_table}")