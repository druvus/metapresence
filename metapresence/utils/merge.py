"""
Utilities for merging abundance outputs with memory optimization.
"""

import os
import logging
import sys
import gc

# Get logger
logger = logging.getLogger('metapresence')


def merge_abundances(abundance_folder, output_table, low_memory=False):
    """
    Merge multiple abundance outputs into a single table with memory optimization.

    Args:
        abundance_folder: Folder containing abundance outputs.
        output_table: Path to the output table file.
        low_memory: Whether to optimize for low memory usage.
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

    # First pass: collect all unique genome names across all files
    if low_memory:
        # Low memory mode: process files sequentially without storing all data in memory
        return merge_abundances_low_memory(abundance_folder, files, output_table)
    else:
        # Standard mode: collect all data in memory
        return merge_abundances_standard(abundance_folder, files, output_table)


def merge_abundances_standard(abundance_folder, files, output_table):
    """
    Merge abundance files using standard memory usage.

    Args:
        abundance_folder: Folder containing abundance outputs.
        files: List of files to process.
        output_table: Path to the output table file.
    """
    printing_dict = {}
    skipped_files = 0

    # Process files in batches if there are many
    batch_size = 100
    file_batches = [files[i:i+batch_size] for i in range(0, len(files), batch_size)]

    for batch_idx, file_batch in enumerate(file_batches):
        logger.debug(f"Processing batch {batch_idx+1}/{len(file_batches)} with {len(file_batch)} files")

        for n, file in enumerate(file_batch, start=batch_idx*batch_size):
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
                            # Fill in gaps for missing files
                            while len(printing_dict[genome]) < n:
                                printing_dict[genome].append('0.0')
                            printing_dict[genome].append(abundance)
                        else:
                            printing_dict[genome] = ['0.0'] * n
                            printing_dict[genome].append(abundance)

                # Ensure all genomes have the same number of values for this file
                for genome in printing_dict:
                    if len(printing_dict[genome]) < n + 1:
                        printing_dict[genome].append('0.0')

                logger.debug(f"Added {genome_count} genomes from {file}")

            except Exception as e:
                logger.error(f"Error processing file {file}: {str(e)}")
                skipped_files += 1

        # Force garbage collection after each batch
        gc.collect()

    if skipped_files > 0:
        logger.warning(f"Skipped {skipped_files} files due to errors")

    if not printing_dict:
        logger.error("No data found in files. Cannot create merged table.")
        sys.exit(1)

    # Get total unique genomes
    total_genomes = len(printing_dict)
    logger.info(f"Found {total_genomes} unique genomes across all files")

    # Ensure all genomes have values for all files
    for genome in printing_dict:
        while len(printing_dict[genome]) < len(files):
            printing_dict[genome].append('0.0')

    # Sort genomes alphabetically
    printing_dict = {genome: printing_dict[genome] for genome in sorted(printing_dict.keys())}

    logger.info(f"Writing merged table to {output_table}")
    try:
        with open(output_table, 'w') as new:
            new.write(f'fasta\t{"\t".join(files)}\n')

            # Write genomes in batches to reduce memory pressure
            genome_batch_size = 1000
            genomes = list(printing_dict.keys())

            for i in range(0, len(genomes), genome_batch_size):
                for genome in genomes[i:i+genome_batch_size]:
                    new.write(f'{genome}\t{"\t".join(printing_dict[genome])}\n')

                # Clear batch data to free memory
                if i + genome_batch_size < len(genomes):
                    batch_genomes = genomes[i:i+genome_batch_size]
                    for g in batch_genomes:
                        del printing_dict[g]
                    gc.collect()

    except Exception as e:
        logger.error(f"Error writing output file: {str(e)}")
        sys.exit(1)

    logger.info(f"Successfully merged {len(files)} files into {output_table}")


def merge_abundances_low_memory(abundance_folder, files, output_table):
    """
    Merge abundance files using minimal memory (streaming approach).

    Args:
        abundance_folder: Folder containing abundance outputs.
        files: List of files to process.
        output_table: Path to the output table file.
    """
    logger.info("Using low memory mode for merging files")

    # First pass: collect all unique genome names
    all_genomes = set()

    for file in files:
        file_path = os.path.join(abundance_folder, file)
        try:
            with open(file_path) as f:
                # Skip header
                f.readline()

                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 1:
                        all_genomes.add(parts[0])
        except Exception as e:
            logger.error(f"Error reading file {file}: {str(e)}")

    logger.info(f"Found {len(all_genomes)} unique genomes across all files")

    if not all_genomes:
        logger.error("No genomes found in any files")
        sys.exit(1)

    # Sort genomes for consistent output
    sorted_genomes = sorted(all_genomes)

    # Now write directly to the output file without storing everything in memory
    with open(output_table, 'w') as out_file:
        # Write header
        out_file.write(f'fasta\t{"\t".join(files)}\n')

        # Process each genome separately
        for genome in sorted_genomes:
            abundances = []

            # Collect abundance values for this genome from each file
            for file in files:
                file_path = os.path.join(abundance_folder, file)
                abundance_found = False

                try:
                    with open(file_path) as f:
                        # Skip header
                        f.readline()

                        for line in f:
                            parts = line.strip().split('\t')
                            if len(parts) >= 2 and parts[0] == genome:
                                abundances.append(parts[1])
                                abundance_found = True
                                break

                except Exception as e:
                    logger.error(f"Error reading file {file} for genome {genome}: {str(e)}")

                # If genome not found in this file, add 0.0
                if not abundance_found:
                    abundances.append('0.0')

            # Write this genome's row to the output file
            out_file.write(f'{genome}\t{"\t".join(abundances)}\n')

            # Clear memory
            abundances = []

            # Force garbage collection occasionally
            if len(sorted_genomes) > 1000 and sorted_genomes.index(genome) % 1000 == 0:
                gc.collect()

    logger.info(f"Successfully merged {len(files)} files into {output_table}")