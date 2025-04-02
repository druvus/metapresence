"""
Utilities for merging abundance outputs.
"""

import os


def merge_abundances(abundance_folder, output_table):
    """
    Merge multiple abundance outputs into a single table.

    Args:
        abundance_folder: Folder containing abundance outputs.
        output_table: Path to the output table file.
    """
    if not os.path.isdir(abundance_folder):
        print(f"Error: {abundance_folder} is not a valid directory")
        return

    files = [x for x in sorted(os.listdir(abundance_folder)) if os.path.isfile(os.path.join(abundance_folder, x))]

    if not files:
        print(f"Error: No files found in {abundance_folder}")
        return

    printing_dict = {}

    for n, file in enumerate(files):
        file_path = os.path.join(abundance_folder, file)
        try:
            with open(file_path) as f:
                # Skip header
                f.readline()

                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) < 2:
                        continue

                    genome = parts[0]
                    abundance = parts[1]

                    if genome in printing_dict:
                        printing_dict[genome].append(abundance)
                    else:
                        printing_dict[genome] = ['0.0'] * n
                        printing_dict[genome].append(abundance)

            # Ensure all genomes have the same number of values
            for genome in printing_dict:
                if len(printing_dict[genome]) < n + 1:
                    printing_dict[genome].append('0.0')
        except Exception as e:
            print(f"Error processing file {file}: {str(e)}")

    if not printing_dict:
        print("No data found in files. The merged table will be empty.")

    # Sort genomes alphabetically
    printing_dict = {genome: printing_dict[genome] for genome in sorted(printing_dict.keys())}

    with open(output_table, 'w') as new:
        new.write(f'fasta\t{"\t".join(files)}\n')

        for genome in printing_dict:
            new.write(f'{genome}\t{"\t".join(printing_dict[genome])}\n')

    print(f"Successfully merged {len(files)} files into {output_table}")