#!/usr/bin/env python3

import sys
import csv

def parse_score_file(input_path, output_path):
    """
    Reads a file of lines starting with 'SCORE:'
    and converts them to a CSV file with the columns:
      binder_aligned_rmsd, pae_binder, pae_interaction, pae_target,
      plddt_binder, plddt_target, plddt_total, target_aligned_rmsd,
      time, description
    """
    # Storage for column names (header) and rows
    header = []
    rows = []

    with open(input_path, 'r') as infile:
        for line in infile:
            # Only process lines that start with 'SCORE:'
            if line.strip().startswith("SCORE:"):
                # Remove the leading 'SCORE:' and extra whitespace
                # Then split by any whitespace
                parts = line.replace("SCORE:", "").split()

                # If this line is the header (check if it has non-numeric fields like 'binder_aligned_rmsd')
                # we store it as the header
                if all(not p.replace('.', '', 1).isdigit() for p in parts):
                    header = parts
                else:
                    # Otherwise, it's a data row
                    rows.append(parts)

    # If no header is found, define a default (ordered) header
    # that matches the columns in your scoreboard example
    if not header:
        header = [
            "binder_aligned_rmsd",
            "pae_binder",
            "pae_interaction",
            "pae_target",
            "plddt_binder",
            "plddt_target",
            "plddt_total",
            "target_aligned_rmsd",
            "time",
            "description"
        ]

    # Write out to CSV
    with open(output_path, 'w', newline='') as outfile:
        writer = csv.writer(outfile)
        writer.writerow(header)
        for row in rows:
            writer.writerow(row)

def main():
    if len(sys.argv) < 3:
        print("Usage: python convert_scores_to_csv.py <input_file.txt> <output_file.csv>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    parse_score_file(input_file, output_file)
    print(f"CSV saved to {output_file}")

if __name__ == "__main__":
    main()

