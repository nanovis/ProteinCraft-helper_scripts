#!/usr/bin/env python3
"""
Script to fix residue types in PDB files based on an outs_tetris.tsv (or similar tab-delimited) mapping.

Usage:
    python fix_residues.py outs_tetris.tsv /path/to/pdb_folder /path/to/output_folder
"""
import os
import csv
import json
import argparse


def parse_tsv(tsv_file):
    """
    Parse the TSV and build a mapping from PDB filename to the residues to fix.
    Assumes a tab-delimited file with a column named exactly 'BondsJson'.
    Returns:
        dict: { filename: { (chain, resnum): new_resname, ... }, ... }
    """
    mapping = {}
    with open(tsv_file, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        # require exact BondsJson column
        if 'BondsJson' not in reader.fieldnames:
            raise KeyError(f"'BondsJson' column not found in {tsv_file}. Available columns: {reader.fieldnames}")

        for row in reader:
            cell = row['BondsJson']
            try:
                bonds = json.loads(cell)
            except json.JSONDecodeError as e:
                print(f"Warning: could not parse JSON in 'BondsJson' for row {row}: {e}")
                continue
            for entry in bonds:
                pdb_file = entry.get('structureFile')
                res1 = entry.get('res1')
                if not pdb_file or not res1:
                    continue
                parts = res1.split(':')
                if len(parts) != 4:
                    print(f"Warning: unexpected res1 format '{res1}' in {pdb_file}")
                    continue
                chain, resnum_str, _, new_resname = parts
                try:
                    resnum = int(resnum_str)
                except ValueError:
                    print(f"Warning: invalid residue number '{resnum_str}' in {res1}")
                    continue
                mapping.setdefault(pdb_file, {})[(chain, resnum)] = new_resname
    return mapping


def fix_pdb_file(input_path, output_path, fixes):
    """
    Read a PDB file, change specified residues, and append REMARK lines.

    Args:
        input_path (str): Path to the original PDB.
        output_path (str): Path where the modified PDB will be saved.
        fixes (dict): {(chain, resnum): new_resname, ...}
    """
    with open(input_path, 'r') as f:
        lines = f.readlines()

    new_lines = []
    for line in lines:
        if line.startswith(('ATOM  ', 'HETATM')):
            chain = line[21]
            try:
                resnum = int(line[22:26].strip())
            except ValueError:
                resnum = None
            if (chain, resnum) in fixes:
                new_name = fixes[(chain, resnum)].rjust(3)
                line = line[:17] + new_name + line[20:]
        new_lines.append(line)

    new_lines.append("\n")
    # Append FIXED remarks at end
    for (chain, resnum), _ in sorted(fixes.items(), key=lambda x: x[0][1]):
        remark = f"REMARK PDBinfo-LABEL:{resnum:>5} FIXED\n"
        new_lines.append(remark)

    with open(output_path, 'w') as f:
        f.writelines(new_lines)


def main():
    parser = argparse.ArgumentParser(
        description="Fix residue types in PDB files based on an outs_tetris.tsv mapping"
    )
    parser.add_argument(
        'tsv_file', help='Path to the tab-delimited file with BondsJson column'
    )
    parser.add_argument(
        'pdb_folder', help='Directory containing the .pdb files'
    )
    parser.add_argument(
        'output_folder', help='Directory to save modified .pdb files'
    )
    parser.add_argument(
        '--min-fixes', type=int, default=2,
        help='Minimum number of fixes required to process a PDB file (default: 2)'
    )
    args = parser.parse_args()

    mapping = parse_tsv(args.tsv_file)
    # Filter mapping to keep only entries with minimum number of fixes
    mapping = {k: v for k, v in mapping.items() if len(v) >= args.min_fixes}
    os.makedirs(args.output_folder, exist_ok=True)

    for pdb_file, fixes in mapping.items():
        # try exact filename first
        input_path = os.path.join(args.pdb_folder, pdb_file)
        output_name = pdb_file
        if not os.path.exists(input_path):
            # fallback: strip everything after first '_dldesign' and add .pdb
            prefix = pdb_file.split('_dldesign')[0]
            alt_name = prefix + '.pdb'
            alt_path = os.path.join(args.pdb_folder, alt_name)
            if os.path.exists(alt_path):
                input_path = alt_path
                output_name = alt_name
            else:
                print(f"Warning: '{pdb_file}' not found and no fallback '{alt_name}', skipping.")
                continue
        output_path = os.path.join(args.output_folder, output_name)
        fix_pdb_file(input_path, output_path, fixes)
        print(f"Processed {output_name}: applied {len(fixes)} residue fix(es).")

if __name__ == '__main__':
    main()
