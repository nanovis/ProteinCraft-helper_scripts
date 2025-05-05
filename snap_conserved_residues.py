#!/usr/bin/env python3
"""
Script to snap conserved residues to backbone positions based on backbone coordinates.

Usage:
    python snap_conserved_residues.py outs_tetris.csv old_RFdiffusion_backbone_folder new_RFdiffusion_backbone_folder
"""
import os
import csv
import json
import argparse
import numpy as np
from collections import defaultdict


def parse_tsv(tsv_file):
    """
    Parse the TSV and build a mapping from PDB filename to the residues to snap.
    Assumes a tab-delimited file with a column named exactly 'BondsJson'.
    Returns:
        dict: { filename: { (chain, resnum): (backbone_coords, resname), ... }, ... }
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
                backbone_coords = entry.get('res1_backbone_coords')
                if not pdb_file or not res1 or not backbone_coords:
                    continue
                parts = res1.split(':')
                if len(parts) != 4:
                    print(f"Warning: unexpected res1 format '{res1}' in {pdb_file}")
                    continue
                chain, resnum_str, _, resname = parts
                try:
                    resnum = int(resnum_str)
                except ValueError:
                    print(f"Warning: invalid residue number '{resnum_str}' in {res1}")
                    continue
                mapping.setdefault(pdb_file, {})[(chain, resnum)] = (backbone_coords, resname)
    return mapping


def calculate_backbone_distance(backbone1, backbone2):
    """
    Calculate the sum of distances between corresponding backbone atoms.
    Args:
        backbone1: List of 4 coordinate arrays [N, CA, C, O]
        backbone2: List of 4 coordinate arrays [N, CA, C, O]
    Returns:
        float: Sum of distances between corresponding atoms
    """
    if not backbone1 or not backbone2:
        return float('inf')
    
    total_dist = 0
    for coords1, coords2 in zip(backbone1, backbone2):
        dist = np.sqrt(np.sum((np.array(coords1) - np.array(coords2)) ** 2))
        total_dist += dist
    return total_dist


def parse_backbone_coords(pdb_file):
    """
    Parse backbone coordinates (N, CA, C, O) from a PDB file.
    Returns a dict mapping (chain, resnum) -> [N_coords, CA_coords, C_coords, O_coords]
    where each coords is [x, y, z]
    """
    backbone_coords = {}
    current_res = None
    current_coords = {}
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if not line.startswith('ATOM'):
                continue
                
            atom_name = line[12:16].strip()
            if atom_name not in ['N', 'CA', 'C', 'O']:
                continue
                
            chain = line[21]
            resnum = int(line[22:26])
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            
            res_key = (chain, resnum)
            if res_key != current_res:
                if current_res is not None and len(current_coords) == 4:
                    backbone_coords[current_res] = [
                        current_coords['N'],
                        current_coords['CA'],
                        current_coords['C'],
                        current_coords['O']
                    ]
                current_res = res_key
                current_coords = {}
            
            current_coords[atom_name] = [x, y, z]
    
    # Don't forget the last residue
    if current_res is not None and len(current_coords) == 4:
        backbone_coords[current_res] = [
            current_coords['N'],
            current_coords['CA'],
            current_coords['C'],
            current_coords['O']
        ]
    
    return backbone_coords


def snap_pdb_file(input_path, output_path, bouquet_mapping, dist_threshold=2.0):
    """
    Read a PDB file, snap specified residues to their nearest backbone positions,
    and append REMARK lines.

    Args:
        input_path (str): Path to the original PDB.
        output_path (str): Path where the modified PDB will be saved.
        bouquet_mapping (dict): { filename: { (chain, resnum): (backbone_coords, resname), ... }, ... }
        dist_threshold (float): Maximum distance threshold for snapping
    """
    # Get the backbone coordinates from the input PDB
    backbone_coords = parse_backbone_coords(input_path)
    
    # Find the best matches for each bouquet residue
    snaps = {}  # {(chain, resnum): (target_chain, target_resnum, resname)}
    for (chain, resnum), (bouquet_backbone, resname) in bouquet_mapping.items():
        min_dist = dist_threshold
        best_match = None
        
        for (target_chain, target_resnum), target_backbone in backbone_coords.items():
            dist = calculate_backbone_distance(bouquet_backbone, target_backbone)
            if dist < min_dist:
                min_dist = dist
                best_match = (target_chain, target_resnum)
        
        if best_match:
            snaps[(chain, resnum)] = (*best_match, resname)
    
    # Read and modify the PDB file
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
            if (chain, resnum) in snaps:
                target_chain, target_resnum, new_name = snaps[(chain, resnum)]
                # Update chain and residue number
                line = line[:21] + target_chain + f"{target_resnum:>4}" + line[26:]
                # Update residue name
                line = line[:17] + new_name.rjust(3) + line[20:]
        new_lines.append(line)

    new_lines.append("\n")
    # Append SNAPPED remarks at end
    for (chain, resnum), (target_chain, target_resnum, resname) in sorted(snaps.items(), key=lambda x: x[0][1]):
        remark = f"REMARK PDBinfo-LABEL:{resnum:>5} SNAPPED to {target_chain}:{target_resnum} as {resname}\n"
        new_lines.append(remark)

    with open(output_path, 'w') as f:
        f.writelines(new_lines)
    
    return len(snaps)


def main():
    parser = argparse.ArgumentParser(
        description="Snap conserved residues to backbone positions based on backbone coordinates"
    )
    parser.add_argument(
        'tsv_file', help='Path to the tab-delimited file with BondsJson column'
    )
    parser.add_argument(
        'old_folder', help='Directory containing the original backbone PDB files'
    )
    parser.add_argument(
        'new_folder', help='Directory to save modified PDB files'
    )
    parser.add_argument(
        '--dist-threshold', type=float, default=2.0,
        help='Maximum distance threshold for snapping (default: 2.0)'
    )
    args = parser.parse_args()

    # Parse the TSV file to get bouquet mappings
    bouquet_mapping = parse_tsv(args.tsv_file)
    os.makedirs(args.new_folder, exist_ok=True)

    # Process each PDB file
    for pdb_file, fixes in bouquet_mapping.items():
        # try exact filename first
        input_path = os.path.join(args.old_folder, pdb_file)
        output_name = pdb_file
        if not os.path.exists(input_path):
            # fallback: try different filename patterns
            # strip everything after first '_dldesign' and add .pdb
            prefix = pdb_file.split('_dldesign')[0]
            alt_name = prefix + '.pdb'
            alt_path = os.path.join(args.old_folder, alt_name)
            
            if not os.path.exists(alt_path):
                print(f"Warning: '{pdb_file}' not found and no fallback '{alt_name}', skipping.")
                continue
            else:
                input_path = alt_path
                output_name = alt_name
        
        output_path = os.path.join(args.new_folder, output_name)
        num_snaps = snap_pdb_file(input_path, output_path, fixes, args.dist_threshold)
        print(f"Processed {output_name}: snapped {num_snaps} residue(s).")


if __name__ == '__main__':
    main()
