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
import sys
from collections import defaultdict, Counter


def parse_tsv(tsv_file):
    """
    Parse the TSV and build a mapping from PDB filename to the residues to snap.
    Assumes a tab-delimited file with a column named exactly 'BondsJson'.
    Returns:
        dict: { filename: { (chain, resnum): (backbone_coords, resname, bond_details), ... }, ... }
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
                
                # Store additional bond details
                bond_details = {
                    'structureFile': pdb_file,
                    'dssp1': entry.get('dssp1'),
                    'resType1': entry.get('resType1'),
                    'interaction': entry.get('interaction'),
                    'res1': res1
                }
                
                mapping.setdefault(pdb_file, {})[(chain, resnum)] = (backbone_coords, resname, bond_details)
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


def calculate_backbone_normal(backbone_coords):
    """
    Calculate the normal vector of the backbone plane (N-CA-C).
    Args:
        backbone_coords: List of 4 coordinate arrays [N, CA, C, O]
    Returns:
        numpy.ndarray: Normal vector of the backbone plane
    """
    if not backbone_coords or len(backbone_coords) < 3:
        return None
    
    N = np.array(backbone_coords[0])
    CA = np.array(backbone_coords[1])
    C = np.array(backbone_coords[2])
    
    # Calculate two vectors in the plane
    v1 = N - CA
    v2 = C - CA
    
    # Calculate normal vector using cross product
    normal = np.cross(v1, v2)
    
    # Normalize the vector
    norm = np.linalg.norm(normal)
    if norm == 0:
        return None
    
    return normal / norm


def calculate_backbone_rotation_diff(backbone1, backbone2):
    """
    Calculate the angle between normal vectors of two backbone planes.
    Args:
        backbone1: List of 4 coordinate arrays [N, CA, C, O]
        backbone2: List of 4 coordinate arrays [N, CA, C, O]
    Returns:
        float: Angle between normal vectors in degrees
    """
    if not backbone1 or not backbone2:
        return float('inf')
    
    normal1 = calculate_backbone_normal(backbone1)
    normal2 = calculate_backbone_normal(backbone2)
    
    if normal1 is None or normal2 is None:
        return float('inf')
    
    # Calculate angle between normal vectors
    cos_angle = np.dot(normal1, normal2)
    cos_angle = np.clip(cos_angle, -1.0, 1.0)  # Ensure value is in valid range
    angle = np.degrees(np.arccos(cos_angle))
    
    return angle


def parse_backbone_coords(pdb_file):
    """
    Parse backbone coordinates (N, CA, C, O) from a PDB file, only for chain A.
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
                
            chain = line[21]
            if chain != 'A':  # Skip if not chain A
                continue
                
            atom_name = line[12:16].strip()
            if atom_name not in ['N', 'CA', 'C', 'O']:
                continue
                
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


def snap_pdb_file(input_path, output_path, bouquet_mapping, dist_threshold=2.0, angle_threshold=6.0, min_snaps=2):
    """
    Read a PDB file, snap specified residues to their nearest backbone positions,
    and append REMARK lines. Only updates the residue name, preserving chain and residue number.

    Args:
        input_path (str): Path to the original PDB.
        output_path (str): Path where the modified PDB will be saved.
        bouquet_mapping (dict): { (chain, resnum): (backbone_coords, resname, bond_details), ... }
        dist_threshold (float): Maximum distance threshold for snapping (default: 2.0)
        angle_threshold (float): Maximum angle difference threshold for snapping in degrees (default: 6.0)
        min_snaps (int): Minimum number of snaps required to write output file
    """
    # Get the backbone coordinates from the input PDB
    backbone_coords = parse_backbone_coords(input_path)
    
    # Find the best matches for each backbone position
    snaps = {}  # {(chain, resnum): (resname, bond_details)}
    snap_details = []  # List to store details for printing
    for (chain, resnum), target_backbone in backbone_coords.items():
        min_dist = dist_threshold
        min_angle_diff = angle_threshold
        best_match = None
        best_details = None
        
        for (_, _), (bouquet_backbone, resname, bond_details) in bouquet_mapping.items():
            # First check distance
            dist = calculate_backbone_distance(bouquet_backbone, target_backbone)
            if dist < min_dist:
                # Only check angle if distance passes
                angle_diff = calculate_backbone_rotation_diff(bouquet_backbone, target_backbone)
                if angle_diff < min_angle_diff:
                    min_dist = dist
                    min_angle_diff = angle_diff
                    best_match = resname
                    best_details = bond_details
        
        if best_match:
            snaps[(chain, resnum)] = (best_match, best_details)
            details = best_details
            action = "Fix" if min_angle_diff < 0.01 else "Snap"
            snap_details.append(
                f"  {action} {best_match} -> {chain}:{resnum} with distance {min_dist:.2f} and angle difference {min_angle_diff:.2f}°\n"
                f"    Structure: {details['structureFile']}\n"
                f"    DSSP: {details['dssp1']}\n"
                f"    ResType: {details['resType1']}\n"
                f"    Interaction: {details['interaction']}\n"
                f"    Res1: {details['res1']}"
            )
    
    # Only proceed if we have enough snaps
    if len(snaps) < min_snaps:
        return len(snaps), snap_details

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
                # Only update the residue name
                new_name = snaps[(chain, resnum)][0]
                line = line[:17] + new_name.rjust(3) + line[20:]
        new_lines.append(line)

    new_lines.append("\n")
    # Append SNAPPED remarks at end
    for (chain, resnum), (resname, _) in sorted(snaps.items(), key=lambda x: x[0][1]):
        remark = f"REMARK PDBinfo-LABEL:{resnum:>5} FIXED\n"
        new_lines.append(remark)

    with open(output_path, 'w') as f:
        f.writelines(new_lines)
    
    return len(snaps), snap_details


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
    parser.add_argument(
        '--angle-threshold', type=float, default=6.0,
        help='Maximum angle difference threshold for snapping in degrees (default: 6.0)'
    )
    parser.add_argument(
        '--min-snaps', type=int, default=2,
        help='Minimum number of snaps required to process a PDB file (default: 2)'
    )
    args = parser.parse_args()

    # Log the full command line
    print(f"\nCommand line: snap_conserved_residues.py {' '.join(sys.argv)}")
    print(f"Working directory: {os.getcwd()}")
    print("-" * 80)

    # Parse the TSV file to get bouquet mappings
    bouquet_mapping = parse_tsv(args.tsv_file)
    os.makedirs(args.new_folder, exist_ok=True)

    # Pool all residue mappings from all files
    pooled_fixes = {}
    for file_mapping in bouquet_mapping.values():
        pooled_fixes.update(file_mapping)

    # Collect statistics about number of snaps
    snap_counts = Counter()

    # Process all PDB files in the old folder
    for filename in os.listdir(args.old_folder):
        if not filename.endswith('.pdb'):
            continue
            
        input_path = os.path.join(args.old_folder, filename)
        output_path = os.path.join(args.new_folder, filename)
        
        # Use the pooled fixes for every file
        num_snaps, snap_details = snap_pdb_file(input_path, output_path, pooled_fixes, args.dist_threshold, args.angle_threshold, args.min_snaps)
        snap_counts[num_snaps] += 1
        
        if num_snaps >= args.min_snaps:
            print(f"Processed {filename}: snapped {num_snaps} residue(s):")
            for detail in snap_details:
                print(detail)
        else:
            print(f"Skipped {filename}: only {num_snaps} snap(s), below minimum threshold of {args.min_snaps}")

    # Print summary table
    print("\nSummary of snap counts:")
    print("Number of snaps | Total PDBs | Processed PDBs")
    print("-" * 50)
    for num_snaps, count in sorted(snap_counts.items()):
        processed = count if num_snaps >= args.min_snaps else 0
        print(f"{num_snaps:^14} | {count:^10} | {processed:^13}")


if __name__ == '__main__':
    main()
