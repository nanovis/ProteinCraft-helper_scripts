#!/usr/bin/env python
import os
import sys

def reindex_chainB_in_pdb(in_pdb_path, out_pdb_path):
    """
    Read a PDB file from in_pdb_path, shift the residue numbering of chain B
    so that it starts at 1, and write the new file to out_pdb_path.
    """
    chainB_offset = None  # Will store (firstBres - 1) once first B residue is encountered
    
    with open(in_pdb_path, 'r') as f_in, open(out_pdb_path, 'w') as f_out:
        for line in f_in:
            record = line[0:6].strip()
            
            # We only modify lines that contain residue numbering: ATOM, HETATM, ANISOU, TER, etc.
            # (PDB columns are 1-based; in Python, slicing is 0-based)
            if record in ("ATOM", "HETATM", "ANISOU", "TER"):
                # Chain ID is in column 22 (index 21 in Python 0-based)
                chain_id = line[21]
                
                # Residue number is in columns 23-26 (indexes 22:26 in Python)
                try:
                    old_resnum_str = line[22:26]
                    old_resnum = int(old_resnum_str)
                except ValueError:
                    # If it doesn't parse, just write line as-is
                    f_out.write(line)
                    continue
                
                # If this is chain B, shift the residue number
                if chain_id == 'B':
                    if chainB_offset is None:
                        # Initialize offset so that the current old_resnum becomes 1
                        chainB_offset = old_resnum - 1
                    new_resnum = old_resnum - chainB_offset
                    
                    # Rebuild the line with the new residue number, preserving PDB format
                    # Residue number goes from columns 23-26, so we do line[:22] + "%4d" + line[26:]
                    new_resnum_str = f"{new_resnum:4d}"
                    new_line = line[:22] + new_resnum_str + line[26:]
                    f_out.write(new_line)
                else:
                    # Keep everything else (other chains, e.g., A) unchanged
                    f_out.write(line)
            else:
                # Non-ATOM/HETATM lines, just copy them
                f_out.write(line)


def main():
    if len(sys.argv) != 3:
        print("Usage: python reindex_chainB.py <input_folder> <output_folder>")
        sys.exit(1)
        
    input_folder = sys.argv[1]
    output_folder = sys.argv[2]

    # Create output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)
    
    # Process each PDB file in the input folder
    for filename in os.listdir(input_folder):
        if filename.lower().endswith(".pdb"):
            in_path = os.path.join(input_folder, filename)
            out_path = os.path.join(output_folder, filename)
            
            reindex_chainB_in_pdb(in_path, out_path)

    print("Reindexing complete. Output written to:", output_folder)


if __name__ == "__main__":
    main()

