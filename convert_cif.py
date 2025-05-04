#!/home/luod/miniconda3/envs/umap/bin/python
"""
convert_cif.py: Use Biopython to read a mmCIF file (even if missing label_asym_id or occupancy) and export a valid mmCIF format.

Usage:
    python convert_cif.py input.cif output.cif

This script will:
 1. Preprocess the CIF text to add missing `label_asym_id` (from auth_asym_id) and `occupancy` columns with defaults.
 2. Parse with MMCIFParser.
 3. Ensure every atom has `label_asym_id` in Atom.xtra.
 4. Write out a proper mmCIF using MMCIFIO.
"""
import sys
import argparse
import tempfile
import os
from Bio.PDB import MMCIFParser, MMCIFIO


def patch_cif(infile_path):
    lines = open(infile_path).read().splitlines(keepends=True)
    out_lines = []
    in_loop = False
    header_lines = []
    data_started = False
    idx_label = None
    idx_auth = None
    idx_occ = None
    insert_label = False
    insert_occ = False

    for line in lines:
        if line.strip() == 'loop_':
            in_loop = True
            header_lines = []
            out_lines.append(line)
            continue

        # Collect atom_site headers
        if in_loop and line.startswith('_atom_site.'):
            header_lines.append(line.rstrip('\n'))
            continue

        # End of headers
        if in_loop and not line.startswith('_atom_site.'):
            in_loop = False
            # Determine insertion indices
            idx_comp = next((i for i,h in enumerate(header_lines)
                             if '.label_comp_id' in h), len(header_lines)-1)
            idx_biso = next((i for i,h in enumerate(header_lines)
                             if '.B_iso_or_equiv' in h), None)
            # Check/preset occupancy
            if not any('.occupancy' in h for h in header_lines):
                insert_occ = True
                idx_occ = (idx_biso + 1) if idx_biso is not None else len(header_lines)
                header_lines.insert(idx_occ, '_atom_site.occupancy')
            # Check/preset label_asym_id
            if not any('.label_asym_id' in h for h in header_lines):
                insert_label = True
                idx_label = idx_comp + 1
                header_lines.insert(idx_label, '_atom_site.label_asym_id')
            # Find auth_asym_id
            idx_auth = next((i for i,h in enumerate(header_lines)
                             if '.auth_asym_id' in h), None)
            # Write headers
            for h in header_lines:
                out_lines.append(h + '\n')
            data_started = True
            # Process this line as first data row
            if line.strip() and not line.startswith('#'):
                out_lines.append(patch_data_line(line, header_lines, idx_auth, idx_label, idx_occ, insert_label, insert_occ))
            else:
                out_lines.append(line)
            continue

        # Data rows
        if data_started:
            if line.strip() == '' or line.startswith('#') or line.startswith('_'):
                out_lines.append(line)
            else:
                out_lines.append(patch_data_line(line, header_lines, idx_auth, idx_label, idx_occ, insert_label, insert_occ))
        else:
            out_lines.append(line)

    # Write to temp file in text mode
    tmp = tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.cif')
    tmp.writelines(out_lines)
    tmp_path = tmp.name
    tmp.close()
    return tmp_path


def patch_data_line(line, headers, idx_auth, idx_label, idx_occ, add_label, add_occ):
    parts = line.strip().split()
    # If already matches headers, return unchanged
    if len(parts) == len(headers):
        return line
    # insert occupancy default if needed
    if add_occ and idx_occ is not None:
        parts.insert(idx_occ, '1.00')
    # insert label_asym_id from auth_asym_id or blank
    if add_label and idx_label is not None:
        auth_val = parts[idx_auth] if idx_auth is not None and idx_auth < len(parts) else ''
        parts.insert(idx_label, auth_val)
    return ' '.join(parts) + '\n'


def main():
    parser = argparse.ArgumentParser(description='Convert CIF to valid mmCIF via Biopython')
    parser.add_argument('infile', help='Input CIF file')
    parser.add_argument('outfile', help='Output mmCIF file')
    args = parser.parse_args()

    # Preprocess and patch missing columns
    patched_path = patch_cif(args.infile)

    # Parse structure from patched file
    cif_parser = MMCIFParser(QUIET=True)
    structure = cif_parser.get_structure('converted', patched_path)

    # Ensure label_asym_id in Atom.xtra for Biopython IO
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if 'label_asym_id' not in atom.xtra:
                        atom.xtra['label_asym_id'] = atom.xtra.get('auth_asym_id', chain.id)

    # Write out mmCIF
    io = MMCIFIO()
    io.set_structure(structure)
    io.save(args.outfile)
    print(f'Wrote fixed mmCIF to {args.outfile}')

    # Clean up
    os.remove(patched_path)

if __name__ == '__main__':
    main()

