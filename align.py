#!/usr/bin/env python3
"""
Usage:
    python align.py <input_folder> <output_folder> <reference_file> <reference_chain> <align_chain>

Aligns the Cα atoms of <align_chain> in each PDB under <input_folder> to the Cα atoms
of <reference_chain> in <reference_file>, and writes the aligned PDBs to <output_folder>.

Requires: BioPython (pip install biopython)
"""

import os
import sys
from Bio.PDB import PDBParser, Superimposer, PDBIO

def get_ca_atoms_by_resid(chain):
    """Return a dict mapping residue IDs to Cα Atom objects for that chain."""
    atoms = {}
    for residue in chain:
        if 'CA' in residue:
            atoms[residue.get_id()] = residue['CA']
    return atoms

def align_structure(mov_struct, ref_atoms, mov_atoms):
    """Superimpose mov_atoms onto ref_atoms, apply to all atoms of mov_struct, and return RMSD."""
    sup = Superimposer()
    sup.set_atoms(ref_atoms, mov_atoms)
    sup.apply(mov_struct.get_atoms())
    return sup.rms

def main(input_folder, output_folder, ref_file, ref_chain_id, mov_chain_id):
    # make sure output dir exists
    os.makedirs(output_folder, exist_ok=True)

    parser = PDBParser(QUIET=True)
    # parse reference
    ref_struct = parser.get_structure('REF', ref_file)
    # grab first model that has the chain
    try:
        ref_chain = next(m for m in ref_struct if ref_chain_id in m)[ref_chain_id]
    except StopIteration:
        sys.exit(f"Error: reference chain {ref_chain_id!r} not found in {ref_file}")

    # dict of Cα atoms in reference chain
    ref_ca = get_ca_atoms_by_resid(ref_chain)

    io = PDBIO()

    for fname in sorted(os.listdir(input_folder)):
        if not fname.lower().endswith('.pdb'):
            continue
        in_path = os.path.join(input_folder, fname)
        mov_struct = parser.get_structure('MOV', in_path)

        # find moving chain
        try:
            mov_chain = next(m for m in mov_struct if mov_chain_id in m)[mov_chain_id]
        except StopIteration:
            print(f"[skip] chain {mov_chain_id!r} not in {fname}")
            continue

        mov_ca = get_ca_atoms_by_resid(mov_chain)

        # find common residues (by resid) to match
        common = sorted(set(ref_ca) & set(mov_ca), key=lambda r: r[1])
        if not common:
            print(f"[skip] no matching residues between ref and {fname}")
            continue

        # build parallel atom lists
        ref_atoms = [ref_ca[r] for r in common]
        mov_atoms = [mov_ca[r] for r in common]

        # align & apply
        rmsd = align_structure(mov_struct, ref_atoms, mov_atoms)

        # save output
        out_path = os.path.join(output_folder, fname)
        io.set_structure(mov_struct)
        io.save(out_path)

        print(f"[ok]  {fname} — RMSD: {rmsd:.3f} Å")

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print(__doc__)
        sys.exit(1)
    _, inp, outp, refp, rc, mc = sys.argv
    main(inp, outp, refp, rc, mc)

