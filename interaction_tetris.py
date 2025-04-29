import os
import glob
from collections import defaultdict
import argparse
import json

# Map 3-letter to 1-letter codes
AA_3TO1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
    'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G',
    'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
    'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

# Choose a fixed order for one-letter amino acids
AA_ORDER = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']

def parse_node_id(node_str):
    """
    Given a NodeId like 'A:50:_:LYS',
      returns (chain_id='A', residue_number='50', residue_name='LYS').
    """
    parts = node_str.split(':')
    chain_id      = parts[0]
    residue_num   = parts[1]
    residue_3name = parts[3]  # e.g. 'LYS'
    return chain_id, residue_num, residue_3name

# This dict maps (b_res_num, b_res_name_3) -> { a_res_name_3: count_of_interactions }
interactions = defaultdict(lambda: defaultdict(int))
# This dict maps (b_res_num, b_res_name_3) -> { a_res_name_3: [(filename, resnum), ...] }
interaction_details = defaultdict(lambda: defaultdict(list))

def main():
    parser = argparse.ArgumentParser(description='Analyze protein interactions from RING output files')
    parser.add_argument('folder_path', help='Path to folder containing RING output files')
    parser.add_argument('--top', type=int, help='Number of top interactions to output (default: all)')
    parser.add_argument('--debug', action='store_true', help='Enable debug output')
    args = parser.parse_args()

    def debug_print(message):
        if args.debug:
            print(message)

    debug_print(f"Processing folder: {args.folder_path}")  # Debug print

    # Find all ring edges files in the specified folder
    ringedges_files = glob.glob(os.path.join(args.folder_path, '*_ringEdges'))
    
    debug_print(f"Found {len(ringedges_files)} ring edges files")  # Debug print
    if not ringedges_files:
        print(f"Warning: No ring edges files found in {args.folder_path}")
        return

    for ringedges_file in ringedges_files:
        debug_print(f"\nProcessing file: {ringedges_file}")  # Debug print
        # Get filename without _ringEdges for the details column
        base_filename = os.path.basename(ringedges_file).replace('_ringEdges', '')
        
        ringnodes_file = ringedges_file.replace('_ringEdges', '_ringNodes')
        if not os.path.exists(ringnodes_file):
            print(f"Warning: No corresponding ringNodes file found for {ringedges_file}, skipping...")
            continue
        
        # Parse the ringNodes file to get DSSP information
        residue_dssp = {}  # Maps (chain, resnum) -> dssp
        with open(ringnodes_file, 'r') as f:
            lines = f.readlines()
            debug_print(f"Found {len(lines)} lines in ringNodes file")  # Debug print
            # Skip header
            for line in lines[1:]:
                fields = line.strip().split('\t')
                if len(fields) < 6:  # Need at least 6 fields to get DSSP
                    continue
                node_id = fields[0]  # e.g. 'A:1:_:GLU'
                chain = fields[1]
                position = fields[2]
                dssp = fields[5]
                
                # Store DSSP for this residue
                try:
                    resnum_int = int(position)
                    residue_dssp[(chain, resnum_int)] = dssp
                except ValueError:
                    continue
        
        debug_print(f"Processed {len(residue_dssp)} residues from ringNodes file")  # Debug print
        
        # Parse the ringEdges file
        with open(ringedges_file, 'r') as f:
            lines = f.readlines()
            debug_print(f"Found {len(lines)} lines in ringEdges file")  # Debug print
            # Skip header
            for line in lines[1:]:
                fields = line.strip().split('\t')
                if len(fields) < 3:
                    continue
                node1       = fields[0]
                interaction = fields[1]
                node2       = fields[2]
                atom1       = fields[5]
                atom2       = fields[6]
                
                # Skip VDW
                if 'VDW' in interaction:
                    continue
                
                chain1, resnum1, res3name1 = parse_node_id(node1)
                chain2, resnum2, res3name2 = parse_node_id(node2)
                
                # Convert to integers if possible
                try:
                    resnum1_int = int(resnum1)
                    resnum2_int = int(resnum2)
                except ValueError:
                    continue
                
                # Get DSSP information
                dssp1 = residue_dssp.get((chain1, resnum1_int), "")
                dssp2 = residue_dssp.get((chain2, resnum2_int), "")
                
                # We only care about chain A <-> chain B
                if chain1 == 'A' and chain2 == 'B':
                    interactions[(resnum2_int, res3name2)][res3name1] += 1
                    interaction_details[(resnum2_int, res3name2)][res3name1].append((
                        base_filename, 
                        dssp1, 
                        interaction,
                        node1,
                        node2,
                        atom1, 
                        atom2))
        
        debug_print(f"Found {len(interactions)} unique B-chain residues with interactions")  # Debug print

    # Summarize the total interactions and pick the top N
    b_info_list = []
    for (b_res_num, b_res_3name), a_counts_3name in interactions.items():
        total = sum(a_counts_3name.values())
        b_info_list.append(((b_res_num, b_res_3name), total))

    debug_print(f"\nTotal unique B-chain residues with interactions: {len(b_info_list)}")  # Debug print

    # Sort by total interactions (descending)
    b_info_list.sort(key=lambda x: x[1], reverse=True)
    
    # If top N is specified, take only those entries
    if args.top is not None:
        b_info_list = b_info_list[:args.top]
    
    # Sort by B residue number (ascending)
    b_info_list.sort(key=lambda x: x[0][0])

    # Print a header row for easier parsing (e.g. to feed into a logo tool)
    # The columns are:
    #  B_Residue  Total   (then freq of each A amino acid in AA_ORDER)  BondsJson
    header_cols = ["B_residue", "Total"] + AA_ORDER + ["BondsJson"]
    print("\t".join(header_cols))

    for ((b_res_num, b_res_3name), total_int) in b_info_list:
        # Build a dict for counts by single-letter code
        a_counts_3name = interactions[(b_res_num, b_res_3name)]
        a_details = interaction_details[(b_res_num, b_res_3name)]
        
        # Convert 3-letter AAs to 1-letter, summing up
        single_letter_counts = defaultdict(int)
        single_letter_details = defaultdict(list)
        for (aa3, count) in a_counts_3name.items():
            aa1 = AA_3TO1.get(aa3, 'X')  # 'X' if unknown
            single_letter_counts[aa1] += count
            single_letter_details[aa1].extend(a_details[aa3])
        
        # Build the row:
        # "B:res_num:1letter:dssp" format
        b_res_1letter = AA_3TO1.get(b_res_3name, 'X')
        b_dssp = residue_dssp.get(("B", b_res_num), "")
        b_label = f"B:{b_res_num}:{b_res_1letter}:{b_dssp}"
        row = [b_label, str(total_int)]
        
        # Frequencies for each standard AA in AA_ORDER
        for aa1 in AA_ORDER:
            c = single_letter_counts[aa1]
            freq = c / total_int if total_int > 0 else 0
            row.append(f"{freq:.3f}")
        
        # First collect all bond JSON objects
        bond_jsons = []
        for aa1 in AA_ORDER:
            if aa1 in single_letter_details:
                # Sort by residue number
                sorted_details = sorted(single_letter_details[aa1], key=lambda x: x[0])  # Sort by filename
                for filename, dssp1, interaction, node1, node2, atom1, atom2 in sorted_details:
                    # Parse residue numbers from node strings
                    chain1, resnum1, res3name1 = parse_node_id(node1)

                    single_bond_json = {
                        "structureFile": filename,
                        "dssp1": dssp1,
                        "resNum1": int(resnum1),
                        "resType1": AA_3TO1.get(res3name1, 'X'),
                        "interaction": interaction,
                        "res1": node1,
                        "res2": node2,
                        "atom1": atom1,
                        "atom2": atom2
                    }
                    bond_jsons.append(single_bond_json)
        
        # Then convert the entire list to a JSON array string
        row.append(json.dumps(bond_jsons))
        
        print("\t".join(row))

if __name__ == '__main__':
    main()

