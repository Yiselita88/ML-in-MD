from Bio import PDB
from Bio.SeqUtils import seq1


def extract_fasta_from_pdb(peptide):
    # Create a PDB parser
    parser = PDB.PDBParser(QUIET=True)

    # Load the PDB file
    structure = parser.get_structure('BRD3', pdb_filename)

    # Initialize an empty sequence
    sequence = ""

    # Iterate through the structure and extract sequence
    for model in structure:
        for chain in model:
            for residue in chain:
                if PDB.is_aa(residue):
                    # Convert three-letter code to one-letter code
                    sequence += seq1(residue.get_resname())

    return sequence


def save_fasta_to_file(fasta_sequence, output_filename):
    # Save the FASTA sequence to a file in the required format
    with open(output_filename, 'w') as output_file:
        output_file.write(">seq\n")
        output_file.write(fasta_sequence)


if __name__ == "__main__":
    # Replace 'your_peptide.pdb' with the actual PDB file name
    pdb_filename = 'BRD3.pdb'

    try:
        fasta_sequence = extract_fasta_from_pdb(pdb_filename)

        # Replace 'output_sequence.fasta' with the desired output file name
        output_filename = 'BRD3.fasta'

        save_fasta_to_file(fasta_sequence, output_filename)
        print(f"Fasta Sequence saved to {output_filename}")
    except Exception as e:
        print(f"Error: {e}")
