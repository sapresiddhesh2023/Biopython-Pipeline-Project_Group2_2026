from Bio import SeqIO 

# Standard 20 amino acids + optional stop (*)
VALID_AA = set("ACDEFGHIKLMNPQRSTVWY*")

input_fasta = "protein_sequence.fasta"
output_fasta = "filtered_protein.fasta"

filtered_records = []

for record in SeqIO.parse(input_fasta, "fasta"):
    seq = str(record.seq).upper()

    # 1. Length filter
    if len(seq) < 30:
        continue

    # 2. Amino acid alphabet validation
    if not set(seq).issubset(VALID_AA):
        print(f"Invalid characters found in {record.id}")
        continue

    filtered_records.append(record)

# Write AFTER loop
SeqIO.write(filtered_records, output_fasta, "fasta")

print(f"Saved {len(filtered_records)} valid protein sequences.")
