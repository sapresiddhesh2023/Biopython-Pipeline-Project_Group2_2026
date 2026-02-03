from Bio import SeqIO 
from Bio.Seq import Seq

VALID_DNA = set("ATGC")
input_fasta = "Otsutsugamushi_Karp_tsa47.fasta" 
output_fasta = "filtered_sequence.fasta"

filtered_records = []
for record in SeqIO.parse(input_fasta, "fasta"): 
    seq = str(record.seq).upper()
    # 1. Length filter 
    if len(seq) < 10: 
        continue

    # 2. Alphabet validation 
    if not set(seq).issubset(VALID_DNA): 
        print(f"Invalid characters in {record.id}") 
        continue

    # 3. Must start with ATG 
    if not seq.startswith("ATG"): 
        continue 
    filtered_records.append(record)

    # Write filtered sequences 
    SeqIO.write(filtered_records, output_fasta, "fasta") 
    print(f"Saved {len(filtered_records)} valid sequences.")

