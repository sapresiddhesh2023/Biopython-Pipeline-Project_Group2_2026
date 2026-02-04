from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

VALID_DNA = set("ATGCN")
VALID_AA = set("ACDEFGHIKLMNPQRSTVWY*")

input_gb = "Otsutsugamushi_Karp_tsa47.gb"
output_fasta = "filtered_proteins.fasta"

filtered_records = []

for record in SeqIO.parse(input_gb, "genbank"):

    for feature in record.features:

        if feature.type == "CDS":

            # -------- DNA extraction --------
            dna_seq = str(feature.extract(record.seq)).upper()

            # ---- DNA QC ----
            if len(dna_seq) < 90:
                continue

            if not set(dna_seq).issubset(VALID_DNA):
                print(f"Invalid nucleotides in {record.id}")
                continue

            if len(dna_seq) % 3 != 0:
                print(f"Frame issue in {record.id}")
                continue

            if not dna_seq.startswith("ATG"):
                print(f"No start codon in {record.id}")
                continue

            # -------- Translate yourself --------
            protein_seq = str(feature.extract(record.seq).translate(table=11, to_stop=False))

            # ---- Protein QC ----
            if not set(protein_seq).issubset(VALID_AA):
                print(f"Invalid AA in {record.id}")
                continue

            if "*" in protein_seq[:-1]:
                print(f"Internal stop codon in {record.id}")
                continue

            # Save
            protein_record = SeqRecord(
                feature.extract(record.seq).translate(table=11),
                id=record.id,
                description="CDS protein (DNA+Protein QC passed)"
            )

            filtered_records.append(protein_record)

SeqIO.write(filtered_records, output_fasta, "fasta")

print(f"Saved {len(filtered_records)} high-quality proteins.")
