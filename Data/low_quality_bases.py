from Bio import SeqIO

valid_aa = set("ACDEFGHIKLMNPQRSTVWY")

record = SeqIO.read("protein_endoprotease.fasta", "fasta")
sequence = str(record.seq)

low_quality_count = sum(1 for aa in sequence if aa not in valid_aa)


print("Protein ID:", record.id)
print("Length:", len(sequence))
print("Low-quality residues:", low_quality_count)

if "*" in sequence[:-1]:
    print("❌ Internal stop codon")




