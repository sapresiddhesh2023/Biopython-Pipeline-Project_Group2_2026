print("Step 1: performed by Siddhesh Uday Sapre")
print("-"*90)
from Bio import SeqIO 
record = SeqIO.read("Otsutsugamushi_Karp_tsa47.fasta", "fasta")
print("The following DNA record has Nucleotide ID:", record.id)
print(record.description)
print("The length of the sequence is:", len(record.seq), "bp")
print("The DNA sequence is:", record.seq)
