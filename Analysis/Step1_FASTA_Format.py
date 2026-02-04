#For this Project, in Step 1 (Biological Sequence Selection/Retrieval),we chose disease related protein.
#A gene from bacteria 'Orientia tsutsugamushi' (OT) which causes human disease 'Scrub typhus' was chosen.
#This gene encodes for a protein of the size 47kDa (molecular weight). It is known as TSA47 or bacterial htrA1.
#Biological role of TSA47 is not yet understood in OT's lif-cycle.
#Abbreviations: TSA: Type-Specific Antigen, htrA: High-Temperature Requirement Protein A
print("Step 1: performed by Siddhesh Uday Sapre")
print("-"*90)
from Bio import SeqIO
record = SeqIO.read("Otsutsugamushi_Karp_tsa47.fasta", "fasta")
print("The following protein record has Nucleotide ID:", record.id)
print(record.description)
print("The length of the sequence is:", len(record.seq), "bp")
print("The protein sequence is:", record.seq)
