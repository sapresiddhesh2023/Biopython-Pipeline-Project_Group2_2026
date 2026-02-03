#For this Project, in Step 1 (Biological Sequence Selection/Retrieval),we chose disease related protein.
#A gene from bacteria 'Orientia tsutsugamushi' (OT) which causes human disease 'Scrub typhus' was chosen.
#This gene encodes for a protein of the size 47kDa (molecular weight). It is known as TSA47 or bacterial htrA1.
#Biological role of TSA47 is not yet understood in OT's lif-cycle.
#Abbreviations: TSA: Type-Specific Antigen, htrA: High-Temperature Requirement Protein A
print("Step 2: performed by Siddhesh Uday Sapre")
print("-"*90)
from Bio import SeqIO
record = SeqIO.read("Otsutsugamushi_Karp_tsa47.gb", "genbank")                   
print("The following DNA record has Nucleotide ID:", record.id)
print(record.description)
print("The details for this GenBank file are:", record.annotations)
print("The number of attributes for this sequence are:", len(record.features))
print("The length of the sequence is:", len(record.seq), "bp")
print("The DNA sequence is", record.seq)
print("-"*30)
for feature in record.features: 
    print(feature.type,feature.location)
    print("-"*30)