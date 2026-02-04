from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio import SeqIO
import ssl
ssl._create_default_https_context = ssl._create_unverified_context


record = SeqIO.read("Kprotein.fasta","fasta")
print(record.id)
print(record.description)
print(len(record))


from Bio.Blast import NCBIXML
with open ("blast_result.xml") as b:
    blast_rec = NCBIXML.read(b)

print(len(blast_rec.alignments)) 

#print("The total aminoacis are:",1401)
#print("The total number of alignments are:",50)
#print("-"*120)

#50 sequences are similar with the query sequence

from Bio.Blast import NCBIXML

with open("blast_result.xml") as result_handle:
    blast_rec = NCBIXML.read(result_handle)

for alignment in blast_rec.alignments:
    print("\nHomolog:", alignment.title)
    print("Length:", alignment.length)

    for hsp in alignment.hsps:
        print(" E-value:", hsp.expect)
        print(" Identity:", hsp.identities, "/", hsp.align_length)
        print(" Query region:", hsp.query_start, "-", hsp.query_end)
        print(" Subject region:", hsp.sbjct_start, "-", hsp.sbjct_end)
