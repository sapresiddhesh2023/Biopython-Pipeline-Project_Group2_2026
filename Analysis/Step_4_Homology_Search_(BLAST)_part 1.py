from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import ssl
from Bio.SeqRecord import SeqRecord
from Bio import Entrez
Entrez.email = "sameenkhaled@gmail.com"

ssl._create_default_https_context = ssl._create_unverified_context

# Full protein sequence
sequence = (
    "MKKAFYLHLIVFALQGISNVHSKSLLNQKALLPQQKSDMHINVNSLSDIVEPLISTVVSIYAVDTNIGISFNNKV"
    "SKYQQEVFLGSGVIIDSSGYIVTNENVIAGAENIKVKLHDGSELIAELVGSDNKINIALLKINSPAALSYATFGD"
    "SNQSRVGDQVIAIGSPFGLRGTVTNGIISSKGRDMGNGIVTDFIQTNAAIHMGSFGGPMFNLEGKIIGINSIHVSY"
    "SGISFAIPSNTVLEAVECLKKGEKIRRGMLNVMLNELTPELNENLGLKQDQNGVLITEVIKEGSAAQCGIAPGDV"
    "ITKFHDKAIKTGRDLQVAVSSTMLNSEREVELLRNGKSMTLKCKIIANKGEDSEQQSNDQSLVVNGVKFVDLTPD"
    "LVKKYNITSANNNGLFVLEVSPNSSWGRYGLKMGLRPRDIILSVKRDDNKKDISVKTLREIVTNIKHNEIFFTVQ"
    "RGDRMLYIALPNINK"
)

print("Performing BLASTP search...")

# Performing BLASTP
result_handle = NCBIWWW.qblast(
    program="blastp",
    database="refseq_protein",
    sequence=sequence,
    hitlist_size=50
)

# Saving results
with open("blast_results.xml", "w") as b:
    b.write(result_handle.read())

print("BLAST completed and saved to blast_results.xml")

