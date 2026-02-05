from Bio.Blast import NCBIWWW
import ssl

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

# Perform BLASTP
result_handle = NCBIWWW.qblast(
    program="blastp",
    database="refseq_protein",
    sequence=sequence,
    hitlist_size=50
)
with open("blast_results.xml", "w") as b:
    b.write(result_handle.read())

print("BLAST completed and saved to blast_results.xml")


