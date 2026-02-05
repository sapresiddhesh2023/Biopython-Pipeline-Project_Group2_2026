from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import Entrez
import ssl
import os
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
ssl._create_default_https_context = ssl._create_unverified_context



from Bio.Blast import NCBIWWW, NCBIXML
from Bio import Entrez, SeqIO
import ssl


#Protein sequence
sequence = (
    "MKKAFYLHLIVFALQGISNVHSKSLLNQKALLPQQKSDMHINVNSLSDIVEPLISTVVSIYAVDTNIGISFNNKV"
    "SKYQQEVFLGSGVIIDSSGYIVTNENVIAGAENIKVKLHDGSELIAELVGSDNKINIALLKINSPAALSYATFGD"
    "SNQSRVGDQVIAIGSPFGLRGTVTNGIISSKGRDMGNGIVTDFIQTNAAIHMGSFGGPMFNLEGKIIGINSIHVSY"
    "SGISFAIPSNTVLEAVECLKKGEKIRRGMLNVMLNELTPELNENLGLKQDQNGVLITEVIKEGSAAQCGIAPGDV"
    "ITKFHDKAIKTGRDLQVAVSSTMLNSEREVELLRNGKSMTLKCKIIANKGEDSEQQSNDQSLVVNGVKFVDLTPD"
    "LVKKYNITSANNNGLFVLEVSPNSSWGRYGLKMGLRPRDIILSVKRDDNKKDISVKTLREIVTNIKHNEIFFTVQ"
    "RGDRMLYIALPNINK"
)

# Step 1: BLAST similarity search
print("Performing BLASTP search...")
result_handle = NCBIWWW.qblast(
    program="blastp",
    database="refseq_protein",
    sequence=sequence,
    hitlist_size=50
)

with open("blast_results.xml", "w") as out:
    out.write(result_handle.read())

print("BLAST completed and saved to blast_results.xml")

# Step 2: Find closest homologs 
closest_hits = []

with open("blast_results.xml") as handle:
    blast_record = NCBIXML.read(handle)

for alignment in blast_record.alignments:
    best_hsp = max(alignment.hsps, key=lambda h: h.bits)
    identity_fraction = best_hsp.identities / best_hsp.align_length
    # Filter for significant homologs
    if best_hsp.expect < 1e-5 and identity_fraction >= 0.35:
        closest_hits.append({
            "title": alignment.title,
            "accession": alignment.accession,
            "identity": identity_fraction,
            "evalue": best_hsp.expect,
            "q_start": best_hsp.query_start,
            "q_end": best_hsp.query_end
        })

print(f"Number of closest homologs found: {len(closest_hits)}")

#Step 3: Sequences of closest homologs
records = []
print("Fetching sequences of closest homologs...")
for hit in closest_hits[:10]:  # top 10 hits
    try:
        handle = Entrez.efetch(
            db="protein",
            id=hit["accession"],
            rettype="fasta",
            retmode="text"
        )
        record = SeqIO.read(handle, "fasta")
        records.append(record)
    except Exception as e:
        print(f"Failed to fetch {hit['accession']}: {e}")

SeqIO.write(records, "homologs.fasta", "fasta")
print("Sequences saved to homologs.fasta")

#Step 4: Detect conserved regions

if len(records) > 1:
    seq_length = min(len(r.seq) for r in records)  # trim to shortest
    conserved_positions = []

    for i in range(seq_length):
        column = [str(r.seq[i]) for r in records]
        most_common = max(set(column), key=column.count)
        conservation = column.count(most_common) / len(records)
        if conservation >= 0.8:  # >=80% conserved
            conserved_positions.append((i + 1, most_common))

    print("Highly conserved positions (1-based index):")
    print(conserved_positions[:20], "...")  # first 20 for brevity
else:
    print("Not enough sequences for conserved region analysis.")

# Step 5: Evolutionary hints (taxonomic distribution)

print("\nEvolutionary hints (species/organism info from top hits):")
for hit in closest_hits[:10]:
    print(f"- {hit['title']}")





