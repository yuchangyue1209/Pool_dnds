# reference/fasta_loader.py
from Bio import SeqIO

def load_reference_fasta(fasta_path):
    ref_seq = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        ref_seq[record.id] = str(record.seq).upper()
    return ref_seq
