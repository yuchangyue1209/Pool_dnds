# annotation/codon_classifier.py
from Bio.Data import CodonTable

standard_table = CodonTable.unambiguous_dna_by_name["Standard"]

def classify_codon_change(ref_codon, alt_codon):
    try:
        ref_aa = standard_table.forward_table[ref_codon]
        alt_aa = standard_table.forward_table[alt_codon]
        return "syn" if ref_aa == alt_aa else "nonsyn"
    except KeyError:
        return "unknown"
