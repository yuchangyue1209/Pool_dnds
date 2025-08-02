# poolseq_dn_ds/parser/sync_parser.py
import pandas as pd

def parse_sync_line(line, ref_base, pop_index):
    fields = line.strip().split('\t')
    chrom, pos, ref = fields[:3]
    base_counts = fields[3:]
    
    # Get counts for the selected population
    counts = base_counts[pop_index]
    A, T, C, G, N, del_ = map(int, counts.split(':'))
    total = A + T + C + G

    if total == 0:
        return None

    freq = {
        'A': A / total,
        'T': T / total,
        'C': C / total,
        'G': G / total
    }

    major = max(freq, key=freq.get)
    minor = min(freq, key=freq.get)

    return {
        'chrom': chrom,
        'pos': int(pos),
        'ref': ref.upper(),
        'freqs': freq,
        'major': major,
        'minor': minor,
        'minor_freq': freq[minor]
    }
