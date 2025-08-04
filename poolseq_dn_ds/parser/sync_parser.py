# poolseq_dn_ds/parser/sync_parser.py
def parse_sync_line(line, pop_index, min_coverage=10):
    fields = line.strip().split('\t')
    chrom, pos, ref = fields[:3]
    base_counts = fields[3:]

    if pop_index >= len(base_counts):
        raise IndexError(f"pop_index {pop_index} out of range: {len(base_counts)} populations available")

    try:
        counts = list(map(int, base_counts[pop_index].split(':')))
        A, T, C, G = counts[:4]
        total = A + T + C + G
    except Exception as e:
        raise ValueError(f"Error parsing base counts: {base_counts[pop_index]} — {str(e)}")

    if total < min_coverage:
        return None

    freqs = {
        'A': A / total,
        'T': T / total,
        'C': C / total,
        'G': G / total
    }

    major = max(freqs, key=freqs.get)
    minor = min(freqs, key=freqs.get)

    return {
        'chrom': chrom,
        'pos': int(pos),
        'ref': ref.upper(),
        'freqs': freqs,
        'major': major,
        'minor': minor,
        'minor_freq': freqs[minor]
    }
