# parser/sync_parser.py
def parse_sync_line(line, min_coverage=10):
    fields = line.strip().split()
    chrom, pos, ref_base = fields[:3]
    pops = fields[3:]
    result = []

    for pop in pops:
        counts = list(map(int, pop.split(':')))
        total = sum(counts[:4])  # only A/T/C/G
        if total < min_coverage:
            result.append(None)
            continue
        freqs = [count / total for count in counts[:4]]
        result.append(freqs)
    return chrom, int(pos), ref_base, result
