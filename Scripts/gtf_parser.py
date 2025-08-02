# parser/gtf_parser.py
def parse_gtf(gtf_path):
    cds_entries = []
    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split('\t')
            if fields[2] != "CDS":
                continue
            chrom, start, end, strand = fields[0], int(fields[3]), int(fields[4]), fields[6]
            gene_id = fields[8].split('gene_id "')[1].split('"')[0]
            cds_entries.append({
                "chrom": chrom, "start": start, "end": end,
                "strand": strand, "gene": gene_id
            })
    return cds_entries
