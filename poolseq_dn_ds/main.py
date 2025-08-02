  # main.py
from parser.sync_parser import parse_sync_line
from parser.gtf_parser import parse_gtf
from annotation.codon_classifier import classify_codon_change
from reference.fasta_loader import load_reference_fasta
from estimator.dn_ds_calculator import estimate_dn_ds

def main(sync_file, gtf_file, fasta_file):
    ref_seqs = load_reference_fasta(fasta_file)
    cds_regions = parse_gtf(gtf_file)
    variants = []

    with open(sync_file) as f:
        for line in f:
            chrom, pos, ref_base, pop_freqs = parse_sync_line(line)
            if not pop_freqs or None in pop_freqs:
                continue
            # TODO: implement mapping SNP to codon
            # → find CDS, extract codon, simulate mutation
            # → classify syn/nonsyn, compute π
            # variants.append({"type": "syn", "pi": 0.001}) or nonsyn

    result = estimate_dn_ds(variants)
    print(result)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--sync", required=True, help="Path to sync file")
    parser.add_argument("--gtf", required=True, help="Path to GTF file")
    parser.add_argument("--fasta", required=True, help="Reference FASTA")
    args = parser.parse_args()

    main(sync_file=args.sync, gtf_file=args.gtf, fasta_file=args.fasta)
