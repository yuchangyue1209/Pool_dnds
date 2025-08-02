# poolseq_dn_ds/main.py

import argparse
from poolseq_dn_ds.parser.sync_parser import parse_sync_line
from poolseq_dn_ds.parser.gtf_parser import parse_gtf
from poolseq_dn_ds.reference.fasta_loader import load_reference_fasta
from poolseq_dn_ds.annotation.codon_classifier import classify_codon_change
from poolseq_dn_ds.estimator.dn_ds_calculator import estimate_dn_ds
from Bio.Seq import Seq

def run():
    parser = argparse.ArgumentParser(description="Estimate dN/dS from Pool-Seq sync + GTF + FASTA")
    parser.add_argument("--sync", required=True, help="Path to .sync file")
    parser.add_argument("--gtf", required=True, help="Path to GTF file")
    parser.add_argument("--fasta", required=True, help="Reference FASTA file")
    parser.add_argument("--output", required=True, help="Output file for dN/dS results")
    args = parser.parse_args()

    print("✅ Loading reference genome...")
    ref_seqs = load_reference_fasta(args.fasta)

    print("✅ Parsing GTF annotation...")
    cds_regions = parse_gtf(args.gtf)

    print("✅ Reading sync file...")
    variants = []
    with open(args.sync) as f:
        for line in f:
            chrom, pos, ref_base, pop_freqs = parse_sync_line(line)
            if not pop_freqs or pop_freqs[0] is None:
                continue

            freq_vec = pop_freqs[0]  # first population only for now
            bases = ['A', 'T', 'C', 'G']
            alt_bases = [b for b, f in zip(bases, freq_vec) if b != ref_base and f > 0.1]
            if not alt_bases:
                continue
            alt_base = alt_bases[0]

            for cds in cds_regions:
                if chrom == cds["chrom"] and cds["start"] <= pos <= cds["end"]:
                    gene = cds["gene"]
                    strand = cds["strand"]
                    cds_offset = pos - cds["start"]
                    codon_index = cds_offset // 3
                    codon_start = cds["start"] + codon_index * 3
                    ref_codon = ref_seqs[chrom][codon_start: codon_start + 3]

                    if strand == "-":
                        ref_codon = str(Seq(ref_codon).reverse_complement())

                    codon_pos = (pos - codon_start) if strand == "+" else (2 - (pos - codon_start))
                    alt_codon = list(ref_codon)
                    if 0 <= codon_pos < 3:
                        alt_codon[codon_pos] = alt_base
                    alt_codon = ''.join(alt_codon)

                    if strand == "-":
                        alt_codon = str(Seq(alt_codon).reverse_complement())

                    change_type = classify_codon_change(ref_codon, alt_codon)
                    if change_type == "unknown":
                        continue

                    pi = 2 * freq_vec[bases.index(alt_base)] * freq_vec[bases.index(ref_base)]
                    variants.append({"type": change_type, "pi": pi, "gene": gene})
                    break

    print("✅ Calculating dN/dS...")
    result = estimate_dn_ds(variants)

    print(f"💾 Writing output to {args.output}")
    with open(args.output, "w") as out:
        out.write("piN\tpiS\tdN/dS\tnonsyn_sites\tsyn_sites\n")
        out.write(f"{result['piN']:.5f}\t{result['piS']:.5f}\t{result['dN/dS']:.5f}\t{result['nonsyn_sites']}\t{result['syn_sites']}\n")

    return 0
