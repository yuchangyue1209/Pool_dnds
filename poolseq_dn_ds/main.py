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
    parser.add_argument("--pop-index", type=int, default=0, help="Population index in sync (0-based)")
    parser.add_argument("--min-coverage", type=int, default=10, help="Minimum coverage threshold")
    args = parser.parse_args()

    print("✅ Loading reference genome...")
    ref_seqs = load_reference_fasta(args.fasta)

    print("✅ Parsing GTF annotation...")
    cds_regions = parse_gtf(args.gtf)

    print("✅ Reading sync file and classifying variants...")
    variants = []
    with open(args.sync) as f:
        for line in f:
            try:
                parsed = parse_sync_line(line, args.pop_index, min_coverage=args.min_coverage)
            except Exception as e:
                continue  # skip malformed lines

            if not parsed:
                continue

            chrom = parsed["chrom"]
            pos = parsed["pos"]
            ref_base = parsed["ref"]
            freqs = parsed["freqs"]
            alt_base = parsed["minor"]
            alt_freq = parsed["minor_freq"]

            if alt_freq < 0.1:
                continue  # skip low-frequency variants

            for cds in cds_regions:
                if chrom == cds["chrom"] and cds["start"] <= pos <= cds["end"]:
                    gene = cds["gene"]
                    strand = cds["strand"]

                    # Find codon start based on strand
                    cds_offset = pos - cds["start"]
                    codon_index = cds_offset // 3
                    codon_start = cds["start"] + codon_index * 3
                    ref_codon = ref_seqs[chrom][codon_start: codon_start + 3]

                    if len(ref_codon) != 3 or "N" in ref_codon:
                        continue

                    if strand == "-":
                        ref_codon = str(Seq(ref_codon).reverse_complement())

                    # Mutate codon
                    codon_pos = (pos - codon_start) if strand == "+" else (2 - (pos - codon_start))
                    if not 0 <= codon_pos < 3:
                        continue
                    alt_codon = list(ref_codon)
                    alt_codon[codon_pos] = alt_base
                    alt_codon = ''.join(alt_codon)

                    if strand == "-":
                        alt_codon = str(Seq(alt_codon).reverse_complement())

                    change_type = classify_codon_change(ref_codon, alt_codon)
                    if change_type == "unknown":
                        continue

                    pi = 2 * freqs[ref_base] * freqs[alt_base]
                    variants.append({"type": change_type, "pi": pi, "gene": gene})
                    break  # only assign one CDS

    print("✅ Calculating dN/dS...")
    result = estimate_dn_ds(variants)

    print(f"💾 Writing output to {args.output}")
    with open(args.output, "w") as out:
        out.write("piN\tpiS\tdN/dS\tnonsyn_sites\tsyn_sites\n")
        out.write(f"{result['piN']:.5f}\t{result['piS']:.5f}\t{result['dN/dS']:.5f}\t{result['nonsyn_sites']}\t{result['syn_sites']}\n")

    return 0

