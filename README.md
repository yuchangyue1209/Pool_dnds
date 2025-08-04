# Pool_dnds

Estimate **dN/dS** from Pool-Seq data using `.sync` files, GTF annotations, and reference FASTA sequences.

This tool enables calculation of non-synonymous (dN) and synonymous (dS) nucleotide diversity from pooled population sequencing data using codon-aware comparisons within annotated coding regions.

---

## Requirements
- Python 3.8+
- Biopython
- pandas

## Contact
Author: Changyue Yu cyu299@wisc.edu
GitHub: yuchangyue1209

## Input Files
1. --sync
A PoPoolation2-style .sync file containing allele counts per population.
Format:
chr  pos  ref  pop1_counts  pop2_counts  ...
Example base count format: 0:0:10:5:0:0 (A:T:C:G:N:del)

2. --gtf
Gene annotation file (GTF) with CDS features, e.g.:
chrI	BED	CDS	878338	879334	.	-	0	gene_id "cox7a1";
chrII	BED	CDS	4892361	4895354	.	-	0	gene_id "cox5a";

Must include:
Chromosome
Feature type = CDS
Start & end positions
Strand (+ or -)
Gene ID in the attribute field

3. --fasta
Reference genome in FASTA format, including all chromosomes present in the .sync file and GTF annotations.

## Usage
### ✅ Example Usage

```bash
python run.py \
  --sync /path/to/your.sync \
  --gtf /path/to/your.gtf \
  --fasta /path/to/reference.fa \
  --output /path/to/output_dnds.tsv \
  --pop-index 0 \
  --min-coverage 10 \
  --min-freq 0.1


## Parameters
  | Argument         | Description                                     |
| ---------------- | ----------------------------------------------- |
| `--sync`         | Path to `.sync` file                            |
| `--gtf`          | Path to gene annotation file (GTF format)       |
| `--fasta`        | Reference genome FASTA file                     |
| `--output`       | Path to output `.tsv` file                      |
| `--pop-index`    | Index of population in sync (default: `0`)      |
| `--min-coverage` | Minimum coverage per site (default: `10`)       |
| `--min-freq`     | Minimum minor allele frequency (default: `0.1`) |

## Output Format

| Column         | Description                                 |
| -------------- | ------------------------------------------- |
| `piN`          | Nucleotide diversity at nonsynonymous sites |
| `piS`          | Nucleotide diversity at synonymous sites    |
| `dN/dS`        | Ratio of πN / πS                            |
| `nonsyn_sites` | Number of nonsynonymous sites               |
| `syn_sites`    | Number of synonymous sites                  |


