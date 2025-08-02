# Pool_dnds
Developed for pooled data dN/dS calculation 
Estimate dN/dS ratios from Pool-Seq data using SNP frequencies and gene annotations.

## Structure

- `parser/`: sync and GTF parsers
- `annotation/`: classify codon changes (syn/nonsyn)
- `reference/`: load reference genome
- `estimator/`: compute dN/dS from SNP-level π values

## Usage

```bash
python main.py --sync your.sync --gtf genes.gtf --fasta genome.fa
