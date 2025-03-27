# Metagenomic Virus Identification Pipeline using Snakemake
You can find laboratory report [here](https://docs.google.com/document/d/1q_JpDl_YqxHrE-WEVlzrE4zzZ50MNCcp/edit?usp=sharing&ouid=110506015076070664163&rtpof=true&sd=true).    
This pipeline is designed to replicate the bioinformatics workflow used in identifying and characterizing novel viruses from sequencing data. Specifically, it follows simplified steps analogous to those used in the original identification of SARS-CoV-2.

## Overview

The pipeline performs the following analyses:
- Genome Assembly Quality Check: Evaluates assembled contigs using QUAST v5.2.0.
- Contig Indexing and Extraction: Uses samtools v1.21 to index and extract specific contigs from assembled sequences.
- Similarity Search: Confirms viral identity by querying contigs against known viruses using BLAST. (website)
- Genome Annotation: Annotates the selected viral genome using Prokka v1.14.6.
- Primer Design (optional): Designs PCR primers targeting conserved viral regions using Primer-BLAST.(website)

## Prerequisites

### Software Requirements
Unix-based OS (Linux/Mac)
Conda environment installed
Snakemake installed (conda install -c bioconda snakemake)


## Input Files Setup

### Directory Structure
```
├── Snakefile                  # Main Snakemake workflow
├── samples.csv                # Metadata for sample files
└──  params.json                # JSON file with pipeline parameters         


```

### samples.csv Format
The `samples.csv` file should contain the following columns:
```csv
sample_id,read_1,read_2,assembly
sample1,/path/to/sample1_R1.fastq,/path/to/sample1_R2.fastq,
sample2,/path/to/sample2_R1.fastq,/path/to/sample2_R2.fastq,
sample3_assembly,,,/path/to/sample3.fasta
```

Notes:
- For raw sequencing data: provide `read_1` and `read_2` paths
- For pre-assembled genomes: add `_assembly` suffix to `sample_id` and provide `assembly` path
- Paths can be relative or absolute

### params.json Configuration
```json
{
    "global_params": {"outdir": "./output", "threads": 8},
    "quast": {"outdir": "./output"},
    "samtools": {"prefix": "faidx", "contig_name": "NODE_1_length_29907_cov_150.822528"},
    "blast": {"db": "nt", "entrez_query": "1900/01/01:2020/01/01[PDAT] AND viruses[ORGN]"},
    "prokka": {"outdir": "./output"}
}

```

## Output Files

The pipeline generates the following directory structure for each sample:

For samples with raw reads:
```
results/
├── {sample}_fastqc/
│   ├── {sample}_1_fastqc.html
│   └── {sample}_2_fastqc.html
├── {sample}_spades/
│   └── scaffolds.fasta
├── {sample}_quast/
│   └── report.html
├── {sample}_prokka/
│   └── annotation.txt
└── {sample}_abricate/
    ├── ncbi_results.txt
    └── resfi_results.txt
```

For pre-assembled samples:
```
results/
├── {sample}_assembly_quast/
│   └── report.html
├── {sample}_assembly_prokka/
    └── annotation.txt

```

## Running the Pipeline

### Basic Usage
```bash
snakemake --cores N --use-conda
```

### Common Options
- `--cores N`: Specify number of cores to use
- `--use-conda`: Activate conda environment management
- `--dryrun`: Show execution plan without running
- `--printshellcmds`: Print shell commands that will be executed

### Examples
1. Dry run to check workflow:
```bash
snakemake -n --use-conda
```

2. Run with 8 cores:
```bash
snakemake --cores 8 --use-conda
```

3. Run specific rule for a sample:
```bash
snakemake --cores 1 --use-conda results/sample1_fastqc/sample1_1_fastqc.html
```


## Notes

- The pipeline automatically detects whether a sample needs assembly (based on presence of raw reads) or starts from a provided assembly
- All paths in the output are determined by the `outdir` parameter in params.json
- Each tool runs in its own conda environment to manage dependencies
