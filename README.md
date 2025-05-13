# NGS Pipeline Project
*A modular workflow for academic NGS analysis (FASTQ → VCF)*  
**Developer**: Aniket Camarushi | UMBC Bioinformatics & Computational Biology '26  

## Overview
This pipeline was developed during my undergraduate research to:
- Process small batches of RNA-seq data (10-15 samples)  
- Demonstrate core NGS workflow principles (STAR alignment → GATK variant calling)  
- Automate QC reporting with MultiQC

## Key Features
| Module       | Purpose                          | Status       |
|--------------|----------------------------------|--------------|
| FastQC       | Raw read quality control         | Stable       |
| STAR         | Genome alignment                 | Functional*  |
| GATK         | Variant calling                  | In refinement|
| MultiQC      | Consolidated QC reporting        | Production-ready |

*\*Some edge cases being optimized for larger datasets*

## Data Availability

### Datasets:

[SRR8668755](https://trace.ncbi.nlm.nih.gov/Traces/?run=SRR8668755)

[SRR8668756](https://trace.ncbi.nlm.nih.gov/Traces?run=SRR8668756)

[SRR8668757](https://trace.ncbi.nlm.nih.gov/Traces?run=SRR8668757)

[SRR8668758](https://trace.ncbi.nlm.nih.gov/Traces?run=SRR8668758)

[SRR8668759](https://trace.ncbi.nlm.nih.gov/Traces?run=SRR8668759)

[SRR8668769](https://trace.ncbi.nlm.nih.gov/Traces?run=SRR8668769)

[SRR8668771](https://trace.ncbi.nlm.nih.gov/Traces?run=SRR8668771)

[SRR8668772](https://trace.ncbi.nlm.nih.gov/Traces?run=SRR8668772)

[SRR8668773](https://trace.ncbi.nlm.nih.gov/Traces?run=SRR8668773)

[SRR8668774](https://trace.ncbi.nlm.nih.gov/Traces?run=SRR8668774)



### Other files required for the pipeline

[Homo_sapiens.GRCh38.cdna.all.fa](https://ftp.ensembl.org/pub/release-114/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz)


Below is a snippet of bash which is required to run before running `run_pipeline.sh`:

```bash
wget ftp://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -P ./Data/
gunzip ./Data/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

wget ftp://ftp.ensembl.org/pub/release-108/gtf/homo_sapiens/Homo_sapiens.GRCh38.108.gtf.gz -P ./Data/
gunzip ./Data/Homo_sapiens.GRCh38.108.gtf.gz

# HapMap (for SNP training)
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz -P ./Data/
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz.tbi -P ./Data/

# Omni (for SNP training)
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz -P ./Data/
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz.tbi -P ./Data/

# 1000G SNPs (for SNP training)
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz -P ./Data/
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi -P ./Data/
```


## Dependencies
```bash
conda create -n rna_seq_pipeline -c bioconda -c conda-forge \
    fastqc \
    kallisto \
    multiqc \
    star \
    samtools \
    gatk4 \
    picard \
    bcftools \
    trimmomatic \
    snpeff \
    openssl=1.0  # Needed for some older GATK tools

conda activate rna_seq_pipeline
```

## How to Use
```bash
# Create 'Data' directory and input all .fastq files with index files
mkdir -p Data

# For academic/learning purposes only
./run_pipeline.sh
```
