# NGS Pipeline Project
*A modular workflow for academic NGS analysis (FASTQ → VCF)*  
**Developer**: Aniket Camarushi | UMBC Bioinformatics '26  
*Please switch to `version_2` branch for code.*

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
