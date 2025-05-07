# NGS Pipeline Project
*A modular workflow for academic NGS analysis (FASTQ → VCF)*  
**Developer**: Aniket Camarushi | UMBC Bioinformatics '26  

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

## How to Use
```bash
# For academic/learning purposes only
./run_pipeline.sh -i /path/to/fastq -o results_dir
```
