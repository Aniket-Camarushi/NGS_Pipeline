# This script checks the quality of our fastq files and performs an alignment to the human cDNA transriptome reference with Kallisto.
# To run this 'shell script' you will need to open your terminal and navigate to the directory where this script resides on your computer.
# This should be the same directory where your fastq files and reference fasta files are found.
# Change permissions on your computer so that you can run a shell script by typing: 'chmod u+x readMapping.sh' (without the quotes) at the terminal prompt
# Then type './readMapping.sh' (without the quotes) at the prompt.
# This will begin the process of running each line of code in the shell script.

# Exit on error
set -e

# First use fastqc to check the quality of our fastq files:
echo "Running FastQC..."
mkdir -p FastQC
fastqc ./Data/*.gz -t 20 -o FastQC/

# Next, we want to build an index from our reference fasta file
# I get my reference mammalian transciptome files from here: https://useast.ensembl.org/info/data/ftp/index.html
echo "Building kallisto index..."
mkdir -p Kallisto
kallisto index -i ./Kallisto/Homo_sapiens.GRCh38.cdna.all_my.index ./Data/Homo_sapiens.GRCh38.cdna.all.fa

# Now map reads to the indexed reference host transcriptome
# Use as meany 'threads' as your machine will allow in order to speed up the read mapping process
# Use 30-60 bootstraps to model technical variance for each sample
# Note that we're also including the '&>' at the end of each line
## this takes the information that would've been printed to our terminal & outputs this in a log file that is saved in your parent directory
echo "Running kallisto quantifications..."

# First the healthy subject (H5)
kallisto quant -i ./Kallisto/Homo_sapiens.GRCh38.cdna.all_my.index  -o ./Kallisto/HS01 -t 26 --single -l 250 -s 30 ./Data/SRR8668755.fastq.gz &> ./Kallisto/HS01.log
kallisto quant -i ./Kallisto/Homo_sapiens.GRCh38.cdna.all_my.index  -o ./Kallisto/HS02 -t 26 --single -l 250 -s 30 ./Data/SRR8668756.fastq.gz &> ./Kallisto/HS02.log
kallisto quant -i ./Kallisto/Homo_sapiens.GRCh38.cdna.all_my.index  -o ./Kallisto/HS03 -t 26 --single -l 250 -s 30 ./Data/SRR8668757.fastq.gz &> ./Kallisto/HS03.log
kallisto quant -i ./Kallisto/Homo_sapiens.GRCh38.cdna.all_my.index  -o ./Kallisto/HS04 -t 26 --single -l 250 -s 30 ./Data/SRR8668758.fastq.gz &> ./Kallisto/HS04.log
kallisto quant -i ./Kallisto/Homo_sapiens.GRCh38.cdna.all_my.index  -o ./Kallisto/HS05 -t 26 --single -l 250 -s 30 ./Data/SRR8668759.fastq.gz &> ./Kallisto/HS05.log

# Then the cutaneous leishmaniasis (CL) patients
kallisto quant -i ./Kallisto/Homo_sapiens.GRCh38.cdna.all_my.index  -o ./Kallisto/CL08 -t 26 --single -l 250 -s 30 ./Data/SRR8668769.fastq.gz &> ./Kallisto/CL08.log
kallisto quant -i ./Kallisto/Homo_sapiens.GRCh38.cdna.all_my.index  -o ./Kallisto/CL10 -t 26 --single -l 250 -s 30 ./Data/SRR8668771.fastq.gz &> ./Kallisto/CL10.log
kallisto quant -i ./Kallisto/Homo_sapiens.GRCh38.cdna.all_my.index  -o ./Kallisto/CL11 -t 26 --single -l 250 -s 30 ./Data/SRR8668772.fastq.gz &> ./Kallisto/CL11.log
kallisto quant -i ./Kallisto/Homo_sapiens.GRCh38.cdna.all_my.index  -o ./Kallisto/CL12 -t 26 --single -l 250 -s 30 ./Data/SRR8668773.fastq.gz &> ./Kallisto/CL12.log
kallisto quant -i ./Kallisto/Homo_sapiens.GRCh38.cdna.all_my.index  -o ./Kallisto/CL13 -t 26 --single -l 250 -s 30 ./Data/SRR8668774.fastq.gz &> ./Kallisto/CL13.log

# summarize fastqc and kallisto mapping results in a single summary html using MultiQC
echo "Running MultiQC summary..."
mkdir -p multiqc
multiqc . -o multiqc/

echo "Finished"

