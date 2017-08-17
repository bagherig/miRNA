# Project Title
Genetic Variation in miRNA Primary Transcripts

# Getting Started
These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

# Prerequisites
The R code requires a specific directory structure. Create an R project (anywhere) from this repository and follow the structure below:

     ─── miRNA                   # The project directory.
         ├── R                   # The folder containing .R files
         ├── miRNA.gff3          # The file containing miRNA genome coordinates.
         └── VCF                 # The folder containing VCF files. All VCF files for a specific population must be
                                   placed in a folder named as the population code. 
               ├── ALL           # The folder containing the VCF data for all populations. Downloaded from 1000 
                                   genomes project.           
               ├── ACB  # The folder containing the VCF files for population 1. Each of these folders must
                                   contain the VCF files for this population seperated by chromosome number (1-22).
                   ├── chr1.vcf.gz
                   ├── chr2.vcf.gz
                   ├── chr3.vcf.gz
                   .
                   .
                   .
                   └── chr22.vcf.gz
               ├── ASW
               ├── BEB
               .
               .
               .
               └── YRI

# Preparing Data
1. Obtain the VCF folder with the correct directory structure from this repository. This folder contains 'samples.txt' for each population, which contains the sample names for that population.
2. Download the original VCF files from ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/, and store them in a folder named 'ALL' under the 'VCF' folder.
3. Run the following command in BASH from the 'VCF' directory. This command subsets the VCF data for one population from the original VCF files. Run this command once for each population (Replace '???' with populations code. ex: 'ACB').

          for file in ALL/*.vcf.gz; do echo "Subsetting $(basename $file)"; bcftools view --min-ac=1 --force-samples -Oz -S ???/samples.txt $file > ???/$(basename $file); done

4. Download miRNA.gff3.zip from this repository and unzip it under the 'miRNA' directory.

You are now ready to run the analysis using main.R.

# Running the analysis
Simply follow main.R to run the analysis.
