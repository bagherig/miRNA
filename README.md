# Project Title
Genetic Variation in miRNA Primary Transcripts

# Getting Started
These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

# Prerequisites
The R code requires a specific directory structure. Create an R project (anywhere) from this repository and follow the structure below:

     ├── miRNA                   # The project directory.
         ├── R                   # The folder containing .R files
         ├── miRNA.gff3          # The file containing miRNA genome coordinates.
         ├── VCF                 # The folder containing VCF files. All VCF files for a specific population must be
                                   placed in a folder named as the population code.                          
               ├── population 1  # The folder containing the VCF files for population 1. Each of these folders must
                                   contain the VCF files for this population seperated by chromosome number (1-22).
                              ├── chr1.vcf.gz
                              ├── chr2.vcf.gz
                              ├── chr3.vcf.gz
                                      .
                                      .
                                      .
                              ├── chr22.vcf.gz
               ├── population 2
               ├── population 3
                        .
                        .
                        .
               ├── population n


# Running the analysis
Simply follow main.R to run the analysis.
