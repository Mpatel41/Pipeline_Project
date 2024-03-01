# Pipeline_Project
This Pipeline is part of Comp 483 Computational Biology Assignment in Spring 2024

# Description
Human herpesvirus 5 is also known as Human cytomegalovirus and is typically abbreviated as HCMV. Cheng et al. 2017 (https://www.ncbi.nlm.nih.gov/pubmed/29158406) sequenced the transcriptomes of HCMV post infection. We would like to compare HCMV transcriptomes 2- and 6-days post-infection (dpi) using Bowtie2, SPAdes, and Blast+). 

# Software 
 - [SRA ToolKit][https://github.com/ncbi/sra-tools]
	- Fastq-Dump
 - [Bowtie-2][https://github.com/BenLangmead/bowtie2]
 - [SPAdes][https://github.com/ablab/spades]
 - [Command Line Blast][https://www.ncbi.nlm.nih.gov/books/NBK279690/]

# Packages 
 - os
 - glob
 - subprocess
 - Bio
 - SeqIO
 - re
 - shutil
 - pandas 
 - sys 
 - argv
# Getting Started 
1. Use the following command to clone this repository to your local device

`git clone https://github.com/Mpatel41/Pipeline_Project.git`

2. Make sure you are in the Pipeline_Project directory using the following command

`cd Pipeline_Project` 

3. 2 Options:
- Run the Script with Test Data (which is provided)
- Run the Script with Actual Data 

4. Option 1 : Run Script with Test Data that is provided - Run the command below 
 `python track2_assesmbly.py test_data`

The test data was derived from the script subset_test.py (no need to run this script as the test data was generated and provided)

5. Option 2: Run Script with Actual Data. Run the following script
 `nohup python transcriptome_files.py &`

Run nohup as this script will take ~10 minutes. After that, run the following command -
 `nohup python track2_assesmbly.py transcriptome_data &`

This script will take ~2-3 hours to run. 

 
