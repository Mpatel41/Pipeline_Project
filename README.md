# Pipeline_Project
This Pipeline is part of Comp 483 Computational Biology Assignment in Spring 2024

# Description
Human herpesvirus 5 is also known as Human cytomegalovirus and is typically abbreviated as HCMV. Cheng et al. 2017 (https://www.ncbi.nlm.nih.gov/pubmed/29158406) sequenced the transcriptomes of HCMV post infection. We would like to compare HCMV transcriptomes 2- and 6-days post-infection (dpi) using bioinformatics tools(Bowtie2, SPAdes, and Blast+). 

# Software 
 - [SRA ToolKit](https://github.com/ncbi/sra-tools)
	- Fastq-Dump
 - [Bowtie-2](https://github.com/BenLangmead/bowtie2)
 - [SPAdes](https://github.com/ablab/spades)
 - [Command Line Blast](https://www.ncbi.nlm.nih.gov/books/NBK279690/)

## Enviornment 
Used Python 


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

```python
git clone https://github.com/Mpatel41/Pipeline_Project.git 
```

2. Make sure you are in the Pipeline_Project directory using the following command

```cd Pipeline_Project ```

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

To view the PipelineProject.log output follow the instructions in the below section 


# Output
 
Within Pipeline_Project there will be a folder Pipeline_Project_Mansi_Patel. You can navigate there using the following command -

`cd Pipeline_Project_Mansi_Patel`

Within this directory you will find all the necessary files such as: 

1. The data used to generate the results (.fastq files)
2. The HCMV dataset that contains files related to HCMV genome to be used for bowtie indexing 
3. The bowtie mapped files (_mapped_ files) and sam files (.sam)
4. Spade_Assesmbly folder containing the scaffolds,contigs, and other data 
5. ncbi_dataset folder which contains the local database and blast results for Betaherpesvirinae subfamily. 
6. ProjectPipeline.log with contains -
 
	1. The number of reads in each transcriptome before and after the Bowtie2 mapping.
	2. SPAdes command used 
	3. number of contigs with a length > 1000 &  length of the assembly for contigs >1000bp
	4. Best alignment (HSP) for any single query-subject pair of sequences
		
		A.Subject accession
		
		B.Percent identity
		
		C.Alignment length
		
		D.Start of alignment in query
		
		E.End of alignment in query
		
		F.Start of alignment in subject
		
		G.End of alignment in subject
		
		H.Bit score
		
		I.E-value
		
		J. Subject Title 

To view this PipelineProject.log, follow the following command -

If you are in the Pipeline_Project_Mansi_Patel directory ignore this command. If you are not in this directory follow the command 

``` cd Pipeline_Project_Mansi_Patel```

Then open the PipelineProject.log file which contains results as indicated above

``` cat PipelineProject.log```
