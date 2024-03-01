#Import the necessary libraries 
import glob
import os

#Create a folder to save the data 
os.mkdir('transcriptome_data')
os.chdir('transcriptome_data/')

#Getting the Transcriptomes from SRA 

#Donor 1 (2dpi)
donor1_2dpi ='wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030'
os.system(donor1_2dpi)

#Donor 1 (6dpi)
donor1_6dpi = 'wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660033/SRR5660033'
os.system(donor1_6dpi)

#Donor 3 (2dpi)
donor3_2dpi = 'wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660044/SRR5660044'
os.system(donor3_2dpi)

#Donor 3 (6dpi)
donor3_6dpi = 'wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660045/SRR5660045'
os.system(donor3_6dpi)

#Get all the SRR files together
donor_pattern = "SRR*"
donor_files = glob.glob(donor_pattern)
donor_files.sort()

#Run Fastq Dump 
for file in donor_files:
  fastq_dump = 'fastq-dump -I --split-files ' + file
  os.system(fastq_dump)
  



