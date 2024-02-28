#Step 1 of Pipeline Project 

import os

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

#Fastq-Dump 
do1_2 = 'fastq-dump -I --split-files SRR5660030'
os.system(do1_2)
do1_6 = 'fastq-dump -I --split-files SRR5660033'
os.system(do1_6)
do3_1= 'fastq-dump -I --split-files SRR5660044'
os.system(do3_1)
do3_6= 'fastq-dump -I --split-files SRR5660045'
os.system(do3_6)

#Subset the Fastq-Dump for test files 
S30_1 = 'head -n 40000 SRR5660030_1.fastq > sample_SRR5660030_1.fastq'
S30_2 = 'head -n 40000 SRR5660030_2.fastq > sample_SRR5660030_2.fastq'
S33_1 = 'head -n 40000 SRR5660033_1.fastq > sample_SRR5660033_1.fastq'
S33_2 = 'head -n 40000 SRR5660033_2.fastq > sample_SRR5660033_2.fastq'

S44_1 = 'head -n 40000 SRR5660044_1.fastq > sample_SRR5660044_1.fastq'
S44_2 = 'head -n 40000 SRR5660044_2.fastq > sample_SRR5660044_2.fastq'

S45_1 = 'head -n 40000 SRR5660045_1.fastq > sample_SRR5660045_1.fastq'
S45_2 = 'head -n 40000 SRR5660045_2.fastq > sample_SRR5660045_2.fastq'