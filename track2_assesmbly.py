#Downloading the Genome HCPV
#Track 2 

#Import the necessary libraries 
import os
import glob
import subprocess
from Bio import SeqIO
import re 
import shutil
import pandas as pd
import sys
from sys import argv

# Get either the actual data or test data based on what was provided 
if len(sys.argv) != 2:
    print("Usage: python file.py <folder1>")
    sys.exit(1)

# Get the folder paths from above
folder = sys.argv[1]

#Create a Working Directory that with contain all the output files 
os.mkdir('Pipeline_Project_Mansi_Patel')

#Moving the Fastq files to our new working directory /Pipeline_Project_Mansi_Patel
original = os.getcwd() 
which_path = os.path.join(original, folder)
new = "Pipeline_Project_Mansi_Patel/"
for filename in os.listdir(which_path):
    if filename.endswith(".fastq"):
        source_file = os.path.join(which_path, filename)
        destination_file = os.path.join(original,new, filename)
        shutil.move(source_file, destination_file)

os.chdir('Pipeline_Project_Mansi_Patel/')

#Import the dataset HCMV, unzip, build index
rename_after = "HCMV_ncbi_dataset"
HCMV = 'datasets download genome accession GCF_000845245.1 --include gff3,rna,cds,protein,genome,seq-report --filename HCMV_datasets.zip'
unzip_HCMV = 'unzip HCMV_datasets.zip'
bowtie_build_index = 'bowtie2-build ' +rename_after+ '/data/GCF_000845245.1/GCF_000845245.1_ViralProj14559_genomic.fna beta_herpesvirus5'

#Run the dataset HCMV, unzip and build index using os.system()
#Rename the ncbi_dataset folder as it saves that as default (avoids trouble later on)
os.system(HCMV)
os.system(unzip_HCMV)
os.rename("ncbi_dataset",rename_after)
os.system(bowtie_build_index)
os.remove('README.md')


#Still in the Pipeline_Project Working directory 
#The fastq files that were copied will be used here 
#Using Glob to get all the .fastq files in a list 
fastq_pattern = "*.fastq"
fastq_files = glob.glob(fastq_pattern)
fastq_files.sort() #sort so its paired _1 and _2 

#Get the filenames (prefix,suffix to match and output)
for item in fastq_files: 
    pre_pattern = r"SRR(\d+)" #extract the SRR name 
    suff_pattern = r"1.fastq" #extract the one with 1.fastq 
    match = re.search(pre_pattern,item) #search the pattern in item 
    suf_match = re.search(suff_pattern,item) #search the pattern in item

    if match:
        prefix = match.group(0) #if matched I want to store ex. SRR5660030 name in the prefix 

#IF it is a _1.fasta file, I want to get the _2.fasta matching pair to use in generating output 
    if suf_match:
        potential_fastq_2 = rf"{prefix}_2.fastq"
        potenital_match = suf_match.group()
        suffix_file = [file for file in fastq_files if potential_fastq_2 in file] #gets the _2.fastq file 

#Bowtie command in command variable & unzip the file 
        if suffix_file:
            fastq_1 = item #_1.fastq file 
            fastq2 = str(suffix_file).strip('[]').strip('\'"') #get the _2.fastq file and format the name
            sam_output = f"{prefix}_HCMV.sam" #sam output format 
            mapped_output = f"{prefix}" #map output format 
            command = 'bowtie2 --quiet -x beta_herpesvirus5 -1 ' +fastq_1+ ' -2 ' +fastq2+ ' -S ' +sam_output+ ' --al-conc-gz ' + mapped_output + '_mapped_%.fq.gz'
            unzipped_1 = 'gunzip ' + mapped_output + '_mapped_1.fq' 
            unzipped_2 = 'gunzip ' + mapped_output + '_mapped_2.fq'

            #Execute the commands for Bowtie & Unzip 
            os.system(command)
            os.system(unzipped_1)
            os.system(unzipped_2)

            #Check how many read pairs before and after filtering 
            before_filtered = f'grep ' + '@' +prefix + ' ' + item + ' | wc -l' 
            after_filtered = f'grep ' + '@' + prefix + ' ' + mapped_output + '_mapped_1.fq | wc -l'
            
            #match the names to write to file 
            if 'SRR5660030' in prefix:
                donor = 'Donor 1 (2dpi)'
            elif 'SRR5660033' in prefix:
                donor = 'Donor 1 (6dpi)'
            elif 'SRR5660044' in prefix:
                donor = 'Donor 3 (2dpi)'
            elif 'SRR5660045' in prefix:
                donor = 'Donor 3 (6dpi)'
            else:
                donor = prefix

            #open log file and write the output of before and after filtered 
            with open('PipelineProject.log','a') as f:
                output = subprocess.getoutput(before_filtered)
                output_after = subprocess.getoutput(after_filtered)
                f.write(donor + ' had ' + str(output)+ ' read pairs before Bowtie2 filtering and ' + str(output_after) + ' read pairs after' + '\n')
                f.write('\n')

#Running Spades 
#Get all the .fq mapping files using glob 
spades_pattern = "*.fq"
spades_files = glob.glob(spades_pattern)
spades_files.sort()

#Assign them so they are paired (sorted so it will be in pairs)
map1_1 = spades_files[0]
map1_2 = spades_files[1]
map2_1 = spades_files[2]
map2_2 = spades_files[3]
map3_1 = spades_files[4]
map3_2 = spades_files[5]
map4_1 = spades_files[6]
map4_2 = spades_files[7]

#Run Spades 
os.mkdir("Spade_Assembly")

#Spades Command & execute it 
spades = 'spades.py -k 77,99,127 -t 2 --only-assembler --pe1-1 ' + map1_1 + ' --pe1-2 ' + map1_2 + ' --pe2-1 ' + map2_1 + ' --pe2-2 ' + map2_2 + ' --pe3-1 ' + map3_1 + ' --pe3-2 ' + map3_2 + ' --pe4-1 ' + map4_1 + ' --pe4-2 ' + map4_2 + ' -o Spade_Assembly/'
os.system(spades)

#append output to log file 
with open('PipelineProject.log','a') as f:
    f.write(spades + '\n')
    f.write('\n')

#Counting the number of contigs >1000
#Go into the Spade_Assembly folder 
os.chdir('Spade_Assembly/')

#The first sequence has the greatest length so I am storing that to use for blast
records = list(SeqIO.parse('contigs.fasta','fasta'))
id = records[0].id 
seq = records[0].seq

#Get the record.id and the 3rd index says the sequence length
#Get Contigs > 1000 and add the sequence length of those 
greater_contigs = 0 
total_length = 0 
for record in records:
    description = record.id
    splitted = description.split('_')
    length = splitted[3] #the length of the sequence 

    if int(length) > 1000:
        greater_contigs += 1
        total_length = total_length + int(length)

#Go back to home directory 
os.chdir('..')

#Write to log file
with open('PipelineProject.log','a') as f: 
    f.write('There are ' + str(greater_contigs) + ' contigs > 1000bp in the assembly. ' + '\n')
    f.write('There are ' + str(total_length) + ' bp in the assembly' + '\n')
    f.write('\n')

#Blast Command (Download the dataset)
subfamily = 'betaherpesvirinae'
download = 'datasets download virus genome taxon ' + subfamily + ' --include genome --filename betaherpes_blast.zip'
download_unzip = 'unzip betaherpes_blast.zip'
os.system(download)
os.system(download_unzip)

#Make local database 
os.chdir('ncbi_dataset/data')
fastq_f = 'genomic.fna'
database = 'makeblastdb -in ' + fastq_f + ' -out ' + subfamily + ' -title ' + subfamily + ' -dbtype nucl'
os.system(database)

#Create a fasta file with the longest contig length
with open('large.fasta', 'w') as output_handle:
    output_handle.write(">" + str(id) + "\n")
    output_handle.write(str(seq) + "\n")
    output_handle.close()

#Using the longest contig sequence we stored in a file run blast 
query_seqfile = 'large.fasta'
output_file = "my_results.txt"
blast_command = 'blastn -query ' +query_seqfile+ ' -db '+subfamily+ ' -out '+output_file + ' -max_hsps 1 -max_target_seqs 10 -outfmt "6 sacc pident length qstart qend sstart send bitscore evalue stitle"'
os.system(blast_command)

#Read the output blast file as a dataframe and add in header 
df = pd.read_csv("my_results.txt", sep='\t', header=None, names=["sacc", "pident", "length", "qstart", "qend", "sstart", "send", "bitscore", "evalue", "stitle"],index_col=False)

#Go back to home directory /Pipeline_Project_Mansi_Patel 
os.chdir('..')
os.chdir('..')

#Append the blast output to log file 
with open('PipelineProject.log','a') as f:
    f.write(str(df.to_string(index=False)))  
    f.close()








