#Import the necessary libraries 
import os 
import glob
import shutil

#Make a folder to store the test subsetted data 
os.mkdir('test_data')
original = os.getcwd() 
#Go into the transcriptome_data folder 
os.chdir('transcriptome_data/')

#Subset the Fastq-Dump for test files 
fastq_pattern = "*.fastq"
fastq_files = glob.glob(fastq_pattern)
fastq_files.sort()

#Subset the data 
for files in fastq_files:
  dump = 'head -n 40000 '+ files + ' > sample_' + files
  os.system(dump)
  
#Move the files to a different folder 
new_original = 'transcriptome_data'
current = os.path.join(original,new_original)
test = 'test_data'
new = os.path.join(original,test)
for filename in os.listdir(current):
    if filename.startswith("sample"):
        source_file = os.path.join(current, filename)
        destination_file = os.path.join(new, filename)
        shutil.move(source_file, destination_file)
  
