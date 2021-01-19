# RPC 170920
# This will merge sample bam files from different lanes into one using samtools
# !/usr/bin/env python
import os
import argparse
import time
import collections
import pandas as pd

parser = argparse.ArgumentParser()
files = parser.add_argument('--bam', required=True)
#files = "S1208-S12.bam S1208-S13.bam S1205-12.bam" this is test input
files2 = files.split(" ")

bam_samples = list()
for f in files2:
    list = f.split("-")
    sample = list[0]
    bam_samples.append(sample)

# Convert dict to df and count number of samples in each
count_dict = collections.Counter(bam_samples)
count_df = pd.DataFrame(count_dict.items())
count_df.columns = ['samples','count']

# Print samples which are no sequenced on multiple lanes
single_samples =  count_df.query('count == 1')
print('\n \n ')
print('The following samples are not sequenced on multiple lanes: \n \n ')
print(single_samples)

# Slice df for samples across multiple lanes to merge
merge_samples =  count_df.query('count > 1')

###############################
## Send samtools job: SBATCH ##
###############################


# For loop submits each sample to be lane merged
for i in merge_samples['samples']:
    sample = str(i)
    script = '"samtools merge ' + sample + '_merged.bam ' + sample + '*bam -@ 5"'
    script2 = '#SBATCH -p compute-16-24 --mem 10G --job-name=bam_merge_' + sample + ' -e %j.err ' + script
    p = subprocess.run(script2, check=True, stdout=subprocess.PIPE, universal_newlines=True, shell=True)
    print(p)
    time.sleep(180)


###############################
## Check the file exists     ##
###############################

# Checks that the file exists (ie samtools output)
def loop():
    for i in merge_samples['samples']:
        sample = str(i)
        file = sample + "_merged.bam"
        if os.path.isfile(file):
            print ("File exist")
        else:
            print('file does not exist yet')
            time.sleep(360)
            loop()


###############################
## Check size of file        ##
###############################

# Check that the file is not still increasing in size (ie being processed)
# Loop will end once the file stops increasing in size. Ie bam merger has finished.

def size():
    for i in merge_samples['samples']:
        sample = str(i)
        file = sample + "_merged.bam"
        size1 = os.stat(file).st_size
        time.sleep(180)
        size2 = os.stat(file).st_size
        if size1 == size2: # compares same files two minutes apart. If they are the same size it is presumed processing has finished.
            print("files are the same size")
            os.system("echo 'size1, size2 " + size1 + size2 + file + "' >> size.log")
        else:
            size()

loop()
size()
time.sleep(500)
