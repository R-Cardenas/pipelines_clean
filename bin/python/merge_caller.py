# RPC 210920
# This will merge sample bam files from different variant callers - SNPs only!
# !/usr/bin/env python
import os
import argparse
import subprocess

files = 'sample1-fam1-GATK.vcf.gz sample1-fam1-freebayes.vcf.gz '


# Split the input from nextflow
parser = argparse.ArgumentParser()
files = parser.add_argument('--bam', required=True)
files2 = files.split(" ")

# Extract samples names from each input and create unique list
bam_samples = list()
for f in files2:
    list = f.split("-")
    sample = list[0]
    bam_samples.append(sample)

unique = set(bam_samples)
print(unique)

# For loop submits each sample to be lane merged
for i in unique:
    list = i.split("-")
    sample = list[0]

    # count how many files to merge
    sample_wild = sample + '*'
    cmd_count = 'ls -l ' + sample_wild + ' | wc -l'
    count = subprocess.run([cmd_count], stdout=subprocess.PIPE, shell = True)
    count_number = str(int(count.stdout))

    script2 = 'bcftools isec -c indels -n +' + count_number + ' -o ' + sample + '.caller.snps.merged.vcf -p ' + sample + ' ' + sample_wild
    print(script2)
    os.system(script2)


