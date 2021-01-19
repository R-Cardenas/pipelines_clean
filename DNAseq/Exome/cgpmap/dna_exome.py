# This will be run by the run_pipeline.py if the master_user_config.yaml
# select dna-exome. This script will select merge or nomerge pipeline
# and the corrrect config files for the corresponding genome files
import os
# RPC 080720
# aft19qdu@uea.ac.uk
import glob
from bin.python.data_yaml import data
import subprocess

################
## check BAIT ##
################

print("Checking for interval files...")
interval_files = glob.glob('input/bait/*.bed')
if interval_files == 0:
    raise SyntaxError("No bed files found in input/bait .. Exome-seq analysis requires these file.")
elif interval_files == 1:
    raise SyntaxError("Only 1 bed file found in input/bait .. Both target and bait files required")
else:
    print("Found the following bed files in input/bait folder:")
    print(interval_files)

#################
## MERGE LANES ##
#################

# Select merge or no merge nextflow
if data['merged_lanes'] == 'no':
    cgpmap_nf = "nextflow run dna-exome-nomerge.nf"
    cmd = "cp DNAseq/Exome/cgpmap/dna-exome-nomerge.nf ."
    p = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, universal_newlines=True, shell=True)
    print(p)

elif data['merged_lanes'] == 'yes':
    cgpmap_nf = "nextflow run dna-exome-merge.nf"
    cmd = "cp DNAseq/Exome/cgpmap/dna-exome-merge_v2.nf ."
    p = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, universal_newlines=True, shell=True)
    print(p)
else:
    raise SyntaxError("dna_exome.py: Incorrect 'merged_lanes' input. Please revise")

###################
## Genome Config ##
###################

# Select config files for hg19 or hg38
if data['genome_assembly'] == 'hg38':
    cgpmap_config = "-c DNAseq/Exome/cgpmap/cgpmap_hg38.config"
elif data['genome_assembly'] == 'hg19':
    cgpmap_config = "-c DNAseq/Exome/cgpmap/cgpmap_hg19.config"
else:
    print('dna_exome.py - line24')
    raise SyntaxError("dna_exome.py: Incorrect 'genome_assembly' input. Please revise")


### Concat nextflow run and config into bash script in repo home Dir

cgpmap_bash = cgpmap_nf + ' ' + cgpmap_config
cgpmap_bash2 = "echo '" + cgpmap_bash + "' >> run_selected_pipeline.sh"
print(cgpmap_bash2)

p = subprocess.run(cgpmap_bash2, check=True, stdout=subprocess.PIPE, universal_newlines=True, shell=True)
print(p)
