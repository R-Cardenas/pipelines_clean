# This will be run by the run_pipeline.py if the master_user_config.yaml
# select dna-exome. This script will select merge or nomerge pipeline
# and the corrrect config files for the corresponding genome files
import os
# RPC 080720
# aft19qdu@uea.ac.uk
from run_pipeline import data

# Select merge or no merge nextflow
if data['merged_lanes'] == 'no':
    cgpmap_nf = "nextflow run DNAseq/Exome/cgpmap/dna-exome-nomerge.nf"
elif data['merged_lanes'] == 'yes':
    cgpmap_nf = "nextflow run DNAseq/Exome/cgpmap/dna-exome-merge.nf"
else:
    print('dna_exome.py - line16')
    raise SyntaxError("dna_exome.py: Incorrect 'merged_lanes' input. Please revise")

# Select config files for hg19 or hg38
if data['genome_assembly'] == 'hg38':
    cgpmap_config = "-C DNAseq/Exome/cgpmap/dna-exome_hg38.config"
elif data['genome_assembly'] == 'hg19':
    cgpmap_config = "-C DNAseq/Exome/cgpmap/dna-exome_hg19.config"
else:
    print('dna_exome.py - line24')
    raise SyntaxError("dna_exome.py: Incorrect 'genome_assembly' input. Please revise")


### Concat nextflow run and config into bash script in repo home Dir

cgpmap_bash = cgpmap_nf + ' ' + cgpmap_config

cgpmap_bash2 = "echo '" + cgpmap_bash + "' >> run_selected_pipeline.sh"
print(cgpmap_bash2)
os.system(cgpmap_bash2)