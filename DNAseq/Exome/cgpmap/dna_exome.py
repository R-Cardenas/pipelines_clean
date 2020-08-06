# This will be run by the run_pipeline.py if the master_user_config.yaml
# select dna-exome. This script will select merge or nomerge pipeline
# and the corrrect config files for the corresponding genome files
import os
# RPC 080720
# aft19qdu@uea.ac.uk
from run_pipeline import data, home_dir

# Select merge or no merge nextflow
if data['merged_lanes'] == 'no':
    cgpmap_nf = "nextflow run dna-exome-nomerge.nf"
    cmd = "cp dna-exome-nomerge.nf " + home_dir
    os.system(cmd)
elif data['merged_lanes'] == 'yes':
    cgpmap_nf = "nextflow run dna-exome-merge.nf"
    cmd = "cp dna-exome-merge.nf " + home_dir
    os.system(cmd)
else:
    print('dna_exome.py - line10')
    raise SyntaxError("dna_exome.py: Incorrect 'merged_lanes' input. Please revise")

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
os.system(cgpmap_bash2)