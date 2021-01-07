# This will be run by the run_pipeline.py if the master_user_config.yaml
# select rna-seq. This script will select merge or nomerge pipeline
# and will select the correct genome
import os
# RPC 080720
# aft19qdu@uea.ac.uk
from run_pipeline import data
import glob


#################
## MERGE LANES ##
#################

# Select merge or no merge nextflow
if data['merged_lanes'] == 'no':
    variant_nf = "nextflow run dna-exome-nomerge.nf -c RNAseq/expression/nextflow.config"
    cmd = "cp RNAseq/expression/nfcore-rnaseq-nomerge.nf ."
    os.system(cmd)
elif data['merged_lanes'] == 'yes':
    variant_nf = "nextflow run dna-exome-merge.nf -c RNAseq/expression/nextflow.config"
    cmd = "cp RNAseq/expression/nfcore-rnaseq-merge.nf ."
    os.system(cmd)
else:
    print('dna_exome.py - line10')
    raise SyntaxError("dna_exome.py: Incorrect 'merged_lanes' input. Please revise")

########################
## Output Dir         ##
########################

output_dir = data['output_dir']
if output_dir == "":
    raise SyntaxError("YAML output_dir value is empty")
else:
    output_dir2 = "--outdir " + output_dir + " \\"
    replace_string = f"sed -i 's/--outdir .*/{output_dir2}/g' nfcore-rnaseq*.nf'" # replace rna-seq nf with YAML input
    os.system(replace_string)

###################
## Genome Config ##
###################

# Select config files for hg19 or hg38
if data['genome_assembly'] == 'hg38':
    rnaseq_genome = "--genome GRCh38 \\"
    replace_string = f"sed -i 's/--genome .*/{rnaseq_genome}/g' nfcore-rnaseq*.nf'"
    os.system(replace_string)

elif data['genome_assembly'] == 'hg19':
    rnaseq_genome = "--genome GRCh37 \\"
    replace_string = f"sed -i 's/--genome .*/{rnaseq_genome}/g' nfcore-rnaseq*.nf'"
    os.system(replace_string)
else:
    print('dna_exome.py - line24')
    raise SyntaxError("dna_exome.py: Incorrect 'genome_assembly' input. Please revise")
