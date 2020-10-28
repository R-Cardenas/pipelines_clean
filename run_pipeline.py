# DESCRIPTION
# RPC 080720
# aft19qdu@uea.ac.uk
import os
import sys
import subprocess
sys.path.append("python-packages/pyyaml")
import yaml
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper

## Open the user master yaml file
with open('master_user_config.yaml') as f:

    data = yaml.load(f,Loader=Loader)

##################
## HOUSEKEEPING ##
##################


# Remove historical bash scripts
rm_existing = "rm -fr run_selected_pipeline.sh"
new_sh = "echo '#!/bin/bash' >> run_selected_pipeline.sh"
os.system(rm_existing) ## should this be bsub? or slurm?
os.system(new_sh)

# Remove historical variables from config files
rm_projectname = f"""for f in $(find . -name '*config'); do sed -i '/env.projectname/d' $f; done"""
rm_build = f"""for f in $(find . -name '*config'); do sed -i '/env.build/d' $f; done"""
os.system(rm_projectname)
os.system(rm_build)

#####################
### PROJECT NAME ####
#####################


projectname = data['projectname'].lower() # extract the projects name
projectname2 = 'env.projectname = "' + projectname + '"'
print("Project name is: " + projectname2)
## Below line requires python >3.6
add_projectname = f"""for f in $(find . -name '*config'); do echo 'env.projectname={projectname2}' >> $f; done""" # add projectname to all config files
os.system(add_projectname)


genome_assembly = data['genome_assembly'].lower() # extract the projects name
genome_assembly2 = 'env.projectname = "' + genome_assembly + '"'
print("Project name is: " + genome_assembly2)
## Below line requires python >3.6
add_build = f"""for f in $(find . -name '*config'); do echo 'env.build={genome_assembly2}' >> $f; done""" # add projectname to all config files
os.system(add_build)



#####################
### EXOME VS WGS ####
#####################

# The next chunk identifies data type and then imports the python script specific for the datatype e.g. exome
if data['samples'].lower() == 'dna-exome':
    from DNAseq.Exome.cgpmap.dna_exome import *
elif data['samples'].lower() == 'dna-wgs':
    from DNAseq.WGS.dna_wgs import *
elif data['samples'].lower() == 'rna-seq':
     print('needs py file to be added')
else:
    raise SyntaxError('incorrect "samples" values input in master_user_config.yaml')

###########################
### SOMATIC VS GERMLINE ###
###########################

if data['samples'].lower() == 'dna-exome' and data['variant'] == 'germline':
    variant_nf = "nextflow run freebayes_individual.nf & \"" \
                 "nextflow run haplotypecaller_individual.nf"
    cmd1 = "cp DNAseq/Exome/germline/freebayes/freebayes_individual.nf ."
    cmd2 = "cp DNAseq/Exome/germline/gatk/haplotypecaller_individual.nf ."
    os.system(cmd1)
    os.system(cmd2)
elif data['samples'].lower() == 'dna-exome' and data['variant'] == 'somatic':
    variant_nf = "nextflow run cgpwxs_v0.1.nf & \"" \
                 "nextflow run mutect2_individual.nf"
    cmd1 = "cp DNAseq/Exome/somatic/cgpwxs_v0.1.nf ."
    cmd2 = "cp DNAseq/Exome/somatic/mutect2_individual.nf ."
    os.system(cmd1)
    os.system(cmd2)
elif data['samples'].lower() == 'rna-seq':
    ### NEEDS inputs..
else:
    raise SyntaxError('incorrect "variant" values input in master_user_config.yaml')


### Concat nextflow run and config into bash script in repo home Dir
cgpmap_bash2 = "echo '" + variant_nf + "' >> run_selected_pipeline.sh"
print(cgpmap_bash2)
os.system(cgpmap_bash2)


###################################
### SOMATIC ONLY - SAMPLE NAMES ###
##################################

if data['matched_identical'].lower() == 'yes':
    print('Tumor/Normal samples haved identical names. Extracting from bam files name...')
    from DNAseq.Exome.somatic.cgpwxs.tumor_normal_var import *
elif data['matched_identical'].lower() == 'no':
    print('Tumor/Normal samples do NOT identical names. Extracting from YAML file')
    from DNAseq.Exome.somatic.cgpwxs.samples_parse import *
else:
    raise SyntaxError('ERROR with sample names. Please check YAML #10 and #11')

############
### RUN ###
###########

### We will run whats in the run_selected_pipeline.sh once it is filled.

run_master = ["bash","run_selected_pipeline.sh"]
output = subprocess.Popen(run_master, stdout=subprocess.PIPE ).communicate()[0]
print(output)

## Remove the nextflow files that were copied to keep clean
cp_nf = 'rm -fr *.nf'
os.system(cp_nf)
