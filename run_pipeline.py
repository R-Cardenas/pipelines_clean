# DESCRIPTION
# RPC 080720
# aft19qdu@uea.ac.uk
import os
import sys
import subprocess
import glob
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

## Remove the nextflow files that may have been copied to baseDir
cp_nf = 'rm -fr *.nf'
os.system(cp_nf)

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
    from RNAseq.rnaseq import *
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
    variant_nf = "nextflow run nfcore-rnaseq-merge.nf"
    cmd1 = "cp RNAseq/nfcore-rnaseq-merge.nf ."
    os.system(cmd1)
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
    pass # is not essential for germline / rna-seq

##########################
### Input & output Dir ### and out dir? - not for rna-seq
########################## is different

# Input fastqs
fastq_input = data['fastq_dir']
remove_inputDir = f"""find . -name "*.nf" -exec sed -i '/params.fq*/d' {{}} \;"""
print(remove_inputDir)
os.system(remove_inputDir)

# Output dir (not for rna-seq - see rnaseq.py)
output_base = data['output_dir']
remove_outputDir = f"""find . -name "*.nf" -exec sed -i '/params.outputdir*/d' {{}} \;"""
print(remove_outputDir)
os.system(remove_outputDir)

# Add input and output dirs back in with YAML input:
input1 = f"""params.fq = \"{fastq_input}\" \n"""
output1 =  f"""params.outputdir = \"{output_base}\" \n"""

nf_files = glob.glob("*.nf")
for f in nf_files:
    file_object = open(f, 'a')
    file_object.write(input1)
    file_object.write(output1)
    file_object.close()

############
### RUN ###
###########

### We will run whats in the run_selected_pipeline.sh once it is filled.

run_master = ["bash","run_selected_pipeline.sh"]
output = subprocess.Popen(run_master, stdout=subprocess.PIPE ).communicate()[0]
print(output)

## Remove the nextflow files that were copied to keep clean
#cp_nf = 'rm -fr *.nf'
#os.system(cp_nf)

exit()
