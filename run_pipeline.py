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

home_dir = os.getcwd()

## Create a bash script that will run all the pipelines selected by the config
rm_existing = "rm -fr run_selected_pipeline.sh"
new_sh = "echo '#!/bin/bash' >> run_selected_pipeline.sh"

os.system(rm_existing) ## should this be bsub? or lsf?
os.system(new_sh)

## Open the user master yaml file
with open('master_user_config.yaml') as f:

    data = yaml.load(f,Loader=Loader)

## Extract data relevant to pipeline one and write the bash script.


projectname = data['projectname'] # extract the projects name
projectname2 = 'env.projectname = "' + projectname + '"'
print("Project name is: " + projectname2)
add_projectname = f"""for f in $(find . -name '*config'); do echo 'env.projectname={projectname2}' >> $f; done""" # add projectname to all config files
os.system(add_projectname)

# The next chunk identifies data type and then imports the python script specific for the datatype e.g. exome

if data['samples'] == 'dna-exome':
    from DNAseq.Exome.cgpmap.dna_exome import *
elif data['samples'] == 'dna-wgs':
    from DNAseq.WGS.dna_wgs import *
elif data['samples'] == 'rna-seq':
     print('needs py file to be added')
else:
    raise SyntaxError('incorrect "samples" values input in master_user_config.yaml')

## Remember to includeconfig user.config..

### We will run whats in the run_selected_pipeline.sh once it is filled.

run_master = ["bash","run_selected_pipeline.sh"]
output = subprocess.Popen(run_master, stdout=subprocess.PIPE ).communicate()[0]
print(output)

## Remove the nextflow files that were copied to keep clean
cp_nf = 'rm -fr *.nf'
os.system(cp_nf)
