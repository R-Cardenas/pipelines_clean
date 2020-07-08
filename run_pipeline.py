# DESCRIPTION
# RPC 080720
# aft19qdu@uea.ac.uk

import yaml
import os

## Create a bash script that will run all the pipelines selected by the config

os.system("rm -fr run_selected_pipeline.sh")
os.system("echo '#!/bin/bash' >> run_selected_pipeline.sh") ## should this be bsub? or lsf?


## Open the user master yaml file
with open('master_user_config.yaml') as f:

    data = yaml.load(f,Loader=yaml.FullLoader)

## Extract data relevant to pipeline one and write the bash script.

# import python scripts based on samples input from YAML

if data['samples'] == 'dna-exome':
    from DNAseq.Exome.cgpmap.dna_exome import *
elif data['samples'] == 'dna-wgs':
    from DNAseq.WGS.dna_wgs import *
elif data['samples'] == 'rna-seq':
     print('needs py file to be added')
else:
    raise SyntaxError('incorrect "samples" values input in master_user_config.yaml')



### We will run whats in the run_selected_pipeline.sh once it is filled.