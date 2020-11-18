# cgpMAP pipeline

This pipeline utilises the sanger cgpMAP pipeline (https://github.com/cancerit/dockstore-cgpmap), in addition to some QC tools. This manual describes how this works and how to use it on its own without the master-YAML.

There are two versions of this pipeline, 'dna-exome-merge' and 'dna-exome-no-merge' - the merging factor pertains to samples run on different lanes that require to be merged once the BAMs have been formed. This results in an extra 'process' step in the 'merge' pipeline that will merge BAMs before the QC steps, but all other processes are identical in both pipelines. Below the different inputs/processes are described, only the most complicated



## Input

Below the params.read1 variable specifies the fastq directory which are read as pairs of files.
The files are piped into two channels where the fastq files are processed separately. The *$baseDir* is specifies the working directory that nextflow is run from - similar to the $PWD linux command line utility.


```go
params.read1 = "$baseDir/input/*{1,2}.fq.gz"

read1_ch = Channel .fromFilePairs( params.read1 )

read1_ch.into { read2_ch; read3_ch }

```

## Trim_galore

This is the first process. The files are copied into the work dir as usually the files are kept on different drives and sym/hardlinks tend to fail. The destination of the files is set by the *storeDir* variable.

```go
stageInMode = 'copy' // trim_galore doesnt like sym/hardlinks.
storeDir "$baseDir/output/cgpMAP/trim_galore"
```

The script in this process uses trim-galore and then renames the files replacing '_val' extension with '.trim.'. This is done as nextflow allows the automatic removal of file extensions using the '.simpleName' or '.baseName' function, this helps prevent the accumulation of unnecessary file extensions. The copied files (ie not the trimmed files) are deleted to preserve space.

```bash
trim_galore --paired --fastqc --illumina \
--basename ${reads[0].simpleName} \
${reads[0]} ${reads[1]}

# rename val with trim
mv ${reads[0].simpleName}_val_1.fq.gz ${reads[0].simpleName}_1.trim.fq.gz
mv ${reads[0].simpleName}_val_2.fq.gz ${reads[0].simpleName}_2.trim.fq.gz

# delete copied files
rm -fr ${reads[0]} # remove the copied files to prevent memory loss
rm -fr ${reads[1]}
```

## FastQ
This process uses the fastq headers to extract platform_unit (PU), library preparation (LB) and lane information to form a YAML file to inform cgpMAP to annontate the BAM correctly.

First extract the header using fqtools - but only using a header containing no 'N' in the sequencing primer field. This was done to ensure the sequencing primer for read1 and read2 are indentical.

The second part is an if statement - if the txt file produced is empty then repeat the extraction, but excluding the regex grep. Due to differences in fastq header this may occur, however for the majority this will work.

(For more information about fastq headers see: https://angus.readthedocs.io/en/2017/Read_group_info.html)

```bash
# Read1
# Extract header
fqtools -d header ${read1} | grep ":[C,A,T,G]*[+][C,A,T,G]" | head -1 > ${read1.simpleName}.txt

### Counts lines in file1 and will repeat if it is empty
words=`wc -l ${read1.simpleName}.txt  | awk '{print \$1}'`;
if [ \$words -eq 0 ]
then
fqtools -d header ${read1} | head -1 > ${read1.simpleName}.txt
else
echo 'alls good!'
fi
```

Next a python script is used to form a YAML file in the correct format for cgpMAP (described here: https://github.com/cancerit/PCAP-core/wiki/File-Formats-groupinfo.yaml)

This script assumes the platform unit is Illumina, and will add the sample name obtained from the fastq filename (E.g. S0102-P-TUMOR-L01-R1.fastq.gz sample name is S0102-P).

```bash
python $baseDir/nextflow_pipelines/bin/python/fastq2config_cgpmap.py \
--fq1 ${read1.simpleName}.txt --fq2 ${read2.simpleName}.txt \
--n1 ${read1} --n2 ${read2} --o ${read1}.yaml

cp *WARNING* $baseDir/logs 2>/dev/null || :
```

An example of the YAML file produced is shown below:

```yaml
SM: S0102-P
READGRPS:
  1.fq:
    PL: ILLUMINA
    LB: HGYYCBBXX_NTTCCT
    PU: HGYYCBBXX.6
  2.fq:
    PL: ILLUMINA
    LB: HGYYCBBXX_NTTCCT
    PU: HGYYCBBXX.6
```
