# cgpMAP pipeline

This pipeline utilises the sanger cgpMAP pipeline (https://github.com/cancerit/dockstore-cgpmap), in addition to some QC tools. This manual describes how this works and how to use it on its own without the master-YAML.

There are two versions of this pipeline, 'dna-exome-merge' and 'dna-exome-no-merge' - the merging factor pertains to samples run on different lanes that require to be merged once the BAMs have been formed. This results in an extra 'process' step in the 'merge' pipeline that will merge BAMs before the QC steps, but all other processes are identical in both pipelines. Below the different inputs/processes are described, only the most complicated



## Input

Below the params.read1 variable specifies the fastq directory which are read as pairs of files.
The files are piped into two channels where the fastq files are processed separately. The *$baseDir* is specifies the working directory that nextflow is run from - similar to the $PWD linux command line utility.



* specify path
```shell
params.read1 = "$baseDir/input/*{1,2}.fq.gz"

read1_ch = Channel .fromFilePairs( params.read1 )

read1_ch.into { read2_ch; read3_ch }

```

## Trim_galore

This is the first process. The files are copied into the work dir as usually the files are kept on different drives and sym/hardlinks tend to fail. The destination of the files is set by the *storeDir* variable.

```shell
stageInMode = 'copy' // trim_galore doesnt like sym/hardlinks.
storeDir "$baseDir/output/cgpMAP/trim_galore"
```

The script in this process uses trim-galore and then renames the files replacing '_val' extension with '.trim.'. This is done as nextflow allows the automatic removal of file extensions using the '.simpleName' or '.baseName' function, this helps prevent the accumulation of unnecessary file extensions. The copied files (ie not the trimmed files) are deleted to preserve space.

```shell
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
