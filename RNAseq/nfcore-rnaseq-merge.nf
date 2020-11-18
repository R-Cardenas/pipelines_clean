/*
 * create a channel for bam files produced by cgpmap_processing pipeline
 */
params.outputdir = "/gpfs/afm/cg_pipelines/Pipelines/Williams_RNASeq_processing"
params.fq = "/gpfs/afm/cg_pipelines/Pipelines/Williams_RNASeq/*{fq,fastq}.gz"
fq_ch = Channel .fromPath( params.fq )

println """\
	\
	\
	\
         ==================================
         NF - core: RNAseq
				 v0.1
         ===================================



         """
         .stripIndent()

process merge_lanes{
	stageInMode = "copy" // trim_galore doesnt like sym/hardlinks.
	storeDir "$baseDir/test_merge"
	input:
	file fastq from fq_ch.collect()
	output:
	file "*merged*" into main_ch
	script:
	"""
	module add python/anaconda/2019.10/3.7

	python $baseDir/bin/python/merge_fastq.py --fq '${fastq}'

	# remove non merged fastq files to conserve space
	find . -type f ! -name '*-merged-*' -delete
	"""
}

process main_nf{
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
	executor 'slurm'
	memory '55 GB'
	storeDir "$baseDir/main_outputs"
	input:
	file fastq from main_ch.collect()
	output:
	file '*'
	script:
	"""
	nextflow run nf-core/rnaseq -resume -profile singularity \
	-c /gpfs/afm/cg_pipelines/Pipelines/singularity/nextflow_configs/UEA.config \
	--reads '*{1,2}.{fastq,fq}.gz' \
	--genome GRCh37 \
	--outdir 'test' \
	--max_memory '55.GB' \
	--saveAlignedIntermediates \
	--saveTrimmed \
	--email_on_fail aft19qdu@uea.ac.uk
	"""
}
