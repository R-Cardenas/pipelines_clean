/*
 * create a channel for bam files produced by cgpmap_processing pipeline
 */
params.outputdir = ""
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
	scratch true
	input:
	file fastq from fq_ch.collect()
	output:
	val "*merged*" into main_ch
	script:
	"""
	module add python/anaconda/2019.10/3.7

	python $baseDir/bin/python/merge_fastq.py --fq '${fastq}'

	# remove non merged fastq files to conserve space
	find . -type f ! -name '*-merged-*' -delete
	"""
}

process main_nf{
	storeDir "$baseDir/test"
	input:
	val fastq from main_ch
	output:
	file "${fastq}.txt"
	script:
	"""
	nextflow run nf-core/rnaseq \
	-c /gpfs/afm/cg_pipelines/Pipelines/singularity/nextflow_configs/rna-seq_v2.6.1.config \
	-profile singularity \
	--reads '*{1,2}.{fastq,fq}.gz' \
	--genome GRCh37 \
	--outdir '${outputdir}' \
	--max_memory '60.GB' \
	--email_on_fail aft19qdu@uea.ac.uk
	"""
}
