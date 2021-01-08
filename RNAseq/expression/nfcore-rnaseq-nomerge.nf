/*
 * create a channel for bam files produced by cgpmap_processing pipeline
 */
params.outputdir = "/gpfs/afm/cg_pipelines/Pipelines/Williams_RNASeq_processing"
params.fq = "$baseDir/input/*{fq,fastq}.gz"
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

// needs to be tested
process main_nf{
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
	executor 'slurm'
	memory '55 GB'
	input:
	file fastq from fq_ch.collect()
	script:
	"""
	nextflow run nf-core/rnaseq -resume -profile singularity \
	-c RNAseq/expression/UEA.config \
	--reads './*{1,2}.{fastq,fq}.gz' \
	--genome GRCh37 \
	--outdir '$baseDir/test' \
	--max_memory '55.GB' \
	--saveAlignedIntermediates \
	--saveTrimmed \
	--email_on_fail aft19qdu@uea.ac.uk
	"""
}
