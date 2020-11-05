/*
 * create a channel for fastq pairss
 */

// Input Reads
params.read1 = "$baseDir/input/*{1,2}.fq.gz"


read1_ch = Channel .fromFilePairs( params.read1 )



println """\
	\
	\
	\
         ==================================
         E X O M E - N F   P I P E L I N E
         N O  M E R G E

         Mapping and Bam Processing
         v0.1
         ===================================



         """
         .stripIndent()


process cgpMAP {


	storeDir "$baseDir/output/cgpMAP/${read1.simpleName}"
  input:
	val read1 from read5_ch
  output:
  file "*.bam" into cgp_ch
  script:
  """

  name=\$(echo '${read2}' | sed -e 's/.*[/]//' -e 's/_.*//')

  ds-cgpmap.pl  \
  -outdir $baseDir/output/cgpMAP/${read1.simpleName} \
  -r $cgpmap_genome \
  -i $cgpmap_index \
  -s \$name \
  -t 5 \
	-g ${read1}.yaml \
  ${read1} ${read2}

	mv $baseDir/output/cgpMAP/${read1.simpleName}/*.bam \
	$baseDir/output/cgpMAP/${read1.simpleName}/${read1.simpleName}.bam

	echo 'fq1: ${read1} fq2: ${read2} bam_name: ${read1.simpleName}' >> $baseDir/${projectname}_cgpmap_samples.log
	ls -l  ${read1} >> $baseDir/logs/symbolic_test_fastq.log
	ls -l  ${read2} >> $baseDir/logs/symbolic_test_fastq.log
  """
}
