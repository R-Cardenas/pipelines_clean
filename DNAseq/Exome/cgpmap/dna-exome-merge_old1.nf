/*
 * create a channel for fastq pairss
 */

// Input Reads

params.

read1_ch = Channel .fromFilePairs( params.fq)
read1_ch.into { read2_ch; read3_ch }

params.csv = "$baseDir/bin/williams_batch2_info.csv"
csv_ch = Channel .fromPath( params.csv )

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


myLongCmdline = "git clone https://github.com/R-Cardenas/nextflow_pipelines.git"
result = myLongCmdline.execute().text


process trim_galore{
	stageInMode = 'copy' // trim_galore doesnt like sym/hardlinks.
  storeDir "$baseDir/output/cgpMAP/trim_galore"
	input:
	tuple val(read2), file(reads) from read2_ch
	output:
	file "${reads[0].simpleName}_1.trim.fq.gz" into (read5_ch, read7_ch)
	file "${reads[0].simpleName}_2.trim.fq.gz" into (read10_ch, read12_ch)
	file("*.html") optional true
	script: {
	"""
	mkdir -p $baseDir/logs

	trim_galore --paired --fastqc --illumina \
	--basename ${reads[0].simpleName} \
	${reads[0]} ${reads[1]}

	# rename val with trim
	mv ${reads[0].simpleName}_val_1.fq.gz ${reads[0].simpleName}_1.trim.fq.gz
	mv ${reads[0].simpleName}_val_2.fq.gz ${reads[0].simpleName}_2.trim.fq.gz

	# delete copied files
	rm -fr ${reads[0]} # remove the copied files to prevent memory loss
	rm -fr ${reads[1]}
	"""
}
}


process fqtools{


  storeDir "$baseDir/output/cgpMAP/trim_galore"
	input:
	file read1 from read7_ch
	file read2 from read12_ch
	output:
	file "${read1}.yaml" into yaml_ch
	file("fqtools_WARNING_?.txt") optional true
	script:
	"""
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

	# Read2
	# Extract header

  fqtools -d header ${read2} | grep ":[C,A,T,G]*[+][C,A,T,G]" | head -1 > ${read2.simpleName}.txt

	### Counts lines in file2 and will repeat if it is empty
	words=`wc -l ${read2.simpleName}.txt  | awk '{print \$1}'`;
	if [ \$words -eq 0 ]
	then
	fqtools -d header ${read2} | head -1 > ${read2.simpleName}.txt
	else
	echo 'alls good!'
	fi


	python $baseDir/nextflow_pipelines/bin/python/fastq2config_cgpmap.py \
	--fq1 ${read1.simpleName}.txt --fq2 ${read2.simpleName}.txt \
	--n1 ${read1} --n2 ${read2} --o ${read1}.yaml

	cp *WARNING* $baseDir/logs 2>/dev/null || :
	"""
}

process cgpMAP {


	storeDir "$baseDir/output/cgpMAP/${read1.simpleName}"
  input:
	val read1 from read5_ch
	val read2 from read10_ch
	val yaml from yaml_ch.collect()
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


process sam_sort {


  storeDir "$baseDir/output/BAM/sorted"
  input:
  file bam from cgp_ch
  output:
  file "${bam.simpleName}.sorted.bam" into bam_merge_ch
  script:
  """
	# create the tmp file as picard can create large temp files

	mkdir -p tmp
  picard SortSam I=${bam} O=${bam.simpleName}.sorted.bam SORT_ORDER=coordinate TMP_DIR=tmp
	TMP_DIR=tmp
	rm -fr tmp
  """
}

// dont forget to add singularity with python3 installed
process bam_merge {
  storeDir "$baseDir/output/BAM/sorted"
  input:
  file bam from bam_merge_ch.collect()
  output:
	file "*_merged.bam" into dup_ch
  script:
  """
	# Create log file
	echo '#RPC bam_merge logs from dna-exome-merge' >> $baseDir/logs/bam_merge_log.txt
	echo '#All sample names on one line are expected to be the same. Otherwise a bug has occured' >> $baseDir/logs/bam_merge_log.txt
	echo 'wild_card samples' >> $baseDir/logs/bam_merge_log.txt

	# Samtools MERGE
	# This will extract the sample name and create a samtools wild card
	for f in \$(ls *.bam | sed -e 's/.*[/]//' -e 's/_.*//' | sort | uniq)
	do
	samtools merge \${f}_merged.bam \${f}*

	printf "`echo \${f}` `echo $(ls \${f}*)`\\n" >> $baseDir/logs/bam_merge_log.txt

	done

  """
}

process picard_pcr_removal {


	storeDir "$baseDir/output/BAM/merged_lanes"
  input:
  file bam from dup_ch.flatten()
  output:
  file "${bam.simpleName}.rmd.bam" into (index1_ch, hs_ch, bam10_ch, bam11_ch, bam12_ch)
	file "${bam.simpleName}.log"
  script:
  """
	mkdir -p tmp
  picard MarkDuplicates I=${bam} O=${bam.simpleName}.rmd.bam M=${bam.simpleName}.log TMP_DIR=tmp
	TMP_DIR=tmp
	rm -fr tmp
  """
}

process bam_index {


	storeDir "$baseDir/output/BAM/merged_lanes"
  input:
  file bam from index1_ch
  output:
  file "${bam}.bai" into (index_3ch, index_4ch)

  script:
  """
	mkdir -p tmp
  picard BuildBamIndex \
	I=${bam} \
	O=${bam}.bai \
	TMP_DIR=tmp
	rm -fr tmp
  """
}


process hybrid_stats {


	storeDir "$baseDir/output/BAM/hybrid_stats"
  input:
  file bam from hs_ch
  output:
  file "${bam.simpleName}_hs_metrics.txt"
  script:
  """

	# create interval files from BED supplied
	picard BedToIntervalList \
	I=${bait_interval} \
	O=${bait_interval}.interval \
	SD=$genome_fasta

	picard BedToIntervalList \
	I=${target_interval} \
	O=${target_interval}.interval \
	SD=$genome_fasta


	# perform hybrid stats
	mkdir -p tmp
  picard CollectHsMetrics I=${bam} O=${bam.simpleName}_hs_metrics.txt \
  R=$genome_fasta \
  BAIT_INTERVALS=${bait_interval}.interval \
  TARGET_INTERVALS=${target_interval}.interval \
	TMP_DIR=tmp
	rm -fr tmp
  """
}

process alignment_stats{


	storeDir "$baseDir/output/BAM/alignment_stats"
	input:
	file bam from bam10_ch
	output:
	file "${bam.simpleName}_align_stats.txt" into verify_ch
	script:
	"""
	mkdir -p tmp
	picard CollectAlignmentSummaryMetrics \
  R=$genome_fasta \
	I=${bam} \
	O=${bam.simpleName}_align_stats.txt \
	TMP_DIR=tmp
	rm -fr tmp
	"""
}

process verifybamid{

	stageInMode = 'copy' // somalier doesnt like sym/hardlinks.
	storeDir "$baseDir/output/BAM/verifyBamID"
	input:
	file bam from bam11_ch
	file idx from index_3ch.collect()
	output:
	file "*.depthRG"
	file "*.depthSM"
	file "*.log"
	file "*.selfRG"
	file "*.selfSM"
	script:
	"""
	verifyBamID --vcf $verifybamid \
	--bam ${bam} \
	--out ${bam.simpleName} \
	--maxDepth 1000 \
	--precise \
	--verbose

	rm -fr *.bam
	"""
}



process somalier{

	stageInMode = 'copy' // somalier doesnt like sym/hardlinks.
  storeDir "$baseDir/output/BAM/somalier"
  input:
  file bam from bam12_ch.collect()
	file idx from index_4ch.collect()
  output:
	file "*.html"
  script:
  """
	mkdir -p bin
	wget -P bin/ https://github.com/brentp/somalier/files/3412456/sites.hg38.vcf.gz

	for f in *.bam; do
    /somalier:v0.2.11/somalier extract -d extracted/ --sites bin/sites.hg38.vcf.gz -f $genome_fasta \$f
	done

	/somalier:v0.2.11/somalier relate --ped $baseDir/bin/project.PED  extracted/*.somalier

	wget -P ancestry_files https://raw.githubusercontent.com/brentp/somalier/master/scripts/ancestry-labels-1kg.tsv

	wget https://zenodo.org/record/3479773/files/1kg.somalier.tar.gz
	tar -xzf 1kg.somalier.tar.gz

	/somalier:v0.2.11/somalier ancestry --labels ancestry_files/ancestry-labels-1kg.tsv \
	1kg-somalier/*.somalier ++ extracted/*.somalier

	rm -fr *.bam
  """
}



workflow.onComplete {
	// create log files and record output
	process finish{
		script:
		""" echo ' Pipeline cgpmap v0.3 completed
		Project: $projectname
		Time: ${nextflow.timestamp}

	 trimmomatic - completed
	 cgpmap - completed
	 (reference = $cgpmap_genome)
	 (index = $cgpmap_index)
	 samsort - completed
	 picard pcr removal - completed
	 picard rename bam - completed
	 samtools bam index - completed
	 picard collect insert size - completed
	 picard collect HS stats - completed
	 (target_intervals = $target_interval)
	 (bait_intervals = $target_interval)
	 merge bams - completed
	 ' >> $baseDir/${projectname}_log.txt

	 mail -s "cgpMAP successful" aft19qdu@uea.ac.uk < $baseDir/${projectname}_log.txt
	 mail -s "cgpMAP info" aft19qdu@uea.ac.uk < $baseDir/${projectname}_cgpmap_samples.log
	 """
	}
}
