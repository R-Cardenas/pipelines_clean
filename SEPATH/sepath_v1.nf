// This script is a rewrite of Abe's SEPATH (https://github.com/UEA-Cancer-Genetics-Lab/sepath_tool_UEA/blob/master/bin/fastq_filtering/all_SEPATH.snake)
//but using Nextflow

bam.params = "$baseDir/input/*.bam"
bam_ch = Channel .fromPath( params.bam )


// ask abe which is the big file
// need to put in contam and kraken db

process bam_filter {
  storeDir: "$baseDir/output/bam_filter"
  input:
  file bam from bam_ch
  output:
  file "*2.fq.gz" into fq1_ch
  file "*1.fq.gz" into fq2_ch
  script:
  """
  $baseDir/bin/bam_filter.py --input_bam $bam
  """
}

process trimmomatic {
  storeDir: "$baseDir/output/trimmomatic"
  input:
  file fq1 from fq1_ch
  file fq2 from fq2_ch
  output:
  file "*1.paired.fq.gz" into fq1pair_ch
  file "*1.unpaired.fq.gz" into fq1unpair_ch
  file "*2.paired.fq.gz" into fq2pair_ch
  file "*2.unpaired.fq.gz" into fq2unpair_ch
  script:
  """
  trimmomatic PE $fq1 $fq2 ${fq1.simpleName}.paired.fq.gz ${fq1.simpleName}.unpaired.fq.gz \
  ${fq2.simpleName}.paired.fq.gz ${fq2.simpleName}.paired.fq.gz \
  SLIDINGWINDOW:4:20 MINLEN:35
  """
}

process bbduk {
  storeDir: "$baseDir/output/bbduk"
  input:
  file fq1 from fq1unpair_ch
  file fq2 from fq2unpair_ch
  file fq11 from fq1pair_ch
  file fq22 from fq2pair_ch
  output:
  file "${fq1.simpleName}_unpaired_cat_bbduk.fq.gz​​" into f1_unpair_ch

  file "${fq11.simpleName}.out.R1​​" into f1_pair_ch
  file "${fq11.simpleName}.out.R2​​" into f2_pair_ch

  script:
  """
  zcat ${fq1} ${fq2} | gzip > ${fq1.simpleName}_unpaired_cat.fq.gz

  # single human depletion
  bbduk.sh \
  in1=${fq1.simpleName}_unpaired_cat.fq.gz \
  out=${fq1.simpleName}_unpaired_cat_bbduk.fq.gz​​ \
  k=30 \
  -Xmx230g \
  ref={​​contaminant_db}​​ \
  mcf=0.5

  # PE bbduk
  bbduk.sh in1=$fq11 in2=$fq22 \
  out1=${fq11.simpleName}.out.R1 \
  out2=${fq11.simpleName}.out.R2 \
  outs=${fq11.simpleName}.out.R3 \
  k=30 \
  -Xmx230g \
  ref={contaminant_db} \
  mcf=0.5

  """
}

// ask abe what is it R1 and R2 we put into kraken
process kraken paired {
  storeDir: "$baseDir/output/kraken"
  input:
  file f1_unpair from f1_unpair_ch

  file f1_pair from f1_pair_ch
  file f2_pair from f2_pair_ch
  output:
  file "${pair.simpleName}.merge.kraken" into thresh_ch
  script:
  """
  #_unpaired
  /gpfs/software/kraken/0.10.6/kraken \
  --preload --db {​​krakendb}​​ --threads 10 --fastq-input --gzip-compressed ${f1_pair}​​ --output ${f1_pair.simpleName}.unpair.kraken

  #paired
  /gpfs/software/kraken/0.10.6/kraken --preload \
  --db {​​krakendb}​​ --threads 10 --fastq-input \
  --paired --gzip-compressed ${f1_pair}​​ ${f2_pair} --output ${f1_pair.simpleName}.pair.kraken​

  cat ${f1_pair.simpleName}.pair.kraken​ ${f1_pair.simpleName}.unpair.kraken > ${pair.simpleName}.merge.kraken

  rm -fr *.unpair.kraken *.pair.kraken
  """
}

// kraken threshold on the concatentated files
// 20% of kmers has to meet the taxonomy
process threshold {
  storeDir: "$baseDir/output/bam_filter"
  input:
  file kraken from thresh_ch
  output:
  file "${kraken.simpleName}.theshold.kraken​​" into report_ch
  script:
  """
  kraken-filter --db {​​krakendb}​​ --threshold 0.2 ${kraken}​​ > ${kraken.simpleName}.theshold.kraken​​
  """
}

//produces readable file
process report {
  storeDir: "$baseDir/output/bam_filter"
  input:
  file kraken from report_ch
  output:
  file "*2.trim.fq.gz" into fq3_ch
  script:
  """
  kraken-report --db {​​krakendb}​​ ${kraken} > ${kraken.simpleName}.report.final
  """
}
