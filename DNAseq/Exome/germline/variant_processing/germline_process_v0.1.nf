params.vcf= "$baseDir/output/VCF_collect/*.vcf.gz"

vcf_ch = Channel. fromPath (params.vcf)
vcf_ch.into { vcf1_ch; vcf2_ch }


// use pipeline bundle 1 has python3 installed
process merge_caller_indels {


  storeDir "$baseDir/output/VCF_collect/caller_merged"
  input:
  file vcf from vcf2_ch.collect()
  output:
  file "*indels.merged.vcf" into indels_filter_ch
  script:
  """
  for file in *.vcf*; do
    bcftools index \$file
  done

  python $baseDir/bin/python/merge_caller_indel.py --bam '$vcf'
  """
}

process indels_filter {


  storeDir "$baseDir/output/VCF_collect/caller_merged"
  input:
  file vcf from indels_filter_ch.flatten()
  output:
  file "${vcf.simpleName}_indels.vcf " into indels_sort_ch
  script:
  """
  bcftools view -V snps ${vcf} > ${vcf.simpleName}_indels.vcf
  """
}

process indels_sort {


  storeDir "$baseDir/output/VCF_collect/caller_merged"
  input:
  file vcf from indels_sort_ch
  output:
  file "${vcf.baseName}_sorted.vcf" into vep3_ch
  script:
  """
  mkdir -p tmp
  bcftools sort -O v -o ${vcf.baseName}_sorted.vcf ${vcf} -T tmp
  rm -fr tmp
  """

}

// dir needs to be changed
// and VEP_fasta added
process VEP2 {


  storeDir "$baseDir/output/VCF_collect/caller_merged/VEP"
  input:
  file vcf from vep3_ch
  output:
  file "${vcf.baseName}_VEP.txt" into vep4_ch
  script:
  """
  /ensembl-vep/vep -i ${vcf} \
  --dir /var/spool/mail/VEP_hg38/.vep \
  -o ${vcf.baseName}_VEP.txt \
  --offline \
  --fasta $VEP_fasta \
  --fork 5 \
  --cache homo_sapiens \
  --sift p \
  --polyphen p \
  --variant_class \
  --af_gnomad \
  --no_stats \
  --tab \
  --show_ref_allele \
  --symbol \
  --verbose \
  --domains \
  --regulatory
  """
}

process vep_header {
  storeDir "$baseDir/output/VCF_collect/caller_merged/VEP"
  input:
  file txt from vep4_ch
  output:
  file "${txt.baseName}_noheader.txt"
  script:
  """
  sed -i 's/#Uploaded_variation/Uploaded_variation/g' ${txt}
  awk '!/\\#/' ${txt} > ${txt.baseName}_noheader.txt
  """
}


// and VEP_fasta added
// use pipeline bundle 1 has python3 installed
process merge_caller_snps {


  storeDir "$baseDir/output/VCF_collect/caller_merged"
  input:
  file vcf from vcf1_ch.collect()
  output:
  file "*snps.merged.vcf" into snps_filter_ch
  script:
  """
  for file in *.vcf*; do
    bcftools index \$file
  done
  
  python $baseDir/bin/python/merge_caller_snps.py --bam '$vcf'
  """
}

process snps_filter {


  storeDir "$baseDir/output/VCF_collect/caller_merged"
  input:
  file vcf from snps_filter_ch.flatten()
  output:
  file "${vcf.simpleName}_snps.vcf" into snps_sort_ch
  script:
  """
  bcftools view -v snps ${vcf} > ${vcf.simpleName}_snps.vcf
  """
}

process snps_sort {


  storeDir "$baseDir/output/VCF_collect/caller_merged"
  input:
  file vcf from snps_sort_ch
  output:
  file "${vcf.baseName}_sorted.vcf" into vep_ch
  script:
  """
  mkdir -p tmp
  bcftools sort -O v -o ${vcf.baseName}_sorted.vcf ${vcf} -T tmp
  rm -fr tmp
  """

}


// dir needs to be changed
// and VEP_fasta added
process VEP3 {


  storeDir "$baseDir/output/VCF_collect/caller_merged/VEP"
  input:
  file vcf from vep_ch
  output:
  file "${vcf.baseName}_VEP.txt" into vep2_ch
  script:
  """
  /ensembl-vep/vep -i ${vcf} \
  --dir /var/spool/mail/VEP_hg38/.vep \
  -o ${vcf.baseName}_VEP.txt \
  --offline \
  --fasta $VEP_fasta \
  --fork 5 \
  --cache homo_sapiens \
  --sift p \
  --polyphen p \
  --variant_class \
  --af_gnomad \
  --no_stats \
  --tab \
  --show_ref_allele \
  --symbol \
  --verbose \
  --domains \
  --regulatory
  """
}

process vep_header2 {
  storeDir "$baseDir/output/VCF_collect/caller_merged/VEP"
  input:
  file txt from vep2_ch
  output:
  file "${txt.baseName}_noheader.txt"
  script:
  """
  sed -i 's/#Uploaded_variation/Uploaded_variation/g' ${txt}
  awk '!/\\#/' ${txt} > ${txt.baseName}_noheader.txt
  """
}
