params.vcf= "$baseDir/output/VCF_collect/*.vcf.gz"

vcf_ch = Channel. fromPath (params.vcf)
vcf_ch.into { vcf1_ch; vcf2_ch }


// use pipeline bundle 1 has python3 installed
process merge_caller_indels {


  storeDir "$baseDir/output/VCF_collect/caller_merged"
  input:
  file vcf from vcf2_ch.collect()
  output:
  file "*caller.merged.indels.vcf.gz" into fam_ch
  script:
  """
  for file in *.vcf*; do
    bcftools index \$file
  done

  python $baseDir/bin/python/merge_caller_indel.py --bam '$vcf'

  bgzip *.vcf
  """
}

process merge_family_indels {


  storeDir "$baseDir/output/VCF_collect/caller_merged"
  input:
  file vcf from fam_ch
  output:
  file "*family.merged.indels.vcf.gz" into indels_filter_ch
  script:
  """

  for file in *.vcf*; do
    bcftools index \$file
  done

  python $baseDir/bin/python/merge_family_indel.py --bam '$vcf'

  bgzip *.vcf
  """
}

process indels_filter {


  storeDir "$baseDir/output/VCF_collect/caller_merged"
  input:
  file vcf from indels_filter_ch.flatten()
  output:
  file "${vcf.simpleName}_indels.vcf.gz" into indels_sort_ch
  script:
  """
  bcftools view -O z -V snps ${vcf} > ${vcf.simpleName}_indels.vcf.gz
  """
}

process indels_sort {


  storeDir "$baseDir/output/VCF_collect/caller_merged"
  input:
  file vcf from indels_sort_ch
  output:
  file "${vcf.baseName}_sorted.vcf.gz" into vep3_ch
  script:
  """
  mkdir -p tmp
  bcftools sort -O z -o ${vcf.baseName}_sorted.vcf.gz ${vcf} -T tmp
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
  --dir $VEP_files \
  -o ${vcf.baseName}_VEP.txt \
  --offline \
  --fasta $genome_fasta \
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
  file "*callermerged.vcf.gz" into snps_fam_ch
  script:
  """
  for file in *.vcf*; do
    bcftools index \$file
  done

  python $baseDir/bin/python/merge_caller_snps.py --bam '$vcf'

  mcp '*/0001.vcf' '#1.vcf'

  bgzip *.vcf

  # The mcp function renames the file based on the path (where the * is )
  # this will take 0001.vcf which uses freebayes vcf info fields.
  """
}

process merge_fam_snps{


  storeDir "$baseDir/output/VCF_collect/caller_merged"
  input:
  file vcf from snps_fam_ch
  output:
  file "*familymerged.vcf.gz" into snps_filter_ch
  script:
  """
  for file in *.vcf.gz; do
    bcftools index \$file
  done

  python $baseDir/bin/python/merge_family_snps.py --bam '$vcf'

  mcp '*/0001.vcf' '#1.vcf'

  bgzip *.vcf

  # The mcp function renames the file based on the path (where the * is )
  # this will take 0001.vcf which uses freebayes vcf info fields.
  """
}



process snps_filter {


  storeDir "$baseDir/output/VCF_collect/caller_merged"
  input:
  file vcf from snps_filter_ch.flatten()
  output:
  file "${vcf.simpleName}_snps.vcf.gz" into snps_sort_ch
  script:
  """
  bcftools view -O z -v snps ${vcf} > ${vcf.simpleName}_snps.vcf.gz
  """
}

process snps_sort {


  storeDir "$baseDir/output/VCF_collect/caller_merged"
  input:
  file vcf from snps_sort_ch
  output:
  file "${vcf.baseName}_sorted.vcf.gz" into vep_ch
  script:
  """
  mkdir -p tmp
  bcftools sort -O z -o ${vcf.baseName}_sorted.vcf.gz ${vcf} -T tmp
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
  --dir $VEP_files \
  -o ${vcf.baseName}_VEP.txt \
  --offline \
  --fasta $genome_fasta \
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
