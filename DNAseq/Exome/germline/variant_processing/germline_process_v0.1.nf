params.vcf= "$baseDir/output/VCF_collect/variant_R_test/*.vcf.gz"

vcf_ch = Channel. fromPath (params.vcf)


// dir needs to be changed
// and VEP_fasta added
process VEP2 {
  errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
  maxRetries 6
  storeDir "$baseDir/output/VCF_collect/variant_R_test/"
  input:
  file vcf from vcf_ch
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

process vep_header {
  storeDir "$baseDir/output/VCF_collect/variant_R_test/"
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
