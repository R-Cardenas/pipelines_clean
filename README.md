# UEA BCRE pipelines

![logo](misc/logo.png)


<br />

<!-- TABLE OF CONTENTS -->
## Table of Contents

<br />

* [Introduction](#Introduction)
* [Quality Control](#Quality-Control)

<br />


## Introduction

<br />

At UEA bob champion genomics we have developed a number of pipelines specialising for processing sequencing data. We have developed a system that allows the easy use of these pipelines by downloading the github repository, moving fastq files into input directory and configuring which analyses to perform by configuring the "XXX.config" file.

<br />

Irrespective of the pipelines chosen, all workflows perform similar analyses:

<br />

![figure-1](misc/figure1.png)

<br />



### Quality Control

The quality control pipeline contains a static backbone that is present with all data types (black). Dependent on the type of



### Quality Control

Current pipelines are able to process the following data input:

DNA-seq:

  - Exome-germline (Freebayes and GATK HaplotypeCaller)
  - Exome-somatic (Sanger-cgpWXS and GATK Mutect2)

  - Whole-genome germline (Freebayes and GATK HaplotypeCaller)
  - Whole-genome somatic (Sanger-cgpWGS and GATK Mutect2)
  - Whole-genome structural variants (Sanger-XXX)

RNA-seq:

  - mRNA hisat2









####DNA:

  Exome:

  +-- Mapping and QC (cgpMAP)

    Germline:
      +-- Freebayes
      +-- GATK HaplotypeCaller

    Somatic:
      +-- Sanger cgpWXs
      +-- GATK Mutects


  Whole Genome:

  +-- Mapping and QC (cgpMAP)

    Germline:
      +-- Freebayes
      +-- GATK HaplotypeCaller

    Somatic:
      +-- Sanger cgpWGS
      +-- GATK Mutects
