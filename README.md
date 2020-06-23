# UEA BCG pipelines

<!-- TABLE OF CONTENTS -->
## Table of Contents

* [Introduction](#Introduction)


## Introduction

At UEA bob champion genomics we have developed a number of pipelines specialising in processing sequencing data. A number of pipelines are available for the processing of both RNA and DNA sequencing experiments - the pipelines can be used in different combinations depending on the nature of the analysis required. Below shows the pipelines that are currently working:


Attempt | #1 | #2 | #3 | #4 | #5 | #6 | #7 | #8 | #9 | #10 | #11
--- | --- | --- | --- |--- |--- |--- |--- |--- |--- |--- |---
Seconds | 301 | 283 | 290 | 286 | 289 | 285 | 287 | 287 | 272 | 276 | 269

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
