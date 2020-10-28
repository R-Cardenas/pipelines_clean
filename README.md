# UEA BCRE pipelines

![logo](misc/logo.png)


<br />

<!-- TABLE OF CONTENTS -->
## Table of Contents

<br />

* [Introduction](#Introduction)
  - [Quality Control](#Quality-Control)
  - [Pipelines](#Pipelines-available)

* [Usage](#Usage)
  - [Software Requirements](#Software-Requirements)
  - [Example](#Example)

<br />


## Introduction

<br />

At UEA bob champion genomics we have developed a number of pipelines specialising for processing sequencing data. We have developed a system that allows the easy use of these pipelines by downloading the github repository, moving fastq files into input directory and configuring which analyses to perform by configuring the "XXX.config" file. The pipeline has been written using the nextflow workflow management.

<br />

The pipelines we have currently are shown below. Click on each individual link to access each README md. For the user README, which explains how to use the pipelines please click here (link).

DNA-seq:

  Mapping:
    - cgpMAP (link)

  Somatic variant discovery Exome:
    - cgpWXS (link)
    - GATK mutect2 (link)

  Somatic variant discovery WGS:
    - cgpWGS (link)
    - GATK mutect2 (link)

  Germline variant discovery:
    - Nextflow Sarek (exome/WGS; link)

RNA-seq:
  - Nextflow RNA-seq pipeline (link)
