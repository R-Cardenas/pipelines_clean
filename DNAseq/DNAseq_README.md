# UEA BCRE pipelines - DNAseq

# Introduction

This README explains how the inputs / outputs are organised and how this connects multiple nextflow pipeline scripts together. For more in-depth details about the individual pipelines please select from the list below:

cgpMAP ([link](Exome/cgpmap/README.md))


# Connecting pipelines

In order to connect the pipelines the outputs from one pipeline mush match what is described in the next pipelines input variable - in both name and location. A graphic is shown below how the pipelines are linked:


![figure-1](../misc/DNAseq_summary.png)
