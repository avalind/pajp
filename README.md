PAJP
=====

Introduction
------------

Snakemake based toolchain for somatic variant calling on tumor/normal pairs of
WES data. Implementing GATK best practices for .bam-file processing and then uses
MuTect and Scalpel to call SNVs and Indels. Furthermore use contest to estimate
cross-sample contamination and generate quality control reports for the bam-files
using QualiMap.

The vcf files of called somatic variants are annotated using Annovar.


The config.py script makes some pretty harsh assumptions on how the fastq-files
are ordered.



TODO
----

* Subworkflows.
* uploading of metadata to dropbox etc.
