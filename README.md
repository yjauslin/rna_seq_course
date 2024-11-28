This repository contains code as well as result files produced during the rna-seq course of the University of Bern.

Description of all the files:
01_run_fastqc: produces fastqc files of all the reads, which were used to check the quality of the reads.
02_run_fastp: cleans up reads using fastp for given read number.
03_submit_fastp: submits all the reads files to clean-up. Only usefull in combination with 02_run_fastp.
04_get_assmbly_annt: downloads reference genome as well as annotation file used for this project from Ensembl.
05_run_hisat2_indexing: produces index-files for the reference genome using hisat2.
06_run_hisat2_mapping: maps the reads to the reference genome using hisat2. Produces .sam-files as result.
07_convert_sam_to_bam: converts the result .sam-files from 06_run_hisat2_mapping into .bam-files.
08_run_samtools_sort: Sorts the bam files using samtools.
09_run_samtools_index: Produces index-files for all the produced .bam files.
10_run_featureCounts: Runs the sorted bam-files through featureCounts in order to analyse differential gene-expression

Contact: yannick.jauslin@students.unibe.ch

