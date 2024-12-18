# RNA-Seq Course

This repository contains code as well as result files produced during the RNA-Seq course of the University of Bern.

### Description of Files:
The following list describes what each file does. Some file numbers have an a and a b version. Here it's important that for this specific project only the b versions should be submitted to the cluster, as those files submit the a files. The a version of the files, however, can be reused for other projects by using commandline arguments for input-file, output-file and if data is pairwise input-file2 inbetween the other two.

- **[01_run_fastqc](01_run_fastqc)**: Produces FastQC files of all the reads, which are used to check the quality of the reads.
  
- **[02a_run_fastp](02a_run_fastp)**: Cleans up reads using fastp for the given read number.
- **[02b_submit_fastp](02b_submit_fastp)**: Submits all the read files for clean-up. Only useful in combination with 02a_run_fastp.

- **[03_run_fastqc_cleaned](03_run_fastqc_cleaned)**: Analyses the cleaned outputs gained through fastp and summarizes the resulting statistics in a multiqc-output
  
- **[04_get_assmbly_annt](04_get_assmbly_annt)**: Downloads the reference genome as well as the annotation file used for this project from Ensembl.
  
- **[05_run_hisat2_indexing](05_run_hisat2_indexing)**: Produces index files for the reference genome using HISAT2.
  
- **[06a_run_hisat2_mapping](06a_run_hisat2_mapping)**: Maps the reads to the reference genome using HISAT2 and produces .sam files as a result.
- **[06b_submit_hisat2_mapping](06b_submit_hisat2_mapping)**: Submits cleaned reads to mapping. Only useful in combination with 06a_run_hisat2_mapping
   
- **[07a_convert_sam_to_bam](07a_convert_sam_to_bam)**: Converts the resulting .sam files from 06a_run_hisat2_mapping into .bam files.
- **[07b_submit_sam_to_bam](07b_submit_sam_to_bam)**: Submits all sam files to conversion to bam. Only useful in combination with 07a_convert_sam_to_bam
  
- **[08a_run_samtools_sort](08a_run_samtools_sort)**: Sorts the BAM files using samtools.
- **[08b_submit_samtools_sort](08b_submit_samtools_sort)**: Submits all BAM files to sorting. Only useful in combination with 08a_run_samtools_sort.
  
- **[09a_run_samtools_index](09a_run_samtools_index)**: Produces index files for all the produced .bam files.
- **[09b_submit_samtools_index](09b_submit_samtools_index)**: Submits sorted BAM files to indexing. Only useful in combination with 09a_run_samtools_index.
  
- **[10_run_featureCounts](10_run_featureCounts)**: Runs the sorted BAM files through featureCounts to analyze differential gene expression.
  
- **[11_run_multiqc](11_run_multiqc)**: Submits the whole project folder to analysis by multiqc. Useful to summarize statistics produced by featureCounts and hisat2 during mapping                         of the reads

### Contact

For any questions, please contact: yannick.jauslin@students.unibe.ch

