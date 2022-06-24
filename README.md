# pipeline_motif_analysis

This pipeline is meant to be ran after the pipeline pipeline_slamdunk_umis.
First, the Rscripts for the LASSO regression have to be ran on their own.
When bed files have been generated in appropriate directories, the
pipeline_motif_analysis can be run.



## The pipeline performs the following:
   * Run STREME on highly stable and lowly stable sequences

   * Run HOMER (findMotifs.pl) on highly stable and lowly stable sequences
   http://homer.ucsd.edu/homer/motif/fasta.html
   * Run FIRE
   https://tavazoielab.c2b2.columbia.edu/FIRE/
   * Convert HOMER and FIRE outputs to MEME motif format
   * Run Tomtom to merge motifs from STREME, HOMER and FIRE together


## Inputs needed      
   1. STREME inputs
   halflife|residual_lowstab|highstab.bed : bed files consisting of 3'UTR
   sequences of transcript with highest or lowest half-lives or residuals.
   2. HOMER inputs
   halflife|residual_lowstab|highstab.bed : bed files consisting of 3'UTR
   sequences of transcript with highest or lowest half-lives or residuals.
   backgroud.bed : bed file consisting of all 3'UTR sequences of transcript
   detected by the pipeline_slamdunk_umis.
   3. fire
   fire.bed : bed file consisting of all 3'UTR sequences with a length >6 and
   <10000
   fire_residual|halflife.txt : transcript id and halflife|residual table


## The pipeline outputs the different files/directories
   * sample_description.tsv file with sample names and time points asked by
       slamdunk_all
   * xxx_processed.fastq files : ouputs of the umi tools extract
   * map directory : contains outputs from slamdunk map
       "xxx_slamdunk_mapped.bam" and  .log
   * filter directory : contains output from
       - slamdunk filter "xxx_slamdunk_mapped_filtered.bam" with index .bai
       and .log
       - umi_tools dedup "xxx_slamdunk_mapped_filtered_dedup.bam" with
       index .bai
       - slamdunk alleyoop summary "xxx_slamdunk_mapped_filtered_summary.tsv"
       and .log, "xxx_slamdunk_mapped_filtered_summary_PCA.txt"


## Requirements

On top of the default CGAT setup, the pipeline requires the following
* Software:
    - python (v3.8.12 with pysam v0.17.0 when built)
    - meme (v5.3.0 when built)
    - HOMER in path (command findMotifs.pl in path)
* R modules:
   - optparse
   - stringr

## Configuration
The pipeline requires a configured :file: `pipeline.yml` file.

Make a directory with your project name.
Configure the pipeline with `python [path_to_repo]/pipeline_motif_analysis.py config`.
A pipeline.log and pipeline.yml file(s) will be added to your new directory.
Modify the pipeline.yml according to your project (specify annotation database and directory, database for uploading the outputs; specify options for Salmon quantification).

## Pipeline use
Run the pipeline with `python [path_to_repo]/pipeline_motif_analysis.py make full -v5`.

For running the pipeline on a large set of samples, submit the pipeline onto the cluster (sharc), using a submit_pipeline custom script.

## Pipeline under construction

So far runs STREME, HOMER, FIRE and generate meme format outputs
