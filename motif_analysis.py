"""===========================
motif_analysis.py
===========================

Overview
========

This pipeline computes conversion Rate from fastq SLAM-seq data

files :file:``pipeline.yml` and :file:`conf.py`.

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general
information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.yml` file.
CGATReport report requires a :file:`conf.py` and optionally a
:file:`cgatreport.ini` file (see :ref:`PipelineReporting`).

Default configuration files can be generated by executing:

   python <srcdir>/motif_analysis.py config

Input files
-----------

Inputs:

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


The pipeline configuration file pipeline.yml.

Requirements
------------

On top of the default CGAT setup, the pipeline requires the following
* Software:
    - python (v3.8.12 with pysam v0.17.0 when built)
    - meme (v5.3.0 when built)
    - HOMER in path (command findMotifs.pl in path)
* R modules:
   - optparse
   - stringr


Pipeline output
===============

Code
====

"""

import sys
import os
import sqlite3
import csv
import re
import glob

from cgatcore import pipeline as P
import cgat.GTF as GTF
import cgatcore.iotools as IOTools
from ruffus import *

#Load config file options
PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])
########################

@transform("*.bed",
           regex("(.+).bed"),
           r"\1.fasta")
def getfasta(infile, outfile):
    '''From most and less stable bed file, get sequences and generate fasta for streme'''
    genome_file = os.path.abspath(os.path.join(PARAMS["genome_dir"], PARAMS["genome"] + ".fa"))
    statement = '''
    bedtools getfasta -bed %(infile)s
    -fi %(genome_file)s
    -fo %(outfile)s
    -nameOnly
    '''
    P.run(statement)

#STREME
@subdivide(getfasta,
           regex("(.+)_highstab.fasta"),
           add_inputs(r"\1_lowstab.fasta"),
           [r"\1_highstab_streme.dir/streme.txt",
            r"\1_lowstab_streme.dir/streme.txt"])
def streme(infiles, outfiles):
    '''STREME motif enrichment motif_analysis
    https://meme-suite.org/meme/doc/streme.html?man_type=web'''
    highstab, lowstab = infiles
    highout, lowout = outfiles
    highout = P.snip(highout, "/streme.txt")
    lowout = P.snip(lowout, "/streme.txt")
    eval = PARAMS["e_value"]
    min = PARAMS["min_motif_width"]
    max = PARAMS["max_motif_width"]
    statement = '''
    streme --p %(highstab)s
    --n %(lowstab)s
    --minw %(min)s --maxw %(max)s --rna
    --thresh %(eval)s
    --evalue %(eval)s
    --patience 0
    --oc %(highout)s
    '''
    statement2 = '''
    streme --p %(lowstab)s
    --n %(highstab)s
    --minw %(min)s --maxw %(max)s --rna
    --thresh %(eval)s
    --evalue %(eval)s
    --patience 0
    --oc %(lowout)s
    '''
    P.run([statement,  statement2],
    job_memory="6G",
    job_threads=1)


#HOMER
@transform(getfasta,
           regex("(.+)_(highstab|lowstab).fasta"),
           add_inputs("background.fasta"),
           r"\1_\2_homer.dir/homerMotifs.all.motifs")
def homer(infiles, outfile):
    '''HOMER motif enrichment analysis
    http://homer.ucsd.edu/homer/motif/fasta.html'''
    stab, bg = infiles
    out = P.snip(outfile, "homerMotifs.all.motifs")
    length = PARAMS["motif_sizes"]
    statement = '''
    findMotifs.pl %(stab)s fasta %(out)s
    -fasta %(bg)s -rna -len %(length)s
    -noknown -noweight -nogo -fdr 100 -p 2
    '''
    P.run(statement ,
    job_memory="8G",
    job_threads=2)

#FIRE
@follows(getfasta)
@transform("fire_*.txt",
       regex("(fire_.+).txt"),
       add_inputs("fire.fasta"),
       r"\1.txt_FIRE/RNA/\1.txt.signif.motifs.rep")
def fire(infiles, outfiles):
    '''FIRE : https://tavazoielab.c2b2.columbia.edu/FIRE/tutorial.html
    When searching for RNA motifs (typically in 3’UTRs), FIRE examines all 16,384 7-mers. For
    each  7-mer,  the  mutual  information  between  its  profile  and  the  expression  profile  is
    evaluated. All 7-mers are then sorted based on their information values and a simple and
    information is not significant, within the sorted list. All 7-mers sorted above these 10 are
    retained  for  further  analysis,  and  are  henceforth  termed  motif  seeds.  Recall,  that  the
    information  associated  with  a  particular  7-mer  is  considered  significant  if  and  only  if  it
    passes  the  randomization  test,  i.e.,  if  it  is  greater  than  all  Nr  random  information  values
    obtained for this 7-mer profile over Nr randomly shuffled expression profiles. To correct
    for  multiple  hypothesis  testing,  Nr  is  set  by  default  to  the  number  of  k-mers  initially
    examined,'''
    expression, fasta = infiles
    expression = os.path.join(os.getcwd()+"/"+expression)
    fasta = os.path.join(os.getcwd()+"/"+fasta)
    statement = '''
    module load bio/fire &&
    fire --expfiles=%(expression)s
    --exptype=continuous
    --fastafile_rna=%(fasta)s
    --nodups=1 --k 8
    --dodna=0 --dodnarna=0
    '''
    P.run(statement,
    job_memory="10G",
    job_threads=1)


@follows(getfasta)
@transform("background.fasta",
           regex("(.+)"),
           r"\1.bg")
def fasta_to_bg(infile, outfile):
    '''Get markov bg for fasta'''
    statement = '''
    fasta-get-markov %(infile)s > %(outfile)s
    '''
    P.run(statement)


@transform(homer,
           regex("(.+)_(highstab|lowstab)_homer.dir/homerMotifs.all.motifs"),
       add_inputs(fasta_to_bg),
       r"\1_\2_homer.dir/homerMotifs.all.motifs.meme")
def homer_to_meme(infiles, outfile):
    '''Convert HOMER output to MEME format for tomtom'''
    homer_file, background = infiles
    out_dir = P.snip(outfile, "/tomtom.tsv")
    script_path = os.path.join((os.path.dirname(__file__)),
                               "Rscripts",
                               "homer2meme.R")
    statement = '''
    Rscript %(script_path)s
    -i %(homer_file)s
    -s %(background)s
    '''


@follows(streme,fire, homer_to_meme)
def full():
    '''Later alligator'''
    pass

P.main()
