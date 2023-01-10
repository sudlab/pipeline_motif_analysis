"""===========================
motif_analysis.py
===========================

Overview
========

This pipeline can be run after the pipeline pipeline_slamdunk_umis but it is
not necessary.
First, the Rscripts for the LASSO regression have to be ran on their own.
When bed files have been generated in appropriate directories, the
pipeline_motif_analysis can be run.
Afterwadrs, different python scripts exist out of the pipeline to get list of wanted
linkers to build libraries.
They are meant to be used in this order:
RE_recombined_sites.py (optional, need specific environment) -> selectLinkers.py ->
createSequencesToOrder.py.

A report of the pipeline can be built using build_report after finishing the
full pipeline.

files :file:``pipeline.yml` and :file:`conf.py`.

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general
information how to use CGAT pipelines.

Run the pipeline with `python [path_to_repo]/motif_analysis.py make full -v5`.

Run the report render (after doing full):
`python [path_to_repo]/motif_analysis.py make build_report -v5`

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
   [X]_lowstab|highstab.bed : bed files consisting of 3'UTR
   sequences of transcript with highest or lowest half-lives or residuals.
2. HOMER inputs
   [X]_lowstab|highstab.bed : bed files consisting of 3'UTR
   sequences of transcript with highest or lowest half-lives or residuals.
   backgroud.bed : bed file consisting of all 3'UTR sequences of transcript
   detected by the pipeline_slamdunk_umis.
3. fire
   fire.bed : bed file consisting of all 3'UTR sequences with a length >6 and
   <10000
   fire_[X].txt : transcript id and ranking value table
were [X] can be any name where you data (ranking value) originate from (halflife, ...).

The pipeline configuration file pipeline.yml.

Requirements
------------

On top of the default CGAT setup, the pipeline requires the following
* Software:
    - python (v3.8.12 with pysam v0.17.0 when built)
    - meme (v5.3.0 when built)
    - HOMER in path (command findMotifs.pl in path)
    - fire in path (https://tavazoielab.c2b2.columbia.edu/FIRE/)
* R modules:
   - Biostrings
   - tidyverse
   - optparse
   - stringr
   - tools
   - universalmotif
   - msa


Pipeline output
===============

Outputs by directories

* [X]_[lowstab|highstab]_streme.dir
- Outputs generate by streme:
    streme.html - an HTML file that provides the results in an interactive,
    human-readable format
    streme.txt - a text file containing the motifs discovered by STREME in
    MEME format
    sequences.tsv - a TSV (tab-separated values) file that lists the true- and
    false-positive sequences identified by STREME for each motif
    streme.xml - an XML file that provides the results in a format designed
    for machine processing
(source: https://meme-suite.org/meme/doc/streme.html)

- tomtom.self
Final output directory, obtained from running tomtom on the ouput file to get
rid of redundant motifs (log of run: streme.txt.tomtom.log)

* [X]_[lowstab|highstab]_homer.dir
- Outputs generated by runing homer findMotifs.pl script
http://homer.ucsd.edu/homer/motif/fasta.html
(NB: randomization folder can be deleted after pipeline as finished running
homer)
homerMotifs.all.motifs contains all discovered motifs by Homer.

- homerMotifs.all.motifs.meme
homerMotifs.all.motifs file in meme format.

- homerMotifs.all.motifs.tomtom.self
Output obtained from running tomtom on the ouput file (homerMotifs.all.motifs)
to get rid of redundant motifs.

* fire_[X].txt.[6-7-8]imer_FIRE
Directories for each kmer size (I have little control over naming of directories,
which explain their weird names.)
- Outputs generated fron FIRE
/RNA directory contains interesting outputs. Fire generates a lot of different
outputs and they don't explain most of them in their tutorial, the important
file are "...signif.motifs.rep" and "....signif.motifs"
(https://tavazoielab.c2b2.columbia.edu/FIRE/tutorial.html)
- ..._[highstab|lowstab].signif.motifs
List of motifs enriched in high or low stability transcripts discovered by fire.

* fire.dir
- [X]_[highstab|lowstab].allkmer.signif.motifs
Merged results from the different kmer sizes
- [highstab|lowstab].allkmer.fireMotifs
Merging of all fire results in either low and high stability motifs
Contains the name of the motif and it's consensus sequence given by fire in the
file "....signif.motifs".
- highstab_in_lowstab.list
List of motifs sequecnces from lowstab also present in lowstab (later filtered)

* final_motifs
- [highstab|lowstab]_final_motifs.list
The most important output.
Table with the name of the motif and its associated consensus sequence.

- [highstab|lowstab]_final_motifs.list.log
Log file giving:
 a. the number of motifs coming from Streme/Homer or Fire,
 b. the number of consensus sequences coming from these motifs,
 c.the number of motifs shared between highstab and lowstab, for
  Homer/Streme
 d. the number of sequences shared between highstab and lowstab for Fire
Once Streme/Homer and Fire results have been merged
 c. number of sequences that contained miRNA seed targets,
 d. number of sequences with one of the most 2 common polyA signals (AAUAAA or AUAAA)
 e. Final number of sequences, once unwanted ones have been removed
(For lowstab the list is not filtered)

- [highstab|lowstab]_final_motifs.matching.mirna.seeds
Table of sequences matching miRNA seed targets
struture: miRNA_name:miRNA_seed_target_sequence   name_matching_motif   motif_sequence

- [X]_merge_homer_streme.meme
Merged motifs from Homer and Streme, similar ones have been clustered using tomtom

- [X]_merge_homer_streme.motifs.tomtom
Output obtained from running tomtom on the merge motifs originating from homer and streme
highstab_merge_homer_streme.meme.log is the tomtom output log.

Once this pipeline has been run, you can then run merge_motifs.Rmd to merge all
motif sequences from the different pipeline runs and gerenate a general list of
highstab and lowstab motifs. Then the python scrypts RE_recombined_sites.py
and/or select_linkers.py can be run to generate linkers.

The report render ouputs in final_motifs pipeline_report.html and associated files.

Code
====

"""

import sys
import os
import sqlite3
import csv
import re
import glob
import pandas
import shutil

from cgatcore import pipeline as P
from cgatcore import experiment as E
from ruffus.combinatorics import *

import cgat.GTF as GTF
import cgatcore.iotools as IOTools
from ruffus import *
from cgat.MEME import MemeMotifFile, MotifCluster
import PipelineTomtom

#Load config file options
PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])
########################

@active_if(PARAMS["from_fasta"])
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
    -s
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
    --evalue
    --patience 0
    --oc %(lowout)s
    '''
    P.run([statement,  statement2],
    job_memory="6G",
    job_threads=1)
    #If Streme didn't find any significant motifs, it didn't ouput it
    for outf in outfiles:
        if os.path.exists(outf) == False:
            IOTools.touch_file(outf)


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
    iter = PARAMS["FDRiteration"]
    statement = '''
    findMotifs.pl %(stab)s fasta %(out)s
    -fasta %(bg)s -rna -len %(length)s
    -noknown -noweight -nogo -fdr %(iter)s -p 8
    '''
    P.run(statement ,
    job_memory="8G",
    job_threads=6)


if (PARAMS["8kmer"] == True):
    end_kmer = 9
else:
    end_kmer = 8

#FIRE
@follows(getfasta)
@subdivide(["fire_*.txt"],
       formatter(),
       add_inputs("fire.fasta"),
       [r"{path[0]}/{basename[0]}.txt.%imer_FIRE/RNA/{basename[0]}.txt.%imer.signif.motifs.rep" % (i,i) for i in range(6,end_kmer)])
       #"{path[0]}/{basename[0]}.txt_FIRE/RNA/{basename[0]}.txt.*.signif.motifs.rep")
       #[r"\1_FIRE/RNA/\1%imer.signif.motifs.rep" % i for i in range(6,8)])
def fire(infiles, outfiles):
    '''FIRE : https://tavazoielab.c2b2.columbia.edu/FIRE/tutorial.html
    When searching for RNA motifs (typically in 3’UTRs), FIRE examines all
    16,384 7-mers ad esempio. For each  k-mer,  the  mutual  information
    between  its  profile  and  the  expression  profile  is evaluated. All
    k-mers are then sorted based on their information values and a simple and
    information is not significant, within the sorted list. All k-mers sorted
    above these 10 are retained  for  further  analysis,  and  are  henceforth
    termed  motif  seeds.  Recall,  that  the information  associated  with
    a  particular  k-mer  is  considered  significant  if  and  only  if  it
    passes  the  randomization  test,  i.e.,  if  it  is  greater  than  all
    Nr  random  information  values obtained for this k-mer profile over Nr
    randomly shuffled expression profiles. To correct for  multiple  hypothesis
    testing,  Nr  is  set  by  default  to  the  number  of  k-mers  initially
    examined'''
    expression, fasta = infiles
    expression = os.path.join(os.getcwd()+"/"+expression)
    fasta = os.path.join(os.getcwd()+"/"+fasta)

    statements = list()
    for kmer in range(6,end_kmer):
        ikmer = str(kmer)+"mer"
        statement = '''
        module load bio/fire &&
        fire --expfiles=%(expression)s
        --exptype=continuous
        --fastafile_rna=%(fasta)s
        --nodups=1 --k=%(kmer)s
        --suffix=%(ikmer)s
        --dodna=0 --dodnarna=0
        --oribiasonly=0
        ''' % locals()
        statements.append(statement)
    P.run(statements,
    job_memory="12G",
    job_threads=1)
    #If Fire didn't find any significant motifs, it didn't ouput it
    for outf in outfiles:
        if os.path.exists(outf) == False:
            IOTools.touch_file(outf)

@follows(getfasta)
@transform("background.fasta",
           regex("(.+)"),
           r"\1.bg")
def fasta_to_bg(infile, outfile):
    '''Get markov bg for fasta'''
    statement = '''
    fasta-get-markov %(infile)s -norc > %(outfile)s
    '''
    P.run(statement)

#HOMER conversion
@transform(homer,
           regex("(.+)_(highstab|lowstab)_homer.dir/homerMotifs.all.motifs"),
           add_inputs(fasta_to_bg),
           r"\1_\2_homer.dir/homerMotifs.all.motifs.meme")
def homer_to_meme(infiles, outfile):
    '''Convert HOMER output to MEME format for tomtom'''
    homer_file, background = infiles
    #out_dir = P.snip(outfile, "/homerMotifs.all.motifs.meme")
    script_path = os.path.join((os.path.dirname(__file__)),
                               "Rscripts",
                               "homer2meme.R")
    fdr_thresh = length = PARAMS["fdr"]
    statement = '''
    Rscript %(script_path)s
    -i %(homer_file)s
    -b %(background)s
    -t %(fdr_thresh)s
    '''
    P.run(statement)

#Fire merge and conversion
@subdivide(fire,
           regex("(.+)txt.([0-9]mer).signif.motifs.rep"),
           [r"\1\2_highstab.signif.motifs", r"\1\2_lowstab.signif.motifs"])
def extractFire(infile, outfiles):
    '''Extract significant motifs from FIRE enriched in top or bottom bin'''
    if (IOTools.is_empty(infile)):
        for outf in outfiles:
            IOTools.touch_file(outf)
    else:
        script_path = os.path.join((os.path.dirname(__file__)),
                                "Rscripts",
                                "extract_fire_motifs.R")
        statement = '''
        Rscript %(script_path)s
        -f %(infile)s
        '''
        P.run(statement)

@follows(mkdir("fire.dir"))
@collate(extractFire,
        regex(r"fire_(.+).txt.(?:[0-9]mer).+([0-9]mer)_(highstab|lowstab).signif.motifs"),
        r"fire.dir/\1_\3.allkmer.signif.motifs")
def mergeFireKmers(infiles, outfile):
    '''Merge the fire kmer results together'''
    input_string = " ".join(infiles)
    statement= '''
    cat %(input_string)s > %(outfile)s
    '''
    P.run(statement)

@collate(mergeFireKmers,
         regex("(fire.dir/).+(highstab.allkmer|lowstab.allkmer).+"),
         r"\1\2.fireMotifs")
def removeRedundantFire(infiles, outfile):
    '''Remove redundants motifs from fire output'''
    motif_filtered = [x for x in infiles if IOTools.is_empty(x) == False]
    if len(motif_filtered) == 0:
        IOTools.touch_file(outfile)
    else:
        input_string = ",".join(motif_filtered)
        script_path = os.path.join((os.path.dirname(__file__)),
                                   "Rscripts",
                                   "removeRedundant.R")
        statement = '''
        Rscript %(script_path)s
        -f %(input_string)s
        -o %(outfile)s
        '''
        P.run(statement)


@merge(removeRedundantFire,
       r"fire.dir/highstab_in_lowstab.list")
def filterHighstabFire(infiles,outfile):
    """"Generate list of fire highstab motifs in lowstab for final merge"""
    highstab = [x for x in infiles if "highstab" in x]
    highstab = "".join(highstab)
    lowstab = [x for x in infiles if "lowstab" in x]
    lowstab = "".join(lowstab)
    if (IOTools.is_empty(highstab) or IOTools.is_empty(lowstab)):
        IOTools.touch_file(outfile)
    else:
        script_path = os.path.join((os.path.dirname(__file__)),
                                   "FilterLists.py")
        statement = """
        python %(script_path)s
        -I %(highstab)s
        -l %(lowstab)s
        -S %(outfile)s
        """
        P.run(statement)


#Merge streme and Homer results, then do Tomtom to remove redundants
@follows(mkdir("final_motifs"))
@collate([streme,homer_to_meme],
         regex(".*(lowstab|highstab).+"),
         add_inputs(fasta_to_bg),
         r"final_motifs/\1_merge_homer_streme.meme")
def tomtom_combine(infiles, outfile):
    '''Merge all motifs together, then run tomtom on merge and eliminate
    redundant motifs to create the final list of motifs'''
    motif_files = [i[0] for i in infiles]
    background = infiles[0][1]
    motif_filtered = []
    for i in motif_files:
        if (i.count("streme.txt") >= 1) and (IOTools.get_num_lines(i) <= 38):
            E.debug("No motifs in ", str(i))
        elif (i.count("homer") >= 1) and (IOTools.get_num_lines(i) <= 9):
            E.debug("No motifs in ", str(i))
        else:
            motif_filtered.append(i)
    if len(motif_filtered) == 0 :
        E.debug("No motifs present for ", str(motif_files))
        IOTools.touch_file(outfile)
        return
    if len(motif_filtered) == 1 :
        motif_filtered = motif_filtered[0]
        E.debug("Only motifs for ", str(motif_filtered))
        statement = '''
        mv %(motif_filtered)s %(outfile)s
        '''
        P.run(statement, job_options = "-P gen2reg -l h_rt=1:00:00")
        return
    input_string = " ".join(motif_filtered)
    temp_file = P.snip(outfile, ".meme")+".temp.merge"
    statement = '''
    meme2meme -bg %(background)s
    %(input_string)s > %(temp_file)s
    '''
    P.run(statement)
    q_val = PARAMS["thresh_merge"]
    tomtom_log = outfile+".log"
    temp_tomtom = P.snip(outfile, ".meme")+".motifs.tomtom"
    statement = '''
    tomtom -verbosity 1 -text -norc -thresh %(q_val)s
    %(temp_file)s %(temp_file)s
    2> %(tomtom_log)s | sed 's/#//' > %(temp_tomtom)s
    '''
    P.run(statement)
    PipelineTomtom.getSeedMotifs(temp_file, temp_tomtom, outfile)
    statement2 = """
    rm %(temp_file)s
    """
    P.run(statement2)


@merge(tomtom_combine,
       r"final_motifs/merge_homer_streme_highstab_vs_lowstab.tomtom")
def TomtomHighVSlow(infiles, outfile):
    '''Scan highstab motifs vs lowstab motifs with tomtom'''
    tomtom_log = outfile+".log"
    E.debug(infiles)
    motif1, motif2 = infiles
    if (IOTools.is_empty(motif1) | IOTools.is_empty(motif2)):
        IOTools.touch_file(outfile)
    else:
        qval = PARAMS["thresh_final"]
        statement = '''
        tomtom -verbosity 1 -text -norc -thresh %(qval)s
        %(motif1)s %(motif2)s
        2> %(tomtom_log)s | sed 's/#//' > %(outfile)s
        '''
        P.run(statement)


@collate([tomtom_combine, removeRedundantFire],
         regex(".+(highstab|lowstab)(?:.+)"),
         add_inputs(TomtomHighVSlow, filterHighstabFire),
         r"final_motifs/\1_final_motifs.list")
def memeToList(infiles, outfile):
    '''Convert all outputs to list of sequences for linker finder'''
    motif_file = [i[0] for i in infiles if "merge_homer_streme.meme" in i[0]]
    motif_file = motif_file[0]
    fire_file = [i[0] for i in infiles if "fire" in i[0]]
    fire_file = fire_file[0]
    high_vs_low_tomtom = [x for x in infiles[0] if "highstab_vs_lowstab.tomtom" in x]
    high_vs_low_tomtom = "".join(high_vs_low_tomtom)
    filter_fire = [x for x in infiles[0] if "fire.dir/highstab_in_lowstab.list" in x]
    filter_fire = "".join(filter_fire)
    mirna_seeds = PARAMS["mirna_seeds_db"]
    filter_seeds = PARAMS["filter_mirna"]
    if PARAMS["filter_mirna"]:
        supp = "-r "+PARAMS["relevant_mirna"]
    else:
        supp = ""
    script_path = os.path.join((os.path.dirname(__file__)),
                               "Rscripts",
                               "merge_filter4mirna_meme2list.R")
    statement = '''
    Rscript %(script_path)s
    -i %(motif_file)s
    -f %(fire_file)s
    -m %(mirna_seeds)s
    -t %(high_vs_low_tomtom)s
    -l %(filter_fire)s
    -o %(outfile)s
    %(supp)s
    '''
    P.run(statement)

@follows(memeToList)
def full():
    '''Later alligator'''
    pass

@follows(memeToList)
@originate("final_motifs/pipeline_report.html")
def renderReport(infile):
    '''build pipeline report'''
    script_path = os.path.join((os.path.dirname(__file__)),
                               "Rscripts",
                               "report_pipeline.Rmd")
    path_directory = os.path.abspath(os.getcwd())
    #Sys.setenv(RSTUDIO_PANDOC="/shared/sudlab1/General/projects/SynthUTR_hepG2_a549/charlotte_cgat/env/bin/pandoc");
    #'Sys.setenv(RSTUDIO_PANDOC="/shared/sudlab1/General/projects/SLAMseq_CHO_Mikayla/env/bin/pandoc");
    statement = '''
    Rscript -e
    'rmarkdown::render(input = "%(script_path)s",
                       knit_root_dir = "%(path_directory)s",
                       output_file = "%(path_directory)s/%(infile)s")'
    '''
    P.run(statement)


@follows(renderReport)
def build_report():
    '''Later alligator'''
    pass


P.main()
