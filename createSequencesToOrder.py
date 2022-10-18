'''
cgat_script_template.py - template for cgat scripts
====================================================

:Author: Charlotte Vandermeulen
:Tags: Python

Purpose
-------

/!\ This script has to be run in a specific environment:
source activate /shared/sudlab1/General/projects/SynthUTR_hepG2_a549/charlotte_cgat/env_IFv2/


From - a linker sequence'
     - a list of motifs and
     - a table of RE and their sites (from RE_recombined_sites)

Checks that the linker combined to the list of motifs doesn't recreate a 
RE site for the 2 REs. If so, it stated which linker(s) are and they are removed
from the output.

The output correspond to sequences to order as:

5' overhang1-motif-overhang2 3' 
and
5' overhang3-RevMotif-overhang4 3'

Usage
-----

/!\ This script has to be run in a specific environment:
source activate /shared/sudlab1/General/projects/SynthUTR_hepG2_a549/charlotte_cgat/env_IFv2/

.. Example use case

Example::

   python createSequencesToOrder.py

'''

import itertools as iOT
from collections import defaultdict
import re
from Bio import Seq


#######START INPUTS###################
Table_RE_sites = "/shared/sudlab1/General/projects/SynthUTR_hepG2_a549/motif_merge_linkers/4linkers/RE_sites_linkers.txt"
out_seq = "/shared/sudlab1/General/projects/SynthUTR_hepG2_a549/motif_merge_linkers/4linkers/linker_motifs_sequences.txt"
motifs = open("/shared/sudlab1/General/projects/SynthUTR_hepG2_a549/motif_merge_linkers/4linkers/highstab.list").readlines()
Linker = "TTCGTC"
####### END INPUTS###################


def extend_ambiguous_dna(seq):
   """return list of all possible sequences given an ambiguous DNA input"""
   d = Seq.IUPAC.IUPACData.ambiguous_dna_values
   possibilities_per_position = map(d.get, seq)
   possible_seq = iOT.product(*possibilities_per_position)
   seqs_as_strings = map("".join, possible_seq)
   return list(seqs_as_strings)

   #  return list(map("".join, iOT.product(*map(d.get, seq))))


def combine_sequences(seq1,seq2):
    """return a list of all possible sequences associations
    inputs are lists"""
    merge = list()

    for i in seq1:
        for j in seq2:
            combine = i + j
            merge.append(combine)
    return list( set(merge) )

file = open(Table_RE_sites, 'r').readlines()[1:]

re_sites_pairs = dict()
for line in file:
    split_entry = line.split("\t")
    split_entry = [s.strip() for s in split_entry]
    re_names =  split_entry[0].split(",")
    re_sites = split_entry[1].split(",")
    re_sites_pairs[(re_names[0],re_names[1])] = (re_sites[0],re_sites[1])
#print(re_sites_pairs)

ambiguous = ("R","Y","M","K","S","W","B","V","D","N","H")
linkers_per_pairs = dict()
for key in re_sites_pairs:
    site1 = re_sites_pairs[key][0]
    site2 = re_sites_pairs[key][1]
    if "^" not in site1 or "^" not in site2:
        #print("The pair: ", key, re_sites_pairs[key],"is an exception and will have to be processed separately")
        continue
    cutindex1 = site1.index("^")
    cutindex2 = site2.index("^")
    #depends where the cut site is for similar sequence
    if (cutindex1 < len(site1)/2 and
        cutindex2 < len(site2)/2):
        similar1 = site1[cutindex1+1:-cutindex1]
        similar2 = site2[cutindex2+1:-cutindex2]
    else:
        similar1 = site1[len(site1)-1-cutindex1:cutindex1]
        similar2 = site2[len(site2)-1-cutindex2:cutindex2]
    #I won't process too complicated cases
    if len(similar1)==0 or len(similar2)==0:
        #print("The pair: ", key, re_sites_pairs[key],"is an exception and will have to be processed separately")
        continue
    if len(similar1) != len(similar2):
        #print("The pair: ", key, re_sites_pairs[key],"is an exception and will have to be processed separately")
        continue
    #keep sequence with no ambiguous
    if similar1 != similar2:
        if (any([letter in similar1 for letter in ambiguous]) and
           all([letter not in similar2 for letter in ambiguous])):
            similarOK = [similar2]
        elif (any([letter in similar2 for letter in ambiguous]) and
           all([letter not in similar1 for letter in ambiguous])):
            similarOK = [similar1]
        else:
            #print("The pair: ", key, re_sites_pairs[key],"is an exception and will have to be processed separately")
            continue
    else:
        similarOK = extend_ambiguous_dna(similar1)
        
    #Get linker sequences
    if (cutindex1 < len(site1)/2 and
        cutindex2 < len(site2)/2):
        site1_end = extend_ambiguous_dna(site1[-cutindex1:])
        linker_site2 = extend_ambiguous_dna(site2[:cutindex2])
        linker_site1 = combine_sequences(similarOK,site1_end)
        linker_sequence = combine_sequences(linker_site2,linker_site1)
    else:
        linker_site1 = extend_ambiguous_dna(site1[cutindex1+1:])
        site2_start = extend_ambiguous_dna(site2[:-cutindex2-1])
        linker_site2 = combine_sequences(site2_start,similarOK)
        linker_sequence = combine_sequences(linker_site2,linker_site1)
    #Create dict
    linkers_per_pairs[key] = linker_sequence

#print(linkers_per_pairs)

possible_re_pairs = set()
for re in linkers_per_pairs:
    #print(l)
    if Linker in linkers_per_pairs[re]:
        possible_re_pairs.add(re)

strip_motifs = list(s.strip() for s in motifs)
list_motifs = set()
for l in strip_motifs:
    if "U" in l:
        T_seq = l.replace("U", "T")
        l = T_seq
        list_motifs.add(l)
    else :
        list_motifs.add(l)

#Check for BsaHi and BstBI sites in motifs 
sequences_RE_sites = set()
for re in possible_re_pairs:
    site1 = re_sites_pairs[re][0]
    site2 = re_sites_pairs[re][1]
    site1 = extend_ambiguous_dna(site1.replace("^", ""))
    site2 = extend_ambiguous_dna(site2.replace("^", ""))
    full = site1 + site2
    sequences_RE_sites.update(full)

filter_motifs = set()
for i in list_motifs:
    sequence = Linker + i + Linker
    for j in sequences_RE_sites:
        if j in sequence:
            filter_motifs.add(i)
print("Motif(s): ", filter_motifs, """create(s) a site for the RE with the
       the linker and ahs been remove from the motifs list""")
list_motifs = list_motifs - filter_motifs

#Write table of sequences
outfile = open(out_seq, "w")
for re in possible_re_pairs:
    outfile.write("Sequences for RE: "+str(re)+"\t"+str(re_sites_pairs[re])+"\n")
    site1 = re_sites_pairs[re][0]
    site2 = re_sites_pairs[re][1]
    cutindex1 = site1.index("^")
    cutindex2 = site2.index("^")
    for m in list_motifs:
        forward = Linker[cutindex2:]+m+Linker[:cutindex2]
        reverse = Seq.Seq(Linker[-cutindex1:]+m+Linker[:-cutindex1]).complement()
        outfile.write(str(forward)+"\t"+str(reverse)+"\n")
outfile.close()
 
