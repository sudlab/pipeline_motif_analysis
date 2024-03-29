'''
RE_recombined_sites - create linkers from RE with cohesive ends
===============================================================

:Author: Charlotte Vandermeulen
:Tags: Python

Purpose
-------

From a list of Restriction Enzymes pairs with cohesive ends, create a list of
possible linkers

Usage
-----

/!\ This script has to be run in a specific environment:
source activate /shared/sudlab1/General/projects/SynthUTR_hepG2_a549/charlotte_cgat/env_IFv2/

Input a dictionnary of possible pairs of REs, in the script, then run the
the script.

   python RE_recombined_sites.py

Output
------
A text file containing a list of linkers based on RE pairs.

'''

import itertools as iOT
from collections import defaultdict
import re
from Bio import Seq

#######START INPUTS###################
name_ouput = "/shared/sudlab1/General/projects/SynthUTR_hepG2_a549/motif_merge_linkers/correct_strandedness/linker/list_linkers.txt"
name_enzyme_RE_site_table = "/shared/sudlab1/General/projects/SynthUTR_hepG2_a549/motif_merge_linkers/correct_strandedness/linker/RE_sites_linkers.txt"
#######Table of RE pairs#######
enzymes_dict = {"Acc65I":["BsiWI", "BsrGI"],
"AgeI":["AvaI", "BspEI", "SgrAI", "NgoMIV"],
"AscI":["AflIII", "MluI"],
"AsiSI":["BsiEI", "PacI"],
"AvrII":["NheI", "SpeI", "XbaI"],
"BamHI":["BclI", "BglII"],
"BclI":["BstYI", "BglII", "MboI"],
"BspEI":["BsrFI", "SgrAI", "AvaI", "XmaI", "BsrFI", "NgoMIV"],
"BssHII":["MluI", "AscI"],
"BstBI":["AccI", "ClaI", "AciI", "BsaHI", "HinP1I", "HpaII", "NarI"],
"ClaI":["AccI", "AciI", "AclI", "BsaHI", "HinP1I", "HpaII", "NarI"],
"MfeI":["ApoI", "EcoRI"],
"NarI":["AccI", "AclI","TaqI", "AciI", "HinP1I", "HpaII"],
"NdeI":["AseI", "BfaI", "Csp6I", "MseI"],
"NgoMIV":["BsaWI", "SgrAI", "AvaI", "XmaI", "BsaWI", "SgrAI"],
"NheI":["SpeI", "StyI"],
"NsiI":["BsiHKAI", "Bsp1286I", "PstI", "SbfI"],
"NspI":["SphI"],
"PciI":["BspHI", "FatI", "NcoI"],
"PpuMI":["RsrII"],
"PspOMI":["EaeI", "EagI", "NotI"],
"PspXI":["SalI"],
"PstI":["BsiHKAI", "Bsp1286I", "SbfI"],
"SacII":["BsiEI"],
"SalI":["XhoI"],
"SbfI":["BsiHKAI", "Bsp1286I"],
"SgrAI":["AvaI", "XmaI"],
"SpeI":["StyI", "XbaI"],
"XbaI":["NheI", "StyI"],
"XmaI":["BsaWI"]}
####### END INPUTS###################

#Create list of pairs
pairs_list = list()
for key in enzymes_dict:
    for value in enzymes_dict[key]:
        pairs_list.append((key,value))
        pairs_list.append((value,key))

#Fetch RE site sequences
file = '/shared/sudlab1/General/mirror/REBASE/db_res_enzymes'
re_db = open(file).readlines()
re_db = [s.strip() for s in re_db]

#Create dict of RE pairs and their sites
re_sites_pairs = dict()
for RE1,RE2 in pairs_list:
    RE_dbname1 = "<1>"+RE1
    RE_dbname2 = "<1>"+RE2
    for i in re_db:
        #Fetch site sequences
        if re.fullmatch(RE_dbname1, i):
            index_site1 = re_db.index(i)
            SITE1 = re_db[index_site1 + 4].split("<5>")[1]
        if re.fullmatch(RE_dbname2, i):
            index_site2 = re_db.index(i)
            SITE2 = re_db[index_site2 + 4].split("<5>")[1]
    re_sites_pairs[(RE1,RE2)] = (SITE1,SITE2)


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



ambiguous = ("R","Y","M","K","S","W","B","V","D","N","H")
linkers_per_pairs = dict()
for key in re_sites_pairs:
    site1 = re_sites_pairs[key][0]
    site2 = re_sites_pairs[key][1]
    if "^" not in site1 or "^" not in site2:
        print("The pair: ", key, re_sites_pairs[key],"is an exception and will have to be processed separately")
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
        print("The pair: ", key, re_sites_pairs[key],"is an exception and will have to be processed separately")
        continue
    if len(similar1) != len(similar2):
        print("The pair: ", key, re_sites_pairs[key],"is an exception and will have to be processed separately")
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
            print("The pair: ", key, re_sites_pairs[key],"is an exception and will have to be processed separately")
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
#Unnest
list_linkers = list(iOT.chain.from_iterable(linkers_per_pairs.values()))
#Replace T for U
linkers = list()
for l in list_linkers:
    if "T" in l:
        U_seq = l.replace("T", "U")
        l = U_seq
        linkers.append(l)
linkers = set(linkers)

#print output
output = open(name_ouput, "w")
for l in linkers:
    output.write(l + '\n')
output.close()

output2 = open(name_enzyme_RE_site_table, "w")
output2.write("RE1 name,RE2 name"+"\t"+"RE1 site, RE2site"+"\n")
for key in re_sites_pairs:
    output2.write(str(key[0])+","+str(key[1])+ "\t"+str(re_sites_pairs[key][0])+","+str(re_sites_pairs[key][1])+"\n")
