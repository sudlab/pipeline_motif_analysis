'''
select_linkers - summarize infomation about a meme motif file
=========================================================

:Author: Charlotte Vandermeulen
:Tags: Python

Purpose
-------

From a list of linkers,
               highstability motifs,
               lowstabibility motifs,
               miRNA seed target sequences;
Process the number of time a linker recreate a lowstabibility motif or a
highstability motif (ecxept itself) or a miRNA seed target sequences when it is
combined with each highstability motifs.

Also does it for all possible 6mer linkers to compare.

Usage
-----

Add inputs to script, then run the the script.

   python select_linkers.py

Output
------

Outputs txt files of each linker and the number of hits they generated
with : - lowstability motifs : RE_linkers_hits_lowstab.txt
       - miRNA seed motifs : RE_linkers_hits_seeds.txt
       - highstability motifs : RE_linkers_hits_highstab.txt
If some linkers didn't create any hits (almost impossible), a txt file listing
these linkers is outputted : RE_linkers_no_hits.txt

Finally, a txt file showing "stats" for the number of hits each 6mer linkers
with the different lists is created : stats_6mer_linkers.txt
'''

import itertools as iOT
from collections import defaultdict
import re
import random

#######START INPUTS#####
high = open("/shared/sudlab1/General/projects/SynthUTR_hepG2_a549/test_linkers/highstab.list").readlines()
low = open("/shared/sudlab1/General/projects/SynthUTR_hepG2_a549/test_linkers/lowstab.list").readlines()
mirseed = open("/shared/sudlab1/General/projects/SynthUTR_hepG2_a549/miRmine/relevant_miRNA_a549_hepg2.seeds.list").readlines()
linkers = open("/shared/sudlab1/General/projects/SynthUTR_hepG2_a549/motif_merge_linkers/linkers_test2/list_linkers.txt").readlines()
out_directory = "/shared/sudlab1/General/projects/SynthUTR_hepG2_a549/motif_merge_linkers/linkers_test2"
####### END INPUTS#####



#Screen wanted linkers
high = set(s.strip() for s in high)
low = set(s.strip() for s in low)
mirseed = set(s.strip() for s in mirseed)
linkers = set(l.strip() for l in linkers)

def checkInMotifs(sequence, motifs):
    '''Check if motifs are in sequence, returns true or false'''
    for i in motifs:
        if i in sequence:
            return True
    return False

#filter linkers that are in motifs low, high, or mir seeds
#Also stretches of 6 same nucleotide
filter_linkers = set()
for l in linkers:
    if str(high).find(l) != -1:
        filter_linkers.add(l)
    if str(low).find(l) != -1:
        filter_linkers.add(l)
    if str(mirseed).find(l) != -1:
        filter_linkers.add(l)
    if l in {"AAAAAAAA","GGGGGGGG", "CCCCCCCC", "UUUUUUUU"} :
        filter_linkers.add(l)
linkers = linkers - (filter_linkers)

print("After filtering, there are ",len(linkers), " linkers.")

#Screen
print("Starting screening of recombined linkers")
highIterator = iOT.product(high, repeat = 2)
MatchLow = defaultdict(int)
MatchSeed = defaultdict(int)
MatchHigh = defaultdict(int)

for i,j in highIterator:
    missingList = high -set([i,j])
    for l in linkers:
        sequence = l + i + l
        if checkInMotifs(sequence, low):
            MatchLow[l] += 1
        if checkInMotifs(sequence, missingList):
            MatchHigh[l] += 1
        if checkInMotifs(sequence, mirseed):
            MatchSeed[l] += 1

NoMatchLowSeed = [l for l in linkers if l not in MatchLow.keys() | MatchSeed.keys()]
print("Finished screening of recombined linkers")

#Create ouputs for linkers
if len(NoMatchLowSeed) != 0:
    print("Following linkers didn't create a match with low and seed motifs:")
    print(NoMatchLowSeed)
    out_OK = open(out_directory+"/RE_linkers_no_hits.txt", "w")
    for l in NoMatchLowSeed:
        out_ok.write(l+"\n")
    out_OK.close()

out_RE_linkers_low = open(out_directory+"/RE_linkers_hits_lowstab.txt", "w")
for l in MatchLow:
    out_RE_linkers_low.write(l+":"+str(MatchLow[l])+"\n")
out_RE_linkers_low.close()

out_RE_linkers_seed = open(out_directory+"/RE_linkers_hits_seeds.txt", "w")
for l in MatchSeed:
    out_RE_linkers_seed.write(l+":"+str(MatchSeed[l])+"\n")
out_RE_linkers_seed.close()

out_RE_linkers_high = open(out_directory+"/RE_linkers_hits_highstab.txt", "w")
for l in MatchHigh:
    out_RE_linkers_high.write(l+":"+str(MatchHigh[l])+"\n")
out_RE_linkers_high.close()


##Test for all possible 6mer linkers
#Redo linkers list
linkers = iOT.product("AUCG", repeat = 6)
linkers = set("".join(l) for l in linkers)
filter_linkers = set()
for l in linkers:
    if str(high).find(l) != -1:
        filter_linkers.add(l)
    if str(low).find(l) != -1:
        filter_linkers.add(l)
    if str(mirseed).find(l) != -1:
        filter_linkers.add(l)
    if l in {"AAAAAA","GGGGGG", "CCCCCC", "UUUUUU"} :
        filter_linkers.add(l)
linkers = linkers - (filter_linkers)
#linkers = random.sample(linkers, 100)

#Screen
print("Starting screening of all possible 6mer linkers")
highIterator = iOT.product(high, repeat = 2)
MatchLow = defaultdict(int)
MatchSeed = defaultdict(int)
MatchHigh = defaultdict(int)
progress = 0
for i,j in highIterator:
    missingList = high -set([i,j])
    for l in linkers:
        sequence = l + i + l
        if checkInMotifs(sequence, low):
            MatchLow[l] += 1
        if checkInMotifs(sequence, missingList):
            MatchHigh[l] += 1
        if checkInMotifs(sequence, mirseed):
            MatchSeed[l] += 1
print("Finished screening of all possible 6mer linkers")
NoMatch = [l for l in linkers if l not in MatchLow.keys() & MatchSeed.keys()]

#Write ouputs
mean_MatchLow = sum(MatchLow.values()) / len(MatchLow)
mean_MatchSeed = sum(MatchSeed.values()) / len(MatchSeed)
mean_MatchHigh = sum(MatchHigh.values()) / len(MatchHigh)

out2 = open(out_directory+"/stats_6mer_linkers.txt", "w")
out2.write("Max number of hits for a lowstab motif (6mer linkers): "+str(max(MatchLow.values()))+"\n")
out2.write("Min number of hits for a lowstab motif (6mer linkers): "+ str(min(MatchLow.values()))+ "\n")
out2.write("Mean number of hits for a lowstab motif (6mer linkers): "+ str(mean_MatchLow)+ "\n")

out2.write("Max number of hits for a seed motif (6mer linkers): "+ str(max(MatchSeed.values()))+ "\n")
out2.write("Min number of hits for a seed motif (6mer linkers): "+ str(min(MatchSeed.values()))+ "\n")
out2.write("Mean number of hits for a seed motif (6mer linkers): "+ str(mean_MatchSeed)+ "\n")

out2.write("Max number of hits for a highstab motif (6mer linkers): "+ str(max(MatchHigh.values()))+ "\n")
out2.write("Min number of hits for a highstab motif (6mer linkers): "+ str(min(MatchHigh.values()))+ "\n")
out2.write("Mean number of hits for a highstab motif (6mer linkers): "+ str(mean_MatchHigh)+ "\n")

out2.close()
