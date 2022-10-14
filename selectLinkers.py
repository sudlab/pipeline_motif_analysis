'''
select_linkers - summarize infomation about a meme motif file
=========================================================

:Author: Charlotte Vandermeulen
:Tags: Python

Purpose
-------

From a list of wanted linkers,
               highstability motifs,
               lowstabibility motifs,
               miRNA seed target sequences;
Process the number of time a linker recreates a lowstabibility motif or a
highstability motif (ecxept itself) or a miRNA seed target sequence, when it is
combined with each highstability motifs.

Also does it for all possible 6mer linkers to compare.

Usage
-----

Add inputs to script, then run the the script.

   python select_linkers.py
          -i highstab_motifs.list
          -l lowstab_motifs.list
          -m mirna_seed_sequences.list
          -r list_linkers.txt
          -o /path/to/output/directory

Output
------

Outputs txt files of each linker and the number of hits they generated
with : - lowstability motifs : RE_linkers_hits_lowstab.txt
       - miRNA seed motifs : RE_linkers_hits_seeds.txt
       - highstability motifs : RE_linkers_hits_highstab.txt
If some linkers didn't create any hits (almost impossible), a txt file listing
these linkers is outputted : RE_linkers_no_hits.txt

A txt file showing "stats" for the number of hits each 6mer linkers
with the different lists is created : stats_6mer_linkers.txt
'''

import sys
import cgatcore.experiment as E
import itertools as iOT
from collections import defaultdict
import re
import random

def checkInMotifs(sequence, motifs):
        '''Check if motifs are in sequence, returns true or false'''
        for i in motifs:
            if i in sequence:
                return True
        return False


def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $1.0$",
                            usage=globals()["__doc__"])

    parser.add_option("-i", "--highstab", dest="highstab", type=str,
                        help="Table of highstab motifs")
    parser.add_option("-l", "--lowstab", dest="lowstab", type=str,
                        help="Table of lowstab motifs to remove from highstab")
    parser.add_option("-m", "--miRNA-seeds", dest="mirseeds", type=str,
                        help="List of miRNA seed targets to avoid")
    parser.add_option("-r", "--RE-linkers", dest="linkers", type=str,
                        help="List of linkers based on RE cohesive sites")
    parser.add_option("-o", "--output-directory", dest="out_dir", type=str,
                        help="Full path directory where outputs are written")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.start(parser, argv=argv)

    high = open(options.highstab).readlines()
    low = open(options.lowstab).readlines()
    mirseed = open(options.mirseeds).readlines()
    linkers = open(options.linkers).readlines()

    #Screen wanted linkers
    high = set(s.strip() for s in high)
    low = set(s.strip() for s in low)
    mirseed = set(s.strip() for s in mirseed)
    linkers = set(l.strip() for l in linkers)

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
        out_OK = open(options.out_dir+"/RE_linkers_no_hits.txt", "w")
        for l in NoMatchLowSeed:
            out_ok.write(l+"\n")
        out_OK.close()

    out_RE_linkers_low = open(options.out_dir+"/RE_linkers_hits_lowstab.txt", "w")
    for l in MatchLow:
        out_RE_linkers_low.write(l+":"+str(MatchLow[l])+"\n")
    out_RE_linkers_low.close()

    out_RE_linkers_seed = open(options.out_dir+"/RE_linkers_hits_seeds.txt", "w")
    for l in MatchSeed:
        out_RE_linkers_seed.write(l+":"+str(MatchSeed[l])+"\n")
    out_RE_linkers_seed.close()

    out_RE_linkers_high = open(options.out_dir+"/RE_linkers_hits_highstab.txt", "w")
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
    #progress = 0
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

    out2 = open(options.out_dir+"/stats_6mer_linkers.txt", "w")
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
    #args.stdout
    #
    # write footer and output benchmark information.
    E.stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
