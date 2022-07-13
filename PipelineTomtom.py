import sys
import os
import sqlite3
import csv
import re
import glob
import pandas

from cgatcore import pipeline as P
from cgatcore import experiment as E
from HomerMemeMotifSubclass import HomerMotif, HomerMotifFile
from cgat.MEME import MotifCluster
import cgatcore.iotools as IOTools



def getSeedMotifs(motif_file, tomtom_file, outfile):

    ungrouped = HomerMotifFile(IOTools.open_file(motif_file))

    if len(ungrouped) == 0:
        with IOTools.open_file(outfile, "w") as outf:
            outf.write(str(ungrouped))
        return

    E.debug("%s: Loaded %i motifs" % (motif_file, len(ungrouped)))
    tomtom = pandas.read_csv(tomtom_file, sep="\t")
    tomtom["Query_ID"] = tomtom["Query_ID"].astype(str)
    tomtom["Target_ID"] = tomtom["Target_ID"].astype(str)

    all_clusters = ungrouped.keys()
    new_index = pandas.MultiIndex.from_product([all_clusters, all_clusters],
                                               names=["#Query_ID", "Target_ID"])
    tomtom = tomtom.set_index(["Query_ID", "Target_ID"])
    tomtom = tomtom.reindex(new_index)
    tomtom = tomtom.sort_index()

    ungrouped.sort("evalue")

    groups = []

    E.debug("%s: Clustering Motifs" % motif_file)
    while len(ungrouped) > 0:
        cur_cluster = MotifCluster(ungrouped.take(0))
        assert len(cur_cluster) == 1
        E.debug("%s: working on cluster %s, %i clusters remaining" %
                (motif_file, cur_cluster.seed.primary_id, len(ungrouped)))

        try:
            seed_distances = tomtom.loc[cur_cluster.seed.primary_id]
        except:
            print(tomtom)
            raise

        close_motifs = seed_distances[seed_distances["E-value"] < 0.05]

        E.debug("%s: Found %i similar motifs" %
                (motif_file, close_motifs.shape[0]))

        for motif in close_motifs.index.values:
            if motif == cur_cluster.seed.primary_id:
                continue

            try:
                cur_cluster.append(ungrouped.take(str(motif)))
            except KeyError:
                pass

        groups.append(cur_cluster)

    E.debug("%s: Got %i clusters, containing a total of %i motifs"
            % (motif_file, len(groups), sum(len(cluster) for cluster in groups)))

    groups.sort(key=lambda cluster: cluster.seed.evalue)

    merged_groups = []

    E.debug("%s: Merging groups with weak similarity" % motif_file)
    while len(groups) > 0:

        cur_cluster = groups.pop(0)

        to_merge = []

        distances = tomtom.loc[cur_cluster.seed.primary_id]
        for other_cluster in groups:

            evals = distances.loc[other_cluster.keys()]["E-value"]

            if (evals < 0.1).all():
                to_merge.append(other_cluster)

        for cluster in to_merge:
            cur_cluster.extend(cluster)
            groups.remove(cluster)

        merged_groups.append(cur_cluster)

    E.debug("%i final clusters found" % len(merged_groups))
    E.debug("%s bulding output" % motif_file)
    for group in merged_groups:
        group.seed.letter_probability_line += " nClustered= %i" % len(group)
        group.seed.letter_probability_line += " totalHits= %i" % sum(
            [motif.nsites for motif in group])

    output = HomerMotifFile(ungrouped)
    output.extend(cluster.seed for cluster in merged_groups)

    E.debug("%s outputting" % motif_file)
    with IOTools.open_file(outfile, "w") as outf:
        outf.write(str(output))
