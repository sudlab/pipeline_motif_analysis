'''
cgat_script_template.py - template for cgat scripts
====================================================

:Author: Charlotte Vqndermeulen
:Tags: Python

Purpose
-------

From list of motif sequences, remove lowstab from highstab

Usage
-----

.. Example use case

Example::

   python FilterLists.py

Type::

   python FilterLists.py --help

for command line help.

Command line options
--------------------

'''

import sys
import cgatcore.experiment as E
import pandas as pd

def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $1.0$",
                            usage=globals()["__doc__"])

    parser.add_option("-l", "--lowstab", dest="lowstab", type=str,
                        help="Table of lowstab motifs to remove from STDIN")
    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.start(parser, argv=argv)

    highstab = pd.read_table(options.stdin, header = None)
    lowstab = pd.read_table(options.lowstab, header = None)[1]

    for m in lowstab[1]:
        inlow_highstab = highstab[highstab[1].str.contains(str(m)) == False]

    inlow_highstab.to_csv(options.stdout, sep ='\t', header = None)
    #args.stdout
    #
    # write footer and output benchmark information.
    E.stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
