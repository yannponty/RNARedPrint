#!/usr/bin/env python3

import argparse
import RNA # ViennaRNA python bindings
import sys
from math import log,exp

def RT():
    return (1.987 / 1000) * 310.15

## convert energy to PF 
def bw(e):
    return exp(-e/RT())

def report_probs(seq,structures):
    MFE = RNA.fold(seq)[1]
    EE  = RNA.pf_fold(seq)[1]

    Zs = [ RNA.energy_of_struct(seq,s) - EE for s in structures ]

    probs = [ bw(z) for z in Zs ]

    probstrings = [ "P{}={:0.02f}".format(i+1, p) for i,p in enumerate(probs) ]

    return "MFE={:0.02f} ".format(MFE) +  "EE={:0.02f} ".format(EE) + " ".join(map(str, probstrings)) + " Psum={:0.02f}".format(sum(probs))

def main():
    parser = argparse.ArgumentParser(description='Design RNA molecules which adopt multiple structural states with specific energies using multi-dimensional Boltzmann sampling.',
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', action='version', version='%(prog)s of RNARedPrint 0.3')
    parser.add_argument("-i", "--input", type=argparse.FileType('r'), help="Read structures from input file. Format must be dot-bracket structures, each per one line with a trailing line containing only a semi-colon.")

    args = parser.parse_args()

    structures = [ s.strip() for s in args.input ]

    print(structures)

    for line in sys.stdin:
        line = line.strip()
        seq = line.split(" ")[0]
        print(line, report_probs(seq,structures))    


if __name__ == "__main__":
    main()

