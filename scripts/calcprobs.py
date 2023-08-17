#!/usr/bin/env python3

import argparse
import RNA # ViennaRNA python bindings
import sys

def report_probs(seq,structures,md):
    fc = RNA.fold_compound(seq, md)
    _, MFE = fc.mfe()
    _, EE  = fc.pf()

    probs = [ fc.pr_structure(s) for s in structures ]

    probstrings = [ "P{}={:0.02f}".format(i+1, p) for i,p in enumerate(probs) ]

    return "MFE={:0.02f} ".format(MFE) +  "EE={:0.02f} ".format(EE) + " ".join(map(str, probstrings)) + " Psum={:0.02f}".format(sum(probs))

def main():
    parser = argparse.ArgumentParser(description='Design RNA molecules which adopt multiple structural states with specific energies using multi-dimensional Boltzmann sampling.',
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', action='version', version='%(prog)s of RNARedPrint 0.3')
    parser.add_argument("-i", "--input", type=argparse.FileType('r'), help="Read structures from input file. Format must be dot-bracket structures, each per one line with a trailing line containing only a semi-colon.")
    parser.add_argument("-T", "--temperature", type=float, default=37.0, help='Temperature of the energy calculations.')

    args = parser.parse_args()

    md = RNA.md()
    md.temperature = args.temperature
    md.uniq_ML = 1

    structures = [ s.strip() for s in args.input ]

    print(structures)

    for line in sys.stdin:
        line = line.strip()
        seq = line.split(" ")[0]
        print(line, report_probs(seq,structures, md))


if __name__ == "__main__":
    main()
