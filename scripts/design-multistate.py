#!/usr/bin/env python3
from __future__ import print_function

import argparse
import sys
import os
import re
import timeit
from collections import Counter
from scipy import stats
import numpy as np

import RNA # ViennaRNA python bindings
from RNARedPrintSampler import RPSampler,gccontent
from Structure import RNAStructure

from itertools import islice
def take(n, iterable):
    "Return first n items of the iterable as a list"
    return list(islice(iterable, n))


def read_input(content):
    '''
    Reads some input and returns all neccessary information in the right container.
    Input is a string, lines separated by linebreaks. Content might be structures,
    a sequence constraint and a start sequence. Additional information for the structural
    states can be provided with any separator ;,: or whitespaces after the structure.

    :param filename: Filename of the file to read
    :return: structures - List of structures in dot-bracket notation
    :return: constraint - Sequence constraint
    :return: sequence - Start sequence
    :return: additions - List of additional information after the structures
    '''
    structures = []
    constraint = ''
    sequence = ''
    additions = []

    lines = content.split("\n")
    for line in lines:
        # if line begins with a semicolon ; stop parsing
        if re.match(re.compile("^\;"), line, flags=0):
            break
        # strip additional information after the structure/sequence string
        m = re.match(re.compile("^([^\s\;\,\:]+)[\s\;\,\:]*(.*)$"), line, flags=0)
        if m:
            line, addition = m.groups()
        if re.match(re.compile("^[\(\)\.\{\}\[\]\<\>\+\&]+$"), line, flags=0):
            structures.append(line.rstrip('\n'))
            additions.append(addition.rstrip('\n'))
        elif re.match(re.compile("^[\ ACGTUWSMKRYBDHVN\&\+]+$", re.IGNORECASE), line, flags=0):
            line = line.replace(" ", "N")
            line = line.upper();
            if re.match(re.compile("^[ACGTU\&\+]+$", re.IGNORECASE), line, flags=0) and sequence == '':
                sequence = line.rstrip('\n')
            elif constraint == '':
                constraint = line.rstrip('\n')
            else:
                raise ValueError('Too many constraints or start sequences: ' + line)

    checklength = len(structures[0])
    if constraint != '' and len(constraint) != checklength:
        raise ValueError('Structures and the sequence constraint must have the same length!')
    elif sequence != '' and len(sequence) != checklength:
        raise ValueError('Structures and the start sequence must have the same length!')

    for s in structures:
        if len(s) != checklength:
            raise ValueError('Structures must all have the same length!')

    return structures, constraint, sequence, additions

def calc_turner_energies(seq,structures,temperature):
    md = RNA.md()
    md.temperature = temperature
    fc = RNA.fold_compound(seq, md)
    return [ fc.eval_structure(s) for s in structures ]


def main():
    parser = argparse.ArgumentParser(description='Design RNA molecules which adopt multiple structural states with equal energies using multi-dimensional Boltzmann sampling.',
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', action='version', version='%(prog)s of RNARedPrint 0.3')
    parser.add_argument("-i", "--input", type=argparse.FileType('r'), default=sys.stdin, help="Read structures from input file. Default: read from stdin. Format must be dot-bracket structures, each per one line with a trailing line containing only a semi-colon.")
    parser.add_argument("-T", "--temperature", type=float, default=37.0, help='Temperature of the energy calculations.')
    parser.add_argument("-n", "--number", type=int, default=1000, help='Number of designs to generate')
    parser.add_argument("-m", "--model", type=str, default='basepairs', help='Model for getting a new sequence: uniform, nussinov, basepairs, stacking')
    parser.add_argument("-g", "--gc", type=float, default=0.5, help='Target GC content.')
    parser.add_argument("-t", "--simple_tolerance", type=float,
            default=0.30, help='Tolerated relative deviation of simple energies.')
    parser.add_argument("--tolerance", type=float, default=1,
            help='Turner energy tolerance of energies (default: 1 kcal/mol, turn off by values <0)')
    parser.add_argument("-c", "--gctolerance", type=float, default=0.1, help='Tolerated relative deviation to target GC content.')
    parser.add_argument("--csv_output", default=False, action='store_true', help='Output csv format (with additional information)')
    parser.add_argument("-d", "--debug", default=False, action='store_true', help='Show debug information of library')
    args = parser.parse_args()

    if args.tolerance < 0:
        args.tolerance = None

    if args.debug:
        print("# Options: number={0:d}, model={1:}, temperature={2:}".format(args.number, args.model, args.temperature))

    data = ''
    for line in args.input:
        data = data + '\n' + line
    (structures, constraint, start_sequence, _) = read_input(data)

    # remove lonely pairs
    structures = [RNAStructure(s).removeLonelyPairs() for s in structures]

    # print input
    if args.debug:
        print("# " + "\n# ".join(structures) + "\n# " + constraint)

    # print header for csv file
    if args.csv_output:
        print(";".join(["sequence",
                    "model",
                    "construction_time",
                    "sample_time"] +
                    ["redprint_energy"]*len(structures) +
                    ["turner_energy"]*len(structures)
                )
        )

    # time the sampling
    time_start = timeit.default_timer()

    nstr = len(structures)
    wastefactor = 20
    stacksize = (wastefactor*args.number)
    if stacksize < 1000:
        stacksize = 1000
    sampler = RPSampler(structures, model=args.model, temperature=args.temperature, stacksize=stacksize, StopConstruct=True, debug=args.debug)
    construction_time = sampler.construction_time
    # get first sequence sample with energies with GC weight optimization and fixed high weights
    initialsample = take(1000,
                         Sample(sampler, nstr, target_energies=None,
                                target_GC=args.gc, args=args))

    # get target energies
    target_energies = getTargetEnergy(structures, initialsample, args)

    target_turner_energies = target_energies[:]

    if args.debug:
        print("# Turner Target Energies are: ", target_energies)
    # get energy offsets
    slope, intercept = getEnergyOffsets(structures, initialsample, args)
    # correct target energies with offsets
    for t in range(0, len(structures)):
        target_energies[t]  = (target_energies[t] - intercept[t]) / slope[t]
    if args.debug:
        print("# Simple Target Energies are: ", target_energies)

    #sampler.weights = [1.0] * nstr
    AdmissibleSample = Sample(sampler, nstr,
            target_energies=target_energies,
            target_turner_energies=target_turner_energies, structures=structures, target_GC=args.gc,
            simple_tolerance = args.simple_tolerance, tolerance = args.tolerance, target_GC_eps=args.gctolerance, args=args)

    admissible_count = 0
    for a in AdmissibleSample:
        admissible_count += 1

        md = RNA.md()
        md.temperature = args.temperature
        fc = RNA.fold_compound(a['seq'], md)
        turner_energies = []
        for s in structures:
            turner_energies.append(fc.eval_structure(s))

        #print('$ simple model: ', a['energies'], ' viennaRNA: ', fc.eval_structure(str(i)))
        if not args.csv_output:
            print(a['seq'],"GC={:.2f}".format(gccontent(a['seq']))," ".join(["E{}={:.2f}".format(i+1,e) for (i,e) in enumerate(turner_energies)]))
        else:
            print(a['seq'],
                    "\"" + args.model + "\"",
                    construction_time,
                    timeit.default_timer() - time_start,
                    ";".join([str(a['energies'][e]) for e in sorted(a['energies'].keys())]),
                    ";".join([str(e) for e in turner_energies]),
                    sep=";")
        # return if large enough
        if args.debug:
            print('# Found addmissible Sequences: ',
                    admissible_count/float(args.number), '%')
        if admissible_count >= args.number:
            break

def getTargetEnergy(structures, sample, args):
    nstr = len(structures)
    # find target energy_step
    energies = np.zeros( (len(sample), nstr) )
    target_energies= list()

    for i, s in enumerate(sample):
        for t in range(0, nstr):
            # calculate offset between turner eos and simple model eos
            energies[i,t] = s['energies'][t]

    mean_target = np.mean(energies)

    return [mean_target] * nstr

def getEnergyOffsets(structures, sample, args):
    nstr = len(structures)
    simple = np.zeros( (len(sample), nstr) )
    turner = np.zeros( (len(sample), nstr) )
    # iterate over sample
    for i, s in enumerate(sample):
        md = RNA.md()
        md.temperature = args.temperature
        fc = RNA.fold_compound(s['seq'], md)
        turner_energies = []
        for structure in structures:
            turner_energies.append(fc.eval_structure(structure))
        # iterate over structures
        for t in range(0, nstr):
            # calculate offset between turner eos and simple model eos
            simple[i,t] = s['energies'][t]
            turner[i,t] = turner_energies[t]
    turner = np.where(turner > 1000, np.nan, turner)
    # get linear regression
    slope = {}
    intercept = {}
    for t in range(0, nstr):
        varx=simple[:,t]
        vary=turner[:,t]
        mask = ~np.isnan(varx) & ~np.isnan(vary)
        slope[t], intercept[t], r_value, p_value, std_err = stats.linregress(varx[mask],vary[mask])

    if args.debug:
        print('# Slopes and intercepts are: ', slope, intercept)
    return slope, intercept

def Sample(sampler, nstr, target_GC, args, *, target_energies = None,
        target_turner_energies = None, structures = None, simple_tolerance = 0.1,
        tolerance = None, target_GC_eps=0.05, maxiterations=20):
    '''If target energies is none, we will sample with the initial weights and just
    adopt the GC weight.

    returns: iterator of dictionarys like:
      {'seq': sequence,
       'energies': [ energy0, ...] ,
       'turner_energies': [ turner_energy0, ]}
    '''
    # weights = [math.exp(1/((args.temperature + 273.15)*0.00198717))] * nstr

    # count empty iterations
    empty_iteration_count = 0
    while empty_iteration_count < maxiterations:
        empty_iteration_count += 1
        # get new sequences
        if args.debug:
            print("# Weights: ", sampler.weights)
            print('# GC weight: ', sampler.gcweight)
        newsample, energies = sampler.dump_new_stack()

        # get average structue energies for newsample
        eos = np.zeros( (len(newsample), nstr) )
        GC_freq = list()
        turner_energies = list()

        for i, s in enumerate(newsample):
            # count GC content
            c = Counter(s)
            sigma = 2.0 # one per nucleotide, laplace
            GC = (c['G'] + c['C'] + sigma) / (len(s) + 2*sigma)
            GC_freq.append(GC)

            # add if it is eps-admissible
            admissible = True
            if not (1-target_GC_eps <= GC/target_GC <= 1+target_GC_eps):
                admissible = False

            if target_energies:
                for t in range(0, nstr):
                    # add to eos np array in any case
                    eos[i,t] = energies[i][t]
                    # check if eps admissible
                    #print(eos[i,t], eos[i,t]/target_energies[t], target_energies[t])
                    if not (1-simple_tolerance <= eos[i,t]/target_energies[t] <= 1+simple_tolerance):
                        admissible = False

            te = None
            if admissible and structures is not None:
                te = calc_turner_energies(s,structures,args.temperature)
            turner_energies.append( te )

            if admissible and tolerance is not None:
                for t in range(0, nstr):
                    if abs( turner_energies[i][t] - target_turner_energies[t] ) > tolerance:
                        admissible = False

            if admissible:
                #print('# is admissible:', eos[i,:], GC)
                yield {'seq': s, 'energies': energies[i],
                        'turner_energies': turner_energies[i]}
                empty_iteration_count = 0

        # update weights
        if target_energies:
            for t in range(0, nstr):
                e_mean = np.mean(eos[:,t])
                if args.debug:
                    print('# Energy mean: ', str(t), e_mean)
                # exp version
                sampler.weights[t] = sampler.weights[t] * (1.1**(e_mean-target_energies[t]))
                # Yann old version without positive e_mean
                #weights[t] = weights[t]*target_energies[t]/e_mean
        # update gcweight
        GC_mean = np.mean(GC_freq)
        sampler.gcweight = sampler.gcweight * target_GC/GC_mean
        if args.debug:
            print('# GC mean: ', GC_mean)

if __name__ == "__main__":
    main()
