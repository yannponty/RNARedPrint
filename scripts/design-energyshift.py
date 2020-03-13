#!/usr/bin/env python3
from __future__ import print_function

import argparse
import sys
import os
import re
import timeit
from scipy import stats
import numpy as np
import RNA # ViennaRNA python bindings

from Structure import RNAStructure
from RNARedPrintSampler import RPSampler,gccontent

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
    parser = argparse.ArgumentParser(description='Design RNA molecules which adopt multiple structural states with specific energies using multi-dimensional Boltzmann sampling.',
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--input", type=argparse.FileType('r'), default="-", help="Read structures from input file. Default: read from stdin. Format must be dot-bracket structures, each per one line with a trailing line containing only a semi-colon.")
    parser.add_argument("-T", "--temperature", type=float, default=37.0, help='Temperature of the energy calculations.')
    parser.add_argument("-n", "--number", type=int, default=1000, help='Number of designs to generate')
    parser.add_argument("-m", "--model", type=str, default='basepairs', help='Model for getting a new sequence: uniform, nussinov, basepairs, stacking')
    parser.add_argument("-e", "--energies", type=str, default='', help='Target Energies for design. String of comma separated float values.')
    parser.add_argument("-g", "--gc", type=float, default=0.5, help='Target GC content.')
    parser.add_argument("-t", "--simple_tolerance", type=float,
            default=0.2, help='Tolerated relative deviation to target energies in the simple model.')
    parser.add_argument("-c", "--gctolerance", type=float, default=0.05, help='Tolerated relative deviation to target GC content.')
    parser.add_argument("--tolerance", type=float, default=1,
            help='Turner energy tolerance to the target energies (default: 1 kcal/mol, turn off by values <0)')
    parser.add_argument("--csv_output", default=False, action='store_true', help='Output csv format (with additional information)')
    parser.add_argument("-d", "--debug", default=False, action='store_true', help='Show debug information of library')
    args = parser.parse_args()

    if args.tolerance < 0:
        args.tolerance = None

    if args.debug:
        print("# Options: number={0:d}, model={1:}, temperature={2:}".format(args.number, args.model, args.temperature))
    # define structures
    data = ''
    for line in sys.stdin:
        data = data + '\n' + line
    (structures, constraint, start_sequence, _) = read_input(data)

    # time the sampling
    start = timeit.default_timer()

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

    # read target energies
    target_energies = []
    if args.energies:
        for i, w in enumerate(args.energies.split(',')):
            target_energies.append(float(w))
    else:
        exit(1)

    target_turner_energies = target_energies[:]

    if args.debug:
        print("# Turner Target Energies are: ", target_turner_energies)
    # get energy offsets
    slope, intercept = getEnergyOffsets(structures, args)
    # correct target energies with offsets
    for t in range(0, len(structures)):
        target_energies[t]  = (target_energies[t] - intercept[t]) / slope[t]

    if args.debug:
        print("# Simple Target Energies are: ", target_energies)

    nstr = len(structures)
    wastefactor = 20
    sampler = RPSampler(structures, model=args.model, weights=([1.0] * nstr), gcweight=1.0,
                temperature=args.temperature, stacksize=(wastefactor*args.number))

    admissible_count = 0
    for a in Sample(sampler, nstr, target_energies,
                target_turner_energies = target_turner_energies,
                structures=structures,
                target_GC=args.gc, simple_tolerance = args.simple_tolerance,
                target_GC_eps=args.gctolerance, args=args):
        admissible_count += 1

        # output sequence
        if not args.csv_output:
            print(a['seq'],"GC={:.2f}".format(gccontent(a['seq']))," ".join(["E{}={:.2f}".format(i+1,e) for (i,e) in
                        enumerate(a['turner_energies'])]))
        else:
            print(a['seq'],
                "\"" + args.model + "\"",
                sampler.construction_time,
                timeit.default_timer() - start, # sample time until now
                ";".join([str(a['energies'][e]) for e in sorted(a['energies'].keys())]),
                ";".join([str(e) for e in a['turner_energies']]), sep=";"
            )

        # terminate if enough samples generated
        if args.debug:
            print('# Found for current Target: ',
                    admissible_count/float(args.number)*100, '%')
        if admissible_count >= args.number:
            break

    if admissible_count < args.number:
        print("#! Terminate unsuccessfully (maximum number of iterations exceeded)")

def getEnergyOffsets(structures, args):
    sampler = RPSampler(structures, model=args.model, temperature=args.temperature, stacksize=1000, StopConstruct=True, debug=args.debug)

    # get new sequecne
    newsample, energies = sampler.dump_new_stack()

    nstr = len(structures)
    simple = np.zeros( (len(newsample), nstr) )
    turner = np.zeros( (len(newsample), nstr) )
    # iterate over sample
    for i, s in enumerate(newsample):
        # get design object
        md = RNA.md()
        md.temperature = args.temperature
        fc = RNA.fold_compound(s, md)
        # iterate over structures
        for t, structure in enumerate(structures):
            # calculate offset between turner eos and simple model eos
            simple[i,t] = energies[i][t]
            turner[i,t] = fc.eval_structure(structure)
    # get linear regression
    slope = {}
    intercept = {}
    r_value = {}
    for t in range(0, nstr):
        #varx=simple[:,t] # > varx[mask]
        #vary=turner[:,t] # > vary[mask]
        #mask = ~np.isnan(varx) & ~np.isnan(vary)
        slope[t], intercept[t], r_value[t], p_value, std_err = stats.linregress(simple[:,t],turner[:,t])

    if args.debug:
        plotRegression(simple, turner, slope, intercept, r_value)
        print('# Slopes, intercepts and r-value are: ', slope, intercept, r_value)
    return slope, intercept

def plotRegression(simple, turner, slope, intercept, r_value):
    import matplotlib.pyplot as plt
    import matplotlib.pylab as pylab
    import matplotlib.cm as cm

    # make histogram plot
    params = {'legend.fontsize': 'x-large',
             'figure.figsize': (7, 5),
             'axes.labelsize': 'x-large',
             'axes.titlesize':'x-large',
             'xtick.labelsize':'x-large',
             'ytick.labelsize':'x-large'}
    pylab.rcParams.update(params)
    #fig, axes = plt.subplots(3, 1, sharex=False, sharey=False)
    fig = plt.figure()
    plt.xlabel("Energy of Structure in Simple Model [kcal/mol]")
    plt.ylabel("Energy of Structure in Turner Model [kcal/mol]")
    fig.legend().set_visible(False)
    color_i = 0.0

    for t in slope.keys():
        plt.plot(simple[:,t], turner[:,t], '.', color=cm.viridis(color_i), label=str(t), alpha=0.5)
        plt.plot(simple[:,t], slope[t] * simple[:,t]+ intercept[t], '-', color=cm.viridis(color_i))
        fig.text(0.16, 0.81-color_i/7, "$R^{2}$ = "+"{0:.3f}".format(r_value[t]), color=cm.viridis(color_i))
        # increment color
        color_i += 0.3
    #plt.legend(loc='upper left')
    fig.savefig('regression.svg', dpi=300)

def Sample(sampler, nstr, target_energies, target_GC, args,
        *,
        target_turner_energies,
        structures,
        simple_tolerance = 0.1, target_GC_eps=0.05, maxiterations=25):
    # weights = [math.exp(1/((args.temperature + 273.15)*0.00198717))] * nstr

    # count iteration
    empty_iteration_count = 0 # count iterations that don't produce anything
    while empty_iteration_count < maxiterations:
        empty_iteration_count += 1
        # get new sequences
        if args.debug:
            print("# Weights: ", sampler.weights)
            print('# GC weight: ', sampler.gcweight)
        newsample, energies = sampler.dump_new_stack()

        # get average structure energies for newsample
        eos = np.zeros( (len(newsample), nstr) )
        GC_freq = list()
        turner_energies = list()

        def admissible(s, GC, eos, turner_energies):
            if not (1-target_GC_eps <= GC/target_GC <= 1+target_GC_eps):
                return False

            # check if eps admissible
            if any(not (1-simple_tolerance <= eos[t]/target_energies[t] <= 1+simple_tolerance)
                    for t in range(0, nstr)):
                return False

            if args.tolerance is not None:
                if any(abs(e-t) > args.tolerance
                       for (e,t) in
                          zip(turner_energies,target_turner_energies)):
                    return False
            return True

        # record features for statistics
        for i, s in enumerate(newsample):
            gc = gccontent(s)
            GC_freq.append(gc)
            for t in range(0, nstr):
                # add to eos np array in any case
                eos[i,t] = energies[i][t]

            turner_energies.append(
                    calc_turner_energies(s,structures,args.temperature) )

            if admissible(newsample, gc, energies[i], turner_energies[i]):
                #print('# is admissible:', eos[i,:], GC)
                yield ({'seq': s, 'energies': energies[i],
                    'turner_energies': turner_energies[i]})
                empty_iteration_count = 0

        # update weights
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
        if args.debug:
            print('# GC mean: ', GC_mean)
        sampler.gcweight = sampler.gcweight * target_GC/GC_mean

if __name__ == "__main__":
    main()
