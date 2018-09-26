#!/usr/bin/env python3
from __future__ import print_function

import argparse
import re
import sys
import os
import timeit
import RNA # ViennaRNA python bindings
from collections import Counter
import numpy as np

from Structure import RNAStructure
from RNARedPrintSampler import RPSampler


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


def main():
    parser = argparse.ArgumentParser(description='Design a multi-stable riboswitch similar using Boltzmann sampling with equal target energies.')
    parser.add_argument("-T", "--temperature", type=float, default=37.0, help='Temperature of the energy calculations.')
    parser.add_argument("-n", "--number", type=int, default=1000, help='Number of designs to generate')
    parser.add_argument("-m", "--model", type=str, default='stacking', help='Model for getting a new sequence: uniform, nussinov, basepairs, stacking')
    parser.add_argument("-g", "--gc", type=float, default=0.5, help='Target GC content.')
    parser.add_argument("-o", "--eps", type=float, default=0.1, help='Offset eps for target energies.')
    parser.add_argument("-p", "--gceps", type=float, default=0.05, help='Offset eps for GC content.')
    parser.add_argument("-d", "--debug", default=False, action='store_true', help='Show debug information of library')
    args = parser.parse_args()

    if (args.debug):
        print("# Options: number={0:d}, model={1:}, temperature={2:}".format(args.number, args.model, args.temperature))

    data = ''
    for line in sys.stdin:
        data = data + '\n' + line
    (structures, constraint, start_sequence, _) = read_input(data)

    # time the sampling
    start = timeit.default_timer()

    # remove lonely pairs
    structures = [RNAStructure(s).removeLonelyPairs() for s in structures]

    # print input
    if (args.debug):
        print("# " + "\n# ".join(structures) + "\n# " + constraint)

    # print header for csv file
    print(";".join(["sequence",
                "model",
                "construction_time",
                "sample_time"] +
                ["redprint_energy"]*len(structures) +
                ["turner_energy"]*len(structures)
            )
    )

    target_energies, offsets, construction_time = getTargetEnergy(structures, args)
    if (args.debug):
        print("# Target Energies are: ", target_energies)
    bs = BalancedSamples(structures, target_energies, offsets, energy_step=args.gc, args=args)

    count = 0
    for b in sorted(bs.keys()):
        if count > args.number:
            break
        count += 1

        md = RNA.md()
        md.temperature = args.temperature
        fc = RNA.fold_compound(bs[b]['seq'], md)
        turner_energies = []
        for s in structures:
            turner_energies.append(fc.eval_structure(s))

        # output sequence
        print(bs[b]['seq'],
            "\"" + args.model + "\"",
            construction_time,
            timeit.default_timer() - start, # sample time until now
            ";".join([str(bs[b]['energies'][e]) for e in sorted(bs[b]['energies'].keys())]),
            ";".join([str(e) for e in turner_energies]), sep=";"
        )

def getTargetEnergy(structures, args):
    sampler = RPSampler(structures, model=args.model, weights=[1]*len(structures), temperature=args.temperature, stacksize=1000, StopConstruct=True, debug=args.debug)

    # get new sequecne
    newsample, energies = sampler.dump_new_stack()
    # get energy offsets
    offsets = getEnergyOffsets(structures, newsample, energies, args)
    # find target energy_step
    phi = 99999.9
    target_energies= []
    for i, s in enumerate(newsample):
        current_phi = getPhi(energies[i], offsets)
        if current_phi < phi:
            phi = current_phi
            target_energies = [np.mean(list(energies[i].values()) + list(offsets.values()))]*len(structures)
            if (args.debug):
                print('# Curren Phi and Simple Target Energies: ', current_phi, energies[i])
    # correct target energies with offsets
    for i, o in offsets.items():
        target_energies[i] -= o

    return target_energies, offsets, sampler.construction_time

def getEnergyOffsets(structures, newsample, energies, args):
    nstr = len(structures)
    offsets = np.zeros( (len(newsample), nstr) )
    # iterate over sample
    for i, s in enumerate(newsample):
        # get ViennaRNA fold compound object
        md = RNA.md()
        md.temperature = args.temperature
        fc = RNA.fold_compound(s, md)
        # iterate over structures
        for structure_i, structure in enumerate(structures):
            # calculate offset between turner eos and simple model eos
            offsets[i,structure_i] = fc.eval_structure(structure) - energies[i][structure_i]
    # calculate mean offsets
    mean_offsets = {}
    for t in range(0, nstr):
        mean_offsets[t] = np.mean(offsets[:,t])
    if (args.debug):
        print('# mean offsets are: ', mean_offsets)
    return mean_offsets

def getPhi(energies, offsets):
    mean_eos = np.mean(list(energies.values()) + list(offsets.values()))
    phi = 0
    for i, eos in energies.items():
        phi += abs(eos + offsets[i] - mean_eos)
    return phi

def BalancedSamples(structures, target_energies, offsets, args, energy_step=0.5):
    BalancedSample = {}
    # construct redprint sampler object
    nstr = len(structures)
    wastefactor = 20
    number = 1000
    sampler = RPSampler(structures, model=args.model, weights=([1.0] * nstr), gcweight=1.0, temperature=args.temperature, stacksize=(wastefactor*number), debug=args.debug)

    for shift in np.arange(0, 9999, (energy_step)):
        te = [x-shift for x in target_energies]
        if (args.debug):
            print("# Current Target energies are: ", te)
        AdmissibleSample = Sample(sampler, nstr, te, target_GC=args.gc, number=number, target_energy_eps = args.eps, target_GC_eps=args.gceps, args=args)
        for s in AdmissibleSample:
            BalancedSample[getPhi(s['energies'], offsets)] = s
        # Stop criterion
        eos = []
        for s in AdmissibleSample:
            eos_mean = np.mean(list(s['energies'].values()) + list(offsets.values()))
            eos.append(eos_mean)
        tartet_energy = target_energies[0]+offsets[0]
        if (args.debug):
            print('# Stop: ', abs(np.mean(eos)-tartet_energy))
            print("# Already found: ", len(BalancedSample)/float(args.number), "%")
        if (abs(abs(np.mean(eos)-tartet_energy)) > energy_step) and len(BalancedSample) > args.number:
            break
    return BalancedSample

def Sample(sampler, nstr, target_energies, target_GC, args, target_energy_eps = 0.10, target_GC_eps=0.05, maxiterations=10, number=1000):
    # weights = [math.exp(1/((args.temperature + 273.15)*0.00198717))] * nstr

    AdmissibleSample = []
    # count iterations
    count = 0
    while count < maxiterations:
        count += 1
        # get new sequences
        if (args.debug):
            print("# Weights: ", sampler.weights)
            print('# GC weight: ', sampler.gcweight)
        newsample, energies = sampler.dump_new_stack()

        # get average structue energies for newsample
        eos = np.zeros( (len(newsample), nstr) )
        GC_freq = []

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
            for t in range(0, nstr):
                # add to eos np array in any case
                eos[i,t] = energies[i][t]
                # check if eps admissible
                #print(eos[i,t], eos[i,t]/target_energies[t], target_energies[t])
                if not (1-target_energy_eps <= eos[i,t]/target_energies[t] <= 1+target_energy_eps):
                    admissible = False
            if admissible:
                #print('# is admissible:', eos[i,:], GC)
                AdmissibleSample.append({'seq': s, 'energies': energies[i]})
        # update weights
        for t in range(0, nstr):
            e_mean = np.mean(eos[:,t])
            if (args.debug):
                print('# Energy mean: ', str(t), e_mean)
            # exp version
            sampler.weights[t] = sampler.weights[t] * (1.1**(e_mean-target_energies[t]))
            # Yann old version without positive e_mean
            #weights[t] = weights[t]*target_energies[t]/e_mean
        # update gcweight
        GC_mean = np.mean(GC_freq)
        if (args.debug):
            print('# GC mean: ', GC_mean)
            print('# Found for current Target: ', len(AdmissibleSample)/float(number), '%')
        sampler.gcweight = sampler.gcweight * target_GC/GC_mean
        # return if large enough
        if len(AdmissibleSample) >= number:
            break
    return AdmissibleSample

if __name__ == "__main__":
    main()
