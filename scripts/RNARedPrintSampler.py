#!/usr/bin/env python3
from __future__ import print_function

'''
    RNARedPrintSampler.py: Class that uses system calls to sample sequences with the RNARedPrint program. It implements the same interface as the DependencyGraph of RNAblueprint. This makes it possible to use it within RNAsketch to optimize sequences.
'''

__author__ = "Stefan Hammer"
__copyright__ = "Copyright 2017"
__version__ = "0.1"
__maintainer__ = "Stefan Hammer"
__email__ = "s.hammer@univie.ac.at"

from distutils.dir_util import copy_tree
import os
import sys
import uuid
import shutil
import math
import re
import timeit
import subprocess as sp
import numpy as np
import scipy.stats as stats

from Structure import RNAStructure

def gccontent(s):
    from collections import Counter
    # count GC content
    c = Counter(s)
    sigma = 2.0 # one per nucleotide, laplace
    return (c['G'] + c['C'] + sigma) / (len(s) + 2*sigma)

class RPSampler(object):
    '''
    RPSampler is a python wrapper for RNARedPrint and implements a similar interface to the RNAblueprint DependencyGraphMT object.

    :param structures: List of Dot-bracket structure strings
    :param constraint: Sequence constraints in IUPACK notation
    :param model: String specifying the model, can be 'uniform',  'nussinov', 'basepairs' or 'stacking'
    :param weights: List of integers specifying the boltzmann weight for the target structures, range in ]0, +infinity]
    :param gcweight: Integer specifying the GC-weight, ]0, +infinity]
    :param temperature: Temperature for the folding predictions (default: 37.0 degree celsius)
    :param stacksize: Size of one sequence sampling batch (default: 1000 sequences)
    :param StopConstruct: Bool specifying if we want to benchmark the construction time. Measured time can be obtained from construction_time after construction (default: False)
    :param RedPrintFolder: Location of the RNARedPrint project folder (default: './RNARedPrint/')
    :param debug: Bool to print debug statements (default: False)
    '''
    def __init__(self, structures, constraint='', model='nussinov', weights=[], gcweight=1, temperature=37.0, stacksize=1000, StopConstruct=False, RedPrintFolder = None, debug=False):
        self._structures = structures
        self._constraint = constraint
        self.model = model
        self._stacksize = stacksize
        self.sample_time = 0
        self.construction_time = 0
        self.gcweight = gcweight
        self._temperature = temperature
        self.weights = weights
        self._debug = debug

        if not RedPrintFolder:
            RedPrintFolder = self._get_path('')
            if not RedPrintFolder:
                RedPrintFolder = os.path.join(os.path.dirname(os.path.abspath(__file__)),"../../RNARedPrint/")
        # copy RNARedprint binary to temp toDirectory for multithread
        self._copy_RNAredprint_folder(RedPrintFolder)

        # call RNAredprint to get construction time
        if StopConstruct:
            self._call_RNAredprint(0)

    def __del__(self):
         shutil.rmtree(self._RedPrintFolder, ignore_errors=True)

    def _get_path(self, subfolder):
        '''
        Try to detect the RNARedPrint project folder location

        :param subfolder: string specifying a subfolder, e.g. 'bin/'
        :return: string with the complete project path
        '''
        if 'REDPRINT' in os.environ:
            return os.path.join(os.environ['REDPRINT'], subfolder)
        else:
            return subfolder

    @property
    def model(self):
        '''
        String specifying the sampling model.

        :return: string, can be 'uniform',  'nussinov', 'basepairs' or 'stacking'
        '''
        return self._model
    @model.setter
    def model(self, value):
        modelarg = None
        if value == 'uniform':
            modelarg = 0
        elif value == 'nussinov':
            modelarg = 1
        elif value == 'basepairs':
            modelarg = 2
        elif value == 'stacking':
            modelarg = 3
        else:
            raise ValueError('Modelarg not set!')
        self._modelarg = modelarg
        self._model = value

    @property
    def weights(self):
        '''
        List of integers specifying the boltzmann weight for the target structures, range in ]0, +infinity]

        '''
        return self._weights
    @weights.setter
    def weights(self, value):
        # standard weights for boltzmann sampling
        if not value:
            value = [math.exp(1/((self._temperature + 273.15)*0.00198717))] * len(self._structures)
        # check amount of weights
        if len(value) != len(self._structures):
            raise ValueError('Amount of weights must be equal to the amount of structures')
        for w in value:
            if w <= 0:
                raise ValueError('Weights must be in range ]0, +infinity]')
        self._weights = value

    @property
    def gcweight(self):
        '''
        Integer specifying the GC-weight, ]0, +infinity]
        '''
        return self._gcweight
    @gcweight.setter
    def gcweight(self, value):
        if value <= 0:
            raise ValueError('GC weight must be in range ]0, +infinity]')
        else:
            self._gcweight = value

    def dump_new_stack(self, stacksize=None):
        '''
        Draws a sample from RNARedPrint and returns the list of generated sequences at once. The variable stacksize specifies the amount of sequences.

        :return: List of RNA sequence strings
        '''
        if not stacksize:
            stacksize = self._stacksize
        newseqs, energies = self._call_RNAredprint(number=stacksize)
        return newseqs, energies

    def _copy_RNAredprint_folder(self, RedPrintFolder):
        # copy subdirectory example
        fromDirectory = RedPrintFolder
        toDirectory = "/tmp/RNARedPrint-" + str(uuid.uuid4())
        self._RedPrintFolder = toDirectory
        copy_tree(fromDirectory, toDirectory)

    def _call_RNAredprint(self, number=1000):
        # resolve PKs
        structuresNoPK = {}
        weights = {}
        for i, s in enumerate(self._structures):
            structuresNoPK[s] = RNAStructure(s).resolvePKs()
            weights[s] = self._weights[i]

        # create structure and weight list
        structures_cmd = []
        weights_cmd = []
        for struct, NoPKs in structuresNoPK.items():
            for NoPK in NoPKs:
                structures_cmd.append('"'+str(NoPK)+'"')
                weights_cmd.append(weights[struct])

        cmd = ['bin/RNARedPrint'] + structures_cmd + ['--num', str(number), '--model', str(self._modelarg), '--gcw', str(self._gcweight), '--weights', ','.join(map(str, weights_cmd))]
        #cmd=['ls']
        cmdstring = " ".join(cmd)
        if (self._debug):
            print("# ", cmdstring)
        #cmd = ['RNAblueprint', '-n', str(number)]
        #stdin = '\n'.join(structures+[constraint])
        #p = sp.Popen(cmd, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE)
        start = timeit.default_timer()
        p = sp.Popen(cmdstring, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE, shell=True, cwd=self._RedPrintFolder )
        #out, err = p.communicate(input=stdin)
        out, err = p.communicate(input=None)
        if err:
            exit(err)
        if number == 0:
            self.construction_time = timeit.default_timer() - start
        else:
            self.sample_time += timeit.default_timer() - start

        if not isinstance(out, str):
            out = out.decode()

        if (self._debug):
            #print("Output:\n", out)
            if err:
                print("Error:\n",err)


        lines = out.splitlines()
        seqpattern = re.compile(r"^[AUGC]+")
        strucpattern = re.compile(r"^[\.\(\)]+$")
        epattern = re.compile(r"E(\d)=(-?[\d]+\.?[\d]*)")

        structures = []
        sequences = []
        energies = []
        for l in lines:
            s = re.search(strucpattern, l, flags=0)
            if s:
                structures.append(s.group(0))
            else:
                s = re.search(seqpattern, l, flags=0)
                if s:
                    sequences.append(s.group(0))
                    struct_energy = {}
                    for e in re.finditer(epattern, l, flags=0):
                        struct_energy[structures[int(e.group(1))-1]] = float(e.group(2))
                    energy = {}
                    for i, s in enumerate(self._structures):
                        energy[i] = 0
                        for sNoPKs in structuresNoPK[s]:
                            energy[i] += struct_energy[sNoPKs]
                    energies.append(energy)
        return sequences, energies
