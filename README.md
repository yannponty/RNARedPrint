# RNARedPrint
Tree-decomposition based dynamic programming algorithm for multiple target RNA design

## Prerequisites
 * g++ 5.1.1+
 * java7 compiler and runtime environment

## Usage
This first prototype only computes the number of sequences compatible with a set of secondary structures.

To count the number of sequences compatible with a bunch of secondary structures
 * ((((....))))
 * ((((..))))..
 * ..((((..))))

first clone the project, move into the directory:
 * cd bin
 * make -C ..
 * RNARedPrint  "((((....))))" "((((..)))).." "..((((..))))"


# Multi-dimensional boltzmann sampling

## Prerequisites

 * ViennaRNA >2.4
 * numpy
 * scipy

If the scripts are located outside of the project folder, please specify a `REDPRINT` environmental variable pointing to the RNARedPrint project folder:
`export REDPRINT='/path/to/RNARedPrint/'`

## Usage

`scripts/design-energyshift.py`: Design RNA molecules which adopt multiple structural states with specific energies using multi-dimensional Boltzmann sampling.
`scripts/design-multistate.py`: Design RNA molecules which adopt multiple structural states with equal energies using multi-dimensional Boltzmann sampling.

 * Input: Secondary RNA structures in dot-bracket notation, every state in one line, last line needs to contain only `;` to stop listening for input.
 * Models: We support different models: basepairs, stacking and uniform.
 * Tolerance: During the sampling procedure we will only collect candidates, which are with the given tolerance range. This applies for the given target energy as well as for the target GC content. Tolerance `t` defines the range `1 + t <= current / target <= 1 + t`.
 * Targets: You can define the target energy (in kcal/mol) as well as the target GC content (in percent).

 For more information run the scripts with `--help`!

## Disclaimer

The scripts provided here only use the `ViennaRNA` package for the evaluation of the full turner model energies and thus no pseudoknotted structures are supported as an input. If you want to design pseudoknotted structures, please use the scripts provided in `RNAsketch`: https://github.com/ViennaRNA/RNAsketch
