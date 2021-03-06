# RNARedPrint
Positive design of RNA sequences with specific energies for multiple target RNA structures

RNARedPrint uses Boltzmann sampling to produce tailored samples, such that highly specific sequences can be effectively generated. The approach is applies dynamic programming over tree decomposition to adapt its efficiency to the complexity of the design targets.

We provide scripts to directly target specific Turner energies (design-energyshift.py) and to generate 
sequences which balance the energies of target structures and are therefore well suited as seed sequences in negative design approaches (design-multistate.py).

Both scripts make use of the work-horse tool RNARedPrint, which implements fixed-parameter tractable Boltzmann sampling based on target structures and weights that control the distribution bias towards each target structure energy and the GC-content. 

## Prerequisites
 * g++ 5.1.1+
 * java7 compiler and runtime environment

## Installation
To install the software, clone the project from github and move into the RNARedPrint directory:

Compile via
```
make 
```

## Usage of the fundamental tool RNARedPrint
We provide convenient high-level scripts for targeting specific energies
and GC content or generating start structures for further multi-target
design optimization using "Multi-dimensional Boltzmann Sampling" below. These scripts are based on the Boltzmann
sampling engine implemented in the program RNARedPrint.

To use this low-level tool RNARedPrint directly, first change to the RNARedPrsnt directory

The following call samples 20 sequences for three target secondory structures with weights 1,2,5 in the stacking energy model.
```
bin/RNARedPrint --model 3 --weights 1,2,5 --gcw 0.5 --num 20 "((((....))))" "((((..)))).." "..((((..))))" 
```

Note that the weights directly control the sample distribution; however, there is no easy relationship between multiple weights and the final target values (energies and GC-content). Therefore, weights are typically inferred automatically as demonstrated by the provided scripts.

For detailed info on the available arameters, please see
```
bin/RNARedPrint --help
```


# Multi-dimensional boltzmann sampling

## Prerequisites

 * ViennaRNA >2.4
 * numpy
 * scipy

If the scripts are located outside of the project folder, please specify a `REDPRINT` environmental variable pointing to the RNARedPrint project folder:
`export REDPRINT='/path/to/RNARedPrint/'`

## Usage of the high-level scripts

`scripts/design-energyshift.py`: Design RNA molecules which adopt multiple structural states with specific energies using multi-dimensional Boltzmann sampling.
`scripts/design-multistate.py`: Design RNA molecules which adopt multiple structural states with equal energies using multi-dimensional Boltzmann sampling.

 * Input: Secondary RNA structures in dot-bracket notation, every state in one line, last line needs to contain only `;` to stop listening for input.
 * Models: We support different models: basepairs, stacking and uniform.
 * Targets: You can define the target energy (in kcal/mol) as well as the target GC content (in percent).
 * Tolerance: During the sampling procedure we collect candidates, which stay in a certain range of tolerance around the target. This applies for the given target energy as well as for the target GC content. For GC content and simple energy, `t` defines the range `1 + t <= current / target <= 1 + t`. For the turner energy targets, the tolerance is given in kcal/mol energy difference.

Example calls from the installation directory:
```
echo -e ".((((((......)))))).((((...((((((...((((...(((.......)))..........(((....))).))))..))))))...))))....\n.((((((......)))))).((((...((((((...((((((.(((.......((........))..)))....)).))))..))))))...))))....\n......((.((((.(((((((.((.((...((((..((.....(((.......))).....))..)))).)).)))))))))..)).)).))........" | \
  scripts/design-energyshift.py --energies="-40,-40,-20" --gc 0.7 --num 20 --tolerance=1

echo -e "(((((((((((((....)))))))))))))\n(((((.....)))))(((((.....)))))" | \
  scripts/design-multistate.py --num 20 --gc 0.6 --tolerance=0.5
```

One can as well call the script from other directories after setting the
environment variable REDPRINT like
```
export REDPRINT=</path/to/redprint/>
```


For more information run the scripts with `--help`!
 
## Disclaimer

The scripts provided here only use the `ViennaRNA` package for the evaluation of the full turner model energies and thus no pseudoknotted structures are supported as an input. If you want to design pseudoknotted structures, please use the scripts provided in `RNAsketch`: https://github.com/ViennaRNA/RNAsketch
