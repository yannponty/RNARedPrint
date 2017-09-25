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
