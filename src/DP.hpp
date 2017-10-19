#ifndef DP_H
#define DP_H

#include <vector>
#include <string>
#include <float.h>
#include <math.h>
#include <cassert>

#include "RNAStructure.hpp"
#include "Nucleotide.hpp"
#include "TreeDecomposition.hpp"
#include "EnergyModels.hpp"

using namespace std;


/**
 * @brief show Outputs to cout the nucleotides associated with a set of positions
 * @param indices Sequence of positions in an RNA sequence
 * @param assignments Sequence of nucleotides associated with the positions
 */
void show(vector<int> indices, vector<Nucleotide> assignments);


/**
 * @brief decode Decodes the encoding of a nucleotide assignment into a vector of nucleotides
 * @param index Encoding of a sequence of nucleotides in base NUM_NUCLEOTIDES
 * @param numIndices Number of positions in the bag (needed for leading zeros in the encoding)
 * @return Sequence of nucleotides associated with the positions in a bag
 *
 * Remark: encode and decode should be reciprocous (providing a suitable numIndices is provided)
 */
vector<Nucleotide> decode(long index, int numIndices);

/**
 * @brief encode Encodes an assignment as an integer
 * @param v Sequence of nucleotides
 * @return Integer representation (encoding in base NUM_NUCLEOTIDES) of the input sequence of nucleotides
 *
 * Remark: encode and decode should be reciprocous (providing a suitable numIndices is provided)
 */
long encode(const vector<Nucleotide> & v);

/**
 * @brief project Projects the assignment for a bag onto one of its descendants
 * @param init
 * @param dest
 * @param val
 * @return
 *
 * Computes the shared indices between the parent and child bags, and only retains the nucleotides assigned to
 * shared indices. The retained nucleotides appear in the same order as in the child bag.
 */
vector<Nucleotide> project(Bag * init, Bag * dest, const vector<Nucleotide> & val);

/**
 * @brief checkBounds Checks that array accesses are within prescribed bounds (0<=i1<b1 and 0<=i2<b2)
 * @param i1 Index for the first dimension
 * @param b1 Upper bound for the first dimension
 * @param i2 Index for the second dimension
 * @param b2 Upper bound for the second dimension
 * @param msg Custom message to be printed in case of out-of-bound access
 */
inline void checkBounds(long i1, long b1, long i2, long b2, string msg);

/**
 * @brief computePartitionFunction Computes and prints to standard output the partition function for a tree decomposition
 * @param td Structure-informed tree decomposition
 */
double ** computePartitionFunction(TreeDecomposition * td);


/**
 * @brief PF Computes and returns the partition function over full sequences
 * @param Z Partition function (previously computed)
 * @param td Joint tree decomposition of structure targets
 */
double PF(double ** Z, TreeDecomposition * td);

/**
 * @brief stochasticSampling Generates a list of sequences from the Boltzmann distribution and stores them into the argument vector
 * @param Z Partition function (previously computed)
 * @param n Size of sequences
 * @param numSamples Expected number of samples
 * @param td Joint tree decomposition of structure targets
 * @param result Set of generated sequences
 */
void stochasticSampling(double ** Z, int n, int numSamples, TreeDecomposition *td, vector<string> & result);

#endif
