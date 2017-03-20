#ifndef _NUCLEOTIDE_H
#define _NUCLEOTIDE_H
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

typedef enum {N_A, N_U, N_G, N_C } Nucleotide;

/*
 * Number of nucleotides
 * Can be reduced (eg to 2) for debugging purpose
 */
#define NUM_NUCLEOTIDES 4

/**
 * @brief operator << Prints a nucleotide
 * @param o Output stream
 * @param n Nucleotide
 * @return Reference to the output stream
 */
ostream& operator<<(ostream& o, Nucleotide n);

/**
 * @brief operator << Prints a sequence of nucleotides
 * @param o Output stream
 * @param v Sequence of nucleotides
 * @return Reference to the output stream
 */
ostream& operator<<(ostream& o, const vector<Nucleotide> & v);
#endif
