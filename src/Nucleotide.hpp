#ifndef _NUCLEOTIDE_H
#define _NUCLEOTIDE_H
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

typedef enum {N_A, N_U, N_G, N_C} Nucleotide;

/*
 * Number of nucleotides
 * Can be reduced (eg to 2) for debugging purpose
 */
#define NUM_NUCLEOTIDES 4

/**
 * Converts a nucleotide to its representation as a character
 * @param n Nucleotide to be converted
 * @return Character representation ('A', 'C', 'G', 'U' & 'X' for unknown)
 */
char nt2char(Nucleotide n);


/**
 * Converts a character to its associated nucleotide
 * @param n Character to be converted
 * @return Associated nucleotide
 */
Nucleotide char2nt(char n);

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
