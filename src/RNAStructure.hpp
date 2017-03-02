#ifndef _RNA_STRUCTURE_H
#define _RNA_STRUCTURE_H

#include<vector>
#include<iostream>
#include <float.h>

#include "Nucleotide.hpp"

using namespace std;
/**
 * @brief The BasePair class encodes a basic base pair
 */
class BasePair{
  public:
    int i;
    int j;
    int id;
    /**
     * @brief BasePair Creates a base pair
     * @param a 5' position of base pair
     * @param b 3' position of base pair
     */
    BasePair(int a, int b, int label=-1);

    /**
     * @brief scoreBasePair Associates a free-energy contribution to a base pair, given a pair of nucleotides
     * @param n1 5' nucleotide
     * @param n2 3' nucleotide
     * @return A free-energy contribution, or +infty if the nucleotides are not allowed to pair
     */
    double scoreBasePair(Nucleotide n1, Nucleotide n2);

};

/**
 * @brief operator << Prints a base pair
 * @param o Output stream
 * @param bp Pointer to a base pair
 * @return Reference to the output stream
 */
ostream& operator<<(ostream& o, BasePair* bp);

/**
 * @brief operator << Prints a sequence of base pairs
 * @param o Output stream
 * @param bp Sequence of pointers to base pairs
 * @return Reference to the output stream
 */
ostream& operator<<(ostream& o, const vector<BasePair*> & v);

/**
 * @brief The Loop class represents a loop in the Turner model, ie a set of indices
 */
class Loop{
  public:
    int i;
    int j;
    vector<BasePair *> basePairs;
    /**
     * @brief Loop Constructs a loop
     * @param a 5' position of the enclosing base pair (-1 if exterior face)
     * @param b 3' position of the enclosing base pair (-1 if exterior face)
     * @param ss Set of positions involved in the loop
     */
    Loop(int a, int b, vector<int> ss);
};

/**
 * @brief The SecondaryStructure class represents an RNA secondary structure, possibly with pseudoknots
 */
class SecondaryStructure{
  private:
    vector<BasePair *> basePairs;
    vector<vector<BasePair*> > pos2BPs;
    vector<Loop*> loops;
    vector<vector<Loop*> > pos2Loops;
    unsigned int length;
    int id;

    /**
     * @brief parseSS Parses a dot-bracket representation and fills internal structures
     * @param s Dot-bracket representation of secondary structure
     */
    void parseSS(string s);

  public:
    /**
     * @brief SecondaryStructure Constructs a secondary structure having given id from a dot-bracket representation
     * @param s Dot-bracket representation of secondary structure
     * @param idd Uniquer identifier for the secondary structure
     */
    SecondaryStructure(string s, int idd);
    
    /**
     * @brief getSS Returns a representation of the secondary structure as an array of base pairing partners
     * @return Representation of the secondary structure as an array of base pairing partners
     *
     * Remark: Should only be called for non-pseudoknotted structures
     */
    vector<int> getSS();
    
    /**
     * @brief getLabelledBPs Returns a set of base pairs decorated by the structure identifier
     * @return Set of base pairs decorated by the structure identifier
     */
    vector<BasePair *> getLabelledBPs();
    
    /**
     * @brief getLength Returns the number of nucleotides involved in the secondary structure
     * @return Number of nucleotides involved in the secondary structure
     */
    int getLength();

    /**
     * @brief getPartners Returns the set of base pairs in which a given position is involved
     * @param i Position
     * @return Set of base pairs associated with a given position
     */
    vector<BasePair*> getPartners(int i);

};

/**
 * @brief saveAsDGF Writes to disk the graph (DGF format) associated with a set of secondary structures
 * @param SSs Sequence of secondary structures
 * @param path Path to a file where the DGF-formatted graph is written
 */
void saveAsDGF(vector<SecondaryStructure*> SSs, string path);

#endif
