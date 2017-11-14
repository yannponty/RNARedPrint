#ifndef _RNA_STRUCTURE_H
#define _RNA_STRUCTURE_H

#include<vector>
#include<iostream>
#include <float.h>

#include "Nucleotide.hpp"

using namespace std;

/**
 * @brief The Loop class represents a loop in the Turner model, ie a set of indices
 */
class Loop{
  public:

    vector<int> indices;
    double weight;

    /**
     * @brief getIndices Returns the set of indices in
     * @return A vector of indices contained in the loop
     */
    vector<int> &
    getIndices(){
        return indices;
    }

    /**
     * @brief getWeight Returns the weight associated to the loop
     * @return The weight of this loop
     */
    double getWeight(){
        return weight;
    }

    /**
     * @brief setWeight Assigns a weight to this loop
     * @param w The new weight of this loop
     */
    double setWeight(double w){
        weight = w;
    }

    /**
     * @brief scoreLoop Associates a free-energy contribution to a set of indices, given a set of nucleotides
     * @param nts List of assigned nucleotides
     * @return A free-energy contribution, or +infty if the nucleotides are not allowed within that structural motif
     */
    virtual double scoreLoop(vector<Nucleotide> n){
        return 0.;
    }
};

/**
 * @brief operator << Prints a loop
 * @param o Output stream
 * @param bp Pointer to a loop
 * @return Reference to the output stream
 */
ostream& operator<<(ostream& o, Loop* l);

/**
 * @brief operator << Prints a sequence of loops
 * @param o Output stream
 * @param bp Sequence of pointers to loops
 * @return Reference to the output stream
 */
ostream& operator<<(ostream& o, const vector<Loop*> & v);



/**
 * @brief The BasePair class encodes a basic base pair
 */
class BasePair: public Loop{
  public:
    int i;
    int j;
    bool isTerminal;

    /**
     * @brief BasePair Creates a base pair
     * @param a 5' position of base pair
     * @param b 3' position of base pair
     */
    BasePair(int a, int b, bool isTerm=false);

    /**
     * @brief scoreLoop Associates a free-energy contribution to a set of indices, given a set of nucleotides
     * @param nts List of assigned nucleotides
     * @return A free-energy contribution, or +infty if the nucleotides are not allowed within that structural motif
     */
    double scoreLoop(vector<Nucleotide> n);

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
 * @brief The BasePair class encodes a basic base pair
 */
class Stack: public Loop{
  public:
    int i5;
    int j5;
    int j3;
    int i3;

    /**
     * @brief BasePair Creates a base pair
     * @param a 5' position of base pair
     * @param b 3' position of base pair
     */
    Stack(int a5, int b5, int b3, int a3);

    /**
     * @brief scoreLoop Associates a free-energy contribution to a set of indices, given a set of nucleotides
     * @param nts List of assigned nucleotides
     * @return A free-energy contribution, or +infty if the nucleotides are not allowed within that structural motif
     */
    double scoreLoop(vector<Nucleotide> n);

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
    bool stackModel;

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
    SecondaryStructure(string s, int idd, bool stack);
    
    void setStacked(bool stack);

    /**
     * @brief checkSequence Verifies that a given sequence is a valid sequence for this secondary structure
     * @param seq Verified sequence
     * @return true if valid, false otherwise
     */
    bool checkSequence(const string & seq);

    bool hasIsolatedBasePair();

    /**
     * @brief getEnergy Returns the energy, in the current model of the argument for this secondary structure
     * @param seq Sequence
     * @return Energy in the current energy model
     */
    double getEnergy(const string & seq);

    /**
     * @brief getSS Returns a representation of the secondary structure as an array of base pairing partners
     * @return Representation of the secondary structure as an array of base pairing partners
     *
     * Remark: Should only be called for non-pseudoknotted structures
     */
    vector<int> getSS();

    /**
     * @brief getLoops Returns a set of Loops
     * @return Set of Loops
     */
    vector<Loop*> getLoops();

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
    //vector<BasePair*> getPartners(int i);

};


ostream & operator<<(ostream & o, SecondaryStructure * ss);

/**
 * @brief saveAsDGF Writes to disk the graph (DGF format) associated with a set of secondary structures
 * @param SSs Sequence of secondary structures
 * @param path Path to a file where the DGF-formatted graph is written
 */
void saveAsDGF(vector<SecondaryStructure*> SSs, string path, int type);

#endif
