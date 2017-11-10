#ifndef TREE_DECOMPOSITION_H
#define TREE_DECOMPOSITION_H

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <stack>
#include <unistd.h>
#include <regex>

#include "Utils.hpp"
#include "RNAStructure.hpp"

using namespace std;

/**
 * @brief The Bag class represents a set of indices, ie positions in the sequences
 * that must be jointly considered in the context of a given tree-decomposition.
 */
class Bag{
  private:
    vector<int> indices;
    vector<int> proper_indices;
    vector<int> proper_parent_indices;
    vector<BasePair*> basepairs;
    vector<Loop*> loops;
    vector<vector<int>> indices2loops;
    int id;

  public:
    vector<Bag*> children;
    Bag* parent;

    /**
     * @brief Bag Constructs a bag with given identifier
     * @param i Unique identifier for the Bag
     */
    Bag(int i);

    /**
     * @brief addIndex Adds an index/position to this bag
     * @param i Index/position to add to the bag
     */
    void addIndex(int i);

    /**
     * @brief getIndices Returns the set of indices/positions for this bag
     * @return
     */
    vector<int> &
    getIndices();

    /**
     * @brief orderIndices Reorders the indices such that the proper index is put at the last position in the index list.
     *
     * Checks that there is exactly one proper index per bag.
     * A proper index is an index which is found in the child bag, but not in the parent bag.
     * Puts it at the end of the index/position list for conveniency of future accesses.
     */
    void orderIndices();

    /**
     * @brief addChild Adds a child to the bag
     * @param b Pointer to a bag to add as a child to this bag
     *
     * Adds a child to the bag, considered as an internal node in the tree-decomposition.
     * The parent of the child bag is set to the current bag.
g     */
    void addChild(Bag * b);

    /**
     * @brief replaceChild Replaces any occurrence of a bag in the list of children bags with another bag
     * @param prev Pointer to current child bag to be replaced
     * @param next Pointer to replacement child bag
     */
    void replaceChild (Bag * prev, Bag * next);

    /**
     * @brief setParent Sets the parent of this bag to a given bag
     * @param b Pointer to bag, set as parent to this bag
     */
    void setParent(Bag * b);

    /**
     * @brief numProper Returns the number of proper indices for this bag
     * @return Number of proper indices, ie indices that are in this bag but not in the parent bag
     */
    int numProper();
    
    void
    precomputeProperIndices();

    /**
     * @brief getProperIndices Returns the list of proper indices for this bag
     * @return List of proper indices for this bag, ie indices that are in this bag but not in the parent bag
     */
    const
    vector<int> &
    getProperIndices();
    
    /**
     * @brief getChildren Returns the list of children bags for this bag
     * @return List of children bags for this bag
     */
    vector<Bag*> getChildren();

    void
    precomputeProperParentIndices();

    /**
     * @brief getProperParentIndices Returns the list of indices that are proper to the parent, ie indices
     * that are in the parent list but not in this child list
     * @return list of indices that are proper to the parent, ie indices that are in the parent list but not in this child list
     */
    const vector<int> &
    getProperParentIndices();

    /**
     * @brief numProperParentIndices Returns the number of indices that are proper to the parent, ie indices
     * that are in the parent list but not in this child list
     * @return Number of indices that are proper to the parent, ie indices
     * that are in the parent list but not in this child list
     */
    int numProperParentIndices();

    /**
     * @brief getLoops Returns the set of loops that are proper to this bag
     * @return Set of loops that are proper to this bag
     *
     * The proper loops for a bag are the base pairs such that: a) All indices of the loop
     * are in the list of indices; b) One of the indices of the loop is proper to the current bag.
     * Note that this definition uniquely defines to which bag a loop must be attributed in a
     * given tree decomposition (assuming that the loop is materialized by some edge in the
     * graph used for the TD).
     */
    const vector<Loop*> & getLoops();

    /**
     * @brief width Returns the width of this bag
     * @return Width, ie number of indices, of this bag
     */
    int width();
    
    /**
     * @brief getId Returns the unique identifier for this bag
     * @return Unique identifier for this bag
     */
    int getId();
    
    /**
     * @brief addBasePair Add a base pair to this bag
     * @param bp Pointer to base pair to be added
     */
    void addBasePair(BasePair * bp);

    /**
     * @brief addLoop Add a loop to this bag
     * @param l Pointer to loop to be added
     */
    void addLoop(Loop *l);


    /**
     * @brief topologicalSort Sorts and returns the bags in the tree initiated at this node
     * @param result List of bags such that the leaves can be found first and, more generally,
     * children can be found before the parents.
     */
    void topologicalSort(vector<Bag *> & result);

    /**
     * @brief scoreBag Associates an overall contribution of the proper nucleotide in the bag
     * @param assignments Nucleotides assigned to each positions of the bag
     * @return The overall free-energy contribution of the (proper indices in) the bag
     */
    double scoreBag(vector<Nucleotide> assignment);

    /**
     * @brief contains Checks that the bag contains the set of indices
     * @param indices Set of indices
     * @return True if the bag contains the set of indices, false otherwise
     */
    double contains(vector<int> indices);

    /**
     * @brief contains Checks that the bag contains a single index
     * @param index A single position
     * @return True if the bag contains the index, false otherwise
     */

    double contains(int index);
};

/**
 * @brief operator << Prints a bag
 * @param o Output stream
 * @param b Pointer to a bag
 * @return Reference to the output stream
 */
ostream& operator<<(ostream& o, Bag * b);

class TreeDecomposition{
  private:
    /**
     * @brief normalize Normalizes the tree decomposition
     *
     * A normalized tree decomposition is a tree decomposition where: a) Each bag has exactly one proper index/position,
     * ie introduces a new index that is not found in its parent; and b) The proper index in each bag can be found at the
     * last position. A normalized tree decomposition is otained by: introducing new intermediate bags (carefully rerouting existing
     * parent/child relationships to address a.) and reorder indices in each bag (b.)
     */
    void normalize();
    
    /**
     * @brief replaceRoot Replaces any instance of a bag with another one, both identified by their id
     * @param from Identifier of the bag to be removed from the roots of the tree-decomposition
     * @param to Identifier of the bag to be added to the roots of the tree-decomposition
     */
    void replaceRoot(int from, int to);

    /**
     * @brief showRec Prints this tree decomposition recursively to standard error
     * @param b Identifier for the top level bag
     * @param depth Depth in tree decomposition of current bag
     */
    void showRec(int b, int depth=0);

    /**
     * @brief addLoopsRec Assigns a set of loops to the descendants of a given bag
     * @param b The ancestor bag
     * @param structures a set of loops
     */
    void addLoopsRec(Bag * b, vector<Loop* > & structures);


  public:
    vector<Bag*> bags;  
    vector<int> roots;
    int tw; // treewidth
    /**
     * @brief TreeDecomposition Constructs an empty tree decomposition
     */
    TreeDecomposition();
    
    /**
     * @brief topologicalSort Builds a topological ordering of bags
     * @return List of topologically-sorted (leaves before internal, children before parent) bags
     */
    vector<Bag*> topologicalSort();
    
    /**
     * @brief loadFromFile Loads the content of a tree-decomposition from a file
     * @param path Path to a file describing a tree-decomposition
     */
    void loadFromFile(string path);
    
    /**
     * @brief addStructure Adds the structural elements of a secondary structure to the suitable bags in the tree-decomposition
     * @param ss Secondary structure to be added
     */
    void addStructure(SecondaryStructure * ss);

    /**
     * @brief getBags Returns the set of bags for this tree decomposition
     * @return List of pointers to bags in this tree decomposition
     */
    const vector<Bag*> &
    getBags();

    /**
     * @brief show Prints this tree decomposition recursively to standard error
     * @param depth Depth in tree decomposition of current bag
     */
    void show(int depth=0);
    
    void reset(){bags.clear(); roots.clear(); tw=0;}
    void copyObj(TreeDecomposition * td){this->bags=td->bags; this->roots=td->roots; this->tw=td->tw;}


    /**
     * @brief addLoops Assigns a set of loops to the nodes of the tree-decomposition
     * @param structures a set of loops
     */
    void addLoops(vector<Loop *> structures);
};

/**
 * @brief The TreeDecompositionFactory class is an abstract class aimed at producing a tree decomposition
 * from a set of secondary structures
 */
class TreeDecompositionFactory{
  public:
    virtual TreeDecomposition* makeTD(vector<SecondaryStructure*>& v) = 0;
};

/**
 * @brief The TDLibFactory class is a concrete implementation of TreeDecompositionFactory based on the
 * TDLib.jar Java implementation of the greedy fill-in heuristics
 */
class TDLibFactory : public TreeDecompositionFactory{
  private:
  	void transFileFormat(string path);
  public:
    TreeDecomposition* makeTD(vector<SecondaryStructure*>& v);
    
};

#endif
