#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stack>
#include <float.h>
#include <math.h>
#include <stdexcept> 


#include "Utils.hpp"
#include "Nucleotide.hpp"
#include "RNAStructure.hpp"
#include "TreeDecomposition.hpp"
#include "DP.hpp"


using namespace std;



int main(int argc, char *argv[])
{
  vector<SecondaryStructure *> structures;
  
  for(int i = 1;i<argc;i++)
  {
    if (argv[i][0]!= '-')
    {
      structures.push_back(new SecondaryStructure(string(argv[i]),structures.size()));
    }
  }
  
  TreeDecompositionFactory * tdFact;

  // To be replaced by alternative implementations
  tdFact = new TDLibFactory();

  TreeDecomposition * td = tdFact->makeTD(structures);
  
  if (td->bags.size()==0)
  {
    cerr << "Error: Tree decomposition software probably crashed!"<<endl;
    return EXIT_FAILURE;
  }
  else
  {

    for (unsigned int i=0;i<structures.size();i++){
      td->addStructure(structures[i]);
    }
    computePartitionFunction(td);
    return EXIT_SUCCESS;
  }
}
