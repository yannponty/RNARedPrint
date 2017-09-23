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

#define DEFAULT_NUM_OPTION 10


#define NUM_OPTION "--num"
#define WEIGHTS_OPTION "--weights"
#define COUNT_OPTION "--count"
#define DEBUG_OPTION "--debug"


void usage(string cmd){
  cerr << "Usage: "<<cmd<<" Struct1 Struct2 ... ["<< NUM_OPTION<< " k]"<<endl;
  cerr << "Generates valid designs for the RNA secondary structures from the weighted distribution"<<endl;
  cerr << "  "<<NUM_OPTION<<" k - Sets number of generated sequences (default "<<DEFAULT_NUM_OPTION<<")"<<endl;
  //cerr << "  "<<WEIGHTS_OPTION<<" w1,w2... - Assigns custom weights to targeted structures (default 1.)"<<endl;
  cerr << "  "<<COUNT_OPTION<<" - Simply compute the partition function and report the result."<<endl;
  exit(2);
}

int main(int argc, char *argv[]){
  vector<SecondaryStructure *> structures;
  vector<double> weights;
  bool count_mode = false;
  unsigned int numSamples = DEFAULT_NUM_OPTION;
  int n = 0;

  for(int i = 1;i<argc;i++){
    if (argv[i][0]!= '-'){
      structures.push_back(new SecondaryStructure(string(argv[i]),structures.size()));
      int nn = structures[structures.size()-1]->getLength();
      if ((n>0)&&(nn!=n)){
          cerr << "Error: Inconsistent length across target structures."<< endl;
          cerr << " (arg #"<<(i)<<"="<< n<<" != arg#"<< (i+1)<<"="<<nn<<")"<<endl;
          usage(string(argv[0]));
      }
      n = max(n,nn);
    }
    else{
      if (string(argv[i])==NUM_OPTION){
        i++;
        numSamples = atoi(argv[i]);
      }
      else if (string(argv[i])==COUNT_OPTION){
        count_mode = true;
      }
      else if (string(argv[i])==DEBUG_OPTION){
        DEBUG = true;
      }
    }
  }
  if (structures.size()==0){
      cerr << "Error: Missing target structure(s)."<< endl;
      usage(string(argv[0]));
  }

  TreeDecompositionFactory * tdFact;

  // To be replaced by alternative implementations
  tdFact = new TDLibFactory();

  TreeDecomposition * td = tdFact->makeTD(structures);
  
  if (td->bags.size()==0){
    cerr << "Error: Tree decomposition software probably crashed!"<<endl;
    return EXIT_FAILURE;
  }
  else{
    for (unsigned int i=0;i<structures.size();i++){
      td->addLoops(structures[i]->getLoops());
    }
    double ** Z = computePartitionFunction(td);
    if (count_mode){
        cout << "Partition function: "<<PF(Z,td)<<endl;
        exit(0);
    }

    vector<string> seqs;
    stochasticSampling(Z,n,numSamples,td,seqs);
    for (unsigned int j=0;j<structures.size();j++){
        cout << structures[j] << endl;
    }
    for (unsigned int i=0;i<seqs.size();i++){
        cout << seqs[i] << endl;
        cout.flush();
        for (unsigned int j=0;j<structures.size();j++){
            if (!structures[j]->checkSequence(seqs[i]))
            {
                cerr << "Error: sequence uncompatible with structure "<<structures[j]<<endl;
                exit(2);
            }
        }
    }
    return EXIT_SUCCESS;
  }
}
