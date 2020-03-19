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
#include <iomanip>

using namespace std;

#include "Utils.hpp"
#include "Nucleotide.hpp"
#include "RNAStructure.hpp"
#include "EnergyModels.hpp"
#include "DP.hpp"
#include "TreeDecomposition.hpp"

#define DEFAULT_NUM_OPTION 10


#define VERSION_OPTION "--version"
#define NUM_OPTION "--num"
#define WEIGHTS_OPTION "--weights"
#define GC_WEIGHT_OPTION "--gcw"
#define COUNT_OPTION "--count"
#define DEBUG_OPTION "--debug"
#define ENERGY_MODEL_OPTION "--model"
#define HELP_OPTION "--help"
#define TEMP_OPTION "-T"
#define PREFIX_OPTION "--prefix"


void usage(string cmd){
  cerr << "Usage: "<<cmd<<" Struct1 Struct2 ... ["<< NUM_OPTION<< " k]"<<endl;
  cerr << "Generates valid designs for the RNA secondary structures from the weighted distribution"<<endl;
  cerr << "------ Mode ------------"<<endl;
  cerr << "  "<<NUM_OPTION<<" k           - Sets number of generated sequences (default "<<DEFAULT_NUM_OPTION<<")"<<endl;
  cerr << "  "<<COUNT_OPTION<<"           - Simply compute the partition function and report the result."<<endl;
  cerr << "------ Options ------------"<<endl;
  cerr << "  "<<WEIGHTS_OPTION<<" w1,w2.. - Assigns custom weights to each targeted structure (default 1. for all)"<<endl;
  cerr << "  "<<GC_WEIGHT_OPTION<<" w       - Assigns custom weight to each occurrence of GC, to control GC% (default 1.)"<<endl;
  //cerr << "  "<<TEMP_OPTION<<" t              - Sets the pseudotemperature (default 37.C)"<<endl;
  cerr << "  "<<ENERGY_MODEL_OPTION<<" m         - Set energy model used for stochastic sampling: "<<endl
       << "        m = "<<COMPATIBLE_BP_MODEL<<": Uniform"<<endl
       << "        m = "<<NUSSINOV_BP_MODEL<<": Nussinov (-3/-2/-1 for GC/AU/GU)"<<endl
       << "        m = "<<FITTED_BP_MODEL<<": Base pair energy model (Default; distinguishs GC/AU/GU, inner/exterior)"<<endl
       << "        m = "<<FITTED_STACKING_PAIRS_MODEL<<": Stacking model (no isolated base-pairs!)"<<endl;
  cerr << "  "<<PREFIX_OPTION  <<"            - Prefix path for locating the TD libraries"<<endl;
  cerr << "  "<<VERSION_OPTION <<"            - Show version and exit"<<endl;
  cerr << "  "<<HELP_OPTION    <<"            - Display help message and exit"<<endl;
  exit(2);
}

void
version() {
  cout << "RNARedPrint 0.3" << endl;
  exit(0);
}

vector<double> parseWeights(string str)
{
    vector<double> res;
    stringstream ss(str);

     double i;

     while (ss >> i)
     {
         res.push_back(i);

         if (ss.peek() == ',')
             ss.ignore();
     }
     return res;
}

// @brief GC-contet
// @param s DNA/RNA sequence string
// @returns GC content
// @note w/ pseudocounts (1 per nucleotide)
double
gccontent(const string s) {

    return (2.0+std::count_if(s.begin(),s.end(),[] (char c) {return c=='C' || c=='G';}))/(4.0+s.length());
}

int main(int argc, char *argv[]){
  vector<SecondaryStructure *> structures;
  vector<double> weights;
  bool count_mode = false;
  unsigned int numSamples = DEFAULT_NUM_OPTION;
  int n = 0;

  string prefix_path="";

  dGModel = FITTED_BP_MODEL; // set default model: base pair model

  for(int i = 1;i<argc;i++){
    if (argv[i][0]!= '-'){
      structures.push_back(new SecondaryStructure(string(argv[i]),structures.size(),false));
      int nn = structures[structures.size()-1]->getLength();
      if ((n>0)&&(nn!=n)){
          cerr << "Error: Inconsistent length across target structures."<< endl;
          cerr << " (arg #"<<(i)<<"="<< n<<" != arg#"<< (i+1)<<"="<<nn<<")"<<endl;
          usage(string(argv[0]));
      }
      n = max(n,nn);
    }
    else{
      if (string(argv[i])==VERSION_OPTION){
        version();
      }
      else if (string(argv[i])==NUM_OPTION){
        i++;
        numSamples = atoi(argv[i]);
      }
      else if (string(argv[i])==COUNT_OPTION){
        count_mode = true;
      }
      else if (string(argv[i])==WEIGHTS_OPTION){
        i++;
        weights = parseWeights(string(argv[i]));
      }
      else if (string(argv[i])==GC_WEIGHT_OPTION){
        i++;
        GCWeight = atof(argv[i]);
      }

      else if (string(argv[i])==ENERGY_MODEL_OPTION){
        i++;
        dGModel = (EnergyModel) atoi(argv[i]);
      }
      else if (string(argv[i])==TEMP_OPTION){
        i++;
        TEMP = atof(argv[i]);
      }
      else if (string(argv[i])==PREFIX_OPTION){
          i++;
          prefix_path = argv[i];
      }
      else if (string(argv[i])==HELP_OPTION){
          usage(argv[0]);
      }
      else if (string(argv[i])==DEBUG_OPTION){
        DEBUG = true;
      }
    }
  }
  initStackingModel();

  if (DEBUG) cerr << "RT: "<<RT<<endl;
  if (structures.size()==0){
      cerr << "Error: Missing target structure(s)."<< endl;
      usage(string(argv[0]));
  }
  // No weight specified explicitly => 1. for everyone
  if (weights.size()==0){
    for (int i=0;i<structures.size();i++)
        weights.push_back(1.);
  }

  if (structures.size()!=weights.size()){
      cerr << "Error: #Specified weights does not match #Structures."<< endl;
      usage(string(argv[0]));
  }

  if (prefix_path=="") {
    if (char *p = getenv("RNAREDPRINT_PATH")) {
        prefix_path = p;
    }
  }


  // To be replaced by alternative implementations
  TreeDecompositionFactory * tdFact = new TDLibFactory(prefix_path);

  bool stackedModel = (dGModel==FITTED_STACKING_PAIRS_MODEL);

  for (unsigned int i=0;i<structures.size();i++){
    structures[i]->setStacked(stackedModel);
    if (stackedModel &&  structures[i]->hasIsolatedBasePair()){
        cerr << "Error: Isolated base-pairs currently unsupported in stacking pair mode."<< endl;
        usage(string(argv[0]));
    }
  }

  TreeDecomposition * td = tdFact->makeTD(structures);

  if (td->bags.size()==0){
    cerr << "Error: Tree decomposition software probably crashed!"<<endl;
    return EXIT_FAILURE;
  }
  else{
      for (unsigned int i=0;i<structures.size();i++){
        vector<Loop*> loops = structures[i]->getLoops();
        for (int j=0;j<loops.size();j++){
            loops[j]->setWeight(weights[i]);
        }
        td->addLoops(loops);
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
        cout << seqs[i];

        cout << fixed << setprecision(2);

        cout << " GC="<<gccontent(seqs[i]);

        for (unsigned int j=0;j<structures.size();j++){
            if (!structures[j]->checkSequence(seqs[i]))
            {
                cerr << "Error: sequence uncompatible with structure "<<structures[j]<<endl;
                exit(2);
            }
            else{
                cout<<" E"<<(j+1)<<"="<<structures[j]->getEnergy(seqs[i]);
            }
        }
        cout  << endl;
    }
    return EXIT_SUCCESS;
  }
}
