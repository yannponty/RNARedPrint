#include "DP.hpp"
#include "Utils.hpp"




void show(vector<int> indices, vector<Nucleotide> assignments){
  assert(indices.size()==assignments.size());\
  cout << "[";
  for(unsigned int i=0;i<indices.size();i++)
  {
    if (i!=0) cout << ",";
    cout << indices[i]<<"->"<<assignments[i];
  }
  cout << "]";
}


vector<Nucleotide> decode(long index, int numIndices){
  vector<Nucleotide> res;
  for(int i=0;i<numIndices;i++){
    Nucleotide n = (Nucleotide) (index % NUM_NUCLEOTIDES);
    res.push_back(n);
    index /= NUM_NUCLEOTIDES;
  }
  return res;
}
long encode(const vector<Nucleotide> & v){
  long res=0;
  for(int i=v.size()-1;i>=0;i--){
    res *= NUM_NUCLEOTIDES;
    res += (int)v[i];
  }
  return res;
}

vector<Nucleotide> project(Bag * init, Bag * dest, const vector<Nucleotide> & val){
  vector<int> parentIndices = init->getIndices();
  vector<int> childIndices = dest->getIndices();
  
  vector<Nucleotide> res;
  for(unsigned int j=0;j<childIndices.size();j++)
  {
    for(unsigned int i=0;i<parentIndices.size();i++)
    {
      if (parentIndices[i]==childIndices[j])
      {
        res.push_back(val[i]);
        break;
      }
    }    
  }
  assert(res.size()==childIndices.size()-1);
  return res;
}

inline void checkBounds(long i1, long b1, long i2, long b2, string msg){
  if (i1<0 || i1>=b1){
    cout <<endl<< msg<<" i1:"<<i1<<" out of range [0,"<<b1-1<< "]"<<endl;
  }
  if (i2<0 || i2>=b2){
    cout <<endl<< msg<<" i2:"<<i2<<" out of range [0,"<<b2-1<< "]"<<endl;
  }
  assert(i1>=0 && i1<b1 && i2>=0 && i2<b2);
}

double **  computePartitionFunction(TreeDecomposition * td){
  vector<Bag*> bags = td->topologicalSort();
  int numBags = bags.size();
  double ** Z = new double*[numBags];
  for (int x=0; x<numBags; x++)
  {
    Bag *b = bags[x];
    int numPar = b->numProperParentIndices();
    long range = (long)pow(NUM_NUCLEOTIDES,numPar);
    if (DEBUG) cout << b << "->"<<numPar<< ","<<range <<endl;
    Z[b->getId()] = new double[range];
  }
  double grandTotal = 1.;
  for (unsigned int x=0; x<bags.size(); x++)
  {
    Bag *b = bags[x];
    if (DEBUG) cout <<(void*) b;
    int id = b->getId();
    vector<Bag*> children = b->getChildren();
    if (DEBUG) cout <<"{"<< b<<"}";
    int numParent = b->numProperParentIndices();
    vector<int> tmp2 = b->getIndices();
    if (DEBUG) cout <<numParent<<endl;
    long numCases = (long)pow(NUM_NUCLEOTIDES,numParent);
    for(long y=0;y<numCases;y++)
    {
      vector<Nucleotide> assignment = decode(y, numParent);      
      if (DEBUG) cout << "  "<<assignment<<endl;
      
      //if (DEBUG) cout << "{" <<0<< "<=" << id<< "<"<<numBags<<"}";
      //if (DEBUG) cout << "{" <<0<< "<=" << y<< "<"<<numCases<<"}";
      
      checkBounds(id, numBags, y, numCases,"Init");
      Z[id][y] = 0.;
      
      for(int n=0;n< NUM_NUCLEOTIDES;n++)
      {
        assignment.push_back((Nucleotide) n);
        if (DEBUG) cout << "    ";
        if (DEBUG) show(tmp2,assignment);
        double localZ = b->scoreBag(assignment);
        if (DEBUG) cout << " *-> "<<localZ<< "("<< b->scoreBag(assignment)<< ")";
        for(unsigned int i=0; i<children.size(); i++){
          Bag * c = children[i];
          int idc = c->getId();
          //if (DEBUG) cout << "[0 i="<< i <<"]"<<(void*)c<<"|"<< ((void*)children[i]) <<"|"<<(void*)((b->getChildren())[i]);
          vector<int> tmp = c->getIndices();
          vector<Nucleotide> v = project(b,c,assignment);
          long numCasesChild = (long)pow(NUM_NUCLEOTIDES,c->numProperParentIndices());
          long yc = encode(v);
          vector<Nucleotide> v2 = decode(yc, c->numProperParentIndices());

          if (DEBUG) cout << ", Z["<<idc<<"]"<< v<<"("<< yc << "): "<<Z[idc][yc];
          //if (DEBUG) cout << "{" <<0<< "<=" << c->getId()<< "<"<<numBags<<"}";
          //if (DEBUG) cout << "{" <<0<< "<=" << yc<< "<"<<numCasesChild<<"}";
          checkBounds(c->getId(), numBags, yc, numCasesChild, "Child");
          localZ *= Z[idc][yc];
          if (DEBUG) cout << "->"<<localZ;
        }
        if (DEBUG) cout <<endl;
        
        //if (DEBUG) cout << "{" <<0<< "<=" << id<< "<"<<numBags<<"}";
        //if (DEBUG) cout << "{" <<0<< "<=" << y<< "<"<<numCases<<"}";
        //checkBounds(id, numBags, y, numCases,"Accumulate");
        Z[id][y] += localZ;
        assignment.pop_back();
      }
      if (b->parent == NULL)
      {
          grandTotal *= Z[id][y];
      }
      if (DEBUG) cout << "    *Total* Z["<<id<<"]"<< assignment<<"("<< encode(assignment) <<"): "<<Z[id][y]<<endl;
    }
    if (DEBUG) cout <<endl;

  }
  if (DEBUG) cout << "#Designs: "<<grandTotal<<endl;
  return Z;
}

double PF(double ** Z, TreeDecomposition * td){
  double grandTotal = 1.;
  vector<Bag*> bags = td->topologicalSort();
  for (unsigned int x=0; x<bags.size(); x++)
  {
    Bag *b = bags[x];
    int id = b->getId();
    int numParent = b->numProperParentIndices();
    long numCases = (long)pow(NUM_NUCLEOTIDES,numParent);
    for(long y=0;y<numCases;y++)
    {
        if (b->parent == NULL)
        {
            grandTotal *= Z[id][y];
        }
    }
  }
  return grandTotal;
}

#include <random>
std::default_random_engine re(std::random_device{}());

double fRand(double fMax)
{
    std::uniform_real_distribution<double> unif(0,fMax);
    return unif(re);
}

void stochasticSamplingRec(Bag * b, long y, TreeDecomposition * td, double ** Z, char * result){
    int id = b->getId();
    vector<Bag*> children = b->getChildren();
    int numParent = b->numProperParentIndices();
    vector<Nucleotide> assignment = decode(y, numParent);
    // Assumption: Each bag has exactly one proper index
    int pi = b->getProperIndices()[0];

    double r = fRand(Z[id][y]);
    if (DEBUG) {cerr << "    r <- "<<r<<"(Max="<<Z[id][y]<<")" <<endl;}

    for(int n=0;n<NUM_NUCLEOTIDES;n++)
    {
        assignment.push_back((Nucleotide) n);
        double localZ = b->scoreBag(assignment);
        if (DEBUG) {  cerr << "    Score bag for "<< nt2char((Nucleotide) n)<<": "<<b->scoreBag(assignment)<<endl;}
        for(unsigned int i=0; i<children.size(); i++){
            Bag * c = children[i];
            int idc = c->getId();
            vector<Nucleotide> v = project(b,c,assignment);
            long yc = encode(v);
            localZ *= Z[idc][yc];
        }
        r -= localZ;
        if (DEBUG) cerr << "    r <- "<<r<<endl;
        assignment.pop_back();
        if (r<0){
            result[pi]=nt2char((Nucleotide) n);
            if (DEBUG) cerr << "    Found!"<<pi<<"="<< ((Nucleotide) n) <<endl;
            for(unsigned int i=0; i<children.size(); i++){
                Bag * c = children[i];
                int idc = c->getId();
                vector<Nucleotide> v = project(b,c,assignment);
                long yc = encode(v);
                vector<Nucleotide> v2 = decode(yc, c->numProperParentIndices());
                //cerr << "    Backtracking on "<<idc<<endl;
                stochasticSamplingRec(c, yc, td, Z, result);
            }
            return;
        }
    }
    assert(false&&"Could not locate suitable value for nucleotide!");
}


void stochasticSampling(double ** Z, int n, int numSamples, TreeDecomposition * td, vector<string> & result)
{
    for (unsigned int i=0;i<numSamples;i++)
    {
        char * seq = new char[n+1];
        seq[n] = '\0';
        for(int i=0;i<n;i++){
            seq[i] = 'X';
        }
        for(unsigned int r=0;r<td->roots.size();r++){
            if (DEBUG) cerr << "  Backtracking from " << td->roots[r] << endl;
            stochasticSamplingRec(td->bags[td->roots[r]], 0, td, Z, seq);
            if (DEBUG) cerr << "    After partial backtrack: " << string(seq) << endl;
        }
        if (DEBUG) cerr << "  Final seq " << string(seq) << endl;
        result.push_back(string(seq));
    }
}
