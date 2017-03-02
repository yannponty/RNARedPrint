#include "DP.hpp"



double scoreBag(Bag * b, vector<Nucleotide> assignments){
  double result = 0.;
  vector<int> pi = b->getProperIndices();
  assert(pi.size()==1);
  vector<BasePair*> vbp = b->getProperBasePairs();
  vector<int> ind = b->getIndices();
  for (unsigned int x=0; x<vbp.size(); x++)
  {
    BasePair * bp =vbp[x];
    if (DEBUG) cout << "("<<bp->i<<","<<bp->j<<")";
    int a =-1;
    int b =-1;
    for(unsigned int i=0;i<ind.size();i++){
      if (ind[i]==bp->i)
      {  a = i;}
      if (ind[i]==bp->j)
      {  b = i;}
    }
    assert((a!=-1)&&(b!=-1));
    if (DEBUG) cout << "("<<assignments[a]<<bp->i<<","<<assignments[b]<<bp->j<<")";
    result += bp->scoreBasePair(assignments[a],assignments[b]);
  }
  return result;
}

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

double BF(double dG){
  return exp(-dG/RT);
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
  for(unsigned int i=0;i<v.size();i++){
    res *= NUM_NUCLEOTIDES;
    res += (int)v[i];
  }
  return res;
}

vector<Nucleotide> project(Bag * init, Bag * dest, const vector<Nucleotide> & val){
  vector<int> parentIndices = dest->getProperParentIndices();
  vector<int> indices = init->getIndices();
  
  vector<Nucleotide> res;
  for(unsigned int j=0;j<indices.size();j++)
  {
    for(unsigned int i=0;i<parentIndices.size();i++)
    {
      if (parentIndices[i]==indices[j])
      {
        res.push_back(val[j]);
      }
    }    
  }
  assert(res.size()==parentIndices.size());
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

void computePartitionFunction(TreeDecomposition * td){
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
      
      if (DEBUG) cout << "{" <<0<< "<=" << id<< "<"<<numBags<<"}";
      if (DEBUG) cout << "{" <<0<< "<=" << y<< "<"<<numCases<<"}";
      
      checkBounds(id, numBags, y, numCases,"Init");
      Z[id][y] = 0.;
      
      for(int n=0;n< NUM_NUCLEOTIDES;n++)
      {
        assignment.push_back((Nucleotide) n);
        if (DEBUG) cout << "    ";
        if (DEBUG) show(tmp2,assignment);
        double localZ = BF(scoreBag(b,assignment));
        if (DEBUG) cout << " *-> "<<localZ;
        for(unsigned int i=0; i<children.size(); i++){
          Bag * c = children[i];
          if (DEBUG) cout << "[0 i="<< i <<"]"<<(void*)c<<"|"<< ((void*)children[i]) <<"|"<<(void*)((b->getChildren())[i]);
          vector<int> tmp = c->getIndices();
          vector<Nucleotide> v = project(b,c,assignment);
          long numCasesChild = (long)pow(NUM_NUCLEOTIDES,c->numProperParentIndices());
          long yc = encode(v);
          if (DEBUG) cout << ", Z["<<c->getId()<<"]"<< v<<"("<< yc << "): "<<Z[c->getId()][yc];
          if (DEBUG) cout << "{" <<0<< "<=" << c->getId()<< "<"<<numBags<<"}";
          if (DEBUG) cout << "{" <<0<< "<=" << yc<< "<"<<numCasesChild<<"}";
          checkBounds(c->getId(), numBags, yc, numCasesChild, "Child");
          localZ *= Z[c->getId()][yc];
          if (DEBUG) cout << "->"<<localZ;
        }
        if (DEBUG) cout <<endl;
        
        if (DEBUG) cout << "{" <<0<< "<=" << id<< "<"<<numBags<<"}";
        if (DEBUG) cout << "{" <<0<< "<=" << y<< "<"<<numCases<<"}";
        checkBounds(id, numBags, y, numCases,"Accumulate");
        Z[id][y] += localZ;
        assignment.pop_back();
      }
      if (DEBUG) cout << "    *Total* Z["<<id<<"]"<< assignment<<"("<< encode(assignment) <<"): "<<Z[id][y]<<endl;
      if (b->parent == NULL)
      {
        if (DEBUG) cout << "{" <<0<< "<=" << id<< "<"<<numBags<<"}";
        if (DEBUG) cout << "{" <<0<< "<=" << y<< "<"<<numCases<<"}";
        checkBounds(id, numBags, y, numCases,"Total");
        grandTotal *= Z[id][y];
      }
    }
    if (DEBUG) cout <<endl;
  }
  cout << "#Designs: "<<grandTotal<<endl;
}
