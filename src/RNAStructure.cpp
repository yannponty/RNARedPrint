#include "RNAStructure.hpp"
#include <cassert>
#include <stack>
#include <fstream>

////// Base pairs elements //////

BasePair::BasePair(int a, int b, int label){
  i=a;
  j=b;
  id=label;
}

double BasePair::scoreBasePair(Nucleotide n1, Nucleotide n2){
  switch(n1){
    case N_A:
      switch(n2){
        case N_U:
          return 0.;
          break;
      case N_A:
      case N_G:
      case N_C:
          break;
      }
      break;
    case N_C:
      switch(n2){
        case N_G:
          return 0.;
          break;
        case N_A:
        case N_U:
        case N_C:
          break;
      }
      break;
    case N_G:
      switch(n2){
        case N_C:
          return 0.;
          break;
        case N_U:
          return 0.;
          break;
        case N_A:
        case N_G:
          break;
      }
      break;
    case N_U:
      switch(n2){
        case N_A:
          return 0.;
          break;
        case N_G:
          return 0.;
          break;
        case N_C:
        case N_U:
          break;
      }
      break;
  }
  return DBL_MAX;
}



ostream& operator<<(ostream& o, BasePair* bp){
  o << "("<<bp->i<<","<<bp->j<<")";
  return o;
}

ostream& operator<<(ostream& o, const vector<BasePair*> & v){
      o << "[";
      for(unsigned int i=0;i<v.size();i++)
      {
        if (i!=0)
          o<<", ";
        o << ""<< v[i];
      }
      o<<"]";
      return o;
    }

///////////////////////////////////


////// Secondary structures ///////

SecondaryStructure::SecondaryStructure(string s, int idd)
{
  length = s.size();
  id = idd;
  parseSS(s);
}

void SecondaryStructure::parseSS(string s){
  stack<int> p;
  for (unsigned int i=0;i<s.size();i++)
  {
    pos2BPs.push_back(vector<BasePair*>());
    if (s[i]=='(')
    {
      p.push(i);
    }
    else if (s[i]==')')
    {
      int j=p.top();
      p.pop();
      BasePair * bp = new BasePair(j,i);
      basePairs.push_back(bp);
      pos2BPs[i].push_back(bp);
      pos2BPs[j].push_back(bp);
    }
  }
}

vector<int> SecondaryStructure::getSS()
{
  vector<int> result;
  for (unsigned int i=0;i<length;i++)
  {
    result.push_back(-1);
  }
  
  for(unsigned int k=0;k<basePairs.size();k++)
  {
    BasePair * bp = basePairs[k];
    result[bp->i] = bp->j;
    result[bp->j] = bp->i;
  }
  
  return result;
}

vector<BasePair *> SecondaryStructure::getLabelledBPs()
{
  vector<BasePair *> result;
  for(unsigned int k=0;k<basePairs.size();k++)
  {
    BasePair * bp = basePairs[k];
    result.push_back(new BasePair(bp->i,bp->j,id));
  }    
  return result;
}


int SecondaryStructure::getLength(){
  return length;
}

vector<BasePair*> SecondaryStructure::getPartners(int i){
  return pos2BPs[i];
}

void saveAsDGF(vector<SecondaryStructure*> SSs, string path)
{
  vector<BasePair *> res;
  int len = -1;
  for(unsigned int i=0;i<SSs.size();i++)
  {
    vector<BasePair *> ss = SSs[i]->getLabelledBPs();
    res.insert(res.end(), ss.begin(), ss.end());
    if (len<0)
    { len = SSs[i]->getLength(); }
    else 
    {  assert(len == SSs[i]->getLength()); }
  }
  ofstream outfile;
  outfile.open(path.c_str());
  for (int i = 0; i < len; i++)
  {
       outfile << "n "<< i << '\n';
  }

  for (unsigned  i = 0; i < res.size(); i++)
  {
    BasePair * bp = res[i];
    outfile << "e "<< bp->i << " " << bp->j << '\n';
  }
  
  outfile.close();  
}

///////////////////////////////////
