#include "RNAStructure.hpp"
#include <cassert>
#include <stack>
#include <fstream>
#include <sstream>
#include "Utils.hpp"





ostream& operator<<(ostream& o, Loop* lp){
  o << "(";
  for(int i=0;i<lp->getIndices().size();i++){
      if (i!=0)
          o<<",";
      o<<lp->getIndices()[i];
  }
  o<<")";
  return o;
}

ostream& operator<<(ostream& o, const vector<Loop*> & v){
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



////// Base pairs elements //////

BasePair::BasePair(int a, int b, bool isTerm, int label){
  i=a;
  j=b;
  id=label;
  isTerminal=isTerm;
  indices.push_back(a);
  indices.push_back(b);
  weight=1.0;
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

double BasePair::scoreLoop(vector<Nucleotide> n){
  assert(n.size()==2);
  Nucleotide n1 = n[0];
  Nucleotide n2 = n[1];
  return scoreBasePair(n1,n2);
}


ostream& operator<<(ostream& o, BasePair* bp){
  o << "("<<bp->i<<","<<bp->j<<")";
  if (bp->isTerminal)
  {
      o<<"{T}";
  }
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

ostream & operator<<(ostream & o, SecondaryStructure * s){
    vector<int> ss = s->getSS();
    for(int i=0;i<ss.size();i++){
        if (ss[i]==-1){
            o << '.';
        }
        else if (ss[i]>i){
            o << '(';
        }
        else{
            o << ')';
        }
    }
    return o;
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

bool SecondaryStructure::checkSequence(const string & seq)
{
  for(unsigned int k=0;k<basePairs.size();k++)
  {
      BasePair * bp = basePairs[k];
      if (bp->scoreBasePair(char2nt(seq[bp->i]), char2nt(seq[bp->j])) == DBL_MAX)
      {
        return false;
      }
  }
  return true;
}



vector<Loop*> SecondaryStructure::getLoops()
{
  vector<Loop *> result;
  for(unsigned int i=0; i<basePairs.size(); i++){
      BasePair * bp = basePairs[i];
      result.push_back(bp);
  }
  return result;
}


int SecondaryStructure::getLength(){
  return length;
}

/*
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
*/

void saveAsDGF(vector<SecondaryStructure*> SSs, string path, int type)
{
  if (DEBUG) cerr << "Saving to DGF: "<<path << endl;
  int len = SSs[0]->getLength();

  int edge = 0;
  int vert = len;

  bool** adj = new bool* [len];
  for (int i=0; i<len; i++){
      adj[i] = new bool[len];
      for (int j=0; j<len; j++){
        adj[i][j] = false;
        }
      if (i>0){
          //adj[i-1][i] = true;
      }
  }

  for (int i=0; i<SSs.size(); i++){
      SecondaryStructure * ss = SSs[i];
      vector<Loop*> loops = ss->getLoops();
      for (int j=0; j<loops.size(); j++){
          Loop* lp = loops[j];
          vector<int> indices = lp->getIndices();
          //cerr << "  " << indices <<endl;
          for (int k=0; k<indices.size(); k++){
              for (int l=0; l<indices.size(); l++){
                if (l>k){
                    adj[indices[k]][indices[l]] = true;
                }
              }
          }
      }
  }

  stringstream links;
  for (int i=0; i<len; i++){
      for (int j=0; j<len; j++){
          if (adj[i][j]){
              edge ++;
              if (type == 1)
                  links << "e "<< (i+1) << " " << (j+1) << endl;
              else if (type == 2)
                  links << (i+1) << " " << (j+1) << endl;
          }
      }
  }

  // Final number of edges and vertices
  stringstream problem;
  if (type == 1)
    problem << "p edge " << vert << " " << edge << endl;
  else if (type == 2)
      problem << "p tw " << vert << " " << edge << endl;

  if (DEBUG) cerr <<"  Problem: "<< problem.str() << endl << "   Edges: "<<links.str()<< endl;


  // Finish up the dgf file
  ofstream outfile;
  outfile.open(path.c_str());
  outfile << problem.str();
  outfile << links.str();
  outfile.close();
  
}
///////////////////////////////////

