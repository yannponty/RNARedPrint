#include "RNAStructure.hpp"
#include <cassert>
#include <stack>
#include <fstream>
#include <sstream>
#include "Utils.hpp"
#include "EnergyModels.hpp"





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

BasePair::BasePair(int a, int b, bool isTerm){
  i=a;
  j=b;
  isTerminal=isTerm;
  indices.push_back(a);
  indices.push_back(b);
  weight=1.0;
}

double BasePair::scoreLoop(vector<Nucleotide> n){
  assert(n.size()==2);
  Nucleotide n1 = n[0];
  Nucleotide n2 = n[1];
  return dGBasePair(n1, n2, isTerminal);
}


ostream& operator<<(ostream& o, BasePair* bp){
  o << "("<<bp->i<<"-"<<bp->j<<")";
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


Stack::Stack(int a5, int b5, int b3, int a3)
{
    i5 = a5;
    j5 = b5;
    j3 = b3;
    i3 = a3;
    indices.push_back(a5);
    indices.push_back(b5);
    indices.push_back(b3);
    indices.push_back(a3);
    weight=1.0;
}

double Stack::scoreLoop(vector<Nucleotide> n)
{
    assert(n.size()==4);
    Nucleotide n5a = n[0];
    Nucleotide n5b = n[1];
    Nucleotide n3b = n[2];
    Nucleotide n3a = n[3];
    return dGStacking(n5a,n5b,n3b,n3a);
}

ostream& operator<<(ostream& o, Stack* bp){
  o << "("<<bp->i5<<"-"<<bp->i3<<","<<bp->j5<<"-"<<bp->j3<<")";
  return o;
}

ostream& operator<<(ostream& o, const vector<Stack*> & v){
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

////// Secondary structures ///////

SecondaryStructure::SecondaryStructure(string s, int idd, bool stack)
{
  length = s.size();
  id = idd;
  stackModel = stack;
  parseSS(s);
}

void SecondaryStructure::setStacked(bool stack){
    stackModel = stack;
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

  // Annotates terminal loops
  vector<int>ss =  getSS();
  for(unsigned int k=0;k<basePairs.size();k++)
  {
    BasePair * bp = basePairs[k];
    int l = ss[bp->i+1];
    if ((l==bp->i) || (l!=bp->j-1))
        bp->isTerminal = true;
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
      if (!isCompatible(char2nt(seq[bp->i]), char2nt(seq[bp->j])))
      {
        return false;
      }
  }
  return true;
}

double SecondaryStructure::getEnergy(const string & seq)
{
    return dGStructure(this,seq);
}

bool hasBasePair(int i, int j, vector<int>ss){
    if ((i<0)||(i>=ss.size()))
        return false;
    if ((j<0)||(j>=ss.size()))
        return false;
    return ss[i] == j;
}

bool SecondaryStructure::hasIsolatedBasePair()
{
    vector<int> ss = getSS();
    for(int i=0;i<ss.size();i++){
        if (ss[i]>i){
            if (!hasBasePair(i-1,ss[i]+1,ss) && !hasBasePair(i+1,ss[i]-1,ss)){
                return true;
            }
        }
    }

    return false;
}


vector<Loop*> SecondaryStructure::getLoops()
{
  vector<Loop *> result;
  if (!stackModel){
      if (DEBUG){
          cerr << "[BPs based loops]"<<endl;
          cerr.flush();
      }
  for(unsigned int i=0; i<basePairs.size(); i++){
      BasePair * bp = basePairs[i];
      result.push_back(bp);
  }
  }
  else{
      if (DEBUG){
          cerr << "[Stacks based loops]"<<endl;
          cerr.flush();
      }
      for(unsigned int i=0; i<basePairs.size(); i++){
          BasePair * bp = basePairs[i];
          if (!bp->isTerminal)
              result.push_back(new Stack(bp->i,bp->i+1,bp->j-1,bp->j));
      }
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

