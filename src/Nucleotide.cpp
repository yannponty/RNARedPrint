#include "Nucleotide.hpp"

ostream& operator<<(ostream& o, Nucleotide n)
{
  switch(n){
    case(N_A):
      o << "A";
    break;
    case(N_C):
      o << "C";
    break;
    case(N_G):
      o << "G";
    break;
    case(N_U):
      o << "U";
    break;
  }
  return o;
}


ostream& operator<<(ostream& o, const vector<Nucleotide> & v){
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
