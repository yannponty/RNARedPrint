#include <assert.h>
#include "Nucleotide.hpp"

char nt2char(Nucleotide n){
    switch(n){
      case(N_A):
        return 'A';
      break;
      case(N_C):
        return 'C';
      break;
      case(N_G):
        return 'G';
      break;
      case(N_U):
        return 'U';
      break;
    }
    assert(false && "Bad Nucleotide");
}

Nucleotide char2nt(char n){
    switch(n){
      case('A'):
        return N_A;
      break;
      case('C'):
        return N_C;
      break;
      case('G'):
        return N_G;
      break;
      case('U'):
        return N_U;
      break;
    }
    assert(false && "Bad character");
}


ostream& operator<<(ostream& o, Nucleotide n)
{
    o<<nt2char(n);
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
