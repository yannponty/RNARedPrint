#include <sstream>


#include "Utils.hpp"
#include "EnergyModels.hpp"

bool DEBUG = 0;

EnergyModel  dGModel;

using namespace std;

vector<string> split(const std::string &s, char delim) {
  vector<string> result;
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        result.push_back(item);
    }
    return result;
}

string pad(const char * c, int a){
  string res = "";
  for (int i=0;i<a;i++){
    res += c;
  }
  return res;
}

ostream& operator<<(ostream& o, const vector<int> & v){
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

vector<int> setSubstract(vector<int> indices,vector<int> parent){
  vector<int> res;
  for(unsigned int i=0;i<indices.size();i++){
    if (find(parent.begin(), parent.end(), indices[i] ) == parent.end())
    {
      res.push_back(indices[i]);
    }
  }
  return res;
}

