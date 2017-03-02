#include "TreeDecomposition.hpp"
#include <cassert>


Bag::Bag(int i){
  id = i;  
  parent = NULL;
  //cout << "[+1']";
}

void Bag::addIndex(int i){
  indices.push_back(i);  
}        

void Bag::addChild(Bag * b){
  children.push_back(b);  
  b->setParent(this);
}        

void Bag::setParent(Bag * b){
  parent = b;  
}

// Checks that there is exactly one proper index per bag, and put it 
// at the end of the list for conveniency of future accesses
void Bag::orderIndices(){
  assert(numProper()==1);
  int proper = getProperIndices()[0];      
  (*remove(indices.begin(), indices.end(), proper)) = proper;
}


int Bag::numProper(){
  return getProperIndices().size();
}

vector<int> Bag::getProperIndices(){
  if (parent==NULL){
    return indices; 
  }
  else{
    vector<int> sub = setSubstract(indices,parent->indices);
    return  sub;
  }
}

void Bag::replaceChild (Bag * prev, Bag * next){
  for(unsigned int i=0;i<children.size();i++)
  {
    Bag * b = children[i];
    if (b->getId()==prev->getId()){
      children[i] = next;
      next->setParent(this);
    }
  }
}

vector<Bag*> Bag::getChildren(){
  return children;
}

vector<int> Bag::getProperParentIndices(){
  if (parent==NULL){
    vector<int> empty;
    return empty; 
  }
  else{
    vector<int> sub = setSubstract(indices,getProperIndices());
    return  sub;
  }
}

int Bag::numProperParentIndices(){
  return getProperParentIndices().size();
}

const vector<BasePair*> & Bag::getProperBasePairs(){
  return basepairs;
}

vector<int> Bag::getIndices(){
  
    return indices; 
}

int Bag::width(){
  return indices.size();  
}        

int Bag::getId()
{  return id;  }

void Bag::addBasePair(BasePair * bp)
{  basepairs.push_back(bp); }

void Bag::topologicalSort(vector<Bag *> & result){
  for(unsigned int i=0;i<children.size();i++)
  {
    children[i]->topologicalSort(result);
  }
  result.push_back(this);
}


double scoreBasePair(BasePair * bp, Nucleotide n1, Nucleotide n2){
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



double Bag::scoreBag(vector<Nucleotide> assignments){
  double result = 0.;
  vector<int> pi = getProperIndices();
  assert(pi.size()==1);
  vector<BasePair*> vbp = getProperBasePairs();
  vector<int> ind = getIndices();
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
    result += scoreBasePair(bp, assignments[a],assignments[b]);
  }
  return result;
}

ostream& operator<<(ostream& o, Bag * b){
  o << b->getId() << "(<-"<<((b->parent==NULL)?-1:b->parent->getId()) <<"): [";
  vector<int> indices = b->getIndices();
  for(unsigned int i=0;i<indices.size();i++)
  {
    if (i!=0)
      o<<", ";
    o << ""<< indices[i];
  }
  vector<BasePair*> vbp = b->getProperBasePairs();
  o<<"] "<<vbp ;
  return o;
}


void TreeDecomposition::normalize(){
  if (DEBUG) cout <<"Before normalization:"<<endl;
  if (DEBUG) show(1);
  //cout <<"A'";
  int numBags = bags.size();
  for(int i=0;i<numBags;i++)
  {
    
    //cout <<endl<<"(1)";
    Bag * b = bags[i];
    //cout <<"(1a)";
    vector<int> properIndices = b->getProperIndices();
    //cout <<"(1b)"<<properIndices;
    // More than a single proper index will mess with the DP, 
    // so we introduce a sequence of intermediate bags 
    if (properIndices.size()>1)
    {
      //cout <<"(2)";
      // Parent bags must be propagated, 
      // and proper indices must be added gradually
      vector<int> parentIndices = b->getProperParentIndices();
      //cout <<"(2b)";
      if (DEBUG) cout <<parentIndices;
      // If N proper indices, create N-1 new bags
      for(unsigned int j=0;j<properIndices.size()-1;j++)
      {
        //cout <<"(3)";
        Bag *newBag= new Bag(bags.size());
        // each bag consist of the parent indices...
        for(unsigned int k=0;k<parentIndices.size();k++)
        { newBag->addIndex(parentIndices[k]); }
        // ... augmented with j+1 proper indices            
        for(unsigned int k=0;k<=j;k++)
        { newBag->addIndex(properIndices[k]); }

        // The first created bag needs to replace the initial bag
        // within the parent children list
        //cout <<"(4)";
        if (j==0)
        {
          //cout <<"(5)";
          newBag->parent = b->parent;
          if (newBag->parent!=NULL)
          {  newBag->parent->replaceChild(b,newBag);  }
          // Current bag was one of the roots... update roots list 
          // to include the first bag of the new list
          else 
          {
            replaceRoot(b->getId(), newBag->getId());
          }
          //cout <<"(6)";
        }
        // otherwise, the bag added at the previous step 
        // is the parent of the current bag
        else
        {  bags[bags.size()-1]->addChild(newBag); }
        //cout <<"(7)";

        // and then add the bag to the list
        bags.push_back(newBag);
        //cout <<"(8)";
      }
      //cout <<"(9)";
      // When all is said and done, the last bag of the 
      // created sequence must point to the complete bag
      bags[bags.size()-1]->addChild(b);
      // replace root if necessary
    }
  }
  if (DEBUG) cout <<"After normalization:"<<endl;
  if (DEBUG) show(1);
  for(unsigned int i=0;i<bags.size();i++)
  {
    bags[i]->orderIndices();
  }
  //cout <<"B"<<endl;
  if (DEBUG) cout <<"After reordering:"<<endl;
  if (DEBUG) show(1);
  
}

void TreeDecomposition::replaceRoot(int from, int to){
for (unsigned int i=0;i<roots.size();i++){
  if (roots[i]==from){
    roots[i] = to;
  }
}
}

TreeDecomposition::TreeDecomposition(){
}

vector<Bag*> TreeDecomposition::topologicalSort(){
  vector<Bag*> res;
  for(unsigned int i=0;i<roots.size();i++){
    Bag * b = bags[roots[i]];
    b->topologicalSort(res);
  }
  return res;
}

void TreeDecomposition::loadFromFile(string path){
  vector<vector<int> > edges;
  
  // Parse tree-decomposition into temporary edge structure
  std::string line;
  std::ifstream infile(path.c_str());
  if (!infile.good()){
    cerr << "Error: Missing Tree Dec. file" <<endl;
    return ;
  }
  while (std::getline(infile, line))
  {
    vector<string> data = split(trim(line),' ');  
    if(data.size()>1)
    {
      string first = data[0];  
      if (first[first.length()-1]==':')
      {
        string sidbag = first.substr(3,first.length()-1);
        int idbag;
        std::istringstream(sidbag) >> idbag;
        Bag * b = new Bag(idbag-1) ;
        for(unsigned int j=1;j<data.size();j++)
        {
          int index;
          std::istringstream(data[j]) >> index;
          b->addIndex(index);
        }
        bags.push_back(b);
        edges.push_back(vector<int>());
      }
      else
      {
        int ib1,ib2;
        string sb1 = data[0].substr(3);
        std::istringstream(sb1) >> ib1;            
        string sb2 = data[1].substr(3);
        std::istringstream(sb2) >> ib2;
        edges[ib1-1].push_back(ib2-1);
        edges[ib2-1].push_back(ib1-1);
      }
    }  
  }
  
  // Find lowest degree bag => root
  bool seen[bags.size()];
  for(unsigned int i=0;i<bags.size();i++){
    seen[i]=false;
  }
  
  int minB;
  do{
    minB = -1;
    for(unsigned int i=0;i<bags.size();i++){
      if (!seen[i])
      {
        Bag * b = bags[i];
        if (minB==-1 || b->width()<bags[minB]->width())
        {
          minB = i;
        }
      }
    }
    if (minB!=-1)
    {
      // This bag was not 'seen' by previous iteration(s)
      // => new connected component starts here! 
      roots.push_back(minB);
      // Add children/parents to bags
      stack<int> p;
      p.push(minB);
      while(!p.empty())
      {
        int b = p.top();
        p.pop();
        if (!seen[b]){
          //cout << b << endl;
          seen[b] = true;
          for(unsigned int i=0;i<edges[b].size();i++){
            int child = edges[b][i];
            //cout << "  -> " << child << endl;
            if (!seen[child])
            {
              bags[b]->addChild(bags[child]);
              bags[child]->setParent(bags[b]);
              p.push(child);
            } 
          }
        }
      }
    }
  }
  while(minB!=-1);
  // Replace bags with >1 proper indices with sequences
  normalize();      
}

void TreeDecomposition::addStructure(SecondaryStructure * ss){
  for(unsigned int i=0;i<bags.size();i++){
    vector<int> pi = bags[i]->getProperIndices();
    assert(pi.size()==1);
    vector<int> all = bags[i]->getIndices();
    int proper = pi[0];
    vector<BasePair*> vbp = ss->getPartners(proper);
    for (unsigned int j=0;j<vbp.size();j++)
    {
      BasePair * bp = vbp[j];
      if ( (find(all.begin(), all.end(), bp->i ) != all.end())
      && (find(all.begin(), all.end(), bp->j ) != all.end())
      && ((proper== bp->i) || (proper== bp->j) ))
      { 
        bags[i]->addBasePair(bp); 
      }
      else
      {
      }
    }
  }
}

vector<Bag*> TreeDecomposition::getBags(){
  return bags;
}

void TreeDecomposition::showRec(int b, int depth){
  Bag * bb = bags[b];
  cout << pad(" ",depth) << bb << endl;
  for (unsigned int i=0;i<bb->children.size();i++)
  {  showRec(bb->children[i]->getId(), depth+1);  }
}

void TreeDecomposition::show(int depth){
  for (unsigned int i=0;i<roots.size();i++)
  { showRec(roots[i],depth);}
}

TreeDecomposition* TDLibFactory::makeTD(vector<SecondaryStructure *>& structures)
{
    string tmpfilein = "./tmp.dgf";
    string tmpfileout = "./tmp.td";

    TreeDecomposition * td = new TreeDecomposition();
    saveAsDGF(structures,tmpfilein);

    // TODO: To be replaced by a modular system for calling TD tools (or, even better a call to an external C++ API )
    string cmd = string("java -cp \"../lib/libtw/libtw.jar;../lib/libtw/\" TD ")+tmpfilein+string(" ")+tmpfileout+" > out 2> outerr";
    system(cmd.c_str());
    new TreeDecomposition();

    td->loadFromFile(tmpfileout);
    remove(tmpfilein.c_str());
    remove(tmpfileout.c_str());

    return td;
}