#include "TreeDecomposition.hpp"
#include <cassert>
#include <math.h>


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

const vector<Loop *> &Bag::getLoops()
{
    return loops;
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
{ basepairs.push_back(bp); }

void Bag::addLoop(Loop * l)
{
    vector<int> map;
    for(int i=0;i<indices.size();i++){
        int pos = -1;
        vector<int> loopIndices = l->getIndices();
        for (int j=0;j<loopIndices.size();j++)
        {
            if (indices[i]==loopIndices[j])
            {
                pos = j;
            }
        }
        map.push_back(pos);
    }
    if (DEBUG) cerr << "Mapping: " << this << " -> "<<l<<" => "<<map<<endl;
    loops.push_back(l);
    indices2loops.push_back(map);
}

void Bag::topologicalSort(vector<Bag *> & result){
  for(unsigned int i=0;i<children.size();i++)
  {
    children[i]->topologicalSort(result);
  }
  result.push_back(this);
}





double Bag::scoreBag(vector<Nucleotide> assignments){
  double result = 1.;
  for (int l=0;l<loops.size();l++){
      Loop * lp = loops[l];
      vector<int> map = indices2loops[l];
      assert(map.size()==assignments.size());
      vector<Nucleotide> loopNTs;
      for(int i=0;i<lp->getIndices().size();i++){
        loopNTs.push_back(N_A);
      }
      for(int i=0;i<map.size();i++){
          if (map[i]!=-1){
              loopNTs[map[i]] = assignments[i];
          }
      }
      double dGLoop = lp->scoreLoop(loopNTs);
      //cerr << endl;
      //cerr << "  Loop:"<<lp<<" Map:"<< indices2loops[l]<<" Bag:"<<this<<" Assign.:"<<assignments <<" Projected:"<<loopNTs<<" score:"<<lp->scoreLoop(loopNTs)<<endl;
      if (DEBUG) cerr<<"{"<< loopNTs <<","<< lp->weight<<","<< dGLoop<< "}";
      if (dGLoop==DBL_MAX) return 0.;
      result *= pow(lp->weight,-lp->scoreLoop(loopNTs));
  }
  return result;
}


double Bag::contains(int index){
    for(int i=0;i<indices.size();i++){
        if (indices[i]==index){
            return true;
        }
    }
    return false;
}

double Bag::contains(vector<int> indices){
    for(int i=0;i<indices.size();i++){
        if (!contains(indices[i])){
            return false;
        }
    }
    return true;
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
  vector<Loop*> vbp = b->getLoops();
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
  string line;
  ifstream infile(path.c_str());
  if(!infile.good()){
  	cerr << "Error: Missing Tree Dec. file" <<endl;
  	return ;
  }
  tw = 0;
  while(getline(infile, line)){
    //  cerr << line << endl;
  	if(line[0] == 'c') continue; // comment line
  	if(line[0] == 's') continue; 
  	// solution line s which contains the string td, followed by number of bags, the width, the vertices of the original input graph
  	
  	vector<string> data = split(trim(line),' ');
  	
  	if(line[0] == 'b'){ // bag line
  		int dataSize = data.size();
  		if (dataSize-3 > tw) tw = dataSize-3;
  		int bagIndex;
  		istringstream(data[1]) >> bagIndex;
  		Bag * b = new Bag(bagIndex-1);  // bag index begin at 0;
  		for(unsigned i=2; i<dataSize; i++){
  			int index;
            istringstream(data[i]) >> index;
            if(index>0)
                b->addIndex(index-1);         // the index in the bag also begins at 0;
  		}
  		bags.push_back(b);
  		edges.push_back(vector<int>());
  	} else {  // all the remaining non-comment lines indicate an edge in the tree decomposition
  		int ib1, ib2;
  		string sb1 = data[0];
  		istringstream(sb1) >> ib1;
  		string sb2 = data[1];
  		istringstream(sb2) >> ib2;
  		edges[ib1-1].push_back(ib2-1);
  		edges[ib2-1].push_back(ib1-1);
  	}
  }

  /*
   * #Patch for unpaired positions#
   * Unpaired positions induce a syntax error, leading to an index being unparsable.
   * Since all unpaired positions are equivalent (note that, in loop-based
   * energy models, loops become cliques and are therefore not alone in their
   * connected components), it does not matter which position is represented by
   * such bags, so we:
   * 1) identify, in an array of boolean, which positions are unpaired, by
   * first computing which positions are paired;
   * 2) assign to each unpaired (ie empty) bag one of the positions left unpaired.
   */
  vector<bool> occupied;
  for (int i=0;i<bags.size();i++)
  {
    Bag * b = bags[i];
    vector<int> indices = b->getIndices();
    for (int j=0;j<indices.size();j++)
    {
      int k = indices[j];
      while(occupied.size()<=k)
      {
          occupied.push_back(false);
      }
      occupied[k] = true;
    }
  }
  int lastEmpty=0;
  for (int i=0;i<bags.size();i++)
  {
      Bag * b = bags[i];
      if ((b->getIndices().size()==0))
      {
        for(;lastEmpty<occupied.size();lastEmpty++)
        {
            if (!occupied[lastEmpty])
            {
                occupied[lastEmpty] = true;
                b->addIndex(lastEmpty);
                occupied.push_back(false);
                break;
            }
        }

      }
  }

  /*
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
  */
  
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


void TreeDecomposition::addLoopsRec(Bag * b, vector<Loop* > & structures){
    vector<int> removedLoops;
    for(int i=0;i<structures.size();i++){
        if (b->contains(structures[i]->getIndices())){
            removedLoops.push_back(i);
        }
    }
    for(int i=removedLoops.size()-1;i>=0;i--){
        int j = removedLoops[i];
        Loop * l = structures[j];
        if (DEBUG) cerr << "  Assigning loop "<<l<<" to "<<b<<endl;
        b->addLoop(l);
        structures.erase(structures.begin() + j);
    }
    vector<Bag*> children = b->getChildren();
    for(int c=0;c<b->children.size();c++){
        addLoopsRec(children[c],structures);
    }
}

void TreeDecomposition::addLoops(vector<Loop* > structures){
    vector<Loop *> backup;
    for (int i=0;i<structures.size();i++){
        backup.push_back(structures[i]);
    }
    for (unsigned int i=0;i<roots.size();i++){
        addLoopsRec(bags[roots[i]],backup);
        if (backup.empty())
        { break; }
    }
    assert(backup.empty());
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

#define USE_LIBTW 1
#define USE_FOX_EPSTEIN 0
#define USE_STRASSER 0
#define USE_GASPERS_ET_AL 0
#define USE_BANACH_ET_AL 0
#define USE_JOGLEKAR_ET_AL 0

TreeDecomposition* TDLibFactory::makeTD(vector<SecondaryStructure *>& structures)
{
    string tmpfilein = "./tmp.dgf";
    string tmpfileout = "./tmp.td";
    string baselib = "../lib";
    string basebin = "../../bin";

    TreeDecomposition * result = new TreeDecomposition();
		
    TreeDecomposition * td = new TreeDecomposition();
    saveAsDGF(structures,tmpfilein,1);


    // TODO: To be replaced by a modular system for calling TD tools (or, even better a call to an external C++ API )
    if (USE_LIBTW)
    {
        // library LIBTW
        string cmd = string("java -cp "+baselib+"/treewidth-java nl.uu.cs.treewidth.TreeDecomposer 1 ./"+tmpfilein + string(" ./") + tmpfileout + string(" tmp.dot >out.tmp "));
        system(cmd.c_str());
        //new TreeDecomposition();
        transFileFormat(tmpfileout);
        td->loadFromFile(tmpfileout);
        result->copyObj(td);
        //remove(tmpfilein.c_str());
        //remove(tmpfileout.c_str());
        //cout << "LIBTW treewidth: " << result->tw << endl;
		td->reset();		
    }

    saveAsDGF(structures,tmpfilein,2);

    //by Fox-Epstein (Brown University)
    if (USE_FOX_EPSTEIN)
    {
        string cmd1 = "timeout --signal=SIGTERM 5s "+baselib+"2016-pace-challenge-master/tw-heuristic " + tmpfilein + string(" > ") + tmpfileout;
		system(cmd1.c_str());
		td->loadFromFile(tmpfileout);
        //cout << "Fox-Epstein --- treewidth: " << td->tw << endl;
        remove(tmpfileout.c_str());
        if((td->tw!=0) &&  (td->tw < result->tw)){
            result->reset();
            result->copyObj(td);
        }
        td->reset();
    }
    
    // by Strasser (Karlsruhe Institute of Technology)
    if (USE_STRASSER)
    {
        string cmd1 = "timeout --signal=SIGTERM 5s "+baselib+"flow-cutter-pace16-master/flow_cutter_pace16 " + tmpfilein + string(" > ") + tmpfileout;
		system(cmd1.c_str());
		td->loadFromFile(tmpfileout);
        //cout << "Strasser --- treewidth: " << td->tw << endl;
		remove(tmpfileout.c_str());
        if((td->tw!=0) &&  (td->tw < result->tw)){
            result->reset();
            result->copyObj(td);
        }
        td->reset();
    }

    
    // by Gaspers, Gudmundsson, Jones, Mestre, Rummele  (UNSW and University of Sidney)
    if (USE_GASPERS_ET_AL)
    {
        string cmd1 = baselib+"/pace2016-master/";
        string cmd2 = "./tw-heuristic < "+basebin+"/" + tmpfilein + string(" > ") + basebin+"/" + tmpfileout;
        string cmd3 = basebin+"/";
        chdir(cmd1.c_str());
        system(cmd2.c_str());
        chdir(cmd3.c_str());
        td->loadFromFile(tmpfileout);
        //cout << "Gaspers, Gudmundsson, Jones, Mestre, Rummele --- treewidth: " << td->tw << endl;
		remove(tmpfileout.c_str());
        if((td->tw!=0) &&  (td->tw < result->tw)){
            result->reset();
            result->copyObj(td);
        }
        td->reset();
    }
    
    
    // by Bannach, Berndt, Ehlers (Luebeck University)
    if (USE_BANACH_ET_AL)
    {
        string cmd1 = baselib+"Jdrasil-master/tw-heuristic < " + tmpfilein + string(" > ") + tmpfileout;
        system(cmd1.c_str());
        td->loadFromFile(tmpfileout);
        //cout << "Bannach, Berndt, Ehlers --- treewidth: " << td->tw << endl;
        remove(tmpfileout.c_str());
        if((td->tw!=0) &&  (td->tw < result->tw)){
            result->reset();
            result->copyObj(td);
        }
        td->reset();
    }
    
    // by Joglekar, Kamble, Pandian (IIT Madras)
    if (USE_JOGLEKAR_ET_AL)
    {
        string cmd1 = "timeout --signal=SIGTERM 10s "+baselib+"/pacechallenge-master/tw-heuristic < " + tmpfilein + string(" > ") + tmpfileout;
        system(cmd1.c_str());
        td->loadFromFile(tmpfileout);
        //cout << "Joglekar, Kamble, Pandian --- treewidth: " << td->tw << endl;
        remove(tmpfileout.c_str());
        if((td->tw!=0) &&  (td->tw < result->tw)){
            result->reset();
            result->copyObj(td);
        }
        td->reset();
    }
    if (DEBUG)
    {
        cerr << "TreeDecomposition: "<<endl;
        result->show();
        //cerr << endl;
    }
    return result;
}


void TDLibFactory::transFileFormat(string path)
{
	string result;
	ifstream input(path.c_str());
	if(input.is_open()){
		string line;
		while(getline(input, line)){
			
			std::regex node("[ \\t]*bag(\\d+) *\\[label=\" *(.+?) *\"\\]$");
			std::smatch fnode;
			std::regex_search(line, fnode, node);
		
			if(fnode.size() > 2){
				result.append("b ");
				result.append(fnode[1]);
				result.push_back(' ');
				string bag_content = fnode[2];
				string token;
				vector<string> data = split(trim(bag_content),' '); 
				for(int i=0; i<data.size(); i++){
					string tmp = data[i];
					if(tmp[0] == '\\'){
						result.append(tmp.substr(2,tmp.size()-2));
						result.push_back(' ');
					} else {
						result.append(tmp);
						result.push_back(' ');
					}
				}
				result.push_back('\n');
			} else {
				std::regex arc("[ \\t]*bag(\\d+) -- bag(\\d+)$");
				std::smatch farc;
				std::regex_search(line, farc, arc);
				if(farc.size() > 2){
					string bag1 = farc[1];
					string bag2 = farc[2];
					result = result + bag1 + " " + bag2 + " \n";
				}
			}
		}
	}
	input.close();

	ofstream output(path.c_str());
	output << "s td " << endl;
	output << result;
	output.close();

}
