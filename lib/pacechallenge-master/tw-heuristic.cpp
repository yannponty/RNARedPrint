#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <map>
#include <set>
#include <queue>
#include <limits>
#include <algorithm> 
#include <stdlib.h>
#include <utility> // pair
#include <sys/time.h>
#include <signal.h>
#define DEBUG 0
using namespace std;
class Set{

	public :

	static set<int> intersect (set<int> s1,set<int> s2)
	{
		set<int> s3,smallSet,largeSet;
		set<int>::iterator sit;
		
		if (s1.size()<s2.size())
		   {
			   smallSet=s1;
			   largeSet=s2;
		   }
		   else
		   {
			   smallSet=s2;
			   largeSet=s1;
		   }
		for (sit=smallSet.begin(); sit!=smallSet.end(); ++sit)  
		     if(largeSet.find(*sit)!=largeSet.end()) 
		        s3.insert(*sit);
		 
		 return s3;  
		
	}
	static set<int> setMinus(set<int> setA,set<int> setB){
		set<int> returnSet;
		set<int>::iterator it1,it2,	t1=setB.begin(), t2=setB.end();
		for(it1 = setA.begin(), it2 = setA.end() ; it1!=it2 ; it1++){
			if(setB.find(*it1)==setB.end()){
				returnSet.insert(*it1);
			}		
		}
		return returnSet;
	}
        static int isSubset (set<int> s1,set<int> s2)
    {
	if(s1.size() > s2.size()) return -1;
	
	unsigned int cnt=0;
	
	set<int> :: iterator sit;
	for(sit=s1.begin();sit!=s1.end();sit++)
	{
		if(s2.find(*sit)==s2.end())
		return -1;      //s1 is not subset
		else cnt++;
	}
	if(s1.size()==s2.size())
	{
	  if(s1.size()==cnt) return 0; //both sets equal
	  else return 1; //s1 is proper subset
	}  
	else return 1; //s1 is proper subset
     }
	static void printSet(set <int> s){
		set<int>::iterator it,end;
		cout<< "NB set: " ;
		for(it = s.begin(), end = s.end() ; it != end ; it++){
			cout<< *it << " ";
		}
		cout << endl;
	}
};

class Graph
{
	public:
		unsigned int n;
		unsigned int m;
		
		unsigned int max_clique_size;
		std::map<int, set <int>  > retMapOfSet;
		vector <int> *adjList;
        set <int>  nodes;
        set < set <int> > returnSoS;
        bool cycle_found;
        std::map<int, int> vertexColour;
        
        //set of Chordless Cycles
        set<vector<int> > C;
        //set of Chordless Paths
        set<vector<int> > T;
        //labelling
        //map<int,int> labels;
	int* labels;
	int no_of_components;
        int* blocked;
	
        Graph()
        {
		}
        Graph(int n , vector <int> * a,set <int>  nodes)
        {	
			
			this->n = n;
			this->nodes=nodes;
			max_clique_size=1;
			blocked=new int[n+1];
			labels=new int[n+1];
			
			adjList = new vector <int> [n+1];
            cycle_found=false;
			for(int i=1; i <= n; i++)
                         {
                         
                          adjList[i]=a[i];
                         }
		}
		void printAdjList(){
			
			for(unsigned int i=1; i <= n; i++){
				 cout  <<i << "\t#";
				for( std::vector<int>::iterator it=adjList[i].begin(), end =adjList[i].end(); it!=end; it++){
					 cout << " "<< *it;
				}
				cout << endl;
                                 
			}
                        
                  }
				
		
    void DFS()
{
	std::set<int>::iterator it,end;
	
	for (it=nodes.begin(),end=nodes.end();it!=end;it++)
			vertexColour[*it]=0;
		if(DEBUG) cout<<"DFS Starts"<<endl;	
		unsigned int i=0;
		std::map<int, int>::iterator iter;
    for (it=nodes.begin(),end=nodes.end();it!=end;it++)	
    	{
         iter=vertexColour.find(*it);
         if(DEBUG) cout<<"vertex"<<*it<<" vcolor "<<iter->second<<endl;
         if( iter->second==0)
         {
             i=i+1;
              if(DEBUG) cout<<"Component  "<<i<<endl;
             DFS_VISIT(*it,i);
	}
	no_of_components=i; 
	}
                          
}
void DFS_VISIT(int v,unsigned int setno)
{
	std::vector<int>::iterator it,end;
	std::map<int, set<int> >::iterator iter;
	set<int> currComponent;
	 if(DEBUG) cout<<v<<endl;
	 if(DEBUG) cout<<"Map size="<<retMapOfSet.size()<<"  setno="<<setno<<endl;
	if (retMapOfSet.size()==setno)
	{
		iter=retMapOfSet.find(setno);
		currComponent=iter->second;
	}
		currComponent.insert(v);
		 if(DEBUG) cout<<"currComponent size="<<currComponent.size()<<endl;
		retMapOfSet[setno]=currComponent;
	
	vertexColour.erase(v);
	vertexColour[v]=1;
	//cout<<"Vertex Color changed:"<<v<<"  "<<vertexColour[v]<<"\n";	
	std::map<int, int>::iterator i;
	for (it=adjList[v].begin(),end=adjList[v].end();it!=end;it++)
                                  {
                                  	i=vertexColour.find(*it);
                                  	if(i->second==2)
									{
									 //cout<<"Cycle Detected: "<<v<<"  "<<i->first<<"\n"	;	
									 cycle_found=true;
								    }
                                  	if(i->second==0)
									  DFS_VISIT(*it,setno); 
									
								  }
	vertexColour.erase(v);
	vertexColour[v]=2;	
	//cout<<"Vertex Color changed:"<<v<<"  "<<vertexColour[v]<<"\n";					  
}	

/*set<int> intersect (set<int> s1,set<int> s2)
{
	set<int> s3,smallSet,largeSet;
	set<int>::iterator sit;
	
	if (s1.size()<s2.size())
	   {
		   smallSet=s1;
		   largeSet=s2;
	   }
	   else
	   {
		   smallSet=s2;
		   largeSet=s1;
	   }
	for (sit=smallSet.begin(); sit!=smallSet.end(); ++sit)  
	     if(largeSet.find(*sit)!=largeSet.end()) 
	        s3.insert(*sit);
	 
	 return s3;  
	
}*/
vector<int> lex_bfs()
{
	if(DEBUG) cout<<"inside lex_bfs\n";
	vector<set<int> > Q;
	//int node_2_Q[n];
	Q.push_back(nodes);
	//int setno=0;
	/*for(unsigned int j=0;j<n;j++)
	node_2_Q[j]=setno;*/
	
	map<int,int> sigma,siginv;
	
	std::vector<int>::iterator it,end;
	int size=nodes.size();
	for(int i=size;i>0;i--)
	{
		if(DEBUG) cout<<"i="<<i<<"\n";
		while(Q.back().size()==0) 
			Q.pop_back();
		set<int> curr_set=Q.back();
		int v=*(curr_set.begin());
		if(DEBUG) cout<<"v="<<v<<"\n";
		//(*(Q.end())).erase(curr_set.begin());
		curr_set.erase(v);
		Q.pop_back();
		Q.push_back(curr_set);
		sigma[i]=v;
		siginv[v]=i;
		if(DEBUG) cout<<"Debug 1\n";
		set<int> curr_adj_ver;
		for (it=adjList[v].begin(),end=adjList[v].end();it!=end;it++)
		{
			if(DEBUG) cout<<"adj v="<<*it<<"\n";
			if(siginv.find(*it)==siginv.end())
				curr_adj_ver.insert(*it);
		}
		unsigned int no_adj_vertices=curr_adj_ver.size();
		 if(DEBUG) cout<<"no_adj_vertices="<<no_adj_vertices<<"\n";
		vector<set<int> > Q1;
		vector<set<int> > ::iterator Q_iter;
		Q1.assign(Q.begin(),Q.end());
		 if(DEBUG) cout<<"Q1.size()="<<Q1.size()<<"\n";
		 Q_iter=Q.begin();
		 int k=0;
		 unsigned int q1Size=Q1.size();
		for(unsigned int j=0;j<q1Size && no_adj_vertices>0;j++)
		{
			 if(DEBUG) cout<<"j="<<j<<"\n";
			set<int> s=Set :: intersect(Q1.at(j),curr_adj_ver);
			if(s.size()!=0)
			{				
				if(k+j+1==Q.size())	 Q.push_back(s);
				else
				{
				//Q_iter++;
					
			    Q.insert(Q.begin()+k+j+1,s);
			    	
				//Q_iter--;
					
			    }		
				//setno++;
				set<int> ::iterator iter;
				for(iter=s.begin();iter!=s.end();iter++)
				{
				   //node_2_Q[*iter]=setno;
				   Q.at(j+k).erase(*iter);
				   no_adj_vertices--;
			   }
			   k++;
			}
			//else Q_iter++;
		}
	}
	vector<int> ordering;
	for( int i=1;i<=size;i++)
	ordering.push_back(sigma.at(i));
	
	return ordering;
}

pair <bool , vector <int> > isChordal()
{
	// cout<<"inside isChordal Function\n";
	vector<int> vec=lex_bfs();
	//vector<int> vec1;
	//vec1.push_back(-1);
	//cout<<"After lex_bfs_comp\n";
	int size=nodes.size();
	if(DEBUG) 
	{
	for( int k=0;k<size;k++)
	cout<<vec.at(k)<<"\n";
        }
	map<int,vector<int> > A;
	map<int,int> siginv;
	for( int i=1;i<=size;i++)
	siginv[vec.at(i-1)]=i;
    std::vector<int>::iterator it,end;
    pair <bool , vector <int> >  retPair;
	//bool flag=true;
	for( int i=0;i<size-1 ;i++)
	{
		//cout<<"i="<<i<<"\n";
		int v=vec.at(i);
		set<int> X;
		for (it=adjList[v].begin(),end=adjList[v].end();it!=end;it++)
		if(siginv[v]<siginv[*it])
		X.insert(*it);
		set<int> :: iterator sit;
		if(X.size()!=0)
		{
			
			int min=numeric_limits<int>::max();
			int u;
			for(sit=X.begin();sit!=X.end();sit++)
			{
			//cout<<"Debug1\n";	
                           if(siginv[*sit]<min)
				{
				min=siginv[*sit];
				u=*sit;
			    }
			}
			X.erase(u);
			if(A.find(u)==A.end())
				A[u].assign(X.begin(),X.end());
			else
			    A[u].insert(A[u].begin(),X.begin(),X.end());
		}
		//cout<<"Debug2\n";
		set<int> s1,s2,s3;
		s1.insert(A[v].begin(),A[v].end());
		s2.insert(adjList[v].begin(),adjList[v].end());
		for(sit=s1.begin();sit!=s1.end();sit++)
		   if(s2.find(*sit)==s2.end())
		  	 s3.insert(*sit);
		//cout<<"Debug3\n";
		//cout<<"s3.size="<<s3.size()<<"\n";
		if(s3.size()!=0) 
		{
			 //cout<<"GRaph Not yet Chordal\n";
				retPair = make_pair (false,vec);
				return retPair;
		
				
		}
		//cout<<"Debug4\n";
	}
	//cout<<"Debug5\n";
	
	if(DEBUG) cout<<"Graph is now Chordal\n";
	retPair = make_pair (true,vec);		
	if(DEBUG) {
	for( int k=0;k<size;k++)
	cout<<vec.at(k)<<" ";
        
	cout<<"\n";
	}
	 //cout<<"Exiting isChordal\n";
	return retPair;
}


void degreeLabelling()
{
	/*map<int,int> color;
	map<int,int> degree;*/
	int *color=new int[n+1];
	int *degree=new int[n+1];
	set<int> ::iterator sit;
	for(sit=nodes.begin();sit!=nodes.end();sit++)
	{
		int curr_vertex=*sit;
		color[curr_vertex]=0;
		degree[curr_vertex]=adjList[curr_vertex].size();
	}
	int mindegree,v=0;
	for(unsigned int i=1;i<=n;i++)
	{
		mindegree=n;
		for(sit=nodes.begin();sit!=nodes.end();sit++)
	{
		int curr_vertex=*sit;
		if(color[curr_vertex]==0 && degree[curr_vertex] <mindegree)
		{
			v=curr_vertex;
			mindegree=degree[curr_vertex];
		}
   }
	labels[v]=i;
	color[v]=1;
	
	vector<int>::iterator vit;
	for(vit=adjList[v].begin();vit!=adjList[v].end();vit++)
	{
		int curr_vertex=*vit;
		if(color[curr_vertex]==0)
		   degree[curr_vertex]=degree[curr_vertex]-1;
	}
	
	}
}
void findTriples()
{
	set<int> ::iterator sit;
	for(sit=nodes.begin();sit!=nodes.end();sit++)
	{
		int u=*sit;
		
		vector<int>::iterator vit1,vit2;
		int x,y;
	    for(vit1=adjList[u].begin();vit1!=adjList[u].end();vit1++)
	{	 x=*vit1;
	         for(vit2=adjList[u].begin();vit2!=adjList[u].end();vit2++)
	         {
				 y=*vit2;
				 //ℓ(u) < ℓ(x) < ℓ(y)
				 if(labels[u]<labels[x] && labels[x]<labels[y] 
					&& find(adjList[x].begin(),adjList[x].end(),y)==adjList[x].end()) //test
				 {					 
					 vector<int> vec;
					 vec.push_back(x);
					 vec.push_back(u);
					 vec.push_back(y);
                  			 T.insert(vec);
					/*else
					   C.insert(vec); */  
				 }
	         }
	}
		
	}
}
void blockNeighbours(int v)
{
	vector<int> adj_v=adjList[v];
	unsigned int size=adj_v.size();
	for(unsigned int i=0;i<size;i++)
	    blocked[adj_v[i]]= blocked[adj_v[i]]+1;
}	

void unblockNeighbours(int v)
{
	vector<int> adj_v=adjList[v];
	unsigned int size=adj_v.size();
	for(unsigned int i=0;i<size;i++)
	    if(blocked[adj_v[i]] >0)
	    blocked[adj_v[i]]= blocked[adj_v[i]]-1;
}	
void CC_Visit(vector<int> p,int l)
{
	int last_p=p[p.size()-1];
	
	blockNeighbours(last_p);	
	vector<int> adj_v=adjList[last_p];
	unsigned int size=adj_v.size();
	for(unsigned int i=0;i<size;i++)
	{
		int curr_neigh=adj_v[i];
		
		if(labels[curr_neigh] > l && blocked[curr_neigh]==1)
		{
			vector<int> pdash=p;
			pdash.push_back(curr_neigh);
			
			if(find(adjList[curr_neigh].begin(),adjList[curr_neigh].end(),p[0])!=adjList[curr_neigh].end())
			   C.insert(pdash);
			else
			    CC_Visit(pdash,l);  
		}
	}
	
	unblockNeighbours(last_p);
}
void ChordlessCycles()
{
	degreeLabelling();
	findTriples();
	set<vector<int> > :: iterator sit;
	
    
    for(unsigned int i=1;i<=n;i++)
        blocked[i]=0;
    
    while(!T.empty())
    {
		vector<int> p=*(T.begin());
		int u=p[1];
		T.erase(T.begin());
		
		blockNeighbours(u);
		CC_Visit(p,labels[u]);
		unblockNeighbours(u);
   }
	if(DEBUG)
	{
	cout<<"Chordless Cycles\n";
	
	for(sit=C.begin();sit!=C.end();sit++)
	{
		vector<int> vec=*sit;
		for(unsigned int i=0;i<vec.size();i++)
		    cout<<vec[i]<<" ";
		cout<<"\n";
	}
	}
	
}
/*void kill_funda_cycles()	
 {
	 set<vector<int> > :: iterator sit;
	 for(sit=C.begin();sit!=C.end();sit++)
	{
		vector<int> vec=*sit;
		int size=vec.size();
                //cout<<"cycle size="<<size<<"\n";
                
                
                int no_of_chords=size-3;
                int interval=2;
                int i,j=0;                           
                while(no_of_chords)
                {
                //cout<<"no_of_chords="<<no_of_chords<<"\n";
                //cout<<"interval size="<<interval<<"\n";
		for( i=j;i<size-interval+1;i=i+interval)
		{
			if(find(adjList[vec[i]].begin(),adjList[vec[i]].end(),vec[(i+interval)%size])==adjList[vec[i]].end())
			{
				adjList[vec[i]].push_back(vec[(i+interval)%size]);
				adjList[vec[(i+interval)%size]].push_back(vec[i]);
				if(DEBUG) cout<<"Added Chord "<<vec[i]<<" "<<vec[(i+interval)%size]<<"\n";
                                
			}
                        no_of_chords--;
                        if(no_of_chords==0) break;
                        
		}
                interval*=2;
                if(no_of_chords!=0)
                   if(i==size) j=0;
                   else
                   {
                   if(find(adjList[vec[i]].begin(),adjList[vec[i]].end(),vec[(i+interval-1)%size])==adjList[vec[i]].end())
			{
				adjList[vec[i]].push_back(vec[(i+interval-1)%size]);
				adjList[vec[(i+interval-1)%size]].push_back(vec[i]);
				if(DEBUG) cout<<"Added Chord "<<vec[i]<<" "<<vec[(i+interval)%size]<<"\n";
                                
			}
                        j=(i+interval-1)%size;
                        no_of_chords--;
                   }
               }
                          
	
 }
}*/
void kill_cycles()	
 {
	 set<vector<int> > :: iterator sit;
	
	int cnt_chords_added=0;
	vector<vector<int> > chords_added;	//pair<>int,int
	 for(sit=C.begin();sit!=C.end();sit++)
	{
		bool flag=true;
		
		vector<int> cycle=*sit;
		cnt_chords_added=chords_added.size();
		//if(cnt_chords_added)
		//{
		for(int i=0;i<cnt_chords_added && flag;i++)
		{
		vector<int> curr_chord=chords_added[i];
		if((find(cycle.begin(),cycle.end(),curr_chord[0])!=cycle.end()) &&
			(find(cycle.begin(),cycle.end(),curr_chord[1])!=cycle.end()))
		  flag=false;
		}
		//}
		if(!flag) continue;
		int size=cycle.size();
                //cout<<"cycle size="<<size<<"\n";
		/*int *adjCount=new int[size];
                for(int i=0;i<size;i++)
                    adjCount[i]=0;

                for(int i=0;i<size;i++)
                   for(int j=i+1;j<size;j++)
                   if(find(adjList[cycle[i]].begin(),adjList[cycle[i]].end(),cycle[j])!=adjList[cycle[i]].end())
                      {
                       adjCount[i]++;
                       adjCount[j]++;
                      }
		int sum=0;
		for(int i=0;i<size;i++)
                    sum+=adjCount[i];
		if(sum>=(2+(size-2)*4)) continue;*/
		int no_of_chords=size-3;
		int num1=rand()%size;
		int v1=cycle[num1];
		set<int> s1,s2,s3;
		s1.insert(cycle.begin(),cycle.end());
		
		while(no_of_chords--)
		{
		//cout<<"no_of_chords="<<no_of_chords<<"\n";
		int num2=rand()%(size-3);
		int v2=cycle[(num1+2+num2)%size];
		//cout<<"Chord to be added "<<v1<<" "<<v2<<"\n";
		
		if(find(adjList[v1].begin(),adjList[v1].end(),v2)==adjList[v1].end())
		{
		adjList[v1].push_back(v2);
		adjList[v2].push_back(v1);
		//cout<<"Added Chord "<<v1<<" "<<v2<<"\n";
		//cnt_chords_added++;
		vector<int> chord;
		chord.push_back(v1);
		chord.push_back(v2);
		chords_added.push_back(chord);
		}		
		v1=v2;
		num1=(num1+2+num2)%size;
		}
		
        }
 }
Graph return_GK(set<int> k)
{
	vector <int> adjList_gk[n+1];
	vector<int>::iterator it,end;
	for(int i=1; i <= (int)n; i++)
      {
			if(k.find(i) != k.end())
			{
				for( it=adjList[i].begin(), end =adjList[i].end(); it!=end; it++)
				{
					//cout<<"return_SGK *it1  "<<*it1<<"\n";
					if(k.find(*it) != k.end())
					  adjList_gk[i].push_back(*it);
				}
				
			}                                 
	}
	Graph gk(n,adjList_gk,k);
	return gk;
}
};

class Tree{
public:
	Tree(){
		numberOfBags = 0;
		treewidth = 0;
	}

	Tree(int n){
		numberOfBags = 0;
		treewidth = 0;
		numberOfVertices = n;
	}


	int numberOfBags;
	int treewidth;
	vector < set<int> > bag;
	vector < pair <int , int> > treeEdge;
	int numberOfVertices;
	
	int getNumberOfBags(){
		return numberOfBags;
	}
	
	vector < set<int> > getBags(){
		return bag;
	}
		
	vector < pair <int , int> >  getTreeEdges(){
		return treeEdge;
	}

	set <int> getHigherNeighbourOf(vector <int> v , int index,  Graph g){
		set <int> retSet, higherVs;		
		
		int vertex = v.at(index);
		int size = v.size();
		for(int i = index + 1  ; i < size  ; i++ ){
			higherVs.insert(v.at(i));
		}
		vector<int>::iterator it,end;
		set <int>  adjSetOfV(g.adjList[vertex].begin() , g.adjList[vertex].end() ) ;
		
		retSet = Set::intersect(adjSetOfV,higherVs);
		return retSet;
	}

	int isInBagAlready(set <int> neighbours){
		int retInt = -1;
		int bagNum = 0;

		vector < set<int> >::iterator it,end;
		it = bag.begin();
		end = bag.end();
		
		for( it=bag.begin(), end = bag.end() ; it != end ; it++, bagNum++ ){
			set <int> newSet = *it ;
			if (newSet == neighbours ){
				retInt = bagNum;
				break;
			}			 
		}
		return retInt;
	}

	int getBagContainingAnyVertexOf( set<int> verSet){
		
		int retBagVal=-1;
		int bagNum = 0;
		 
		vector < set<int> >::iterator it,end;
		it  = bag.begin();
		end = bag.end();

		for( it=bag.begin(), end = bag.end() ; it != end ; it++, bagNum++ ){
			set <int> newSet = *it ;
			
			if(Set::intersect(verSet,newSet).size() == verSet.size()){
				retBagVal = bagNum;
				break;
			}
		}
		return retBagVal;
	}


	void computeTreeDecomposition( vector<int> v , Graph g){
		vector<int> :: iterator vit,end;
		int vertex,i;

		i = v.size()-1;
		vertex = v.at(i);
		
		// create bag for the initial 
		set <int> s ;
		s.insert(vertex);
		bag.push_back(s);
		numberOfBags++;
		--i;

		for( ; i >= 0 ; i-- ){
			//cout << "for " <<  i << endl;
			vertex = v.at(i);
			set <int> higherNeighbours = getHigherNeighbourOf(v,i,g);
			int bagNum = isInBagAlready(higherNeighbours);
			///<< vertex << "th isInBagAlready? " << bagNum << endl;
			///Set::printSet(higherNeighbours);
			// if heigher neighbours are in some bag add the vertex to it
			if(bagNum >= 0){
				///cout<< "same bag " << bagNum << endl;
				(bag.at(bagNum)).insert(vertex);
			} // otherwise create a new bag out of v + neighbours and link to any of the node containing it
			else{								
				// add the tree edge.				
				int OldBagNum = getBagContainingAnyVertexOf(higherNeighbours);
				///cout<< i<<" -itr " << OldBagNum << endl;
				higherNeighbours.insert(vertex);
				treeEdge.push_back(make_pair(OldBagNum , numberOfBags));
				bag.push_back(higherNeighbours) ;
				numberOfBags++;
			}

		}
		//printAns();
	}
	
	
	void computeDefaultDecomposition(Graph g){		
		numberOfBags=1;
		bag.push_back(g.nodes) ;
	}

	int getTreeWidthPlusOne(){
		
		//if(treewidth != 0)
		//	return treewidth;
		
		int retVal = numeric_limits<int>::min();
		vector < set<int> >:: iterator it, end;
		for( it=bag.begin() ,end = bag.end(); it != end ; it++){
			int val = (*it).size();
			retVal = max(retVal,val );
		}
		return retVal;
	}
	void printAns(){
		vector < set<int> >:: iterator it, end;
		//cout<<"c Before s td\n";
		cout << "s td " << numberOfBags << " " << getTreeWidthPlusOne() <<" " <<numberOfVertices << endl;
		//cout<<"c After s td\n";
		int i=1;
		//cout<<"c Before bag loop\n";
		for( it=bag.begin(), end = bag.end() ; it != end ; it++, i++){
			//cout<<"c Inside bag loop\n";
			set < int > inBag = *it;
			set < int > ::iterator sit, send;
			cout << "b " << i ; 
			for( sit=inBag.begin(),send = inBag.end() ; sit != send ; sit++){
				 cout << " " << *sit ;
			}
			cout << endl;
		}
		//cout<<"c Outside bag loop\n";
		vector< pair <int , int>  >::iterator it1, vend;
		for( it1 = treeEdge.begin() , vend = treeEdge.end() ; it1 != vend ; it1++ ){
			
			cout<< (*it1).first+1 << " " << (*it1).second+1 << endl;
		}
		//cout<<"c Outside Edge loop\n";
	}
	void computeBigTree(vector <Tree*> smallTrees){
		vector <Tree*> :: iterator it, end;
		
		//dummy bag to connect all the trees
		set <int> emptySet;
		bag.push_back(emptySet);
		numberOfBags=1;
		
		for( it = smallTrees.begin(), end = smallTrees.end(); it != end ; it++ ) {
			Tree *tree = *it;
			int treesize = tree->getNumberOfBags();
			vector < pair <int , int> > edges =  tree->getTreeEdges();
			vector < pair <int , int> > ::iterator vec_it = edges.begin(), vec_end = edges.end();
			
			treeEdge.push_back(make_pair(0, numberOfBags));
			for(  ; vec_it != vec_end ; vec_it++ ){
				pair <int,int> v = *vec_it;
				treeEdge.push_back(make_pair( v.first + numberOfBags, v.second + numberOfBags));
			}
			
			numberOfBags+=treesize;
			vector < set<int> > bags =  tree->getBags();
			vector < set<int> > ::iterator set_it = bags.begin(), set_end = bags.end();
			for(  ; set_it != set_end ; set_it++ )
				bag.push_back(*set_it);
		} 
	}
};
Tree *mintree,*minComptree;
int mintw=numeric_limits<int>::max();
void handleSignal(int signalNum)
{
    /*if (signalNum == SIGINT){
		//int max_bag_size= 2;
		if(DEBUG) cout << "Received SIGINT!\n";
		struct timeval tv;
		gettimeofday(&tv, NULL);
		unsigned long long msSinceEpoch = (unsigned long long)(tv.tv_sec) * 1000
		   + (unsigned long long)(tv.tv_usec) / 1000;
		cout << "c status " << mintw+1 << " " << msSinceEpoch << endl;
		
	}*/
    if (signalNum == SIGINT || signalNum == SIGTERM)
    {
        if(DEBUG) cout << "Received SIGTERM!\n";
	//if(mintree)
        mintree->printAns();
        //else cout<<"No Output\n";
	//cout << "c Printing Trivial Decomposition as timeout occurred before computing tree decomposition\n" ;
        exit(1);
    }
        
}	
int main()
{
int total_nodes, total_edges;
	int n,m;
	
	string lineStr , keyword;
	string problem;
	int node1, node2;
        signal(SIGTERM, handleSignal);
	signal(SIGINT, handleSignal);
	// To scan throught all the comments 
	do{
		getline(std::cin, lineStr);
			
	}while(lineStr[0]!='p');
	
	istringstream words(lineStr);
	words >>  keyword; //  p
	words >>  problem; // tw or fvs
	words >> total_nodes;
	n = total_nodes;	
	words >> total_edges;
	m = total_edges;
	
	if (n == 0){
		cout << "s td 0 0 0" << endl;
		exit(0);
	}
	
	int t = m;
	vector <int> adjList[n+1];
	set <int>  input_nodes;
	for(int j=1;j<=n;j++) 
		input_nodes.insert(j);
	if(DEBUG) cout << "n=" << n << " m=" << m << endl;
	//reading m lines 
	while(t--){		
		
		// I am assuming no "c" in between those m lines
		
		getline(std::cin, lineStr);
		if(lineStr[0]=='c'){
			t++;
			continue;
			
		}	
		istringstream edges(lineStr);
		edges >> node1;
		edges >> node2;
		if(DEBUG) cout << node1 << " " << node2 << endl;
		if (node1<=n && node2<=n){
			//input_nodes.insert(node1);
			//input_nodes.insert(node2);
			std::vector<int>::iterator tmp_it1,tmp_it2;
			tmp_it1 = find (adjList[node1].begin(), adjList[node1].end(),  node2);
			tmp_it2 = find (adjList[node2].begin(), adjList[node2].end(),  node1);
			if( tmp_it1==adjList[node1].end() && tmp_it2==adjList[node2].end()){
				adjList[node1].push_back(node2);
				adjList[node2].push_back(node1);
			}
		}		
	}
	
	Graph g(n,adjList,input_nodes);
	if(DEBUG) 
        {
	cout<<"Adjacency List\n";
	g.printAdjList();
	}
	vector<int> vOfNodes;
	vOfNodes.assign(input_nodes.begin(),input_nodes.end());
	mintree  = new Tree (n);
	mintree->computeDefaultDecomposition(g);
	g.DFS();
	int k;
	if(!g.cycle_found)
	k=1;
    	else 
	k=1000;
	/*g.ChordlessCycles();
	g.kill_cycles();
	if(g.isChordal().first)
	cout<<"Graph is now chordal\n";
	else
         {
	cout<<"Not yet Chordal\n";
        Graph g1(g.n,g.adjList,g.nodes);
        g1.ChordlessCycles();
         }*/
	 // cout<<"c no_of_components="<<g.no_of_components<<"\n";
	 // cout<<"c retMapOfSet.size()="<<g.retMapOfSet.size()<<"\n";
	pair <bool , vector <int> > op,minop;
	int tw;
	if(n==m && g.no_of_components==1)
	{
	if(DEBUG) cout<<"Graph has a Cycle. Tree-width of graph is "<<2<<"\n";
	do{
	g.ChordlessCycles();
	//cout<<"Cycles Computed\n";
	g.kill_cycles();
	//cout<<"Cycles Killed\n";
	g.C.clear();
        g.T.clear(); 
	op=g.isChordal();
        }while(!op.first);
	Tree *tree  = new Tree (n);
	tree->computeTreeDecomposition(op.second,g);
	tree->printAns();
        return 0;
        }

//	Graph* g1;
         if(g.no_of_components > 1)
	{
	k=50;
	vector<Tree *> component_trees;
	int mintw1=numeric_limits<int>::max();
	for(int i=1;i<=g.no_of_components;i++)
	{
	//cout<<"Component "<<i<<"\n";
	mintw1=numeric_limits<int>::max();
	set<int> curr_comp_nodes=g.retMapOfSet[i];
	Graph curr_comp=g.return_GK(curr_comp_nodes);
	int curr_nodes_size=curr_comp.nodes.size();
	//mintree=new Tree(curr_nodes_size);
		for(int j=0;j<k;j++)
		{
			Graph g2(curr_comp.n,curr_comp.adjList,curr_comp.nodes);
			 //cout<<"Iteration "<<j<<"\n";
         		do
        		{
         			//cout<<"Inside IsChordal Loop \n";
				g2.ChordlessCycles();
				//cout<<"Cycles Computed\n";
	 			g2.kill_cycles(); 
				//cout<<"Cycles Killed\n";
         			g2.C.clear();
         			g2.T.clear();  
        			op=g2.isChordal();
        		}while(!op.first);

        		Tree *tree  = new Tree (curr_nodes_size);
			tree->computeTreeDecomposition(op.second,g2);
			tw=tree->getTreeWidthPlusOne();
		//cout<<"Tree width="<<tw<<"\n";
			if(tw<mintw1)
			{
				mintw1=tw;
				minComptree=tree;
				minop=op;
				//g1=new Graph(g2.n,g2.adjList,g2.nodes);
			}
        		if(DEBUG) cout<<"Iteration "<<i<<" Tree width="<<mintw<<"\n";
		}
		//mintree->printAns();
		if(i==1) 
			mintw=mintw1;
		else if(mintw1 >mintw) 
			mintw=mintw1;
		component_trees.push_back(minComptree);
	}
//	int no_of_trees=component_trees.size();
	
	
	//~ for(int i=0;i<no_of_trees;i++)
		//~ component_trees[i]->printAns();
	
	Tree *bigTree = new Tree(n);
	bigTree->computeBigTree(component_trees);	
	bigTree->printAns();
	
	return 0;
	}
	
	
	
	
	
	for(int i=0;i<k;i++)
	{
	Graph g2(g.n,g.adjList,g.nodes);
	if(DEBUG) cout<<"Iteration "<<i<<"\n";
         do
        {
         g2.ChordlessCycles();
	//cout<<"Cycles Computed\n";
	 g2.kill_cycles(); 
	//cout<<"Cycles Killed\n";
         g2.C.clear();
         g2.T.clear();  
         op=g2.isChordal();
        }while(!op.first);

        Tree *tree  = new Tree (n);
	tree->computeTreeDecomposition(op.second,g2);
	tw=tree->getTreeWidthPlusOne();
	//cout<<"Tree width="<<tw<<"\n";
	if(tw<mintw)
	{
	mintw=tw;
	mintree=tree;
	minop=op;
//	g1=new Graph(g2.n,g2.adjList,g2.nodes);
	struct timeval tv;
	gettimeofday(&tv, NULL);
	unsigned long long msSinceEpoch = (unsigned long long)(tv.tv_sec) * 1000
		   + (unsigned long long)(tv.tv_usec) / 1000;
	cout << "c status " << mintw << " " << msSinceEpoch << endl;
	}
        if(DEBUG) cout<<"Iteration "<<i<<" Tree width="<<mintw<<"\n";
	}
	if(DEBUG) cout<<"Optimal Tree width="<<mintw<<"\n";
	//mintree->computeTreeDecomposition(minop.second,*g1);
	mintree->printAns();
	return 0;
}
