////////////////////////////////////////////////////////
// hvb++ 1.0
// Copyleft: Javier Rodríguez Laguna
// 080725-101112-150125-150603
// Graph class, either directed or undirected
#ifndef GRAPH_HEADER
#define GRAPH_HEADER
#include"Common.h"
#include"Matrix.h"
#include"Text.h"

class Graph
{
public:
     long N; // sites number
     long Nl; // number of links
     List *V;  // neighborhood structure
     List I1, I2; // starting and ending point of each link
     bool directed;  // if true, directed graph
     
     Graph();
     Graph(long N);
     Graph(const Graph &);
     ~Graph();
     void Start();
     void Create(long n);
     void Destroy();
     void Set_Directed();
     void Set_Undirected(); 

     long Add_Site(long n=0); // n is the stem of the site
     void Add_Sites(long n); // create n new sites
     bool Add_Link(long, long);
     void Remove_Site(long); 
     void Remove_Sites(const List &);
     void Remove_Link(long, long);
     void Remove_Link(long);
     void Clear(); // delete all links, but retain sites
     void Update_Index(); // Update the link index, neighborhood is canonical

     long Degree(long p) const; // number of neighbours of site p
     long Neighbour(long p, long k) const; // k-th neighbour of site p
     List Neighbours(long p) const;
     bool Is_Link(long, long) const;
     long Link_Index(long,long) const;
     void Link_Sites(long&,long&,long) const;

     Matrix Adjacency_Matrix() const;
     List   Connected_Component(long) const;
     List   Find_Path(long s1, long s2) const;
     List   Find_Path(const Vector &V, long s1, long s2) const;
     long   Distance(long s1, long s2) const;
     Graph  Minimum_Spanning_Tree() const;

     Graph& operator=(const Graph &);
     long   operator()(long) const;
     long&  operator()(long,long);
     long   operator()(long,long) const;

     void   Write() const;
     bool   Save(const char*) const;
     bool   Load(const char*);
};

	 
void Copy(Graph&, const Graph&);

//////////////////////////////////////////////////////////////////////
// A few concrete graphs
//////////////////////////////////////////////////////////////////////

Graph Linear_Graph(long N);
Graph Linear_Graph_PBC(long N);
Graph Square_Graph(long lx, long ly);
Graph Square_Graph_PBC(long lx, long ly);
Graph Complete_Graph(long N);
Graph Matrix_2_Graph(const Matrix &M);
Graph Points_2_Graph(const Matrix &M, double dcutoff);

/////////////////////////////////////////////////////////////////
// A few handy routines
/////////////////////////////////////////////////////////////////

Graph  Remove_Site(const Graph &G, long n);
Graph  Remove_Sites(const Graph &G, const List &);

#endif





