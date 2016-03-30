////////////////////////////////////////////////////////
// Hvb
// Copyleft: Javier Rodríguez Laguna
// 080725-100418-150125-150606
// Graph structure both directed and undirected

#include"Graph.h"

Graph::Graph()
{
     Start();
}

Graph::Graph(long n)
{
     Start();
     Create(n);
}

Graph::Graph(const Graph &G)
{
     Start();
     Copy(*this,G);
}

Graph::~Graph()
{
     Destroy();
}

void Graph::Start()
{
     N=0; Nl=0;
     V=(List*)NULL; 
     I1.Start();
     I2.Start();
     directed=false;
}

void Graph::Create(long n)
{
     Destroy();
     Add_Sites(n);
}

void Graph::Destroy()
{
     if (!N) return;
     for (long i=1;i<=N;i++)
	  V[i].Destroy();
     free(V);
     I1.Destroy();
     I2.Destroy();
     V=(List*)NULL;
     directed=false;
     N=0; Nl=0;
}

void Graph::Set_Directed()
{
     directed=true;
     Update_Index();
}

void Graph::Set_Undirected()
{
     directed=false;
     Update_Index();
}

// stem is 0 by default
long Graph::Add_Site(long stem)
{
     N++;
     if (N==1)
	  V=(List*)malloc(2*sizeof(List));
     else
	  V=(List*)realloc(V,(N+1)*sizeof(List));
   
     V[N].Start();
     if (stem!=0 && stem<=N-1)
	  Add_Link(stem,N);
     return N;
}

void Graph::Add_Sites(long n)
{
     V=(List*)realloc(V,(N+n+1)*sizeof(List));
     for (long i=1;i<=n;i++)
	  V[N+i].Start();
     N+=n;
}

bool Graph::Add_Link(long s1, long s2)
{
     if (s1==s2) return false;
     if (s1<1 || s1>N || s2<1 || s2>N) return false;
     if (Is_Link(s1,s2)) return false;
     V[s1].Append(s2);
     if (!directed) V[s2].Append(s1);
     I1 &= s1; // append site s1 to the I1 list
     I2 &= s2;
     Nl++;
     return true;
}

void Graph::Remove_Site(long s)
{
     if (s<1 || s>N) Error("Trying to remove non-existent site\n");
     // First, detach site from the rest
     for (long k=Degree(s);k>=1;k--) 
	  Remove_Link(s,Neighbour(s,k));
     // second, regroup the rest of the info
     V[s].Destroy();
     for (long i=s;i<N;i++)
     {
	  V[i].D=V[i+1].D;
	  V[i].N=V[i+1].N;
     }
     V=(List*)realloc(V,N*sizeof(List));
}

void Graph::Remove_Sites(const List &L)
{
     for (long k=1;k<=L.N;k++)
	  Remove_Site(L(k));
}

void Graph::Remove_Link(long s1, long s2)
{
     if (!Is_Link(s1,s2)) return;
     V[s1].Substract(s2);
     if (!directed)
	  V[s2].Substract(s1);
     long index=Link_Index(s1,s2);
     I1.Remove(index);
     I2.Remove(index);
     Nl--;
}

void Graph::Remove_Link(long l)
{
     long s1, s2;
     Link_Sites(s1,s2,l);
     Remove_Link(s1,s2);
}


void Graph::Clear()
{
     if (!N) return;
     I1.Destroy(); 
     I2.Destroy(); 
     for (long i=1;i<=N;i++)
	  V[i].Destroy();
     Nl=0;
}

// Re-build the link indices, assuming the neighborhood structure V
// and knowing whether the graph is directed or not
void Graph::Update_Index()
{
     I1.Destroy();
     I2.Destroy();
     for (long i=1;i<=N;i++)
	  for (long k=1;k<=V[i].N;k++)
	  {
	       long j=V[i](k); // so i-j is a pair
	       if (directed || (i<j) )
	       {
		    I1.Append(i);
		    I2.Append(j);
	       }
	  }
     Nl=I1.N;
}

long Graph::Degree(long p) const // number of neighbours of site p
{
     return V[p].N;
}

long Graph::Neighbour(long p, long k) const // k-th neighbour of site p
{
     return V[p](k);
}

List Graph::Neighbours(long p) const
{
     return V[p];
}
     
bool Graph::Is_Link(long s1, long s2) const
// returns true if there is a link from s1 to s2.
{
     if (!N || !Nl) return false;
     if (!V[s1].Find(s2)) return false;
     return true;
}
	  
long Graph::Link_Index(long s1, long s2) const
// Get the index of the link corresponding to sites s1->s2
// Returns 0 if sites are not linked
{
     if (s1<=0 || s1>N || s2<=0 || s2>N) return 0;
     if (!Is_Link(s1,s2)) return 0;
     
     long index=0;
     if (directed)
	  for (long l=1;l<=Nl;l++)
	       if (I1(l)==s1 && I2(l)==s2) 
	       {
		    index=l;
		    break;
	       }
     if (!directed)
	  for (long l=1;l<=Nl;l++)
	       if ( (I1(l)==s1 && I2(l)==s2) ||  (I1(l)==s2 && I2(l)==s1))
	       {
		    index=l;
		    break;
	       }
     if (!index) Error("Link not found! Weird!\n");
     return index;    
}

// returns s1 and s2 as the sites corresponding to link number "index"
void Graph::Link_Sites(long &s1, long &s2, long index) const
{
     s1=I1(index);
     s2=I2(index);
}

Matrix Graph::Adjacency_Matrix() const
{
     Matrix A(N);
     for (long i=1;i<=N;i++)
 	  for (long j=1;j<=N;j++)
	       if (Is_Link(i,j)) A(i,j)=1.0;
     return A;
}

List Graph::Connected_Component(long p) const
{
     List R;  // result
     List A(1); // list of elements to add
     A(1)=p;
     bool cont=true;
     while(cont)
     {
	  R.Append(A);
	  List A_New;
	  for (long i=1;i<=A.N;i++)
	       A_New.Append(Neighbours(A(i))); 
          // all neighbours of A go to A_New
	  A_New.Uniquify();
	  A_New.Substract(R);
	  if (!A_New.N) cont=false;
	  A=A_New;
     }
     return R;
}

// Implementation of Dijkstra's algorithm
List Graph::Find_Path(long s1, long s2) const
{
     List D(N), Prev(N); // distances from s1
     for (long i=1;i<=N;i++)
     {
	  D(i)=LONG_MAX; // infinity
	  Prev(i)=0; // undefined
     }
     D(s1)=0; // s1 to s1 is zero...
     List S; // empty list, sites whose D(i) is OK
     List Q=List_Range(1,N); // all sites, opposite of S
     while(Q.N) // while Q is not empty
     {
	  // Extract the site with the minimum D within Q
	  long imin=0, dmin=LONG_MAX;
	  for (long i=1;i<=Q.N;i++)
	       if (D(Q(i))<dmin) { imin=i; dmin=D(Q(i)); }
	  if (imin==0) printf("Help!\n");
	  long u=Q(imin);
	  S.Append(u);
	  Q.Substract(u);
	  for (long k=1;k<=Degree(u);k++)
	  {
	       long v=Neighbour(u,k);
	       if (D(v)>D(u)+1) // if links are weighted, this should be modif
	       {
		    D(v)=D(u)+1;
		    Prev(v)=u;
	       }
	  }
	  if (u==s2) break;
     }
     long p=s2;
     List P;
     while(Prev(p))
     {
	  P.Append(p);
	  p=Prev(p);
     }
     P.Append(p);
     return P;
}

long Graph::Distance(long s1, long s2) const
{
     List P=Find_Path(s1,s2);
     return P.N;
}

// Implementation of Dijkstra's algorithm
// V is the vector of distances for each link
List Graph::Find_Path(const Vector &V, long s1, long s2) const
{
     List P;
     Vector D(N); 
     List Prev(N); // distances from s1 
     for (long i=1;i<=N;i++)
     {
	  D(i)=1e20; // infinity
	  Prev(i)=0; // undefined
     }
     D(s1)=0.0; // s1 to s1 is zero...
     List S; // empty list, sites whose D(i) is OK
     List Q=List_Range(1,N); // all sites, opposite of S
     while(Q.N) // while Q is not empty
     {
	  // Extract the site with the minimum D within Q
	  long imin=0; double dmin=1e20;
	  for (long i=1;i<=Q.N;i++)
	       if (D(Q(i))<dmin) { imin=i; dmin=D(Q(i)); }
	  if (imin==0) printf("Help!\n");
	  long u=Q(imin);
	  S.Append(u);
	  Q.Substract(u);
	  for (long k=1;k<=Degree(u);k++)
	  {
	       long v=Neighbour(u,k);
	       long jx=Link_Index(u,v);
	       if (D(v)>D(u)+V(jx)) 
	       {
		    D(v)=D(u)+V(jx);
		    Prev(v)=u;
	       }
	  }
	  if (u==s2) break;
     }
     long p=s2;
     while(Prev(p))
     {
	  P.Append(p);
	  p=Prev(p);
     }
     P.Append(p);
     return P;
}

// Prim's (aka DJP) algorithm seems suitable
Graph Graph::Minimum_Spanning_Tree() const
{
     Graph G2(*this);
     G2.Clear(); // remove all links
   
     List Tree;
     List Links=List_Range(1,Nl);
     long p=0;
     do
     {
	  p++;
     }while(!V[p].N);
     Tree.Append(p);
     long s1, s2;
     long ausgang=0;
     long treen_old=0;
     do
     {
	  long i=1;
	  do
	  {
	       Link_Sites(s1,s2,Links(i));
	       // If the link connects the graph and the rest?
	       if ((Tree.Find(s1) && !Tree.Find(s2)) ||
		   (!Tree.Find(s1) && Tree.Find(s2)))
	       {
		    G2.Add_Link(s1,s2);
		    if (Tree.Find(s1)) 
			 Tree.Append(s2);
		    else
			 Tree.Append(s1);
		    Links.Substract(Links(i));
	       }
	       i++;
	  }while(i<=Links.N);
	  // Check if all remaining links are internal to the tree
	  long suma=0;
	  for (long i=1;i<=Links.N;i++)
	  {
	       Link_Sites(s1,s2,Links(i));
	       if (Tree.Find(s1) && Tree.Find(s2)) suma++;
	  }
	  if (suma==Links.N) ausgang=1;
	  else if (Tree.N==treen_old)
	  {
	       printf("Graph was not connected, "
		      "Spanning tree applies to one component.\n");
	       ausgang=1;
	  }
	  treen_old=Tree.N;
     }while(!ausgang);
     return G2;
}

Graph& Graph::operator=(const Graph &G)
{
     if (&G==this) return *this;
     Copy(*this,G);
     return *this;
}
     
long Graph::operator()(long n) const
{
     return V[n].N;
}

long Graph::operator()(long n, long m) const
{
     if (m>V[n].N) return 0;
     return V[n](m);
}

long& Graph::operator()(long n, long m)
{
     if (m>V[n].N) 
	  Error("Error accessing Graph data!\n");
     return V[n](m);
}

void Graph::Write() const 
{
     if (directed) printf("Directed ");
     printf("Graph with %ld sites and %ld links\n",N,Nl);
     for (long i=1;i<=N;i++)
     {
	 printf("%ld (%ld): ",i,V[i].N);
	 V[i].Write();
     }
     printf("\n");
     printf("Links: \n");
     long s1, s2;
     for (long l=1;l<=Nl;l++)
     {
	  Link_Sites(s1,s2,l);
	  if (!directed)
	       printf("Link %ld: (%ld-%ld)\n",l,s1,s2);
	  else printf("Link %ld: %ld->%ld\n",l,s1,s2);
     }
     printf("\n");
}
     
bool Graph::Save(const char *name) const
{
     FILE *fich=fopen(name,"wt");
     if (!fich) return false;
     if (directed) fprintf(fich,"Directed ");
     fprintf(fich,"Graph with %ld sites and %ld links\n",N,Nl);

     for (long i=1;i<=N;i++)
     {
	  fprintf(fich,"%ld (%ld): ",i,V[i].N);
	  for (long j=1;j<=V[i].N;j++)
	       fprintf(fich,"%ld ",V[i](j));
	  fprintf(fich,"\n");
     }
     fprintf(fich,"\n");
     fprintf(fich,"Links: \n");
     long s1, s2;
     for (long l=1;l<=Nl;l++)
     {
	  Link_Sites(s1,s2,l);
	  fprintf(fich,"Link %ld: (%ld-%ld)\n",l,s1,s2);
     }
     fprintf(fich,"\n");
     fclose(fich);
     return true;
}

bool Graph::Load(const char *name)
{
     FILE *fich=fopen(name,"rt");
     long n, nl;
     Text Z;
     Destroy();
     Z.Get_Line(fich);
     if (Z.Is_Here("Directed"))
     {
	  directed=true;
	  sscanf(Z.D,"Directed Graph with %ld sites and %ld links\n",
		 &n,&nl);
     }
     else 
     {
	  directed=false;
	  sscanf(Z.D,"Graph with %ld sites and %ld links\n",&n,&nl);
     }
     Add_Sites(n);
     for (long i=1;i<=n;i++)
     {
	  long i2, d;
	  if (!fscanf(fich,"%ld (%ld):",&i2,&d))
	       return false;
	  if (i2!=i) return false;
	  V[i].Create(d);
	  for (long j=1;j<=d;j++)
	       if (!fscanf(fich,"%ld",&V[i](j)))
		    return false;
     }
     Update_Index();
     fclose(fich);
     return true;
}

/////////////////////////////////////////////////////////////////////
// External functions
/////////////////////////////////////////////////////////////////////

void Copy(Graph &G, const Graph &G2)
{
     if (G.N!=0) G.Destroy();
     G.Create(G2.N);
     G.Nl=G2.Nl;
     for (long i=1;i<=G.N;i++)
	  Copy(G.V[i],G2.V[i]);
     Copy(G.I1,G2.I1);
     Copy(G.I2,G2.I2);
}

//////////////////////////////////////////////////////////////////////
// A few concrete graphs
//////////////////////////////////////////////////////////////////////

Graph Linear_Graph(long N)
{
     Graph G;
     G.Add_Site();
     for (long i=2;i<=N;i++)
	 G.Add_Site(i-1);
//     G.Update_Index();
     return G;
}

Graph Linear_Graph_PBC(long N)
{
     Graph G=Linear_Graph(N);
     G.Add_Link(1,N);
//     G.Update_Index();
     return G;
}

Graph Square_Graph(long lx, long ly)
{
     long N=lx*ly;
     Graph G(N);
    
     // graph index: (x-1)*Ly+y+1
     for (long x=1;x<=lx;x++)
         for (long y=1;y<=ly;y++)
	 {
	      if (x>1) G.Add_Link((x-1)*ly+y,(x-2)*ly+y);
	      if (y>1) G.Add_Link((x-1)*ly+y,(x-1)*ly+y-1);
	 }
//     G.Update_Index();
     return G;
}

Graph Square_Graph_PBC(long lx, long ly)
{
     long N=lx*ly;
     Graph G(N);

     // graph index: (x-1)*Ly+y+1
     for (long x=1;x<=lx;x++)
         for (long y=1;y<=ly;y++)
	 {
	      if (x>1) G.Add_Link((x-1)*ly+y,(x-2)*ly+y);
	      else G.Add_Link((x-1)*ly+y,(ly-1)*ly+y);
	      
	      if (y>1) G.Add_Link((x-1)*ly+y,(x-1)*ly+y-1);
	      else G.Add_Link((x-1)*ly+y,(x-1)*ly+ly);
	     
	 }
//     G.Update_Index();
     return G;
}

Graph Complete_Graph(long N)
{
     Graph G(N);
     for (long i=1;i<N;i++)
	  for (long j=i+1;j<=N;j++)
	       G.Add_Link(i,j);
//     G.Update_Index();
     return G;
}

// Get the Graph representing the non-zero entries of M
Graph Matrix_2_Graph(const Matrix &M)
{
     if (M.N1!=M.N2) Error("Need square Matrix for Matrix_To_Graph\n");
     long N=M.N1;
     Graph G(N);
     for (long i=1;i<=N;i++)
	  for (long j=1;j<=N;j++)
	       if (M(i,j)!=0.0) G.Add_Link(i,j);
//     G.Update_Index();
     return G;
}

// Get the Graph from the set of Points given in the columns of M, with cutoff
Graph Points_2_Graph(const Matrix &M, double dcutoff)
{
     long N=M.N2;
     Graph G(N);
     Matrix D(N);
     for (long i=1;i<=N;i++)
     {
	  Vector Vi=M.Col(i);
	  for (long j=i+1;j<=N;j++)
	  {
	       Vector Vj=M.Col(j);
	       double dist=(Vi-Vj).Norm();
	       if (dist<=dcutoff) G.Add_Link(i,j);
	  }
     }
//     G.Update_Index();
     return G;
}

Graph Remove_Site(const Graph &G,long n)
{
     Graph G2(G);
     G2.Remove_Site(n);
     return G2;
}

Graph Remove_Sites(const Graph &G, const List &L)
{
     Graph G2(G);
     G2.Remove_Sites(L);
     return G2;
}

	  












