// MANYBODY - manybody calculations, complex numbers, high spin, no bmatrices
// Javier Rodriuez-aguna
// 100615-150127
#include"Manybody.h"

// return spin 1/2 matrix, 0, +1 or -1
MatrixC Spin_Op(int s) // for a single 1/2
{
     return Spin_Op(2,s);
}

// give multiplicity and component (0,+1,-1)
MatrixC Spin_Op(long mult, int s)
{
     double spin=mult/2.0-0.5;
     MatrixC S(mult);
     if (s==0)
     {
	  for (long i=1;i<=mult;i++)
	       S(i,i)=i-spin-1.0;
	  return S;
     }
     for (long i=1;i<mult;i++)
     {
	  double ms=i-spin-1.0;
	  double factor=sqrt(spin*(spin+1.0) - ms*(ms+1.0));
	  S(i+1,i)=factor;     
     }
     if (s==-1) S.Herm();
     return S;
}

MatrixC Sz_Op(long mult)
{
     return Spin_Op(mult,0);
}

MatrixC Sx_Op(long mult)
{
     return 0.5*(Spin_Op(mult,1)+Spin_Op(mult,-1));
}

MatrixC Sy_Op(long mult)
{
     return -0.5*M_I*(Spin_Op(mult,1)-Spin_Op(mult,-1));
}

MatrixC C_Op(int s) //  particle operator; 0: number, 1: creator, -1: annih
{
     MatrixC B(2);
     switch(s)
     {
     case 0: B(2,2)=1.0; break;
     case -1: B(1,2)=1.0; break;
     case 1: B(2,1)=1.0; break;
     }
     return B;
}

// for 2-states per site
MatrixC Site_Op(const MatrixC &B, long i, long N)
{
     if (i==1) return Tens_Prod_Unit(B,Pow_2(N-1),Right);
     if (i==N) return Tens_Prod_Unit(B,Pow_2(N-1),Left);
     MatrixC Rparc, R;
     Tens_Prod_Unit(Rparc,B,Pow_2(N-i),Right);
     Tens_Prod_Unit(R,Rparc,Pow_2(i-1),Left);
     return R;
}

// for general number of states per site
// Ldim is a list with the dimensionality of each site
// B is a matrix with the correct dimensionality, of course!
MatrixC Site_Op(const MatrixC &B, long i, const List &Ldim)
{
     long Ntot=Ldim.Prod();
     long Np=Ntot/Ldim(i);
     long N=Ldim.N;
     if (i==1) return Tens_Prod_Unit(B,Np,Right);
     if (i==N) return Tens_Prod_Unit(B,Np,Left);
     long Np1=Ldim.Prod(1,i-1);
     long Np2=Ldim.Prod(i+1,N);
     MatrixC Rparc, R;
     Tens_Prod_Unit(Rparc,B,Np2,Right);
     Tens_Prod_Unit(R,Rparc,Np1,Left);
     return R;
}

// Multiply If*..*If*c*I*..*I (tensor!)
// Not very efficient, but it works
MatrixC Ferm_Site_Op(const MatrixC &c, long i, long L)
{
     static MatrixC If;
     If.Create(2); If.Zero(); If(1,1)=1.0; If(2,2)=-1.0;
     MatrixC Left(1); Left(1,1)=1.0; 
     for (long j=1;j<i;j++)
	  Left=Tens_Prod(Left,If);
     MatrixC C=Tens_Prod(Left,c);
     return Tens_Prod_Unit(C,1<<(L-i),Right);
}

// Multiply If*..*If*c*I*..*I (tensor!)
MatrixC Ferm_Site_Op_2(const MatrixC &A, long i, long L)
{
     if (i==1) return Tens_Prod_Unit(A,1<<(L-1),Right); // no difference!
     // First, build the diagonal of (-1)^{Nf} on the left side
     long Nleft=Pow_2(i-1);
     VectorC D(Nleft);
     for (long j=0;j<Nleft;j++)
	  D(j+1)=Mop(Count_Ones(j));
     MatrixC Rtemp;
     if (i==L) Copy(Rtemp,A);
     else Tens_Prod_Unit(Rtemp,A,Pow_2(L-i),Right);
     MatrixC R;
     Tens_Prod_Diag(R,Rtemp,D,Left);
     return R;
}

// Spin 1/2 system, trace Rho on the sites in list L
MatrixC Trace_On(const MatrixC &Rho, const List &Lsites)
{
     long Lt=Log_2(Rho.N1); 
     List Ldim(Lt); Ldim.Set(2);
     return Trace_On(Rho,Lsites,Ldim);
}

// General spin system, trace Rho on the sites in list L1, they have
// dimensions in Ldim
MatrixC Trace_On(const MatrixC &Rho, const List &Ls1, const List &Ldim)
{
     long l=Ldim.N; // l1=Ls1.N; // l2=l-l1;
     // Find the dimension of the resulting Hilbert space
     List Ldim1=Combine(Ls1,Ldim);
     List Ls=List_Range(1,l);
     List Ls2=Substract(Ls,Ls1);
     List Ldim2=Combine(Ls2,Ldim);
     long N1=Ldim1.Prod();
     long N=Ldim.Prod(); long N2=N/N1;

     MatrixC R(N1);
     for (long i1=1;i1<=N1;i1++)
	  for (long j1=1;j1<=N1;j1++)
	       for (long k2=1;k2<=N2;k2++)
	       {
		    List I1=Int_2_List(i1-1,Ldim1);
		    List J1=Int_2_List(j1-1,Ldim1);
		    List K2=Int_2_List(k2-1,Ldim2);
		    List I(l), J(l);
		    for (long p=1;p<=l;p++)
			 if (Ls1.Find(p)) 
			 {
			      I(p)=I1(Ls1.Find(p));
			      J(p)=J1(Ls1.Find(p));
			 }
			 else 
			 {
			      I(p)=K2(Ls2.Find(p));
			      J(p)=K2(Ls2.Find(p));
			 }
		    long i=List_2_Int(I,Ldim)+1;
		    long j=List_2_Int(J,Ldim)+1;
		    R(i1,j1)+=Rho(i,j);
	       }
     return R;
}

double Shannon(const Vector &V)
{
     double S=0.0;
     for (long i=1;i<=V.N;i++)
	  S+=(V(i)<=0.0 ? 0.0 : -V(i)*log(V(i)));
     return S; 
}

double Renyi(const Vector &V, double alpha)
{
     double S=0.0;
     for (long i=1;i<=V.N;i++)
	  S+=pow(V(i),alpha);
     return (1.0/(1.0-alpha)) * log(S);
}

double Von_Neumann(const MatrixC &Rho)
{
     Vector Prob; Rho.Spectrum(Prob);
     return Shannon(Prob);
}

double Renyi(const MatrixC &Rho, double alpha)
{
     Vector Prob; Rho.Spectrum(Prob);
     return Renyi(Prob,alpha);
}

// returns operator Ntot
MatrixC Total_Number_Op(long L)
{
     long Ntot=Pow_2(L);
     MatrixC Nop(Ntot);
     for (long i=0;i<Ntot;i++)
	  Nop(i+1,i+1)=Count_Ones(i);
     return Nop;
}

// ITF hamiltonian with graph G, couplings J and field Gamma
// s is the spin multiplicity
MatrixC ITF_Ham(const Graph &G, const Vector &J, double Gamma, long s)
{
     long N=G.N; long Nl=G.Nl;
     MatrixC Ham;
     MatrixC Sz=Spin_Op(s,0);
     MatrixC Sx=Sx_Op(s);

     List Ldim(N); Ldim.Set(s);
     for (long k=1;k<=Nl;k++)
     {
	  long s1, s2;
	  G.Link_Sites(s1,s2,k);
	  MatrixC Sz1=Site_Op(Sz,s1,Ldim);
	  MatrixC Sz2=Site_Op(Sz,s2,Ldim);
	  MatrixC Szz=Sz1*Sz2;
	  Ham-=J(k)*Szz;
     }
     for (long i=1;i<=N;i++)
     {
	  MatrixC SX=Site_Op(Sx,i,Ldim);
	  Ham-=Gamma*SX;
     }
     return Ham;
}

// Heisenberg hamiltonian with graph G and couplings J
// s is the spin multiplicity
MatrixC Heisenberg_Ham(const Graph &G, const Vector &J, long s)
{
     long N=G.N; long Nl=G.Nl;
     MatrixC Ham;
     MatrixC Sz=Spin_Op(s,0);
     MatrixC Sp=Spin_Op(s,+1);
     MatrixC Sm=Spin_Op(s,-1);
     List Ldim(N); Ldim.Set(s);
     for (long k=1;k<=Nl;k++)
     {
	  long s1, s2;
	  G.Link_Sites(s1,s2,k);
	  MatrixC Sz1=Site_Op(Sz,s1,Ldim);
	  MatrixC Sz2=Site_Op(Sz,s2,Ldim);
	  MatrixC Szz=Sz1*Sz2;
	  Sz1.Destroy(); Sz2.Destroy();
	  MatrixC Sp1=Site_Op(Sp,s1,Ldim);
	  MatrixC Sm2=Site_Op(Sm,s2,Ldim);
	  MatrixC Spm=Sp1*Sm2;
	  Sp1.Destroy(); Sm2.Destroy();
	  Spm+=Herm(Spm);
	  Ham+=J(k)*( Szz + 0.5*( Spm ));
     }
     return Ham;
}
