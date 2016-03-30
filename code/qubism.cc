// Functions needed for qubistic plots
// Javier Rodríguez-Laguna and Piotr Migdał
// 111107-160230
#include"qubism.h"

// Standard Qubistic Plot
// {00-01;10-11}
MatrixC Standard_2D_Plot(const VectorC &Psi, long L)
{
     List Ls(L/2); Ls.Set(4); // A list with 
     long Lt=Pow_2(L/2);
     MatrixC Im(Lt);
     for (long i=1;i<=Psi.N;i++)
     {
	  List D=Int_2_List(i-1,Ls);; // number in base-4
	  long x=0, y=0, l=Lt/2;
	  for (long k=1;k<=D.N;k++)
	  {
	       if (D(k)==1) x+=l;
	       if (D(k)==2) y+=l;
	       if (D(k)==3) { x+=l; y+=l;}
	       l/=2;
	  }
	  Im(x+1,y+1)=Psi(i);
     }
     return Im;
} 

// Alternate Qubistic Plot
// {00-01;11-10}
MatrixC Alt_2D_Plot(const VectorC &Psi, long L)
{
     List Ls(L/2); Ls.Set(4);
     long Lt=1<<(L/2);
     MatrixC Im(Lt);
     for (long i=1;i<=Psi.N;i++)
     {
	  List D=Int_2_List(i-1,Ls);; // number in base-4
	  long x=0, y=0, l=Lt/2;
	  for (long k=1;k<=D.N;k++)
	  {
	       if (D(k)==1) x+=l;
	       if (D(k)==2) { x+=l; y+=l;}
	       if (D(k)==3) y+=l;
	       l/=2;
	  }
	  Im(x+1,y+1)=Psi(i);
     }
     return Im;
} 

// Spin-1 Standard Qubistic Plot
MatrixC Spin1_2D_Plot(const VectorC &P, long L)
{
     long ltot=(long)pow(3,L/2);
     long Ltot=ltot*ltot;
     MatrixC M(ltot);
     List Ls(L); Ls.Set(3);
     for (long i=1;i<=Ltot;i++)
     {
	  List Li=Int_2_List(i-1,Ls);
	  long x=0, y=0, l=ltot/3;
	  for (long k=1;k<=L;k+=2)
	  {
	       y+=Li(k)*l;
	       x+=Li(k+1)*l;
	       l/=3;
	  }
	  M(x+1,y+1)=P(i);
     }
     return M;
}

// Auxiliary plotting routines for EX
// Plot a pxp grid on the current window
void Plot_Grid(long p)
{
     static long xsize=EX_CW->width;
     double ps=(double)xsize/(double)p;
     for (long i=0;i<=p;i++)
     {
	  long ips=(long)round(i*ps);
	  if (ips<0) ips=0;
	  if (ips>=xsize) ips=xsize-1;
	  EX_Line(ips,0,ips,xsize-1);
	  EX_Line(0,ips,xsize-1,ips);
     }
}

// Plot auxiliary lines for the qubistic plot
void Plot_Lines(long s)  // s:2, qubits; s:3, qutrits
{
     long color1=EX_Alloc_RGB_Color(0.4,0.4,0.4); 
     long color2=EX_Alloc_RGB_Color(0.7,0.7,0.7);

     EX_Set_Color(color2);
     Plot_Grid(Sqr(s));

     EX_Set_Color(color1);
     Plot_Grid(s);
}
