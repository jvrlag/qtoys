// Optimization routines
// 100317-150126
#include"Optimize.h"

double Func_RR::operator()(double x) const
{
     return Func(x,p);
}

double Func_VR::operator()(const Vector &X) const
{
     return Func(X,p);
}

double Func_Grad(Vector &G, const Vector &X, Func_VR &F)
{
     double f=F(X);
     double dx=1e-6; // CAUTION!! This should be an open parameter
     Vector X2=X;
     for (long i=1;i<=X.N;i++)
     {
	  double x=X2(i);
          X2(i)=x+dx;
          double f2=F(X2);
	  G(i)=(f2-f)/dx;
          X2(i)=x;
     }
     return f;
}

////////////////////////////////////////////////////////////////////////
// Another "metafunction"
// 1D-mensionalization of a vector function along a line
// Needed to convert a vector function to 1d function for line search
typedef struct 
{
     Func_VR Function;
     Vector X0;
     Vector G;
}Func_1D_Params_Type;

double Func_1D_Line_Optimize(double x, void *Q)
{
     Func_1D_Params_Type Qex=*((Func_1D_Params_Type*)Q);
     return Qex.Function(Qex.X0+x*Qex.G);
}
//////////////////////////////////////////////////////////////////////////

double Line_Optimize(Vector &X, const Vector &G, Func_VR &F, double tol)
{
     // First, prepare parameters for Func1D
     Func_1D_Params_Type Qex;
     Qex.Function=F;
     Qex.X0=X;
     Qex.G=G;
     Func_RR F1D(Func_1D_Line_Optimize,&Qex);

     double a=0.0,c=min(0.01,1.0/G.Norm()),b=c/2.0;
     Bracket_Minimum(a,b,c,F1D);
     double x;
     double f=Brent_Optimize(x,a,b,c,F1D,tol);

     X=X+x*G;
     return f;
}

#define shift(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
void Bracket_Minimum(double &a, double &b, double &c, Func_RR &F)
{
     static double golden=(1.0+sqrt(5.0))/2.0; 
     static double glimit=100.0;
     static double epsilon=1e-20;
     
     double ulim,u,r,q,fu;

     double fa=F(a);
     double fb=F(b);
     if (fb > fa) 
     {
	  Swap(a,b);
	  Swap(fa,fb);
     }
     c=b+golden*(b-a);
     double fc=F(c);
     while (fb > fc) 
     {
	  r=(b-a)*(fb-fc);
	  q=(b-c)*(fb-fa);
	  u=b-((b-c)*q-(b-a)*r)/
	       (2.0*Sign(Max(fabs(q-r),epsilon),q-r));
	  ulim=b+glimit*(c-b);
	  if ((b-u)*(u-c) > 0.0) 
	  {
	       fu=F(u);
	       if (fu < fc) 
	       {
		    a=b;
		    b=u;
		    fa=fb;
		    fb=fu;
		    return;
	       } 
	       else if (fu > fb) 
	       {
		    c=u;
		    fc=fu;
		    return;
	       }
	       u=c+golden*(c-b);
	       fu=F(u);
	  } 
	  else if ((c-u)*(u-ulim) > 0.0) 
	  {
	       fu=F(u);
	       if (fu < fc) 
	       {
		    shift(b,c,u,c+golden*(c-b));
		    shift(fb,fc,fu,F(u));
	       }
	  } else if ((u-ulim)*(ulim-c) >= 0.0) 
	  {
	       u=ulim;
	       fu=F(u);
	  } else 
	  {
	       u=c+golden*(c-b);
	       fu=F(u);
	  }
	  shift(a,b,c,u);
	  shift(fa,fb,fc,fu);
     }

}

double Brent_Optimize(double &x, double ax, double bx, double cx, 
		      Func_RR &F, double tol)
{
     static double cgold=(3.0-sqrt(5.0))/2.0;
     static long itmax= 200;
     static double epsilon=1e-10;     

     double a,b,etemp,fu,f3,f2,fx,tol1,tol2,u,x3,x2,xm;
     double d=0.0,e=0.0;
     
     a=(ax < cx ? ax : cx);
     b=(ax > cx ? ax : cx);
     x=x2=x3=bx;
     f2=f3=fx=F(x);
     long iter=0;
     do
     {
	  xm=0.5*(a+b);
	  tol2=2.0*(tol1=tol*fabs(x)+epsilon);
	  if (fabs(x-xm) <= (tol2-0.5*(b-a))) 			
	       return fx;
	  
	  if (fabs(e) > tol1) 
	  {
	       double r=(x-x2)*(fx-f3);
	       double q=(x-x3)*(fx-f2);
	       double p=(x-x3)*q-(x-x2)*r;
	       q=2.0*(q-r);
	       if (q > 0.0) p = -p;
	       q=fabs(q);
	       etemp=e;
	       e=d;
	       if (fabs(p) >= fabs(0.5*q*etemp) || 
		   p <= q*(a-x) || p >= q*(b-x))
		    d=cgold*(e=(x >= xm ? a-x : b-x));
	       else 
	       {
		    d=p/q;
		    u=x+d;
		    if (u-a < tol2 || b-u < tol2)
			 d=Sign(tol1,xm-x);
	       }
	  } 
	  else 
	  {
	       d=cgold*(e=(x >= xm ? a-x : b-x));
	  }
	  u=(fabs(d) >= tol1 ? x+d : x+Sign(tol1,d));
	  fu=F(u);
	  if (fu <= fx) 
	  {
	       if (u >= x) a=x; else b=x;
	       shift(x3,x2,x,u)
	       shift(f3,f2,fx,fu)		       
	  } 
	  else 
	  {
	       if (u < x) a=u; else b=u;
	       if (fu <= f2 || x2 == x) 
	       {
		    x3=x2;
		    x2=u;
		    f3=f2;
		    f2=fu;
	       } 
	       else if (fu <= f3 || x3 == x || x3 == x2) 
	       {
		    x3=u;
		    f3=fu;
	       }
	  }
	  iter++;
     }while(iter<itmax);
     Error_Flag(Error_Mat); // THIS ERROR!?!?!?!?!
     return fx;
}	
#undef shift

double SD_Optimize(Vector &X,double(*f)(const Vector&,void*),void *p, 
		   double tol)
{
     Func_VR F(f,p);
     return SD_Optimize(X,F,tol);
}

double SD_Optimize(Vector &X, Func_VR &F, double tol)
{
     static long itmax=300;
     static double epsilon=1e-12;
     long N=X.N;
     
     Vector G(N);
     
     double fp=Func_Grad(G,X,F);
     for (long its=1;its<=itmax;its++) 
     {
	  double fret=Line_Optimize(X,G,F,tol); 
	  if (2.0*fabs(fret-fp) <= tol*(fabs(fret)+fabs(fp)+epsilon)) 
	       return fret;
	  fp=Func_Grad(G,X,F);
     }
     Error("Too many iterations in SD_Optimize\n");
     return fp;
}

double CG_Optimize(Vector &X,double(*f)(const Vector&,void*),void *p, 
		   double tol)
{
     Func_VR F(f,p);
     return CG_Optimize(X,F,tol);
}

double CG_Optimize(Vector &X, Func_VR &F, double tol)
{
     static long itmax=300;
     static double epsilon=1e-10;
     long N=X.N;
     
     double gg,fp,dgg,fret;
     
     Vector G(N);
     Vector H(N);
     Vector Xi(N);
     
     fp=Func_Grad(Xi,X,F);
     G=(-1.0)*Xi;
     Xi=H=G;
     for (long its=1;its<=itmax;its++) 
     {
	  fret=Line_Optimize(X,Xi,F,tol); 
	  if (2.0*fabs(fret-fp) <= tol*(fabs(fret)+fabs(fp)+epsilon)) 
	       return fret;
	  fp=Func_Grad(Xi,X,F);
	  gg=Dot(G,G);
	  dgg=Dot(Xi+G,Xi);
	  if (gg==0.0) return fp;
	  double gam=dgg/gg;
	  G=-Xi;
	  Xi=H=G+gam*H;
     }
     Error("Too many iterations in CG_Optimize\n");
     return fp;
}

double Powell_Optimize(Vector &X,double(*f)(const Vector&,void*),void *p, 
		   double tol)
{
     Func_VR F(f,p);
     return Powell_Optimize(X,F,tol);
}

double Powell_Optimize(Vector &P, Func_VR &F, double tol)
{
     static long itmax=300;
     long ibig; double del, fp, fptt, t;
     long N=P.N;

     Vector Pt(P), Ptt(N), Xit(N);
     Matrix Xi=Unit(N);
     double fret=F(P);
     for (long iter=1;;iter++)
     {
	  fp=fret;
	  ibig=0;
	  del=0.0;
	  for (long i=1;i<=N;i++) // for all directions in the set
	  {
	       Xit=Xi.Col(i);
	       double fptt=fret;
	       fret=Line_Optimize(P,Xit,F,tol);
	       if (fabs(fptt-fret)>del)
	       {
		    del=fabs(fptt-fret);
		    ibig=i;
	       }
	  }
	  if (2.0*fabs(fp-fret) <= tol*(fabs(fp)+fabs(fret)))
	  {
	       return fret;
	  }
	  if (iter==itmax) 
	  {
	       // Error("Powell_Optimize, too many iter\n");
	       return fret;
	  }
	  Ptt=2.0*P-Pt;
	  Xit=P-Pt;
	  Pt=P;
	  fptt=F(Ptt);
	  if (fptt < fp)
	  {
	       t=2.0*(fp-2.0*fret+fptt)*Sqr(fp-fret-del)-del*Sqr(fp-fptt);
	       if (t<0.0)
	       {
		    fret=Line_Optimize(P,Xit,F,tol);
		    Xi.Set_Col(Xi.Col(N),ibig);
		    Xi.Set_Col(Xit,N);
	       }
	  }
     }
}


double Annealing_Optimize(Vector &X,double(*f)(const Vector&,void*),void *p, 
     Annealing_Params *Q)
{
     Func_VR F(f,p);
     return Annealing_Optimize(X,F,Q);
}

double Annealing_Optimize(Vector &X, Func_VR &F, Annealing_Params *Q)
{
     double Beta_Min, Beta_Max, R_Beta, A0; long N_Times;
     if (!Q)
     {
	  Beta_Min=1.0;
	  Beta_Max=1000.0;
	  R_Beta=1.02;
	  N_Times=1000;
	  A0=0.001;
     }
     else
     {	  Beta_Min=Q->Beta_Min;
	  Beta_Max=Q->Beta_Max;
	  R_Beta=Q->R_Beta;
	  N_Times=Q->N_Times;
	  A0=Q->A0;
     }

     long N=X.N;
     
     double f=F(X);
     double fbest=f;
     Vector Xbest=X;

     for (double beta=Beta_Min;beta<=Beta_Max;beta*=R_Beta)
     {
	  long naccepted=0;
	  for (long nt=1;nt<=N_Times;nt++)
	  {
	       Vector X2(X);
	       long idx=Rand_I(1,N);
	       X2(idx)+=A0*Rand(-1.0,1.0);
	       double f2=F(X2);
	       if (f2<fbest)
	       {
		    fbest=f2;
		    Xbest=X2;
	       }
	       if (f2<f || Rand()<exp(-beta*(f2-f)))
	       {
		    X=X2;
		    f=f2;
		    naccepted++;
	       }
	       if (nt == N_Times/10) // make it a parameter?
	       {
		    double acceptance_ratio=(double)naccepted/(double)(nt);
		    if (acceptance_ratio<0.4) A0/=1.2;
		    if (acceptance_ratio>0.6) A0*=1.2;
	       }
	  }
     }
     X=Xbest;
     return fbest;
}
