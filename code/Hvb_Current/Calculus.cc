// Calculus functions implementations
// 150822
// Integration / Interpolation / ODE
#include"Calculus.h"

// Sum function F, using trapezoid rule, with n rectangles
double Trapezoid_Sum(Func_RR &F, double a, double b, long n)
{
     if (n == 1) 
	  return (0.5*(b-a)*(F(a) + F(b)));
     else 
     {
	  long Nr=Pow_2(n-2);
	  double dx=(b-a)/(double)Nr;
	  double x=a+0.5*dx;
	  double sum=0.0;
	  for (long j=1;j<=Nr;j++,x+=dx) 
	       sum += F(x);
	  return (b-a)*sum/(double)Nr;
     }
}

double Integrate(Func_RR &F, double a, double b, double tol)
{
     static long jmax=40;
     if (fabs(a-b)<1e-20) return 0.0;

     double sum=Trapezoid_Sum(F,a,b,1);
     double sum_old=sum;
     double S_old=sum;
     for (long j=2;j<=jmax;j++) 
     {
	  sum=0.5*(sum+Trapezoid_Sum(F,a,b,j));
	  double S=(4.0*sum-sum_old)/3.0;
	  if (fabs(S-S_old) < tol*fabs(S_old)) 
	       return S;
	  S_old=S;
	  sum_old=sum;
     }
     printf("Integrator has failed, too low dx is needed!\n");
     return 0.0;
}

double Integrate(double(*f)(double,void*), void* P,
		 double a, double b, double tol)
{
     Func_RR F(f,P);
     return Integrate(F,a,b,tol);
}

double Integrate_To_Inf(Func_RR &F, double a, double tol1, double tol2)
{
     static long jmax=40;
     static double fx=1.5;
     double l=1.0, x0=a;
     double S_old=0.0, S=0.0;
     for (int j=1;j<=jmax;j++)
     {
	  S+=Integrate(F,x0,x0+l,tol1);
	  if (fabs(S-S_old) < tol2*fabs(S)) return S;
	  S_old=S;
	  x0=x0+l;
	  l*=fx;
     }
     printf("Integrator has failed, too long L is needed!\n");
     return 0.0;
}

double Integrate_To_Inf(double(*f)(double,void*),void* P,
			double a, double tol1, double tol2)
{
     Func_RR F(f,P);
     return Integrate_To_Inf(F,a,tol1,tol2);
}


ODE_Solver::ODE_Solver()
{
     P=NULL; t=0.0; fn=NULL;
}

ODE_Solver::ODE_Solver(Func_ODE f, const Vector &x, void *p, double tt, double dtt)
{
     Initialize(f,x,p,tt,dtt);
}

ODE_Solver::~ODE_Solver()
{ }

void ODE_Solver::Initialize(Func_ODE f, const Vector &x, void *p,
			    double tt, double dtt)
{
     fn=f; X=x; P=p; t=tt; dt=dtt; Ndim=X.N; Nstep=0;
     F1.Create(Ndim);
     F2.Create(Ndim);
     F3.Create(Ndim);
     F4.Create(Ndim);
     Acum.Create(Ndim);
}

void ODE_Solver::Runge_Kutta()
{
     fn(F1,X,P,t);
     Acum=X+(dt/2.0)*F1;
     
     fn(F2,Acum,P,t+dt/2.0);
     Acum=X+(dt/2.0)*F2;

     fn(F3,Acum,P,t+dt/2.0);
     Acum=X+dt*F3;

     fn(F4,Acum,P,t+dt);
     X+=(dt/6.0)*(F1+2.0*F2+2.0*F3+F4);

     t+=dt;
     Nstep++;
}

// Function given by table M=X|Y, approximate F(x)
double Interpolate(double x, const Matrix &M)
{
     // xi=M(i,1); yi=M(i,2)
     Vector X;
     X.N=M.N1;
     X.D=M.D;
     long i=Search(x,X);
     // i is the highest index s.t. x>X(i)
     X.D=NULL;
     if (i==0) return M(1,2);
     if (i==M.N1) return M(M.N1,2);
     double dx=M(i+1,1)-M(i,1);
     double wr=(x-M(i,1))/dx;
     return wr*M(i+1,2)+(1-wr)*M(i,2);
     
}

// Interpolation function in Func_RR form, so it can be used in Optimize, etc.
double Interpolate(double x, void *P)
{
     Matrix *M=(Matrix*)P;
     return Interpolate(x,*M);
}

// Return highest index i such that x>X(i), or 0 if none
// Assume X is sorted in increasing order
long Search(double x, const Vector &X)
{
     if (x<X(1)) return 0;
     long N=X.N;
     if (x>X(N)) return N;
     long imin=1, imax=N;
     long itry;
     do
     {
	  if ((imax-imin)<=1) return imin;
	  itry=(imin+imax)/2;
	  if (x>X(itry)) imin=itry;
	  else imax=itry;
     }while(imax-imin>1);
     return imin;
}
