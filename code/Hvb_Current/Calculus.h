// Calculus related routines
// Integration / ODE / Interpolation
// 150822
#ifndef CALCULUS_H
#define CALCULUS_H
#include"Matrix.h"

#ifndef OPTIMIZE_H
// Function R->R
class Func_RR
{
public:
     double(*Func)(double x, void *);
     void *p;
     Func_RR() { Func=NULL; p=NULL;}
     Func_RR(double(*f)(double, void*), void *param): Func(f), p(param) {}
     ~Func_RR() {}
     double operator()(double x) const;
};
#endif 

double Trapezoid_Sum(Func_RR &F, double a, double b, long n);

double Integrate(Func_RR &F, double a, double b, double tol);
double Integrate(double(*f)(double,void*),void*,double a, double b, double tol);

double Integrate_To_Inf(Func_RR &F, double a, double tol1, double tol2);
double Integrate_To_Inf(double(*f)(double,void*),void*,
			double a, double tol1, double tol2);

typedef void(*Func_ODE)(Vector&,const Vector&,void*,double);

class ODE_Solver
{
public:
     Func_ODE fn;
     void *P;
     double t;
     double dt;
     long Ndim;
     long Nstep;
     Vector X;
     Vector F1,F2,F3,F4,Acum;
     ODE_Solver();
     ODE_Solver(Func_ODE f, const Vector&,void*,double,double);
     ~ODE_Solver();
     void Initialize(Func_ODE f, const Vector &x, void *p, double tt, double dtt);
     void Runge_Kutta();
};
     
// Using this interface we can use 
double Interpolate(double x, void *);
double Interpolate(double x, const Matrix &M);

long Search(double x, const Vector &X);

#endif
