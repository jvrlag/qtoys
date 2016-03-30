// Optimization routines
// JaviRL, 100317-150126
#ifndef OPTIMIZE_H
#define OPTIMIZE_H

#include"Matrix.h"

#ifndef CALCULUS_H
#define CALCULUS_H
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

// Function R^n->R
class Func_VR
{
public:
     double(*Func)(const Vector &X, void *);
     void *p;
     Func_VR() { Func=NULL; p=NULL;}
     Func_VR(double(*f)(const Vector &, void*), void *param): 
	  Func(f), p(param) {}
     ~Func_VR() {}
     double operator()(const Vector & x) const;
};

// Returns the value of a function and its gradient
// Caution: G must be allocated!
double Func_Grad(Vector &G, const Vector &X, Func_VR &F);

double Line_Optimize(Vector &X, const Vector &G, Func_VR &F, double tol);

void Bracket_Minimum(double &a, double &b, double &c, Func_RR &F);

double Brent_Optimize(double &x, double a, double b, double c, 
		      Func_RR &F, double tol);


double SD_Optimize(Vector &X, Func_VR &F, double tol);
double SD_Optimize(Vector &X,double(*f)(const Vector&,void*), void*, 
		   double tol);

double CG_Optimize(Vector &X, Func_VR &F, double tol);
double CG_Optimize(Vector &X,double(*f)(const Vector&,void*), void*,
		   double tol);

double Powell_Optimize(Vector &P, Func_VR &F, double tol);
double Powell_Optimize(Vector &P,double(*f)(const Vector&,void*), void*,
		       double tol);

typedef struct 
{
     double Beta_Min;
     double Beta_Max;
     double R_Beta;
     long   N_Times;
     double A0;
} Annealing_Params;

double Annealing_Optimize(Vector &X, Func_VR &F, Annealing_Params *A);
double Annealing_Optimize(Vector &X,double(*f)(const Vector&,void*), void*,
			  Annealing_Params *A);

#endif
