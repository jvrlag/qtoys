// MANYBODY - manybody calculations 
// Javier Rodríguez-Laguna
// 100614-110509-130108-150127
#ifndef MANYBODY_H
#define MANYBODY_H

#include"MatrixC.h"
#include"Graph.h"

MatrixC Spin_Op(int ms); // for a single 1/2; 0: z, 1: p, -1: m.
MatrixC Spin_Op(long mult, int ms); // multiplicity and ms (-1,0,1)
MatrixC Sz_Op(long mult); // given the multiplicity, i.e.: 2 for spin 1/2
MatrixC Sx_Op(long mult);
MatrixC Sy_Op(long mult);
MatrixC C_Op(int s); // hard-core operator; 0: number, 1: creator, -1: annih

MatrixC Site_Op(const MatrixC &B, long i, long N);
MatrixC Site_Op(const MatrixC &B, long i, const List &Ldim);
MatrixC Ferm_Site_Op(const MatrixC &C, long i, long L); // fermionic
MatrixC Ferm_Site_Op_2(const MatrixC &C, long i, long L); // fermionic

MatrixC Trace_On(const MatrixC &Rho, const List &Lsites);
MatrixC Trace_On(const MatrixC &Rho, const List &Ls1, const List &Ldim);

double Shannon(const Vector &V);
double Renyi(const Vector &V, double alpha);
double Von_Neumann(const MatrixC &Rho);
double Renyi(const MatrixC &Rho, double alpha);

MatrixC Total_Number_Op(long L);
MatrixC Heisenberg_Hamiltonian(const Graph &G, const Vector &J, long s);
MatrixC ITF_Hamiltonian(const Graph &G, const Vector &J, double Gamma, long s);

#endif
