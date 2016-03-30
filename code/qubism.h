// Qubism header file
// Javier Rodríguez-Laguna and Piotr Migdał
// 111107-160330
#ifndef QUBISM_H
#define QUBISM_H

#include"MatrixC.h"
#include"EX.h"

// Standard Qubistic Plot
// {00-01;10-11}
MatrixC Standard_2D_Plot(const VectorC &Psi, long L);

// Alternate Qubistic Plot
// {00-01;11-10}
MatrixC Alt_2D_Plot(const VectorC &Psi, long L);

// Spin-1 Standard Qubistic Plot
MatrixC Spin1_2D_Plot(const VectorC &P, long L);

// Auxiliary plotting routines for EX
// Plot a pxp grid on the current window
void Plot_Grid(long p);

// Plot auxiliary lines for the qubistic plot
// s:2, qubits; s:3, qutrits
void Plot_Lines(long s=2);  

#endif
