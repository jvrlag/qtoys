// Toy code, classical simulation of an adiabatic quantum computation
// 160321, Javier Rodríguez-Laguna, UNED.
// http://mononoke.fisfun.uned.es/qcomp

#include"Manybody.h" // Part of the HVB distribution

int main()
{
     long N=8; // number of qubits, increase carefully!
     Rand_Open(time(0)); // start the random number generator
     
     Matrix J(N); // the matrix with the spin-glass couplings
     for (long i=1;i<=N;i++)
	  for (long j=i+1;j<=N;j++)
	       J(i,j)=J(j,i)=Rand(-1,1);

     // Now, build the target Hamiltonian, H_I
     // Site_Op(A,i,N) builds the 2^N x 2^N matrix in the lexicographic
     // basis, in which operator A acts on site "i" of N.
     MatrixC Hising;
     MatrixC Sz(2);
     Sz(1,1)=-1.0;
     Sz(2,2)=1.0;
     for (long i=1;i<=N;i++)
	  for (long j=i+1;j<=N;j++)
	       Hising+=J(i,j)*Site_Op(Sz,i,N)*Site_Op(Sz,j,N);
     // In order to break the symmetry between the two 
     // possible solutions, we force the first site to point "down"
     Hising+=Site_Op(Sz,1,N);
     
     // Now, build the "easy" Hamiltonian, H_X
     MatrixC HX;
     MatrixC Sx(2);
     Sx(1,2)=Sx(2,1)=1.0;
     for (long i=1;i<=N;i++)
	  HX+=-Site_Op(Sx,i,N);
     
     // Now, the adiabatic path, with parameter "s"
     double ds=0.1;
     printf("# AQC simulation spin-glass system, N=%ld qubits\n",N);
     printf("# s   Gap\n");
     double time=0.0;
     MatrixC H, Basis; Vector Eigen;
     for (double s=0.0;s<=1.0;s+=ds)
     {
	  // Build the full Hamiltonian
	  H=s*Hising + (1-s)*HX;
	  // Diagonalize it: find the eigenstates and eigenvalues
	  H.Diagonalize(Basis,Eigen); 
	  // Find the energy gap:
	  double gap=Eigen(2)-Eigen(1);
	  printf("%10g %16.12g \n",s,gap);
	  time+=ds/Sqr(gap);
     }
     printf("# AQC time estimate: %g\n", time);

     // Print the solution!
     printf("# Solution: ");
     VectorC GS=Basis.Col(1); // the ground state is the first column
     for (long i=1;i<=N;i++)
     {
	  // Obtain the matrix representation of Sz-i
	  MatrixC Szi=Site_Op(Sz,i,N);
	  // Obtain the expected value of Sz-i in the GS
	  double szi=real(Szi.Elem(GS,GS));
	  printf("%g ",szi);
     }
     printf("\n");
	  
	  
}


