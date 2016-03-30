////////////////////////////////////////////////////////
// Complex, complex vector and matrix classes
// 080912-101022-150123
// JaviRL
// Part of hvb
#ifndef CMATRIX_HEADER
#define CMATRIX_HEADER
#include "Matrix.h"

///////////////////////////////////////////////////////////////////////////
// Class CVector
///////////////////////////////////////////////////////////////////////////

class VectorC
{
public:
     long N; // length of the vector 
     cmplx *D; // the pointer itself

     VectorC();     
     VectorC(long);
     VectorC(const VectorC &);
     VectorC(const Vector &);
     VectorC(cmplx*, long);
     ~VectorC();
     void Start();
     void Create(long n);
     void Load(cmplx *, long);
     void Load_Copy(cmplx *, long);
     void Transfer(VectorC &);
     void Destroy();

     void Zero();
     void Set(cmplx);
     void Set_Part(const VectorC &, long);
     bool Normalize();
     void Part(long, long);
     void Reverse();
     void Append(const cmplx); 
     void Append(const VectorC &);
     void Insert(cmplx, long i);
     void Insert(const VectorC &, long i);
     void Remove(long i, long j=0);
     void Re();
     void Im();
     void Abs();
     void Conj();
     void Sqr(); // Caution! It's not square, it's squared modulus!

     void Write(int prec=0) const;
     void Write_Col() const;
     bool Save_Binary(const char *name) const;
     bool Load_Binary(const char *name);
     bool Save_Binary(FILE *fich) const;
     bool Load_Binary(FILE *fich);
     bool Save(const char *name) const;
     bool Load(const char *name);
     bool Save(FILE *fich) const;
     bool Load(FILE *fich);

     bool Is_Zero(double tolerance) const;
     double Norm() const;

     cmplx operator() (long) const;
     cmplx& operator() (long);
     VectorC& operator=(const VectorC&);
     VectorC& operator=(const Vector&);
     void operator+=(const VectorC&);
     void operator-=(const VectorC&);
     void operator*=(const cmplx);
     void operator/=(const cmplx);
     void operator*=(const double);
     void operator/=(const double);

     void operator&=(const VectorC &);
     void operator&=(cmplx);
}; 

void Copy(VectorC& B, const VectorC& A);
void Copy(VectorC& B, const Vector& A);
VectorC operator+(const VectorC &V, const VectorC &W);
VectorC operator+(const VectorC &V, const Vector &W);
VectorC operator-(const VectorC &, const VectorC &);
VectorC operator-(const VectorC &V);
VectorC operator*(cmplx, const VectorC &);
VectorC operator*(const VectorC &, cmplx);
VectorC operator*(double, const VectorC &);
VectorC operator*(const VectorC &, double);
VectorC operator*(cmplx, const Vector &);
VectorC operator/(const VectorC &, cmplx);

VectorC operator&(const VectorC &L1, const VectorC &L2);
VectorC operator&(cmplx, const VectorC &L);
VectorC operator&(const VectorC &L, cmplx);

cmplx   Dot(const VectorC &, const VectorC &);
VectorC Tens_Prod(const VectorC &, const VectorC &);
void    Tens_Prod(VectorC &, const VectorC &, const VectorC &);
void    Zaxpy(VectorC &R, const VectorC &V, cmplx alpha);

VectorC Normalize(const VectorC &V);
VectorC Part(const VectorC &V,long,long);
VectorC Part(const VectorC &V,const List &L);
VectorC Reverse(const VectorC &V);
VectorC Set_Part(const VectorC &V, const VectorC &W, long);
Vector  Sqr(const VectorC &V); 
Vector  Abs(const VectorC &V);
Vector  Re(const VectorC &V);
Vector  Im(const VectorC &V);
VectorC Conj(const VectorC &V);
VectorC Cmplx(const Vector &V);
VectorC Cmplx(const Vector &V1, const Vector &V2);
VectorC Constant(cmplx z, long N);

/////////////////////////////////////////////////////////////////////////
// Class matrix:
/////////////////////////////////////////////////////////////////////////

class MatrixC
{
public:
     cmplx *D;
     long N1,N2;

     MatrixC();
     MatrixC(long n, long m=0);
     MatrixC(const MatrixC &);
     MatrixC(const Matrix &);
     ~MatrixC();
     void Start();
     void Create(long n, long=0);
     void Load(cmplx*, long, long=0);
     void Load_Copy(cmplx*, long, long=0);
     void Transfer(MatrixC &M);
     void Destroy();
     void Resize(long,long=0);

     void Zero();
     void Unit();
     void Set(cmplx);
     cmplx& Elem(long, long);
     void Set_Col(const VectorC &,long);
     void Set_Row(const VectorC &, long);      
     void Set_Diag(const VectorC &);      
     void Append_Col(const VectorC &);
     void Append_Col(const MatrixC &);
     void Append_Row(const VectorC &);
     void Append_Row(const MatrixC &);
     void Insert_Col(const VectorC &, long);
     void Insert_Row(const VectorC &, long);
     void Remove_Col(long);
     void Remove_Row(long);
     void Swap_Cols(long,long);
     void Swap_Rows(long,long);
     void Permute_Cols(const List &L);
     void Permute_Rows(const List &L);
     void Sort_Cols(const Vector &V);
     void T();
     void Herm(); // hermitian conjugate
     void Part(long, long, long, long);
     void Add_Part(const MatrixC &,long, long);
     void Set_Part(const MatrixC &,long, long);

     void Change_Basis(const MatrixC &);
     bool Gram_Schmidt();
     bool Invert();
     bool LU_Decomp(int *);          
     void Re();
     void Im();
     void Abs();
     void Sqr();
     void Conj();

     bool Is_Zero(double tolerance) const;
     cmplx Elem(long, long) const;
     VectorC Col(long n) const;
     VectorC Row(long n) const;
     VectorC Diag() const;
     void Col(VectorC &, long) const;
     void Row(VectorC &, long) const;

     cmplx Elem(const VectorC &, const VectorC &) const;
     cmplx Elem(const MatrixC &, long, long) const;
     bool  Solve(VectorC &) const;
     bool  Solve(MatrixC &) const;
     cmplx Det() const; 
     cmplx Trace() const; 
     void  Diagonalize(MatrixC &,Vector &) const; 
     void  Spectrum(Vector &) const;
     void  Tridiagonalize(MatrixC &, VectorC &, VectorC &) const;
     void  NH_Diagonalize(MatrixC&, VectorC &) const; 
     void  SVD(MatrixC &, MatrixC &, Vector &) const;

     void Write(int prec=0) const;
     bool Save_Binary(const char *name) const;
     bool Load_Binary(const char *name);
     bool Save_Binary(FILE *fich) const;
     bool Load_Binary(FILE *fich);
     bool Save(const char *name) const;
     bool Load(const char *name);
     bool Save(FILE *fich) const;
     bool Load(FILE *fich);

     cmplx& operator()(long, long);
     cmplx operator() (long, long) const;
     MatrixC& operator=(const MatrixC&);
     MatrixC& operator=(const Matrix&);
     void operator+=(const MatrixC&);
     void operator+=(const Matrix&);
     void operator-=(const MatrixC&);
     void operator-=(const Matrix&);
     void operator*=(cmplx);  
     void operator*=(double);
     void operator*=(const MatrixC&);
     void operator*=(const Matrix&);

     void operator&=(const VectorC &);
     void operator&=(const MatrixC &);
     void operator|=(const VectorC &);
     void operator|=(const MatrixC &);
};

void Copy(MatrixC& B, const MatrixC& A);
void Copy(MatrixC& B, const Matrix& A); 
MatrixC Unit_C(long,long=0);
MatrixC Zero_C(long,long=0);
MatrixC Diag(const VectorC &V);
MatrixC Diag_C(const Vector &V);

MatrixC operator-(const MatrixC &);
MatrixC operator+(const MatrixC &, const MatrixC &);
MatrixC operator+(const MatrixC &, const Matrix &);
MatrixC operator-(const MatrixC &, const MatrixC &);
MatrixC operator*(cmplx K,  const MatrixC &A);
MatrixC operator*(double K, const MatrixC &A);
VectorC operator*(const MatrixC &A, const VectorC &V);
VectorC operator*(const MatrixC &A, const Vector &V);
MatrixC operator*(const MatrixC &A, const MatrixC &B);
MatrixC operator*(const MatrixC &A, const Matrix &B);

MatrixC operator|(const VectorC &, const VectorC &);
MatrixC operator|(const MatrixC &, const VectorC &);
MatrixC operator|(const VectorC &, const MatrixC &);
MatrixC operator|(const MatrixC &, const MatrixC &);
MatrixC operator&(const MatrixC &, const VectorC &);
MatrixC operator&(const VectorC &, const MatrixC &);
MatrixC operator&(const MatrixC &, const MatrixC &);

void Multiply(VectorC &, const MatrixC &M, const VectorC &V);
void Multiply(MatrixC &, const MatrixC &M1, const MatrixC &M2);
void Multiply_Add(MatrixC &, const MatrixC &M1, const MatrixC &M2,
		  cmplx alpha, cmplx beta, bool T1, bool T2);

MatrixC Tens_Prod(const MatrixC &, const MatrixC &);
void    Tens_Prod(MatrixC &, const MatrixC &, const MatrixC &);
MatrixC Tens_Prod_Unit(const MatrixC &, long, Side);
void    Tens_Prod_Unit(MatrixC &, const MatrixC &, long, Side);
MatrixC Tens_Prod_Diag(const MatrixC &, const VectorC &, Side s);
void    Tens_Prod_Diag(MatrixC &, const MatrixC &, const VectorC &, Side s);

void    Ket_Bra(MatrixC &, const VectorC &, const VectorC &);
MatrixC Ket_Bra(const VectorC &, const VectorC &);
MatrixC Projector(const VectorC &);

void    Sort(Vector &, MatrixC &);
MatrixC Change_Basis(const MatrixC &, const MatrixC &B);
MatrixC T(const MatrixC &);
MatrixC Herm(const MatrixC &);
MatrixC Part(const MatrixC &, long, long, long, long);
void    Part(MatrixC &, const MatrixC &, long, long, long, long);
MatrixC Part(const MatrixC &A, const List &L1, const List &L2);
MatrixC Invert(const MatrixC &M);
VectorC Solve(const MatrixC &A, const VectorC &b);
cmplx   Det(const MatrixC &M);
cmplx   Trace(const MatrixC &M);
double  Norm(const MatrixC &M);

Matrix  Re(const MatrixC &M);
Matrix  Im(const MatrixC &M);
Matrix  Abs(const MatrixC &M);
Matrix  Sqr(const MatrixC &M);
MatrixC Conj(const MatrixC &M);
MatrixC Cmplx(const Matrix &M);


/////////////////////////////////////////////////////////////
/// BLAS-LAPACK HEADERS
extern "C"{
     // V <- alpha *V
     void zscal_(long *N, cmplx *alpha, cmplx* V, long *ix);
     // V <- V + W
     void zaxpy_(long *N, cmplx *alpha, cmplx *X, long *ix, cmplx *Y, long *iy);
     // Dot product
     cmplx zdotc_(long *N, cmplx *X, long *ix, cmplx *Y, long *iy);
     // Matrix-vector product
     void zgemv_(char *, long* n1, long* n2, cmplx* alpha, cmplx* A, long *lda,
		 cmplx* X, long *incx, cmplx* beta, cmplx* Y, long* incy);
     // Matrix-matrix product
     void zgemm_(char*, char*, long* n1, long* n2, long*k, cmplx* alpha, 
		 cmplx* A, long* lda, cmplx* B, long *ldb,
		 cmplx* beta, cmplx* C, long *ldc);
     // LU decomposition
     void zgetrf_(long*, long*, cmplx*, long*, int*, long*);
     // Linear equations solving
     void zgesv_(long*, long*, cmplx*, long*, long*, cmplx*, long*, long*);
     // Tridiagonal Matrix diagonalization
     void zsteqr_(char* compz, long* n, double* d, double* e,
		  cmplx* z, long* ldz, cmplx* work,long* info);
     // Full Diagonalization, expert driver
     void zheevx_(char* jobz, char* range, char* uplo, long* n, cmplx* A,
		  long* lda, double* vl, double* vu,
		  long* il, long* iu, double* abstol, long* M, double* W,
		  cmplx* Z, long* ldz, cmplx* work,
		  long* lwork, double* rwork, long* iwork,
		  long* ifail, long* info);
     // Non-symmetric full diagonalization, non-expert driver
     void zgeev_(char* jobl, char* jobr, long* N, cmplx* A, long* lda,
		 cmplx* w, cmplx* vl, long* ldvl, cmplx* vr, long* ldvr, 
		 cmplx* work, long* lwork, double* rwork, long* info);
     // Reduce to tri-diagonal form, basis is in "strange form"
     void zhetrd_(char* uplo, long* N, cmplx* A, long* lda, cmplx* D, cmplx* E,
		  cmplx* tau, cmplx* work, long* lwork, long* info);
     // Compute the basis in "normal form" for tri-diagonal reduction
     void zungtr_(char* uplo, long* N, cmplx* A, long* lda, cmplx* tau,
		  cmplx* work, long* lwork, long* info);
     // SVD
     void zgesvd_(char* jobu, char* jobvt, long *M, long *N, cmplx *A,
		  long *lda, double *S, cmplx *U, long *ldu, cmplx *Vt,
		  long *ldvt, cmplx *work, long *lwork, double *rwork,
		  long *info);
}

#endif



