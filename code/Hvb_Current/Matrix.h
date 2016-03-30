////////////////////////////////////////////////////////
// hvb++ 1.0
// C++ library for linear algebra, JaviRL, started Nov 1999.
// Copyleft: Javier Rodríguez Laguna
// 1999-2003-2005-2006-2008-2015
#ifndef MATRIX_HEADER
#define MATRIX_HEADER
#include "Common.h"


// If you want to trace your errors easily, write this in your program
// #define DEBUG

#ifdef DEBUG
extern long Num_Vec;
extern long Num_Mat;
extern long Total_Mem;
extern long Max_Total_Mem;
void Mem_Control(long nv, long nm, long tm);
void Mem_Status();
#endif

///////////////////////////////////////////////////////////////////////////
// Class vector:
///////////////////////////////////////////////////////////////////////////

class Vector
{
public:
     long N; // length of the vector 
     double *D; // the pointer itself

     Vector();     
     Vector(long);
     Vector(const Vector &);
     Vector(double *, long);
     ~Vector();
     void Start();
     void Create(long n);
     void Load(double *, long);
     void Load_Copy(double *, long);
     void Transfer(Vector &);
     void Destroy();
     
     void Zero();
     void Set(double);
     void Set_Part(const Vector &, long);
     bool Normalize();
     void Part(long, long);
     void Reverse();
     void Append(double); 
     void Append(const Vector &);
     void Insert(double, long i);
     void Insert(const Vector &V, long i);
     void Remove(long i, long j=0);
     void Sqr();
     void Abs();
     void Sort(int p=1);
     void Rand(double,double);
     void Rand_Gaussian(double,double);

     void Write(long prec=0) const;
     void Write_Col() const;
     bool Save_Binary(const char *name) const;
     bool Load_Binary(const char *name);
     bool Save_Binary(FILE *fich) const;
     bool Load_Binary(FILE *fich);
     bool Save(const char *name) const;
     bool Load(const char *name);
     bool Save(FILE *fich) const;
     bool Load(FILE *fich);

     double Min() const;
     double Min(long &) const;
     double Max() const;
     double Max(long &) const;
     double Sum() const;
     double Sum(long i1, long i2) const;
     double Prod() const;
     double Prod(long i1, long i2) const;
     double Average() const;
     double Deviation() const;
     double Variance() const;
     bool   Is_Zero(double tolerance) const;
     double Norm() const;

     double  operator() (long) const;
     double& operator() (long);
     Vector& operator=(const Vector&);
     void    operator+=(const Vector&);
     void    operator-=(const Vector&);
     void    operator*=(double);
     void    operator/=(double);
     void    operator+=(double);
     void    operator-=(double);

     void operator&=(const Vector &);
     void operator&=(double);
}; 

void   Copy(Vector& B, const Vector& A);
Vector operator-(const Vector &);
Vector operator+(const Vector &, const Vector &);
Vector operator-(const Vector &, const Vector &);
Vector operator+(const Vector &V, double p); 
Vector operator+(double p, const Vector &V); 
Vector operator-(const Vector &V, double p);
Vector operator-(double p, const Vector &V);
Vector operator*(double, const Vector &);
Vector operator*(const Vector &, double);
Vector operator/(const Vector &, double);

Vector operator&(const Vector &L1, const Vector &L2);
Vector operator&(double, const Vector &L);
Vector operator&(const Vector &L, double);

double Dot(const Vector &, const Vector &);
Vector Cross(const Vector &, const Vector &);
Vector Tens_Prod(const Vector &, const Vector &);
void   Tens_Prod(Vector &, const Vector &, const Vector &);
void   Daxpy(Vector &, const Vector &, double);
Vector Elem_Mult(const Vector &, const Vector &);
void   Elem_Mult(Vector &, const Vector &, const Vector &);

Vector Normalize(const Vector &);
Vector Part(const Vector &,long, long);
Vector Part(const Vector &, const List &);
Vector Reverse(const Vector &);
Vector Set_Part(const Vector &, const Vector&, long);
Vector Sqr(const Vector &);
Vector Abs(const Vector &);

double Norm(const Vector &);
double Min(const Vector &);
double Max(const Vector &);
double Sum(const Vector &);
double Average(const Vector &);
double Deviation(const Vector &);
double Variance(const Vector &);

Vector Range(double x0, double x1, long N);
Vector Constant(double x, long N);
Vector To_Vector(const List &L);
List   To_List(const Vector &V);

/////////////////////////////////////////////////////////////////////////
// Class matrix:
/////////////////////////////////////////////////////////////////////////

class Matrix
{
public:
     double *D;
     long N1,N2;

     Matrix();
     Matrix(long n, long m=0);
     Matrix(const Matrix &);
     ~Matrix();
     void Start();
     void Create(long n, long=0);
     void Load(double *, long, long=0);
     void Load_Copy(double *, long, long=0);
     void Transfer(Matrix &M);
     void Destroy();
     void Resize(long,long=0);

     void Zero();
     void Unit();
     void Set(double);
     double& Elem(long, long);
     void Set_Col(const Vector &, long);
     void Set_Row(const Vector &, long);
     void Set_Diag(const Vector &);
     void Append_Col(const Vector &);
     void Append_Col(const Matrix &);
     void Append_Row(const Vector &);
     void Append_Row(const Matrix &);
     void Insert_Col(const Vector &, long);
     void Insert_Row(const Vector &, long);
     void Remove_Col(long);
     void Remove_Row(long);
     void Swap_Cols(long,long);
     void Swap_Rows(long,long);
     void Permute_Cols(const List &L);
     void Permute_Rows(const List &L);
     void Sort_Cols();
     void Sort_Cols(const Vector &V);
     void T();
     void Part(long, long, long, long);
     void Add_Part(const Matrix &, long, long);
     void Set_Part(const Matrix &, long, long);
     void Sqr();
     void Rand(double,double);
     void Rand_Gaussian(double,double);

     void Change_Basis(const Matrix &);
     bool Gram_Schmidt();
     bool Invert();
     bool LU_Decomp(int *);

     bool   Is_Zero(double tolerance) const;
     double Elem(long, long) const;
     Vector Col(long n) const;
     Vector Row(long n) const;
     Vector Diag() const;
     void   Col(Vector &, long) const;
     void   Row(Vector &, long) const;
     double Max() const;
     double Max(long &, long &) const;
     double Min() const;
     double Min(long &, long &) const;

     double Elem(const Vector &, const Vector &) const;
     double Elem(const Matrix &, long, long) const;
     bool   Solve(Vector &) const;
     bool   Solve(Matrix &) const;
     double Det() const; 
     double Trace() const; 
     void   Diagonalize(Matrix &,Vector &) const; 
     void   Spectrum(Vector &) const;
     void   Tridiagonalize(Matrix &, Vector &, Vector &) const;
     void   NS_Diagonalize(Matrix&, Matrix&, Vector &, Vector &) const; 

     void Write(int prec=0) const;
     bool Save_Binary(const char *name) const;
     bool Load_Binary(const char *name);
     bool Save_Binary(FILE *fich) const;
     bool Load_Binary(FILE *fich);
     bool Save(const char *name) const;
     bool Save(const char *name, const char *comment) const;
     bool Load(const char *name);
     bool Save(FILE *fich) const;
     bool Save(FILE *fich, const char *comment) const;
     bool Load(FILE *fich);

     double& operator()(long, long);
     double operator() (long, long) const;
     Matrix& operator=(const Matrix&);
     void operator+=(const Matrix&);
     void operator-=(const Matrix&);
     void operator*=(const Matrix&);
     void operator+=(double);
     void operator-=(double);
     void operator*=(double);  
     void operator/=(double);

     void operator&=(const Vector &);
     void operator&=(const Matrix &);
     void operator|=(const Vector &);
     void operator|=(const Matrix &);
};

void   Copy(Matrix& B, const Matrix& A);
Matrix Zero(long,long=0);
Matrix Unit(long,long=0);
Matrix Diag(const Vector &E);
void Write(const Matrix &M);
Matrix Constant(double x, long, long);

Matrix operator-(const Matrix &);
Matrix operator+(const Matrix &, const Matrix &);
Matrix operator-(const Matrix &, const Matrix &);
Matrix operator*(double K, const Matrix &A);
Matrix operator*(const Matrix &A, double K);
Matrix operator/(const Matrix &A, double K);
Vector operator*(const Matrix &A, const Vector &V);
Matrix operator*(const Matrix &A, const Matrix &B);

Matrix operator|(const Vector &, const Vector &);
Matrix operator|(const Matrix &, const Vector &);
Matrix operator|(const Vector &, const Matrix &);
Matrix operator|(const Matrix &, const Matrix &);
Matrix operator&(const Matrix &, const Vector &);
Matrix operator&(const Vector &, const Matrix &);
Matrix operator&(const Matrix &, const Matrix &);

void Multiply(Vector &, const Matrix &M, const Vector &V);
void Multiply(Matrix &, const Matrix &M1, const Matrix &M2);
void Multiply_Add(Matrix &, const Matrix &M1, const Matrix &M2,
		  double alpha, double beta, bool T1, bool T2);

Matrix Tens_Prod(const Matrix &, const Matrix &);
void   Tens_Prod(Matrix &, const Matrix &, const Matrix &);
Matrix Tens_Prod_Unit(const Matrix &, long, Side);
void   Tens_Prod_Unit(Matrix &, const Matrix &, long, Side);
Matrix Tens_Prod_Diag(const Matrix &, const Vector &, Side s);
void   Tens_Prod_Diag(Matrix &, const Matrix &, const Vector &, Side s);

void   Elem_Mult(Matrix &, const Matrix &, const Matrix &);
Matrix Elem_Mult(const Matrix &, const Matrix &);
void   Ket_Bra(Matrix &, const Vector &, const Vector &);
Matrix Ket_Bra(const Vector &, const Vector &);
Matrix Projector(const Vector &);

// Sort V and swap columns of M accordingly
void   Sort(Vector &V, Matrix &M);
Matrix Change_Basis(const Matrix &M, const Matrix &B);
Matrix T(const Matrix &M);
Matrix Part(const Matrix &M, long,long,long,long);
void   Part(Matrix &R,const Matrix &M,long,long,long,long);
Matrix Part(const Matrix &R, const List &L1, const List &L2);
Matrix Invert(const Matrix &M);
Vector Solve(const Matrix &A, const Vector &b);
double Trace(const Matrix &M);
double Det(const Matrix &M);
double Norm(const Matrix &M);
Matrix Sqr(const Matrix &M);

Matrix To_Matrix(const Table &T);
Table  To_Table(const Matrix &M);

void   Trid_Spectrum(Vector &D, Vector &S);
void   Trid_Diagonalize(Matrix &B, Vector &D, Vector &S);

/////////////////////////////////////////////////////////////
/// BLAS-LAPACK Headers
////////////////////////////////////////////////////////////
extern "C"{
     // obtain machine parameters
     double dlamch_(char *c);
     // V <- alpha *V
     void dscal_(long *N, double *alpha, double* V, long *ix);
     // V <- V + W
     void daxpy_(long *N, double *alpha, double *X, long *ix, double *Y, 
		 long *iy);
     // Dot product
     double ddot_(long *N, double *X, long *ix, double *Y, long *iy);
     // Matrix-vector product
     void dgemv_(char *, long* n1, long* n2, double* alpha, double* A, 
		 long *lda, double* X, long *incx, double* beta, double* Y, 
		 long* incy);
     // Matrix-matrix product
     void dgemm_(char*, char*, long* n1, long* n2, long*k, double* alpha, 
		 double* A, long* lda, double* B, long *ldb,
		 double* beta, double* C, long *ldc);
     // LU decomposition
     void dgetrf_(long*, long*, double*, long*, int*, long*);
     // Linear equations solving
     void dgesv_(long*, long*, double*, long*, long*, double*, long*, long*);
     // Tridiagonal Matrix diagonalization
     void dsteqr_(char* compz, long* n, double* d, double* e,
		  double* z, long* ldz, double* work,long* info);
     // Full Diagonalization, expert driver
     void dsyevx_(char*,char*,char*,long*,double*,long*,double*,double*,
		  long*,long*,double*,long*,double*,double*,long*,double*,
		  long*,long*,long*,long*);
     // Non-symmetric full diagonalization, non-expert driver
     void dgeev_(char* jobl, char* jobr, long* N, double* A, long* lda,
		 double* wr, double* wi, double* vl, long* ldvl,
		 double* vr, long* ldvr, double* work, long* lwork, 
		 long* info);
     // Reduce to tri-diagonal form, basis is in "strange form"
     void dsytrd_(char* uplo, long* N, double* A, long* lda, double* D, 
		  double* E, double* tau, double* work, long* lwork, 
		  long* info);
     // Compute the basis in "normal form" for tri-diagonal reduction
     void dorgtr_(char* uplo, long* N, double* A, long* lda, double* tau,
		  double* work, long* lwork, long* info);
}

#endif



