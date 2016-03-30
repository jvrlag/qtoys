////////////////////////////////////////////////////////
// hvb++ 
// Copyleft: Javier Rodríguez Laguna
// 060224-080725-101022-121004-150123
#include "Matrix.h"
#include "Text.h"

#ifdef DEBUG
long Num_Vec=0;
long Num_Mat=0;
long Total_Mem=0;
long Max_Total_Mem=0;

void Mem_Control(long nv, long nm, long tm)
{
     Num_Vec+=nv;
     Num_Mat+=nm;
     Total_Mem+=tm;
     if (Total_Mem>Max_Total_Mem) Max_Total_Mem=Total_Mem;
}  

void Mem_Status()
{
     // To catch any possible leaks
     printf("** Memory Status:\n");
     printf("Maximum number of doubles alloc'd: %ld\n",Max_Total_Mem);
     printf("Total doubles allocated: %ld\n",Total_Mem);
     printf("Number of stored matrices: %ld\n", Num_Mat);          
     printf("Number of stored vectors: %ld\n\n",Num_Vec);
}
#endif

/////////////////////////////////////////////////////////////////
// Vector creators & memory handling
/////////////////////////////////////////////////////////////////

Vector::Vector(): N(0), D(NULL) {} 

Vector::Vector(long n)
{
      Start();
      Create(n);
      Zero();
}

Vector::Vector(const Vector &V)
{
     Start();
     Create(V.N);
     if (N) memcpy(D,V.D,(N+1)*sizeof(double));
}

Vector::Vector(double *data, long n)
{
     Start();
     Create(n);
     memcpy(D,data,(n+1)*sizeof(double));
}

Vector::~Vector() { Destroy(); }

// Call Start if you assume that the Vector does not exist. 
// Useful when we're initializing an array of them.
void Vector::Start()
{
     N=0; 
     D=(double*)NULL;
}

// Call Create if it is sure that the Vector existed beforehand.
// It is destroyed, then allocated.
void Vector::Create(long n)
{
     Destroy();
     N=n;
     if (!N) { D=NULL; return;}
     D=(double*)malloc((n+1)*sizeof(double));    
     
#ifdef DEBUG 
     if (!D) Error("Error allocating vector."); 
     Mem_Control(1,0,N);
#endif
}

// Load a double* into a Vector
void Vector::Load(double *d, long n)
{
     Destroy();
     N=n;
     D=d;
#ifdef DEBUG
     Mem_Control(1,0,N);
#endif
}

// Copy a double* into a Vector
void Vector::Load_Copy(double *d1, long n)
{
     double *d2=(double*)malloc((n+1)*sizeof(double));
     memcpy(d2,d1,(n+1)*sizeof(double));
     Load(d2,n);
}

// The data of V go to our vector, V is left empty.
void Vector::Transfer(Vector &V)
{
     Destroy();
     D=V.D;
     N=V.N;
     V.Start();
}

void Vector::Destroy()
{
     if (N) 
     { 	
	  free(D);
#ifdef DEBUG
	  Mem_Control(-1,0,-N);
#endif
	  N=0;
	  D=NULL;
     }
}

/////////////////////////////////////////////////////////////////
// Vector Transformations: Element handling
/////////////////////////////////////////////////////////////////

void Vector::Zero()
{
     if (N) memset(D,0,(N+1)*sizeof(double));
}

// Set the value to a given one.
void Vector::Set(double x)
{
     if (N) for (long i=1;i<=N;i++) D[i]=x;
}

void Vector::Set_Part(const Vector& V, long n)
{
#ifdef DEBUG
     if (n+V.N-1>N)
	  Error("Set_Part of a vector which is too big.");
#endif
     memcpy(D+n,V.D+1,V.N*sizeof(double));
}

// return false if normalization was not possible!
bool Vector::Normalize()
{
     double norm=Norm();
     if (norm==0.0) 
     {
	  Error_Flag(Error_Mat);
	  return false;
     }
     (*this)/=norm;
     return true;
}

void Vector::Part(long n1, long n2)
{
#ifdef DEBUG
     if (n1<1 || n2>N) 
	  Error("The part can't be larger than the whole!");
#endif
     Vector R(n2-n1+1);
     memcpy(R.D+1,D+n1,(n2-n1+1)*sizeof(double));
     Transfer(R);
}

void Vector::Reverse()
{
     for (long i=1;i<=N/2;i++)
	  Swap(D[i],D[N+1-i]);
}

void Vector::Append(double x) 
{
     if (!N) Create(1);
     else
     {
	  N++;
	  D=(double*)realloc(D,(N+1)*sizeof(double));
#ifdef DEBUG
	  Mem_Control(0,0,1);
#endif
     }
     D[N]=x;     
}

void Vector::Append(const Vector &V)
{
     long nold=N;
     if (!N) Create (V.N);
     else
     {
          N+=V.N;
          D=(double*)realloc(D,(N+1)*sizeof(double));
#ifdef DEBUG
	  Mem_Control(0,0,V.N);
#endif
     }
     for (long i=1;i<=V.N;i++)
          D[nold+i]=V(i);
}

// Insert x at position i, Vector size increases by 1!
void Vector::Insert(double x, long i)
{
     N++;
     D=(double*)realloc(D,(N+1)*sizeof(double));
#ifdef DEBUG
     Mem_Control(0,0,1);
#endif
     memmove(D+i+1,D+i,(N-i)*sizeof(double));
     D[i]=x;
}

// Insert V at position i, Vector size increases by V.N
void Vector::Insert(const Vector &V, long i)
{
     N+=V.N;
     D=(double*)realloc(D,(N+1)*sizeof(double));
#ifdef DEBUG
     Mem_Control(0,0,V.N);
#endif
     memmove(D+i+V.N,D+i,(N-i)*sizeof(double));
     memcpy(D+i+1,V.D+1,V.N*sizeof(double));
}

void Vector::Remove(long i, long j) // Remove from i to j, both included
{
     long n0=i, nf=(j?j:i);
     long n=nf-n0+1;
#ifdef DEBUG
     Mem_Control(0,0,-n);
#endif
     if (n0==1 && nf==N) { Destroy(); return; }
     memmove(D+n0,D+nf+1,(N-nf)*sizeof(double));
     N-=n;
     D=(double*)realloc(D,(N+1)*sizeof(double));
}

void Vector::Sqr()
{
     for (long i=1;i<=N;i++)
	  D[i]=::Sqr(D[i]);
}

void Vector::Abs()
{
     for (long i=1;i<=N;i++)
	  D[i]=abs(D[i]);
}

void Vector::Sort(int p)
{     
     if (p<0)
	  qsort(D+1,N,sizeof(double),Sort_Double_21);
     else
	  qsort(D+1,N,sizeof(double),Sort_Double_12);
}

void Vector::Rand(double a, double b)
{
#ifdef DEBUG
     if (!Rand_Started) Error("Using Rand without Rand_Open\n");
#endif
     for (long i=1;i<=N;i++)
	  D[i]=::Rand(a,b);
}

void Vector::Rand_Gaussian(double a, double b)
{
#ifdef DEBUG
     if (!Rand_Started) Error("Using Rand without Rand_Open\n");
#endif
     for (long i=1;i<=N;i++)
	  D[i]=::Rand_Gaussian(a,b);
}

/////////////////////////////////////////////////////////////////
// Vector I/O routines
/////////////////////////////////////////////////////////////////

void Vector::Write(long prec) const
{
     char form[40];
     if (!prec) sprintf(form,"%%12.8g ");
     else sprintf(form,"%%%ld.%ldg ",prec+4,prec);
     for (long i=1;i<=N;i++)
	  printf(form,D[i]);
     printf("\n\n");
}

void Vector::Write_Col() const
{
     for (long i=1;i<=N;i++)
	  printf("%16.12g \n",D[i]);
     printf("\n");
}

bool Vector::Save_Binary(FILE *fich) const
{
     int ausgang;
     ausgang=fwrite(&N,sizeof(long),1,fich);
     if (ausgang!=1) return false;
     if (!N) return true;
     ausgang=fwrite(D,sizeof(double),N+1,fich);
     if (ausgang!=N+1) return false;
     return true;
}

bool Vector::Save_Binary(const char *name) const
{
     FILE *fich=fopen(name,"wb");
     if (!fich) return false;
     bool status=Save_Binary(fich);
     fclose(fich);
     return status;
}

bool Vector::Load_Binary(FILE *fich)
{
     int ausgang;
     ausgang=fread(&N,sizeof(long),1,fich);
     if (ausgang!=1) return false;
     if (!N) { D=(double*)NULL; return true; }
     Create(N);
     Zero();
     ausgang=fread(D,sizeof(double),N+1,fich);
     if (ausgang!=N+1) return false;
     return true;
}

bool Vector::Load_Binary(const char *name)
{
     FILE *fich=fopen(name,"rb");
     if (!fich) return false;
     bool status=Load_Binary(fich);
     fclose(fich);
     return status;
}

bool Vector::Save(const char *name) const
{
     FILE *fich=fopen(name,"wt");
     if (!fich) return false;
     bool status=Save(fich);
     fclose(fich);
     return status;
}

bool Vector::Save(FILE *fich) const
{
     if (fprintf(fich,"# %ld\n",N)<0) return false;
     for (long i=1;i<=N;i++)
	  if(fprintf(fich,"%20.16g\n",D[i])<0) return false;
     return true;
}

bool Vector::Load(const char *name) 
{
     FILE *fich=fopen(name,"rt");
     if (!fich) return false;
     bool status=Load(fich);
     fclose(fich);
     return status;
}

bool Vector::Load(FILE *fich)
{
     if (!fich) return false;
     long n; double x;
     if (!fscanf(fich,"# %ld\n",&n)) return false;
     Create(n);
     Zero();
     for (long i=1;i<=n;i++)
     {
	  if (!fscanf(fich,"%lg\n",&x)) return false;
	  D[i]=x;
     }
     return true;
}

/////////////////////////////////////////////////////////////////
// Vector: Functions
/////////////////////////////////////////////////////////////////

// Minimum, don't care for location
double Vector::Min() const
{
     double E0=D[1];
     for (long i=2;i<=N;i++)
	  if (D[i]<E0) E0=D[i];
     return E0;
}

// Minimum, return location in imin
double Vector::Min(long &imin) const
{
     double E0=D[1];
     imin=1;
     for (long i=2;i<=N;i++)
	  if (D[i]<E0) { E0=D[i]; imin=i; }
     return E0;
}

// Maximum, don't care for location
double Vector::Max() const
{
     double E0=D[1];
     for (long i=2;i<=N;i++)
	  if (D[i]>E0) E0=D[i];
     return E0;
}

// Maximum, return the location in imax
double Vector::Max(long &imax) const
{
     double E0=D[1];
     imax=1;
     for (long i=2;i<=N;i++)
	  if (D[i]>E0) { E0=D[i]; imax=i; }
     return E0;
}

double Vector::Sum() const
{
     return Sum(1,N);
}

// Sum from i1 to i2
double Vector::Sum(long i1, long i2) const
{
#ifdef DEBUG
     if (!N) return 0.0;
     if (i1<=0 || i1>i2 || i2>N) Error("Incorrect sum limits\n");
#endif
     double sum=D[i1];
     for (long i=i1+1;i<=i2;i++)
	  sum+=D[i];
     return sum;
}

double Vector::Prod() const
{
     return Prod(1,N);
}

double Vector::Prod(long i1, long i2) const
{
#ifdef DEBUG
     if (!N) return 0.0;
     if (i1<=0 || i1>i2 || i2>N) Error("Incorrect prod limits\n");
#endif
     double prod=D[i1];
     for (long i=i1+1;i<=i2;i++)
	  prod*=D[i];
     return prod;
}

double Vector::Average() const
{
     return Sum()/(double)N;
}

double Vector::Deviation() const
{
     return sqrt(Variance());
}

double Vector::Variance() const
{
     if (N==1) return 0.0;
     double sumsq=0.0;
     double aver=Average();
     for (long i=1;i<=N;i++)
          sumsq+=::Sqr(D[i]-aver);
     sumsq/=(double)N;
     return sumsq;
}

// Find if the vector is zero within a given tolerance
bool Vector::Is_Zero(double tolerance) const
{
     for (long i=1;i<=N;i++)
	  if (fabs(D[i])>tolerance) return false;
     return true;
}

double Vector::Norm() const
{
     long ix=1, n=N;
     double norm=ddot_(&n,D+1,&ix,D+1,&ix);
     return sqrt(norm);
}

/////////////////////////////////////////////////////////////////
// Vector: Overloaded Operators
/////////////////////////////////////////////////////////////////

double Vector::operator() (long n) const
{
#ifdef DEBUG
     if (n<0 || n>N) Error("Error getting vector comp."); 
#endif
    return D[n];
}

double& Vector::operator() (long n)
{
#ifdef DEBUG
    if (n<0 || n>N) Error("Error putting vector comp."); 
#endif
    return D[n];
}

Vector& Vector::operator=(const Vector& W)
{
     if (this==&W) return *this;
     Copy(*this,W);
     return (*this);
}

void Vector::operator+=(const Vector& W)
{
     if (!N) { Create(W.N); Zero(); }
#ifdef DEBUG
      if (N!=W.N) Error("Incompatible sizes in vector +=");
#endif
      double alpha=1.0;
      long ix=1;
      daxpy_(&N,&alpha,W.D+1,&ix,D+1,&ix);
}

void Vector::operator-=(const Vector& W)
{
     if (!N) { Create(W.N); Zero(); }
#ifdef DEBUG
     if (N!=W.N) Error("Incompatible sizes in vector -=");
#endif
     double alpha=-1.0;
     long ix=1;
     daxpy_(&N,&alpha,W.D+1,&ix,D+1,&ix);
}

void Vector::operator*=(double x)
{
#ifdef DEBUG
    if (!N) return;
#endif
    long ix=1;
    double cosa=x;
    dscal_(&N,&cosa,D+1,&ix);
}

void Vector::operator/=(double x)
{
#ifdef DEBUG
     if (!N) return;
#endif
     long ix=1;
     double cosa=1.0/x;
     dscal_(&N,&cosa,D+1,&ix);
}

void Vector::operator+=(double x)
{
#ifdef DEBUG
     if (!N) return;
#endif
     for (long i=1;i<=N;i++)
	  D[i]+=x;
}

void Vector::operator-=(double x)
{
#ifdef DEBUG
     if (!N) return;
#endif
     for (long i=1;i<=N;i++)
	  D[i]-=x;
}

void Vector::operator&=(double p)
{
     Append(p);
}

void Vector::operator&=(const Vector &L)
{
     Append(L);
}

/////////////////////////////////////////////////////////////////
// Vector: External functions
/////////////////////////////////////////////////////////////////

void Copy(Vector& B, const Vector& A)
{
     if (!A.N) return;
     B.Create(A.N);
     memcpy(B.D,A.D,(A.N+1)*sizeof(double));
}

Vector operator-(const Vector &M)
{
     Vector R(M);
     R*=-1.0;
     return R;
}

Vector operator+(const Vector &A, const Vector &B)
{
#ifdef DEBUG
     if (A.N!=B.N) 
	  Error("Adding vectors with different dimensions.");
#endif
     Vector R(A);
     R+=B;
     return R;
}

// It is correct like this, although it doesn't seem to be!!
// This thing takes its arguments opposite: B-A!!
Vector operator-(const Vector &A, const Vector &B)
{
#ifdef DEBUG
     if (A.N!=B.N) 
	  Error("Substracting vectors with different dimensions.");
#endif
     Vector R(A);
     R-=B;
     return R;
}

Vector operator+(const Vector &V, double K)
{
     return V+Constant(K,V.N);
}

Vector operator+(double K, const Vector &V)
{
     return Constant(K,V.N)+V;
}

Vector operator-(const Vector &V, double K)
{
     return V-Constant(K,V.N);
}

Vector operator-(double K, const Vector &V)
{
     return Constant(K,V.N)-V;
}

Vector operator*(double x, const Vector& V)
{
     Vector R(V);
     R*=x;
     return R;
}

Vector operator*(const Vector& V, double x)
{
     Vector R(V);
     R*=x;
     return R;
}

Vector operator/(const Vector& V, double x)
{
     Vector R(V);
     R/=x;
     return R;
}

Vector operator&(const Vector &L1, const Vector &L2)
{
     Vector L(L1);
     L.Append(L2);
     return L;
}

Vector operator&(double p, const Vector &L2)
{
     Vector L(1); L(1)=p; 
     L.Append(L2);
     return L;
}

Vector operator&(const Vector &L2, double p)
{
     Vector L(L2);
     L.Append(p);
     return L;
}

double Dot(const Vector& V1, const Vector& V2)
{
#ifdef DEBUG
     if (V1.N!=V2.N) Error("Dot product of vectors of diferent dim.");
#endif
     long ix=1, n=V1.N;
     return ddot_(&n,V1.D+1,&ix,V2.D+1,&ix);
}

Vector Cross(const Vector& V, const Vector& W)
{
#ifdef DEBUG
     if (V.N!=3 || W.N!=3) Error ("Cross product: vector dim must be 3.");
#endif
     Vector R(3);
     R(1)=V(2)*W(3)-V(3)*W(2);
     R(2)=V(3)*W(1)-V(1)*W(3);
     R(3)=V(1)*W(2)-V(2)*W(1);
     return R;
}

Vector Tens_Prod(const Vector &V, const Vector &W)
{
     long N=V.N*W.N;
     Vector R(N);
     long k=1;
     for (long i=1;i<=V.N;i++)
	  for (long j=1;j<=W.N;j++)
	  {
	       R.D[k]=V.D[i]*W.D[j];
	       k++;
	  }
     return R;
}

void Tens_Prod(Vector &R, const Vector &V, const Vector &W)
{
     long N=V.N*W.N;
     R.Create(N);
     R.Zero();
     long k=1;
     for (long i=1;i<=V.N;i++)
	  for (long j=1;j<=W.N;j++)
	  {
	       R.D[k]=V.D[i]*W.D[j];
	       k++;
	  }      
}

// R <- alpha * V + R
void Daxpy(Vector &R, const Vector &V, double alpha)
{
     long ix=1, n=R.N;
     daxpy_(&n,&alpha,V.D+1,&ix,R.D+1,&ix);     
}

Vector Elem_Mult(const Vector &A, const Vector &B)
{
     Vector R(A.N);
     Elem_Mult(R,A,B);
     return R;
}

void Elem_Mult(Vector &V, const Vector &A, const Vector &B)
{
#ifdef DEBUG
     if (A.N!=B.N) Error("Different dimensions in Elem_Mult\n");
#endif
     long N=A.N;
     if (V.N!=N) V.Create(N);
     for (long i=1;i<=N;i++)
	  V(i)=A(i)*B(i);
}

Vector Normalize(const Vector &V)
{
     Vector R(V);
     R.Normalize();
     return R;
}

Vector Part(const Vector &V, long n1, long n2) 
{
     Vector R(V);
     R.Part(n1,n2);
     return R;
}

Vector Part(const Vector &V, const List &L)
{
     Vector W(L.N);
     for (long i=1;i<=L.N;i++)
	  W(i)=V(L(i));
     return W;
}

Vector Reverse(const Vector &V)
{
     Vector R(V);
     R.Reverse();
     return R;
}

Vector Set_Part(const Vector &V, const Vector &W, long n)
{
     Vector R(V);
     R.Set_Part(W,n);
     return R;
}

Vector Sqr(const Vector &V)
{
     Vector R(V);
     R.Sqr();
     return R;
}

Vector Abs(const Vector &V)
{
     Vector R(V);
     R.Abs();
     return R;
}

double Norm(const Vector &V) 
{
     return V.Norm();
}

double Min(const Vector &V)
{
     return V.Min();
}

double Max(const Vector &V)
{
     return V.Max();
}

double Sum(const Vector &V)
{
     return V.Sum();
}

double Average(const Vector &V)
{
     return V.Average();
}

double Deviation(const Vector &V)
{
     return V.Deviation();
}

double Variance(const Vector &V)
{
     return V.Variance();
}

Vector Range(double x0, double x1, long N)
{
     if (N<2) Error("Range must take at least two values\n");
     Vector R(N);
     double h=(x1-x0)/(double)(N-1);
     for (long i=1;i<=N;i++)
	  R(i)=x0+(i-1)*h;
     return R;
}

Vector Constant(double x, long N)
{
     Vector R(N);
     R.Set(x);
     return R;
}

Vector To_Vector(const List &L)
{
     Vector R(L.N);
     for (long i=1;i<=L.N;i++)
	  R(i)=(double)L(i);
     return R;
}

List   To_List(const Vector &V)
{
     List L(V.N);
     for (long i=1;i<=V.N;i++)
	  L(i)=(long)round(V(i));
     return L;
}

/////////////////////////////////////////////////////////////////
//// MATRIX
/////////////////////////////////////////////////////////////////

Matrix::Matrix(long n1, long n2) : N1(n1), N2(n2)
{
     Start();
     Create(n1,n2);
     Zero();
}

Matrix::Matrix(const Matrix &M)
{
     Start();
     Copy(*this,M);
}

Matrix::Matrix()
{
     Start();
}

Matrix::~Matrix()
{
     Destroy();
}

void Matrix::Start()
{
     N1=N2=0;
     D=(double*)NULL;
}

void Matrix::Create(long n1, long n2)
{
     Destroy();
     N1=n1; N2=n2;
     if (!N2) N2=N1;
     if (!N1)
     {
	  D=NULL;
	  return;
     }
     D=(double*)malloc((N1*N2+1)*sizeof(double));
#ifdef DEBUG
     if (!D) Error ("Error allocating matrix.");
     Mem_Control(0,1,N1*N2);
#endif
}

// CAUTION: the data must be stored in columns
void Matrix::Load(double* d, long n1, long n2)
{
     Destroy();
     N1=n1;
     N2=n2; if (!N2) N2=N1;
     D=d;
#ifdef DEBUG
     Mem_Control(0,1,N1*N2);
#endif
}

// CAUTION: the data must be stored in columns
void Matrix::Load_Copy(double *d1, long n1, long n2)
{
     if (!n2) n2=n1;
     double *d2=(double*)malloc((n1*n2+1)*sizeof(double));
     memcpy(d2,d1,(n1*n2+1)*sizeof(double));
     Load(d2,n1,n2);
}

void Matrix::Transfer(Matrix &M)
{
     Destroy();
     D=M.D;
     N1=M.N1;
     N2=M.N2;
     M.Start();
}

void Matrix::Destroy()
{
     if (N1) free(D);
     D=NULL;
#ifdef DEBUG
     if (N1*N2) Mem_Control(0,-1,N1*N2);
#endif
     N1=N2=0;
}

void Matrix::Resize(long n1, long n2) 
{
     if (!n2) n2=n1;
     if (n1==N1 && n2==N2) return; 
     long nr1=::Min(N1,n1), nr2=::Min(N2,n2);
     long oldN1=N1, oldN2=N2;
     
     double *D2=(double*)malloc((N1*N2+1)*sizeof(double));
     memcpy(D2+1,D+1,N1*N2*sizeof(double));
     
     Create(n1,n2);
     if (n1>oldN1 || n2>oldN2)
	  Zero(); // In case the new matrix is bigger
     for (long i=1;i<=nr2;i++)
	  memcpy(D+N1*(i-1)+1,D2+oldN1*(i-1)+1,nr1*sizeof(double));
     free(D2);
}

///////////////////////////////////////////////////
// Matrix Transformations: Element Handling
///////////////////////////////////////////////////

void Matrix::Zero()
{
     memset(D+1,0,N1*N2*sizeof(double));
}

void Matrix::Unit()
{
     Zero();
     long n=::Min(N1,N2);
     for (long i=1;i<=n;i++)
	  Elem(i,i)=1.0;
}

void Matrix::Set(double x)
{
     for (long i=1;i<=N1;i++)
	  for (long j=1;j<=N2;j++)
	       Elem(i,j)=x;
}

double& Matrix::Elem(long i, long j)
{
#ifdef DEBUG
     if (i<1 || i>N1 || j<1 || j>N2) 
	  Error("Error accessing matrix element.");
#endif
     return D[N1*(j-1)+i];
}

void Matrix::Set_Col(const Vector &V, long n)
{
#ifdef DEBUG
     if ((V.N!=N1) || (n>N2)) Error("Set_Col is impossible.");
#endif
     memcpy(D+(n-1)*N1+1,V.D+1,V.N*sizeof(double));
}

void Matrix::Set_Row(const Vector &V, long n)
{
     for (long i=1;i<=N2;i++)
	  Elem(n,i)=V(i);
}

void Matrix::Set_Diag(const Vector &V)
{
#ifdef DEBUG
     if (V.N>::Min(N1,N2)) Error("Error in Set_Diag.");
#endif
     for (long i=1;i<=V.N;i++)
	  Elem(i,i)=V(i);
}

void Matrix::Append_Col(const Vector &V)
{
     if (N1)
     	  Resize(N1,N2+1);
     else
     	  Create(V.N,1);
     Set_Col(V,N2);
}

void Matrix::Append_Col(const Matrix &T)
{
     if (N1)
     	  Resize(N1,N2+T.N2);
     else
     {
	  Copy(*this,T);
	  return;
     }
     Set_Part(T,1,N2-T.N2+1);
}

void Matrix::Append_Row(const Vector &V)
{
     if (N2)
     	  Resize(N1+1,N2);
     else
     	  Create(1,V.N);
     Set_Row(V,N1);
}

void Matrix::Append_Row(const Matrix &T)
{
     if (N2)
     	  Resize(N1+T.N1,N2);
     else
     {
	  Copy(*this,T);
	  return;
     }
     Set_Part(T,N1-T.N1+1,1);
}

void Matrix::Insert_Col(const Vector &V, long p)
{
#ifdef DEBUG
     if (V.N!=N1) Error("Insert_Col with wrong dimensions\n");
#endif
     if (!N1)
     {
	  Create(V.N,1);
	  Set_Col(V,1);
     }
     Resize(N1,N2+1);
     for (long i=N2;i>=p+1;i--)
	  memmove(D+1+(i-1)*N1,D+1+(i-2)*N1,N1*sizeof(double));
     Set_Col(V,p);
}

void Matrix::Insert_Row(const Vector &V, long p)
{
#ifdef DEBUG
     if (V.N!=N2) Error("Insert_Row with wrong dimensions\n");
#endif
     if (!N1)
     {
	  Create(1,V.N);
	  Set_Row(V,1);
     }
     Resize(N1+1,N2);
     for (long i=N1;i>=p+1;i--)
	  for (long j=1;j<=N2;j++)
	       Elem(i,j)=Elem(i-1,j);
     Set_Row(V,p);
}

void Matrix::Remove_Col(long p)
{
#ifdef DEBUG
     if (p>N2) Error("Remove non-existent col\n");
#endif
     if (N2==1) { Destroy(); return; }
     for (long i=p+1;i<=N2;i++)
	  memcpy(D+(i-2)*N1+1,D+(i-1)*N1+1,N1*sizeof(double));
     Resize(N1,N2-1);
}

void Matrix::Remove_Row(long p)
{
#ifdef DEBUG
     if (p>N1) Error("Remove non-existent row\n");
#endif
     if (N1==1) { Destroy(); return; }
     Matrix R(N1-1,N2);
     for (long j=1;j<=N2;j++)
     {
	  for (long i=1;i<=p-1;i++)
	       R(i,j)=Elem(i,j);
	  for (long i=p+1;i<=N1;i++)
	       R(i-1,j)=Elem(i,j);
     }
     Transfer(R);	  
}

void Matrix::Swap_Cols(long k1,long k2)
{
     double acum;
     for(long i=1;i<=N1;i++)
     {
	  acum=Elem(i,k1);
	  Elem(i,k1)=Elem(i,k2);
	  Elem(i,k2)=acum;
     }
}

void Matrix::Swap_Rows(long k1,long k2)
{
    double acum;
    for(long i=1;i<=N2;i++)
    {
	 acum=Elem(k1,i);
	 Elem(k1,i)=Elem(k2,i);
	 Elem(k2,i)=acum;
    }
}

void Matrix::Permute_Cols(const List &P)
{
     Matrix A(N1,N2);
     for (long i=1;i<=N2;i++)
	  A.Set_Col(Col(P(i)),i);
     Transfer(A);
}

void Matrix::Permute_Rows(const List &P)
{
     Matrix A(N1,N2);
     for (long i=1;i<=N1;i++)
	  A.Set_Row(Row(P(i)),i);
     Transfer(A);
}

// Sort the columns wrt the first element of each column
void Matrix::Sort_Cols()
{
     qsort(D+1,N2,N1*sizeof(double),Sort_Double_12);
}

// Sort the columns wrt the elements of this vector 
// Dirty trick: create a larger matrix, order by first element, remove
void Matrix::Sort_Cols(const Vector &V) // CAUTION: V is invariant!
{
     Matrix A(N1+1,N2);
     A.Set_Row(V,1);
     A.Set_Part(*this,2,1);
     A.Sort_Cols();
     A.Part(2,1,N1+1,N2);
     Transfer(A);
}

void Matrix::T()
{
     long n1=N1, n2=N2;
     if (n1!=n2)
     {
	  long N=::Max(n1,n2);
	  Resize(N);
     }
     for (long i=1;i<=N1;i++)
	  for (long j=i+1;j<=N2;j++)
	       Swap(Elem(i,j),Elem(j,i));
     if (n1!=n2)
	  Part(1,1,n2,n1);
}

void Matrix::Part(long n10, long n20, long n1f, long n2f)
{
     Matrix R;
     ::Part(R,*this,n10,n20,n1f,n2f);
     Transfer(R);
}

void Matrix::Add_Part(const Matrix &M, long n1, long n2)
{
#ifdef DEBUG
     if (n1+M.N1>N1) Error("Incompatible dimensions in Add_Part");
     if (n2+M.N2>N2) Error("Incompatible dimensions in Add_Part");
#endif
     long n=M.N1, ix=1;
     double alpha=1.0;
     for (long i=1;i<=M.N2;i++)
	  daxpy_(&n,&alpha,M.D+M.N1*(i-1)+1,&ix,D+N1*(n1+i-2)+n2,&ix);
}

void Matrix::Set_Part(const Matrix &M, long i, long j)
{
#ifdef DEBUG
     if (i+M.N1-1>N1) Error("Incompatible dimensions in Set_Part.");
     if (j+M.N2-1>N2) Error("Incompatible dimensions in Set_Part.");
#endif
     long m=M.N2;
     for (long k=1;k<=m;k++)
	  memcpy(D+N1*(j+k-2)+i,M.D+M.N1*(k-1)+1,M.N1*sizeof(double));
}

void Matrix::Sqr()
{
     for (long i=1;i<=N1;i++)
	  for (long j=1;j<=N2;j++)
	       Elem(i,j)=::Sqr(Elem(i,j));
}

void Matrix::Rand(double a, double b)
{
#ifdef DEBUG
     if (!Rand_Started) Error("Using Rand without Rand_Open\n");
#endif
     for (long i=1;i<=N1;i++)
	  for (long j=1;j<=N2;j++)
	       Elem(i,j)=::Rand(a,b);
}

void Matrix::Rand_Gaussian(double a, double b)
{
#ifdef DEBUG
     if (!Rand_Started) Error("Using Rand without Rand_Open\n");
#endif
     for (long i=1;i<=N1;i++)
	  for (long j=1;j<=N2;j++)
	       Elem(i,j)=::Rand_Gaussian(a,b);
}


/////////////////////////////////////////////////////
// Matrix Transformations: Linear Algebra
/////////////////////////////////////////////////////

void Matrix::Change_Basis(const Matrix &B)
{
     *this=::Change_Basis(*this,B);
}

bool Matrix::Gram_Schmidt()
{
     Vector V, W;
     double dotprod;
     
     for (long k=1;k<=N2;k++)
     {
	  Col(V,k);
	  for (long j=1;j<=k-1;j++)
	  {
	       Col(W,j);
	       dotprod=-Dot(V,W);
	       Daxpy(V,W,dotprod);
	  }
	  bool success=V.Normalize();
	  if (!success) 
	  {
	       Error_Flag(Error_Mat);
	       return false;
	  }
	  Set_Col(V,k);
     }
     return true;
}

bool Matrix::Invert()
{
#ifdef DEBUG
     if (N1!=N2) Error("Can't invert non-square matrix!");
#endif
     Matrix R(N1);
     R.Unit();
     bool error=Solve(R);
     Transfer(R);
     return error;
}

// what is the return value in this function?
bool Matrix::LU_Decomp(int *I) 
{
     long M=N1;
     long N=N2;
     double *A=D+1;
     long lda=M;
     long info=0;
     dgetrf_(&M,&N,A,&lda,I,&info);
     if (info)
     {
	  // printf("Info: %ld\n",info);
	  Error_Flag(Error_Mat);
	  return false;
     }
     return true;
}

//////////////////////////////////////////////////////////////
// Matrix Functions: Element Handling
//////////////////////////////////////////////////////////////

bool Matrix::Is_Zero(double tolerance) const
{
     for (long i=1;i<=N1*N2;i++)
	  if (fabs(D[i])>tolerance) return false;
     return true;
}

double Matrix::Elem(long i, long j) const
{
#ifdef DEBUG
     if (i<1 || i>N1 || j<1 || j>N2) 
	  Error("Error accessing matrix element.");
#endif
     return D[N1*(j-1)+i];
}

Vector Matrix::Col(long n) const
{
     Vector R;
     Col(R,n);
     return R;
}

Vector Matrix::Row(long n) const
{
     Vector R;
     Row(R,n);
     return R;
}

Vector Matrix::Diag() const
{
     long M=::Min(N1,N2);
     Vector V(M);
     for (long i=1;i<=M;i++)
	  V(i)=Elem(i,i);
     return V;
}

void Matrix::Col(Vector &R, long n) const
{
     double *d=(double*)malloc((N1+1)*sizeof(double));
     memcpy(d+1,D+1+(n-1)*N1,N1*sizeof(double));
     R.Load(d,N1);
}

void Matrix::Row(Vector &R, long n) const
{
     R.Create(N2);
     R.Zero();
     for (long i=1;i<=N2;i++)
	  R(i)=Elem(n,i);
}

double Matrix::Max() const
{
     double X=-1e20;
     for (long i=1;i<=N1;i++)
	  for (long j=1;j<=N2;j++)
	  {
	       double x=Elem(i,j);
	       if (x>X) X=x;
	  }
     return X;
}

double Matrix::Max(long &im, long &jm) const
{
     double X=-1e20;
     for (long i=1;i<=N1;i++)
	  for (long j=1;j<=N2;j++)
	  {
	       double x=Elem(i,j);
	       if (x>X) 
	       {
		    X=x;
		    im=i; jm=j;
	       }
	  }
     return X;
}

double Matrix::Min() const
{
     double X=1e20;
     for (long i=1;i<=N1;i++)
	  for (long j=1;j<=N2;j++)
	  {
	       double x=Elem(i,j);
	       if (x<X) X=x;
	  }
     return X;
}

double Matrix::Min(long &im, long &jm) const
{
     double X=1e20;
     for (long i=1;i<=N1;i++)
	  for (long j=1;j<=N2;j++)
	  {
	       double x=Elem(i,j);
	       if (x<X) 
	       {
		    X=x;
		    im=i; jm=j;
	       }
	  }
     return X;
}

/////////////////////////////////////////////////////////
// Matrix Functions: Linear Algebra
/////////////////////////////////////////////////////////

// Isn't there a BLAS routine to do this?
double Matrix::Elem(const Vector &V1, const Vector &V2) const
{
     return Dot(V1,(*this)*V2);
}

// Take the matrix element between 2 columns of matrix M
double Matrix::Elem(const Matrix &M, long c1, long c2) const
{
#ifdef DEBUG
     if (N1!=N2) Error("Can't take matrix elem of non-square matrix.");
     if (N2!=M.N1) Error("Incompatible dimensions in Elem.");
#endif
     double alpha=1.0, beta=0.0;
     long incx=1, n1=N1, n2=N2;
     char c='N';
     double *d=(double*)malloc(N1*sizeof(double));
     dgemv_(&c, &n1, &n2, &alpha, D+1, &n1, M.D+N1*(c2-1)+1, &incx, 
	    &beta, d, &incx);
     long ix=1, n=M.N1;
     double value=ddot_(&n,M.D+N1*(c1-1)+1,&ix,d,&ix);     
     free(d);
     return value;
}

// OPTIMAL ROUTINE
bool Matrix::Solve(Vector &R) const
{
#ifdef DEBUG
     if (R.N!=N1) Error("Incompatible dimensions in Solve routine.");
#endif
     long N=N1; // number of linear equations
     long M=1; // number of RHS's
     long info=0;
     double *d=(double*)malloc((N1*N2+1)*sizeof(double));
     memcpy(d,D,(N1*N2+1)*sizeof(double));
     long *ip=(long*)malloc(N*sizeof(long));
     dgesv_(&N, &M, d+1, &N, ip, R.D+1, &N, &info);
     free(d); free(ip);
     if (info) 
     {
	  Error_Flag(Error_Mat);
	  return false;
     }
     return true;
}

// OPTIMAL ROUTINE
// Solve Ax=b for all columns of matrix R, put results as columns of R.
bool Matrix::Solve(Matrix &R) const
{
#ifdef DEBUG
     if (N1!=R.N1) Error("Incompatible dimensions in solver.");
#endif
     long N=N1; // number of linear equations
     long M=R.N2; // number of RHS's
     long info=0;
     double *d=(double*)malloc((N1*N2+1)*sizeof(double));
     memcpy(d,D,(N1*N2+1)*sizeof(double));
     long *ip=(long*)malloc(N*sizeof(long));
     dgesv_(&N, &M, d+1, &N, ip, R.D+1, &N, &info);
     free(ip); free(d);
     if (info) 
     {
	  Error_Flag(Error_Mat);
	  return false;
     }
     return true;
}

double Matrix::Det() const
{
#ifdef DEBUG
     if (!N1) return 0.0; // Non-existent matrix, zero det
     if (N1!=N2) Error("Trying to get det of non-square matrix.");
#endif
     int *I=(int*)malloc(N1*sizeof(int));
     Matrix B(*this);
     bool success=B.LU_Decomp(I);
     if (!success)
     {
	  Error_Flag(Error_Mat);
	  return 0.0;
     }
     double det=1.0;
     for (long i=1;i<=N1;i++)
     	  if (i!=I[i-1]) det*=-1.0; // Now it has the right sign
     for (long i=1;i<=N1;i++)
	  det*=B(i,i);
     free(I);
     return det;
}

double Matrix::Trace() const
{
     double result=0.0;
     long N=::Min(N1,N2);
     for (long i=1;i<=N;i++)
	  result+=Elem(i,i);
     return result;
}

// Full diagonalization: all eigenvalues, all eigenvectors
// Symmetric matrix.
// NOT OPTIMAL: matrix is preserved 
void Matrix::Diagonalize(Matrix &B, Vector &E) const
{
#ifdef DEBUG
     if (N1!=N2) Error("Trying to diagonalize non-square matrix.");
#endif
     char jobz='V', range='A', uplo='U';
     long N=N1;
     double *A=(double*)malloc(N*N*sizeof(double));
     memcpy(A,D+1,N*N*sizeof(double));
     long lda=N;
     double vl, vu;
     long il, iu;
     char askmachine='S';
     double abstol=2.0*dlamch_(&askmachine);
     long m; // total number of eigenvalues found
     E.Create(N); double *W=E.D+1;
     B.Create(N); double *Z=B.D+1;
     long ldz=N;
     double *work=(double*)malloc(8*N*sizeof(double));
     long lwork=8*N;
     long *iwork=(long*)malloc(5*N*sizeof(long));
     long *ifail=(long*)malloc(N*sizeof(long));
     long info=0;
     dsyevx_(&jobz,&range,&uplo,&N,A,&lda,&vl,&vu,&il,&iu,&abstol,
	     &m,W,Z,&ldz,work,&lwork,iwork,ifail,&info);
     if (info) Error_Flag(Error_Mat);
     free(iwork);
     free(ifail);
     free(work);
     free(A);
}

void Matrix::Spectrum(Vector &E) const
{
#ifdef DEBUG
     if (N1!=N2) Error("Trying to diagonalize non-square matrix.");
#endif
     char jobz='N', range='A', uplo='U';
     long N=N1;
     double *A=(double*)malloc(N*N*sizeof(double));
     memcpy(A,D+1,N*N*sizeof(double));
     long lda=N;
     double vl, vu;
     long il, iu;
     char askmachine='S';
     double abstol=2.0*dlamch_(&askmachine);
     long m; // total number of eigenvalues found
     E.Create(N); double *W=E.D+1;
     //B.Create(N); double *Z=B.D+1;
     double *Z=NULL;
     long ldz=N;
     double *work=(double*)malloc(8*N*sizeof(double));
     long lwork=8*N;
     long *iwork=(long*)malloc(5*N*sizeof(long));
     long *ifail=(long*)malloc(N*sizeof(long));
     long info=0;
     dsyevx_(&jobz,&range,&uplo,&N,A,&lda,&vl,&vu,&il,&iu,&abstol,
	     &m,W,Z,&ldz,work,&lwork,iwork,ifail,&info);
     free(A); free(work); free(iwork); free(ifail);
     if (info) Error_Flag(Error_Mat);
}

void Matrix::Tridiagonalize(Matrix &B, Vector &Diag, Vector &S) const
{
     // Needs dsytrd y dorgtr
     // Reduce to tri-diagonal form, basis is in "strange form"
     char uplo='U';
     long N=N1;
     double *A=(double*)malloc((N*N+1)*sizeof(double));
     memcpy(A+1,D+1,N*N*sizeof(double));
     long lda=N;
     Diag.Create(N); double* d=Diag.D+1;
     S.Create(N-1); double* e=S.D+1;
     double* tau=(double*)malloc(N*sizeof(double));
     long lwork=N*N;
     double* work=(double*)malloc(lwork*sizeof(double));
     long info=0;
     dsytrd_(&uplo, &N, A+1, &lda, d, e, tau, work, &lwork, &info);
     if (info) Error("Error during tridiagonalization, stage 1.");
     // Compute the basis in "normal form" for tri-diagonal reduction
     dorgtr_(&uplo, &N, A+1, &lda, tau, work, &lwork, &info);
     B.Create(N);
     B.Load(A,N,N);
     free(work); free(tau);
}

void Matrix::NS_Diagonalize(Matrix &BL, Matrix &BR, 
			    Vector &ER, Vector &EI) const
{
#ifdef DEBUG
     if (N1!=N2) Error("Can't Ns_Diagonalize non-square matrix.");
#endif
     char jobl='V';
     char jobr='V';
     long N=N1;
     double *A=(double*)malloc(N*N*sizeof(double));
     memcpy(A,D+1,N*N*sizeof(double));
     long lda=N;
     ER.Create(N); EI.Create(N);
     BL.Create(N); BR.Create(N);
     double* wr=ER.D+1;
     double* wi=EI.D+1;
     double* vl=BL.D+1;
     long ldvl=N;
     double* vr=BR.D+1;
     long ldvr=N;
     long lwork=6*N;
     double* work=(double*)malloc(lwork*sizeof(double));
     long info=0;
     dgeev_(&jobl,&jobr,&N,A,&lda,
	    wr, wi, vl, &ldvl, vr, &ldvr, 
	    work, &lwork, &info);
     if (info) Error_Flag(Error_Mat);
     free(A); free(work); 
}

/////////////////////////////////////////////////
// Matrix I/O
/////////////////////////////////////////////////

void Matrix::Write(int prec) const
{
     char form[40]="%12.8g ";
     if (prec) sprintf(form,"%%%d.%dg ",prec+6,prec);
     for (long i=1;i<=N1;i++)
     {
         for (long j=1;j<=N2;j++)
	 {
	      double x=Elem(i,j);
	      printf(form,x);
	 }
         printf("\n");
     }
     printf("\n");
}

bool Matrix::Save_Binary(FILE *fich) const
{
     int nwrite=fwrite(&N1,sizeof(long),1,fich);
     if (nwrite!=1) { Error_Flag(Error_IO); return false; }
     nwrite=fwrite(&N2,sizeof(long),1,fich);
     if (nwrite!=1) { Error_Flag(Error_IO); return false; }
     if (!(N1*N2)) return true;
     nwrite=fwrite(D+1,sizeof(double),N1*N2,fich);
     if (nwrite!=N1*N2) { Error_Flag(Error_IO); return false;}
     if (ferror(fich)) return false;
     return true;
}

bool Matrix::Save_Binary(const char *s) const
{
     FILE *fich=fopen(s,"wb");
     if (!fich) return false;
     bool status=Save_Binary(fich);
     fclose(fich);
     return status;
}

bool Matrix::Load_Binary(FILE *fich)
{
     int nread=fread(&N1,sizeof(long),1,fich);
     if (nread!=1) 
     { 
	  Error_Flag(Error_IO); return false; 
     }
     nread=fread(&N2,sizeof(long),1,fich);
     if (nread!=1) 
     { 
	  Error_Flag(Error_IO); return false; 
     }
     if (!(N1*N2)) 
     { 
	  D=(double*)NULL; 
	  return true;
     }
     Create(N1,N2);
     Zero();
     nread=fread(D+1,sizeof(double),N1*N2,fich);
     if ((nread!=N1*N2) || (ferror(fich)))
     { 
	  Error_Flag(Error_IO); 
	  return false; 
     }
     return true;
}

bool Matrix::Load_Binary(const char *s) 
{
     FILE *fich=fopen(s,"rb");
     if (!fich) return false;
     bool status=Load_Binary(fich);
     fclose(fich);
     return status;
}

bool Matrix::Save(const char *name) const
{
     FILE *fich=fopen(name,"wt");
     if (!fich) return false;
     Save(fich);
     fclose(fich);
     return true;
}

bool Matrix::Save(const char *name, const char *comment) const
{
     FILE *fich=fopen(name,"wt");
     if (!fich) return false;
     Save(fich,comment);
     fclose(fich);
     return true;
}

bool Matrix::Load(const char *name)
{
     FILE *fich=fopen(name,"rt");
     if (!fich) return false;
     bool status=Load(fich);
     fclose(fich);
     return status;
}

bool Matrix::Save(FILE *fich) const
{
     if (fprintf(fich,"# %ld %ld\n",N1,N2)<0)
     {
	  Error_Flag(Error_IO); return false; 
     }
     for (long i=1;i<=N1;i++)
     {
	  for (long j=1;j<=N2;j++)
	       fprintf(fich,"%20.16g ",Elem(i,j));
	  fprintf(fich,"\n");
     }
     if (ferror(fich)) 
     {
	  Error_Flag(Error_IO); return false;
     }
     return true;
}

bool Matrix::Save(FILE *fich, const char *comment) const
{
     if (fprintf(fich,"# %ld %ld\n",N1,N2)<0)
     {
	  Error_Flag(Error_IO); return false; 
     }
     fprintf(fich,"# %s\n",comment);
     for (long i=1;i<=N1;i++)
     {
	  for (long j=1;j<=N2;j++)
	       fprintf(fich,"%20.16g ",Elem(i,j));
	  fprintf(fich,"\n");
     }
     if (ferror(fich)) 
     {
	  Error_Flag(Error_IO); return false;
     }
     return true;
}

bool Matrix::Load(FILE *fich)
{
     Text Z;
     long nchar=Z.Get_Line(fich);
     if (nchar<=0) { Error_Flag(Error_IO); return false; }
     if (!Z.Is_Here("#")) { Error_Flag(Error_IO); return false; }
     long n1=Z.Get_Int(2), n2=Z.Get_Int(3);
     Create(n1,n2);
     Zero();
     while(Z.Is_Here("#"))
     {
	  nchar=Z.Get_Line(fich);
	  if (nchar<=0) { Error_Flag(Error_IO); return false; }
     }
     for (long i=1;i<N1;i++)
     {
	  Vector V=Z.To_Vector();
	  Set_Row(V,i);
	  nchar=Z.Get_Line(fich);
	  if (nchar<=0) { Error_Flag(Error_IO); return false; }
     }
     Vector V=Z.To_Vector();
     Set_Row(V,N1);
     return true;
}

/////////////////////////////////////////////////////////////////
// Matrix: overloaded operators
/////////////////////////////////////////////////////////////////

double& Matrix::operator()(long i, long j)
{
#ifdef DEBUG
     if (i<1 || i>N1 || j<1 || j>N2) 
	  Error("Error accessing matrix element.");
#endif
     return D[N1*(j-1)+i];
}

double Matrix::operator() (long i, long j) const
{
#ifdef DEBUG
     if (i<1 || i>N1 || j<1 || j>N2) 
	  Error("Error accessing matrix element.");
#endif
     return D[N1*(j-1)+i];
}

Matrix& Matrix::operator=(const Matrix& M)
{
     if (this==&M) return *this;
     Copy(*this,M);
     return (*this);
}

void Matrix::operator+=(const Matrix& M)
{
     if (!N1) { Create(M.N1,M.N2); Zero(); }
#ifdef DEBUG
     if (N1!=M.N1 || N2!=M.N2) 
	  Error("Adding matrices with different dims.");
#endif
     double alpha=1.0;
     long ix=1, n=N1*N2;
     daxpy_(&n,&alpha,M.D+1,&ix,D+1,&ix);
}

void Matrix::operator-=(const Matrix& M)
{
     if (!N1) { Create(M.N1,M.N2); Zero(); }
#ifdef DEBUG
     if (N1!=M.N1 || N2!=M.N2) 
	  Error("Adding matrices with different dims.");
#endif
     double alpha=-1.0;
     long ix=1, n=N1*N2;
     daxpy_(&n,&alpha,M.D+1,&ix,D+1,&ix);
}

void Matrix::operator*=(const Matrix& M)
{
#ifdef DEBUG
     if (N2!=M.N1)
	  Error("Incompatible dimensions in *=");
     if (!N1 || !N2 || !M.N1 || !M.N2)
	  Error("Zero matrix in *=\n");
#endif
     double alpha=1.0, beta=0.0;
     char ca='N', cb='N';
     long n1=N1, n2=M.N2, k=M.N1;
     double *D2=(double*)malloc((n1*n2+1)*sizeof(double));
     dgemm_(&ca,&cb,&n1,&n2,&k,&alpha,D+1,&n1,M.D+1,&k,&beta,D2+1,&n1);
     Load(D2,n1,n2);     
}

// Add x to all entries in the matrix
void Matrix::operator+=(double x)  
{
#ifdef DEBUG
    if (!N1) return;
#endif
    for (long i=1;i<=N1;i++)
	 for (long j=1;j<=N2;j++)
	      Elem(i,j)+=x;
}

void Matrix::operator-=(double x)  
{
#ifdef DEBUG
    if (!N1) return;
#endif
    for (long i=1;i<=N1;i++)
	 for (long j=1;j<=N2;j++)
	      Elem(i,j)-=x;
}

void Matrix::operator*=(double x)  
{
#ifdef DEBUG
    if (!N1) return;
#endif
    long ix=1, n=N1*N2;
    double cosa=x;
    dscal_(&n,&cosa,D+1,&ix);
}

void Matrix::operator/=(double x)  
{
#ifdef DEBUG
    if (!N1) return;
#endif
    long ix=1, n=N1*N2;
    double cosa=1./x;
    dscal_(&n,&cosa,D+1,&ix);
}

void Matrix::operator&=(const Vector &L)
{
     Append_Row(L);
}

void Matrix::operator&=(const Matrix &T)
{
     Append_Row(T);
}

void Matrix::operator|=(const Vector &L)
{
     Append_Col(L);
}

void Matrix::operator|=(const Matrix &T)
{
     Append_Col(T);
}


/////////////////////////////////////////////////////////////////
// Matrix: External Functions
/////////////////////////////////////////////////////////////////

void Copy(Matrix& B, const Matrix& A) // B <- A raw and strict copy. 
{
     if (!A.N1) return;
     B.Create(A.N1,A.N2);
     memcpy(B.D+1,A.D+1,A.N1*A.N2*sizeof(double));
}

Matrix Zero(long N1, long N2)
{
     Matrix Z(N1,N2);
     return Z;
}

Matrix Unit(long N1, long N2)
{
     Matrix U(N1,N2);
     U.Unit();
     return U;
}

Matrix Diag(const Vector &E)
{
     long N=E.N;
     Matrix R(N);
     for (long i=1;i<=N;i++)
	  R(i,i)=E(i);
     return R;
}

Matrix Constant(double x, long n1, long n2)
{
     Matrix A(n1,n2);
     A.Set(x);
     return A;		  
}

void Write(const Matrix &M)
{
     M.Write();
}

Matrix operator-(const Matrix &M)
{
     Matrix R(M);
     R*=-1.0;
     return R;
}

Matrix operator+(const Matrix &A, const Matrix &B)
{
#ifdef DEBUG
     if (A.N1!=B.N1 || A.N2!=B.N2) 
	  Error("Adding matrices with different dimensions.");
#endif
     Matrix R(A.N1,A.N2);
     for (long i=1;i<=A.N1;i++)
	  for (long j=1;j<=A.N2;j++)
	       R(i,j)=A(i,j)+B(i,j);
     return R;
}

Matrix operator-(const Matrix &A, const Matrix &B)
{
     Matrix R(A);
     R-=B;
     return R;
}

Matrix operator*(double K, const Matrix &A)
{
     Matrix R(A);
     R*=K;
     return R;
}

Matrix operator*(const Matrix &A, double K)
{
     Matrix R(A);
     R*=K;
     return R;
}

Matrix operator/(const Matrix &A, double K)
{
     Matrix R(A);
     R/=K;
     return R;
}

Vector operator*(const Matrix &A, const Vector &V)
{
#ifdef DEBUG
     if (A.N2!=V.N) Error("Wrong dims in Matrix-Vector product.");
#endif
     Vector R;
     Multiply(R,A,V);
     return R;
}

Matrix operator*(const Matrix &A, const Matrix &B)
{
#ifdef DEBUG
     if (A.N2!=B.N1) Error("Wrong dims in Matrix-Matrix product.");
#endif
     Matrix R;
     Multiply(R,A,B);
     return R;
}

Matrix operator|(const Vector &A, const Vector &B)
{
     Matrix R;
     R.Append_Col(A);
     R.Append_Col(B);
     return R;
}

Matrix operator|(const Matrix &M, const Vector &B)
{
     Matrix R(M);
     R.Append_Col(B);
     return R;
}

Matrix operator|(const Vector &B, const Matrix &M)
{
     Matrix R(M);
     R.Insert_Col(B,1);
     return R;
}

Matrix operator|(const Matrix &M1, const Matrix &M2)
{
     Matrix R(M1);
     R.Append_Col(M2);
     return R;
}

Matrix operator&(const Matrix &M, const Vector &B)
{
     Matrix R(M);
     R.Append_Row(B);
     return R;
}

Matrix operator&(const Vector &B, const Matrix &M)
{
     Matrix R(M);
     R.Append_Row(B);
     return R;
}

Matrix operator&(const Matrix &M1, const Matrix &M2)
{
     Matrix R(M1);
     R.Append_Row(M2);
     return R;
}

// Optimal routine
void Multiply(Vector &R, const Matrix &M, const Vector &V)
{
#ifdef DEBUG
     if (M.N2!=V.N) Error("Wrong dimensions in Multiply\n");
#endif
     double alpha=1.0, beta=0.0;
     long incx=1, n1=M.N1, n2=M.N2;
     char c='N';
     R.Create(n1);
     R.Zero();
     dgemv_(&c, &n1, &n2, &alpha, M.D+1, &n1, V.D+1, &incx, 
	    &beta, R.D+1, &incx);
}

// Optimal routine
void Multiply(Matrix &R, const Matrix &M1, const Matrix &M2)
{

#ifdef DEBUG
     if (M1.N2!=M2.N1) Error("Incompatible dimensions in Multiply\n");
#endif
     if (M1.N1==0 || M1.N2==0 || M2.N1==0 || M2.N2==0)
     {
	  R.Destroy();
	  return;
     }
     double alpha=1.0, beta=0.0;
     char ca='N', cb='N';
     long n1=M1.N1, n2=M2.N2, k=M2.N1;
     double *D2=(double*)malloc((n1*n2+1)*sizeof(double));
     dgemm_(&ca,&cb,&n1,&n2,&k,&alpha,M1.D+1,&n1,M2.D+1,&k,&beta,D2+1,&n1);
     R.Load(D2,n1,n2);     
}

// The most general Matrix-Matrix product routine
// R <- alpha*M1*M2 + beta*R
// M1 is transposed if T1=true
// M2 is transposed if T2=true
void Multiply_Add(Matrix &R, const Matrix &M1, const Matrix &M2,
		  double alpha, double beta, bool T1, bool T2)
{
     long n1=(T1 ? M1.N2 : M1.N1), 
	  n2=(T2 ? M2.N1 : M2.N2), 
	  k=(T1 ? M1.N1 : M1.N2);
#ifdef DEBUG
     long kk=(T2 ? M2.N2 : M2.N1);
     if (k!=kk) Error("Wrong dimensions in Multiply_Add");
     if (M1.N1==0 || M1.N2==0 || M2.N1==0 || M2.N2==0)
	  Error("Zero matrix in Multiply_Add\n");
#endif
     char ca=(T1==true ? 'T' : 'N');
     char cb=(T2==true ? 'T' : 'N');
     long N1=M1.N1;
     long N2=M2.N1;
     dgemm_(&ca,&cb,&n1,&n2,&k,&alpha,M1.D+1,&N1,M2.D+1,&N2,&beta,R.D+1,&n1);
}

Matrix Tens_Prod(const Matrix &A, const Matrix &B)
{
     Matrix R;
     Tens_Prod(R,A,B);
     return R;
}

void Tens_Prod(Matrix &R, const Matrix &A, const Matrix &B)
{
     long nA1=A.N1, nA2=A.N2, nB1=B.N1, nB2=B.N2;
     R.Create(nA1*nB1,nA2*nB2);
     R.Zero();
     for (long i1=1;i1<=nA1;i1++)
	for (long i2=1;i2<=nA2;i2++)
	   for (long j1=1;j1<=nB1;j1++)
	      for (long j2=1;j2<=nB2;j2++)
		 R( (i1-1)*nB1 + j1, (i2-1)*nB2 + j2 )=
		      A(i1,i2)*B(j1,j2);
}

Matrix Tens_Prod_Unit(const Matrix &A, long m, Side s)
{
     Matrix R;
     Tens_Prod_Unit(R,A,m,s);
     return R;
}

void Tens_Prod_Unit(Matrix &R, const Matrix &A, long m, Side s)
{
     long n1=A.N1, n2=A.N2;
     R.Create(n1*m,n2*m);
     R.Zero();
     if (s==Right)
     {
	  for (long i1=1;i1<=n1;i1++)
	       for (long i2=1;i2<=n2;i2++)
		    for (long j=1;j<=m;j++)
			 R((i1-1)*m+j,(i2-1)*m+j)=A(i1,i2);
     }
     else
     {
	  for (long i=1;i<=m;i++)
	       for (long j1=1;j1<=n1;j1++)
		    for (long j2=1;j2<=n2;j2++)
			 R((i-1)*n1+j1,(i-1)*n2+j2)=A(j1,j2);
     }
     return;
}

Matrix Tens_Prod_Diag(const Matrix &A, const Vector &D, Side s)
{
     Matrix R;
     Tens_Prod_Diag(R,A,D,s);
     return R;
}

void Tens_Prod_Diag(Matrix &R, const Matrix &A, const Vector &D, Side s)
{
     long n1=A.N1, n2=A.N2, m=D.N;
     R.Create(n1*m,n2*m);
     R.Zero();
     if (s==Right)
     {
	  for (long i1=1;i1<=n1;i1++)
	       for (long i2=1;i2<=n2;i2++)
		    for (long j=1;j<=m;j++)
			 R((i1-1)*m+j,(i2-1)*m+j)=A(i1,i2)*D(j);
     }
     else
     {
	  for (long i=1;i<=m;i++)
	       for (long j1=1;j1<=n1;j1++)
		    for (long j2=1;j2<=n2;j2++)
			 R((i-1)*n1+j1,(i-1)*n2+j2)=D(i)*A(j1,j2);
     }
     return;
}

void Elem_Mult(Matrix &R, const Matrix &A, const Matrix &B)
{
#ifdef DEBUG
     if (A.N1!=B.N1 || A.N2!=B.N2)
	  Error("Different dimensions in Elem_Mult\n");
#endif
     long N1=A.N1, N2=A.N2;
     if (R.N1!=N1 || R.N2!=N2)
	  R.Create(N1,N2);
     for (long i=1;i<=N1;i++)
	  for (long j=1;j<=N2;j++)
	       R(i,j)=A(i,j)*B(i,j);
}

Matrix Elem_Mult(const Matrix &A, const Matrix &B)
{
     Matrix R(A.N1,A.N2);
     Elem_Mult(R,A,B);
     return R;
}

void Ket_Bra(Matrix &R, const Vector &V, const Vector &W)
{
     long N1=V.N, N2=W.N;
     R.Create(N1,N2);
     for (long i=1;i<=N1;i++)
	  for (long j=1;j<=N2;j++)
	       R(i,j)=V(i)*W(j);
}

Matrix Ket_Bra(const Vector &V, const Vector &W)
{
     long N1=V.N, N2=W.N;
     Matrix R(N1,N2);
     Ket_Bra(R,V,W);
     return R;
}

Matrix Projector(const Vector &V)
{
     return Ket_Bra(V,V);
}

void Sort_Cols(Vector &V, Matrix &M) 
{
     Matrix A(M);
     A.Insert_Row(V,1);
     A.Sort_Cols();
     Copy(V,A.Row(1));
     Copy(M,Part(A,2,1,A.N1+1,A.N2));
}

// Orthogonal basis change
Matrix Change_Basis(const Matrix &M, const Matrix &B)
{
     long ncolB=B.N2; // number of columns in B
     Matrix R(ncolB);
     for (long i=1;i<=ncolB;i++)
	  for (long j=1;j<=ncolB;j++)
	       R(i,j)=M.Elem(B,i,j);
     return R;
}

Matrix T(const Matrix &M)
{
     Matrix Res(M);
     Res.T();
     return Res;
}

Matrix Part(const Matrix &M, long i0,long j0,long i1,long j1)
{
     Matrix R;
     Part(R,M,i0,j0,i1,j1);
     return R;
}

void Part(Matrix &R, const Matrix &M, long n10, long n20, long n1f, long n2f)
{
#ifdef DEBUG
     if (n1f<n10 || n2f<n20) Error("End before start in Part.\n");
     if (n10<1 || n20<1) Error("Wrong boundaries in Part (begin).\n");
     if (n1f>M.N1 || n2f>M.N2) Error("Wrong boundaries in Part (end).\n");
#endif
     
     long m1=n1f-n10+1;
     long m2=n2f-n20+1;
     R.Create(m1,m2);
     R.Zero();
     for (long i=0;i<m2;i++)
	  memcpy(R.D+m1*i+1,
		 M.D+M.N1*(n20+i-1)+n10,
		 m1*sizeof(double));
}

Matrix Part(const Matrix &R, const List &L1, const List &L2)
{
     Matrix B(L1.N,L2.N);
     for (long i=1;i<=L1.N;i++)
	  for (long j=1;j<=L2.N;j++)
	       B(i,j)=R(L1(i),L2(j));
     return B;
}

Matrix Invert(const Matrix &M)
{
     Matrix R(M);
     R.Invert();
     return R;
}

Vector Solve(const Matrix &M, const Vector &b)
{
     Vector B(b);
     M.Solve(B);
     return B;
}

double Trace(const Matrix &M)
{
     return M.Trace();
}

double Det(const Matrix &M)
{
     return M.Det();
}

double Norm(const Matrix &M)
{
#ifdef DEBUG
     if (!M.N1 || !M.N2) Error("Norm of non-existent matrix\n");
#endif
     double norm=0.0;
     for (long i=1;i<=M.N1;i++)
	  for( long j=1;j<=M.N2;j++)
	       norm+=::Sqr(M(i,j));
     return sqrt(norm);
}

Matrix Sqr(const Matrix &M)
{
     Matrix A(M.N1,M.N2);
     for (long i=1;i<=M.N1;i++)
	  for (long j=1;j<=M.N2;j++)
	       A(i,j)=Sqr(M(i,j));
     return A;
}

Matrix To_Matrix(const Table &T)
{
     Matrix M(T.N1,T.N2);
     for (long i=1;i<=M.N1;i++)
	  for (long j=1;j<=M.N2;j++)
	       M(i,j)=(double)T(i,j);
     return M;
}

Table To_Table(const Matrix &M)
{
     Table T(M.N1,M.N2);
     for (long i=1;i<=T.N1;i++)
	  for (long j=1;j<=T.N2;j++)
	       T(i,j)=(long)round(M(i,j));
     return T;
}

void Trid_Spectrum(Vector &D, Vector &S)
{
     char compz='N';
     long N=D.N;
     double *z=(double*)NULL;
     long ldz=1;
     double *work=(double*)NULL;
     long info;
     dsteqr_(&compz, &N, D.D+1, S.D+1, z, &ldz, work, &info);
     free(work);
}

// Matrix B should be the Unit Matrix
void Trid_Diagonalize(Matrix &B, Vector &D, Vector &S)
{
     char compz='V';
     long N=D.N;
     double *z=B.D+1;
     long ldz=N;
     double *work=(double*)malloc(2*(N-1)*sizeof(double)); 
     long info;
     dsteqr_(&compz,&N,D.D+1,S.D+1,z,&ldz,work,&info);
     free(work);
}
