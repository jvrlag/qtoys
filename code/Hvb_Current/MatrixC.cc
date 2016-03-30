////////////////////////////////////////////////////////
// hvb++ 1.0, CMATRIX
// Copyleft: Javier Rodriguez Laguna
// 080725-121211-150125
#include "MatrixC.h"

///////////////////////////////////////////////////////////////
/// Class VectorC
///////////////////////////////////////////////////////////////

VectorC::VectorC(): N(0), D(NULL) {} 

VectorC::VectorC(long n)
{
      Start();
      Create(n);
      Zero();
}

VectorC::VectorC(const VectorC &V)
{
     Start();
     Create(V.N);
     if (N) memcpy(D,V.D,(N+1)*sizeof(cmplx));
}

VectorC::VectorC(const Vector &V)
{
     Start();
     Create(V.N);
     if (N) for (long i=1;i<=N;i++) D[i]=(cmplx)V(i);
}

VectorC::VectorC(cmplx *data, long n)
{
     Start();
     Create(n);
     memcpy(D,data,(n+1)*sizeof(cmplx));
}

VectorC::~VectorC() { Destroy(); }

void VectorC::Start()
{
     N=0;
     D=(cmplx*)NULL;
}

void VectorC::Create(long n)
{
     Destroy();
     N=n;
     if (!N) { D=NULL; return;}
     D=(cmplx*)malloc((n+1)*sizeof(cmplx));    
     
#ifdef DEBUG 
     if (!D) Error("Error allocating vector."); 
     if (N) Mem_Control(1,0,2*N);
#endif
}

void VectorC::Load(cmplx *d, long n)
{
     Destroy();
     N=n;
     D=d;
#ifdef DEBUG
     if (N) Mem_Control(1,0,2*N);
#endif
}

void VectorC::Load_Copy(cmplx *d1, long n)
{
     cmplx *d2=(cmplx*)malloc((n+1)*sizeof(cmplx));
     memcpy(d2,d1,(n+1)*sizeof(cmplx));
     Load(d2,n);
}

// Acts as if V goes into our VectorC
void VectorC::Transfer(VectorC &V)
{
     Destroy();
     D=V.D;
     N=V.N;
     V.Start();
}

void VectorC::Destroy()
{
     if (N) 
     { 	
	  free(D);
#ifdef DEBUG
	  Mem_Control(-1,0,-2*N);
#endif
	  N=0;
	  D=(cmplx*)NULL;
     }
}

void VectorC::Zero()
{
     if (N) memset(D,0,(N+1)*sizeof(cmplx));
}

void VectorC::Set(cmplx x)
{
     if (N) for (long i=1;i<=N;i++) D[i]=x;
}

void VectorC::Set_Part(const VectorC& V, long n)
{
#ifdef DEBUG
     if (n+V.N-1>N)
	  Error("Set_Part of a cvector which is too big.");
#endif
     memcpy(D+n,V.D+1,V.N*sizeof(cmplx));
}

// return false if normalization was not possible!
bool VectorC::Normalize()
{
     cmplx norm=Norm();
     if (norm==0.0) return false;
     (*this)/=norm;
     return true;
}

void VectorC::Part(long n1, long n2)
{
#ifdef DEBUG
     if (n1<1 || n2>N) 
	  Error("The part can't be larger than the whole!");
#endif
     VectorC R(n2-n1+1);
     memcpy(R.D+1,D+n1,(n2-n1+1)*sizeof(cmplx));
     Transfer(R);
}

void VectorC::Reverse()
{
     for (long i=1;i<=N/2;i++)
	  Swap(D[i],D[N+1-i]);
}

void VectorC::Append(const cmplx x) 
{
     if (!N) Create(1);
     else
     {
	  N++;
	  D=(cmplx*)realloc(D,(N+1)*sizeof(cmplx));
#ifdef DEBUG
	  Mem_Control(0,0,2);
#endif
     }
     D[N]=x;     
}

void VectorC::Append(const VectorC &V)
{
     long nold=N;
     if (!N) Create (V.N);
     else
     {
          N+=V.N;
          D=(cmplx*)realloc(D,(N+1)*sizeof(cmplx));
#ifdef DEBUG
	  Mem_Control(0,0,2*V.N);
#endif
     }
     for (long i=1;i<=V.N;i++)
          D[nold+i]=V(i);
}

// insert x at position i, vector increases size
void VectorC::Insert(const cmplx x, long i)
{
     N++;
     D=(cmplx*)realloc(D,(N+1)*sizeof(cmplx));
#ifdef DEBUG
     Mem_Control(0,0,2);
#endif
     memmove(D+i+1,D+i,(N-i)*sizeof(cmplx));
     D[i]=x;
}

// Insert V at position i, Vector size increases by V.N
void VectorC::Insert(const VectorC &V, long i)
{
     N+=V.N;
     D=(cmplx*)realloc(D,(N+1)*sizeof(cmplx));
#ifdef DEBUG
     Mem_Control(0,0,2*V.N);
#endif
     memmove(D+i+V.N,D+i,(N-i)*sizeof(cmplx));
     memcpy(D+i+1,V.D+1,V.N*sizeof(cmplx));
}

void VectorC::Remove(long i, long j) // Remove from i to j, both included
{
     long n0=i, nf=(j?j:i);
     long n=nf-n0+1;
#ifdef DEBUG
     Mem_Control(0,0,-2*n);
#endif
     if (n0==1 && nf==N) { Destroy(); return; }
     if (nf!=N)
	  memmove(D+n0,D+nf+1,(N-nf)*sizeof(cmplx));
     N-=n;
     D=(cmplx*)realloc(D,(N+1)*sizeof(cmplx));
}

void VectorC::Re()
{
     for (long i=1;i<=N;i++)
	  D[i]=real(D[i]);
}

void VectorC::Im()
{
     for (long i=1;i<=N;i++)
	  D[i]=imag(D[i]);
}

void VectorC::Abs()
{
     for (long i=1;i<=N;i++)
	  D[i]=abs(D[i]);
}

void VectorC::Conj()
{
     for (long i=1;i<=N;i++)
	  D[i]=conj(D[i]);
}

void VectorC::Sqr()
{
     for (long i=1;i<=N;i++)
	  D[i]=norm(D[i]);
}

void VectorC::Write(int prec) const
{    char form[40];
     if (!prec) sprintf(form,"(%%12.8g,%%12.8g) ");
     else sprintf(form,"(%%%d.%dg,%%%d.%dg) ",prec+4,prec,prec+4,prec);
     for (long i=1;i<=N;i++)
	  printf(form,real(D[i]),imag(D[i]));
     printf("\n\n");
}

void VectorC::Write_Col() const
{
     for (long i=1;i<=N;i++)
	  printf(" %16.12g %16.12g\n",real(D[i]),imag(D[i]));
     printf("\n");
}

bool VectorC::Save_Binary(FILE *fich) const
{
     int nread;
     nread=fwrite(&N,sizeof(long),1,fich);
     if (nread!=1) return false;
     if (!N) return true;
     nread=fwrite(D+1,sizeof(cmplx),N,fich);
     if (nread!=N) return false;
     return true;
}

bool VectorC::Save_Binary(const char *name) const
{
     FILE *fich=fopen(name,"wb");
     if (!fich) return false;
     bool status=Save_Binary(fich);
     fclose(fich);
     return status;
}

bool VectorC::Load_Binary(FILE *fich)
{
     int nread;
     nread=fread(&N,sizeof(long),1,fich);
     if (nread!=1) return false;
     if (!N) { D=(cmplx*)NULL; return true; }
     Create(N);
     nread=fread(D+1,sizeof(cmplx),N,fich);
     if (nread!=N) return false;
     return true;
}

bool VectorC::Load_Binary(const char *name)
{
     FILE *fich=fopen(name,"rb");
     if (!fich) return false;
     bool status=Load_Binary(fich);
     fclose(fich);
     return status;
}

bool VectorC::Save(const char *name) const
{
     FILE *fich=fopen(name,"wt");
     if (!fich) return false;
     bool status=Save(fich);
     fclose(fich);
     return status;
}

bool VectorC::Save(FILE *fich) const
{
     if (fprintf(fich,"# %ld\n",N)<0) return false;
     for (long i=1;i<=N;i++)
	  if (fprintf(fich,"(%16.12g,%16.12g)\n",real(D[i]),imag(D[i]))<0)
	       return false;
     return true;
}

bool VectorC::Load(const char *name) 
{
     FILE *fich=fopen(name,"rt");
     if (!fich) return false;
     int status=Load(fich);
     fclose(fich);
     return status;
}

bool VectorC::Load(FILE *fich)
{
     if (!fich) return false;
     long n; double x, y;
     if (!fscanf(fich,"# %ld\n",&n))
	  return false;
     Create(n);
     Zero();
     for (long i=1;i<=n;i++)
     {
	  if (!fscanf(fich,"(%lg,%lg)\n",&x,&y))
	       return false;
	  cmplx z(x,y);
	  D[i]=z;
     }
     return true;
}

bool VectorC::Is_Zero(double tolerance) const
{
     for (long i=1;i<=N;i++)
	  if (abs(D[i])>tolerance) return false;
     return true;
}

double VectorC::Norm() const
{
     long ix=1, n=N;
     cmplx norm=zdotc_(&n,D+1,&ix,D+1,&ix);
     return sqrt(real(norm));
}

cmplx VectorC::operator() (long n) const
{
#ifdef DEBUG
     if (n<0 || n>N) Error("Error getting Cvector comp."); 
#endif
    return D[n];
}

cmplx& VectorC::operator() (long n)
{
#ifdef DEBUG
    if (n<0 || n>N) Error("Error putting vector comp."); 
#endif
    return D[n];
}

VectorC& VectorC::operator=(const VectorC& W)
{
     if (this==&W) return *this;
     Copy(*this,W);
     return (*this);
}

VectorC& VectorC::operator=(const Vector& W)
{
     Copy(*this,W);
     return (*this);
}

void VectorC::operator+=(const VectorC& W)
{
     if (!N) { Create(W.N); Zero(); }
#ifdef DEBUG
      if (N!=W.N) Error("Incompatible sizes in vector +=");
#endif
      cmplx alpha=1.0;
      long ix=1;
      zaxpy_(&N,&alpha,W.D+1,&ix,D+1,&ix);
}

void VectorC::operator-=(const VectorC& W)
{
     if (!N) { Create(W.N); Zero(); }
#ifdef DEBUG
     if (N!=W.N) Error("Incompatible sizes in vector -=");
#endif
     cmplx alpha=-1.0;
     long ix=1;
     zaxpy_(&N,&alpha,W.D+1,&ix,D+1,&ix);
}

void VectorC::operator*=(const cmplx x)
{
    if (!N) return;
    long ix=1;
    cmplx cosa=x;
    zscal_(&N,&cosa,D+1,&ix);
}

void VectorC::operator*=(const double x)
{
     *this*=(cmplx)x;
}

void VectorC::operator/=(const cmplx x)
{
     if (!N) return;
     long ix=1;
     cmplx cosa=1.0/x;
     zscal_(&N,&cosa,D+1,&ix);
}

void VectorC::operator&=(cmplx p)
{
     Append(p);
}

void VectorC::operator&=(const VectorC &L)
{
     Append(L);
}

///////////////////////////////////////////////////////
// VectorC: external functions
///////////////////////////////////////////////////////

void Copy(VectorC& B, const VectorC& A)
{
     B.Destroy(); 
     if (!A.N) return;
     B.Create(A.N);
     memcpy(B.D,A.D,(A.N+1)*sizeof(cmplx));
}

void Copy(VectorC& B, const Vector& A)
{
     B.Destroy(); 
     if (!A.N) return;
     B.Create(A.N);
     for (long i=1;i<=A.N;i++)
	  B(i)=(cmplx)A(i);
}

VectorC operator+(const VectorC &V, const VectorC &W)
{
#ifdef DEBUG
     if (V.N!=W.N) Error("Inconsistent dimensions in vector +\n");
#endif
     VectorC R(V.N);
     for (long i=1;i<=V.N;i++)
	  R(i)=V(i)+W(i);
     return R;
}

VectorC operator+(const VectorC &V, const Vector &W)
{
#ifdef DEBUG
     if (V.N!=W.N) Error("Inconsistent dimensions in vector +\n");
#endif
     VectorC R(V.N);
     for (long i=1;i<=V.N;i++)
	  R(i)=V(i)+(cmplx)W(i);
     return R;
}

VectorC operator-(const VectorC &V, const VectorC &W)
{
#ifdef DEBUG
     if (V.N!=W.N) Error("Inconsistent dimensions in vector +\n");
#endif
     VectorC R(V.N);
     for (long i=1;i<=V.N;i++)
	  R(i)=W(i)-V(i);
     return R;
}

VectorC operator-(const VectorC &V)
{
     VectorC R(V.N);
     return (-1.0)*R;
}

VectorC operator*(cmplx x, const VectorC& V)
{
     VectorC R(V);
     R*=x;
     return R;
}

VectorC operator*(const VectorC& V, cmplx x)
{
     VectorC R(V);
     R*=x;
     return R;
}

VectorC operator*(double x, const VectorC& V)
{
     VectorC R(V);
     R*=x;
     return R;
}

VectorC operator*(const VectorC& V, double x)
{
     VectorC R(V);
     R*=x;
     return R;
}

VectorC operator*(cmplx z, const Vector &V)
{
     VectorC R(V);
     R*=z;
     return R;
}

VectorC operator/(const VectorC& V, cmplx x)
{
     VectorC R(V);
     R/=x;
     return R;
}

VectorC operator&(const VectorC &L1, const VectorC &L2)
{
     VectorC L(L1);
     L.Append(L2);
     return L;
}

VectorC operator&(cmplx p, const VectorC &L2)
{
     VectorC L(1); L(1)=p; 
     L.Append(L2);
     return L;
}

VectorC operator&(const VectorC &L2, cmplx p)
{
     VectorC L(L2);
     L.Append(p);
     return L;
}

cmplx Dot(const VectorC& V1, const VectorC& V2)
{
#ifdef DEBUG
     if (V1.N!=V2.N) Error("Dot product of vectors of diferent dim.");
#endif
     long ix=1, n=V1.N;
     return zdotc_(&n,V1.D+1,&ix,V2.D+1,&ix);
}

VectorC Tens_Prod(const VectorC &V, const VectorC &W)
{
     long N=V.N*W.N;
     VectorC R(N);
     long k=1;
     for (long i=1;i<=V.N;i++)
	  for (long j=1;j<=W.N;j++)
	  {
	       R.D[k]=V.D[i]*W.D[j];
	       k++;
	  }
     return R;
}

void Tens_Prod(VectorC &R, const VectorC &V, const VectorC &W)
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
void Zaxpy(VectorC &R, const VectorC &V, cmplx alpha)
{
     long ix=1, n=R.N;
     zaxpy_(&n,&alpha,V.D+1,&ix,R.D+1,&ix);     
}

VectorC Normalize(const VectorC &V)
{
     VectorC R(V);
     R.Normalize();
     return R;
}

VectorC Part(const VectorC &V, long n1, long n2)
{
     VectorC R(V);
     R.Part(n1,n2);
     return R;
}

VectorC Part(const VectorC &V, const List &L)
{
     VectorC W(L.N);
     for (long i=1;i<=L.N;i++)
	  W(i)=V(L(i));
     return W;
}

VectorC Reverse(const VectorC &V)
{
     VectorC R(V);
     R.Reverse();
     return R;
}

VectorC Insert(const VectorC &V, const VectorC &W, long n)
{
     VectorC R(V);
     R.Insert(W,n);
     return R;
}

Vector Sqr(const VectorC &V) 
{
     Vector R(V.N);
     for (long i=1;i<=V.N;i++)
	  R(i)=norm(V(i));
     return R;
}

Vector Abs(const VectorC &V) 
{
     Vector R(V.N);
     for (long i=1;i<=V.N;i++)
	  R(i)=abs(V(i));
     return R;
}

Vector Re(const VectorC &V) 
{
     Vector R(V.N);
     for (long i=1;i<=V.N;i++)
	  R(i)=real(V(i));
     return R;
}

Vector Im(const VectorC &V) 
{
     Vector R(V.N);
     for (long i=1;i<=V.N;i++)
	  R(i)=imag(V(i));
     return R;
}

VectorC Conj(const VectorC &V)
{
     VectorC R(V);
     R.Conj();
     return R;
}

VectorC Cmplx(const Vector &V)
{
     VectorC R(V.N);
     for (long i=1;i<=V.N;i++)
	  R(i)=(cmplx)V(i);
     return R;
}

VectorC Cmplx(const Vector &V1, const Vector &V2)
{
     VectorC R(V1.N);
     for (long i=1;i<=V1.N;i++)
	  R(i)=(cmplx)V1(i)+M_I*V2(i);
     return R;

}

VectorC Constant(cmplx z, long N)
{
     VectorC R(N);
     for (long i=1;i<=N;i++)
	  R(i)=z;
     return R;
}

///////////////////////////////////////////////////////////////
// MatrixC
///////////////////////////////////////////////////////////////

MatrixC::MatrixC()
{
     Start();
}

MatrixC::MatrixC(long n1, long n2) : N1(n1), N2(n2)
{
     Start();
     Create(n1,n2);
     Zero();
}

MatrixC::MatrixC(const MatrixC & M)
{
     Start();
     Copy(*this,M);
}

MatrixC::MatrixC(const Matrix &M)
{
     Start();
     Copy(*this,M);
}

MatrixC::~MatrixC()
{
     Destroy();
}

void MatrixC::Start()
{
     N1=N2=0;
     D=(cmplx*)NULL;
}

void MatrixC::Create(long n1, long n2) // n2=0
{
     Destroy();
     N1=n1; N2=n2;
     if (!N2) N2=N1;
     if (!N1) return;
     D=(cmplx*)malloc((N1*N2+1)*sizeof(cmplx));
#ifdef DEBUG
     if (!D) Error ("Error allocating matrix.");
     Mem_Control(0,1,2*N1*N2);
#endif
}

// CAUTION: the data must be stored in columns
void MatrixC::Load(cmplx* d, long n1, long n2)
{
     Destroy();
     N1=n1;
     N2=n2; if (!N2) N2=N1;
     D=d;
#ifdef DEBUG
     if (N1*N2) Mem_Control(0,1,2*N1*N2);
#endif
}

// CAUTION: the data must be stored in columns
void MatrixC::Load_Copy(cmplx *d1, long n1, long n2)
{
     if (!n2) n2=n1;
     cmplx *d2=(cmplx*)malloc((n1*n2+1)*sizeof(cmplx));
     memcpy(d2,d1,(n1*n2+1)*sizeof(cmplx));
     Load(d2,n1,n2);
}

void MatrixC::Transfer(MatrixC &M)
{
     Destroy();
     D=M.D;
     N1=M.N1;
     N2=M.N2;
     M.Start();
}

void MatrixC::Destroy()
{
     if (N1) free(D);
     D=(cmplx*)NULL;
#ifdef DEBUG
     if (N1*N2) Mem_Control(0,-1,-2*N1*N2);
#endif
     N1=N2=0;
}

void MatrixC::Resize(long n1, long n2) // m=0
{
     if (!n2) n2=n1;
     if (n1==N1 && n2==N2) return; 
     long nr1=::Min(N1,n1), nr2=::Min(N2,n2);
     long oldN1=N1, oldN2=N2;
     
     cmplx *D2=(cmplx*)malloc((N1*N2+1)*sizeof(cmplx));
     memcpy(D2+1,D+1,N1*N2*sizeof(cmplx));
     
     Create(n1,n2);
     if (n1>oldN1 || n2>oldN2)
	  Zero(); // In case the new matrix is bigger
     for (long i=1;i<=nr2;i++)
	  memcpy(D+N1*(i-1)+1,D2+oldN1*(i-1)+1,nr1*sizeof(cmplx));
     free(D2);
}

void MatrixC::Zero()
{
     memset(D+1,0,N1*N2*sizeof(cmplx));
}

void MatrixC::Unit()
{
     Zero();
     long n=::Min(N1,N2);
     for (long i=1;i<=n;i++)
	  Elem(i,i)=1.0;
}

void MatrixC::Set(cmplx x)
{
     for (long i=1;i<=N1;i++)
	  for (long j=1;j<=N2;j++)
	       Elem(i,j)=x;
}

cmplx& MatrixC::Elem(long i, long j)
{
#ifdef DEBUG
     if (i<1 || i>N1 || j<1 || j>N2) 
	  Error("Error accessing matrix element.");
#endif
     return D[N1*(j-1)+i];
}

void MatrixC::Set_Col(const VectorC &V, long n)
{
#ifdef DEBUG
     if ((V.N!=N1) || (n>N2)) Error("Put_Col is impossible.");
#endif
     memcpy(D+(n-1)*N1+1,V.D+1,V.N*sizeof(cmplx));
}

void MatrixC::Set_Row(const VectorC &V, long n)
{
     for (long i=1;i<=N2;i++)
	  Elem(n,i)=V(i);
}

void MatrixC::Set_Diag(const VectorC &V)
{
#ifdef DEBUG
     if (V.N>::Min(N1,N2)) Error("Error in Set_Diag.");
#endif
     for (long i=1;i<=V.N;i++)
	  Elem(i,i)=V(i);
}

void MatrixC::Append_Col(const VectorC &V)
{
     if (N1)
     	  Resize(N1,N2+1);
     else
     	  Create(V.N,1);
     Set_Col(V,N2);
}

void MatrixC::Append_Col(const MatrixC &T)
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

void MatrixC::Append_Row(const VectorC &V)
{
     if (N2)
     	  Resize(N1+1,N2);
     else
     	  Create(1,V.N);
     Set_Row(V,N1);
}

void MatrixC::Append_Row(const MatrixC &T)
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

void MatrixC::Insert_Col(const VectorC &V, long p)
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
	  memmove(D+1+(i-1)*N1,D+1+(i-2)*N1,N1*sizeof(cmplx));
     Set_Col(V,p);
}

void MatrixC::Insert_Row(const VectorC &V, long p)
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

void MatrixC::Remove_Col(long p)
{
#ifdef DEBUG
     if (p>N2) Error("Remove non-existent Col\n");
#endif
     if (N2==1) { Destroy(); return; }
     for (long i=p+1;i<=N2;i++)
	  memcpy(D+(i-2)*N1+1,D+(i-1)*N1+1,N1*sizeof(cmplx));
     Resize(N1,N2-1);
}

void MatrixC::Remove_Row(long p)
{
#ifdef DEBUG
     if (p>N1) Error("Remove non-existent row\n");
#endif
     if (N1==1) { Destroy(); return; }
     MatrixC R(N1-1,N2);
     for (long j=1;j<=N2;j++)
     {
	  for (long i=1;i<=p-1;i++)
	       R(i,j)=Elem(i,j);
	  for (long i=p+1;i<=N1;i++)
	       R(i-1,j)=Elem(i,j);
     }
     Transfer(R);	  
}

void MatrixC::Swap_Cols(long k1,long k2)
{
     cmplx acum;
     for(long i=1;i<=N1;i++)
     {
	  acum=Elem(i,k1);
	  Elem(i,k1)=Elem(i,k2);
	  Elem(i,k2)=acum;
     }
}

void MatrixC::Swap_Rows(long k1,long k2)
{
    cmplx acum;
    for(long i=1;i<=N2;i++)
    {
	 acum=Elem(k1,i);
	 Elem(k1,i)=Elem(k2,i);
	 Elem(k2,i)=acum;
    }
}

void MatrixC::Permute_Cols(const List &P)
{
     MatrixC A(N1,N2);
     for (long i=1;i<=N2;i++)
	  A.Set_Col(Col(P(i)),i);
     Transfer(A);
}

void MatrixC::Permute_Rows(const List &P)
{
     MatrixC A(N1,N2);
     for (long i=1;i<=N1;i++)
	  A.Set_Row(Row(P(i)),i);
     Transfer(A);
}

// Very dirty trick: read a Matrix with the real and imaginary
// parts, Sort it, read it back as MatrixC
void MatrixC::Sort_Cols(const Vector &V) // CAUTION: V is invariant!
{
     Matrix A;
     A.D=(double*)D;
     A.N1=2*N1;
     A.N2=N2;
     A.Insert_Row(V,1);
     A.Sort_Cols(V);
     A.Remove_Row(1);
     D=(cmplx*)A.D;
     A.N1=A.N2=0;
}

void MatrixC::T()
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

void MatrixC::Herm()
{
     T();
     Conj();
}

void MatrixC::Part(long n10, long n20, long n1f, long n2f)
{
     MatrixC R;
     ::Part(R,*this,n10,n20,n1f,n2f);
     Transfer(R);
}

void MatrixC::Add_Part(const MatrixC &M, long n1, long n2)
{
#ifdef DEBUG
     if (n1+M.N1>N1) Error("Incompatible dimensions in Add.");
     if (n2+M.N2>N2) Error("Incompatible dimensions in Add.");
#endif
     long m=n2-n1+1, n=M.N1, ix=1;
     cmplx alpha=1.0;
     for (long i=1;i<=m;i++)
	  zaxpy_(&n,&alpha,D+N1*(n1+i-1)+1,&ix,M.D+M.N1*(i-1)+1,&ix);
}

void MatrixC::Set_Part(const MatrixC &M, long i, long j)
{
#ifdef DEBUG
     if (i+M.N1-1>N1) Error("Incompatible dimensions in Set_Part.");
     if (j+M.N2-1>N2) Error("Incompatible dimensions in Set_Part.");
#endif
     long m=M.N2;
     for (long k=1;k<=m;k++)
	  memcpy(D+N1*(j+k-2)+i,M.D+M.N1*(k-1)+1,M.N1*sizeof(cmplx));
}

void MatrixC::Change_Basis(const MatrixC &B)
{
     (*this)=::Change_Basis(*this,B);
}

bool MatrixC::Gram_Schmidt()
{
     VectorC V, W;
     cmplx dotprod;
     
     for (long k=1;k<=N2;k++)
     {
	  Col(V,k);
	  for (long j=1;j<=k-1;j++)
	  {
	       Col(W,j);
	       dotprod=-Dot(W,V);
	       Zaxpy(V,W,dotprod);
	  }
	  bool success=V.Normalize();
	  if (!success) return false;
	  Set_Col(V,k);
     }
     return true;
}

bool MatrixC::Invert()
{
#ifdef DEBUG
     if (N1!=N2) Error("Can't invert non-square matrix!");
#endif
     MatrixC R(N1);
     R.Unit();
     bool error=Solve(R);
     Transfer(R);
     return error;
}

// what is the return value in this function?
bool MatrixC::LU_Decomp(int *I) 
{
     long M=N1;
     long N=N2;
     cmplx *A=D+1;
     long lda=M;
     long info=0;
     zgetrf_(&M,&N,A,&lda,I,&info);
     if (info)
     {
	  Error_Flag(Error_Mat);
	  return false;
     }
     return info;
}

void MatrixC::Re()
{
     for (long i=1;i<=N1;i++)
	  for (long j=1;j<=N2;j++)
	       Elem(i,j)=real(Elem(i,j));
}

void MatrixC::Im()
{
     for (long i=1;i<=N1;i++)
	  for (long j=1;j<=N2;j++)
	       Elem(i,j)=imag(Elem(i,j));
}

void MatrixC::Abs()
{
     for (long i=1;i<=N1;i++)
	  for (long j=1;j<=N2;j++)
	       Elem(i,j)=abs(Elem(i,j));
}

void MatrixC::Sqr()
{
     for (long i=1;i<=N1;i++)
	  for (long j=1;j<=N2;j++)
	       Elem(i,j)=norm(Elem(i,j));
}

void MatrixC::Conj()
{
     for (long i=1;i<=N1;i++)
	  for (long j=1;j<=N2;j++)
	       Elem(i,j)=conj(Elem(i,j));
}

bool MatrixC::Is_Zero(double tolerance) const
{
     for (long i=1;i<=N1*N2;i++)
	  if (abs(D[i])>tolerance) return false;
     return true;
}

cmplx MatrixC::Elem(long i, long j) const
{
#ifdef DEBUG
     if (i<1 || i>N1 || j<1 || j>N2) 
	  Error("Error accessing matrix element.");
#endif
     return D[N1*(j-1)+i];
}

VectorC MatrixC::Col(long n) const
{
     VectorC R;
     Col(R,n);
     return R;
}

VectorC MatrixC::Row(long n) const
{
     VectorC R;
     Row(R,n);
     return R;
}

VectorC MatrixC::Diag() const
{
     long M=::Min(N1,N2);
     VectorC V(M);
     for (long i=1;i<=M;i++)
	  V(i)=Elem(i,i);
     return V;
}

void MatrixC::Col(VectorC &R, long n) const
{
     cmplx *d=(cmplx*)malloc((N1+1)*sizeof(cmplx));
     memcpy(d+1,D+1+(n-1)*N1,N1*sizeof(cmplx));
     R.Load(d,N1);
}

void MatrixC::Row(VectorC &R, long n) const
{
     R.Create(N2);
     for (long i=1;i<=N2;i++)
	  R(i)=Elem(n,i);
}

// Isn't there a BLAS routine to do this?
cmplx MatrixC::Elem(const VectorC &V1, const VectorC &V2) const
{
     VectorC W;
     Multiply(W,(*this),V2);
     return Dot(V1,W);
}

// Take the matrix element between 2 columns of matrix M
cmplx MatrixC::Elem(const MatrixC &M, long c1, long c2) const
{
#ifdef DEBUG
     if (N1!=N2) Error("Can't take matrix elem of non-square matrix.");
     if (N2!=M.N1) Error("Incompatible dimensions in Elem.");
#endif
     cmplx alpha=1.0, beta=0.0;
     long incx=1, n1=N1, n2=N2;
     char c='N';
     cmplx *d=(cmplx*)malloc(N1*sizeof(cmplx));
     zgemv_(&c, &n1, &n2, &alpha, D+1, &n1, M.D+N1*(c2-1)+1, &incx, 
	    &beta, d, &incx);
     long ix=1, n=M.N1;
     cmplx value=zdotc_(&n,M.D+N1*(c1-1)+1,&ix,d,&ix);     
     free(d);
     return value;
}

// OPTIMAL ROUTINE!!!
bool MatrixC::Solve(VectorC &R) const
{
#ifdef DEBUG
     if (R.N!=N1) Error("Incompatible dimensions in Solve routine.");
#endif
     long N=N1; // number of linear equations
     long M=1; // number of RHS's
     long info=0;
     cmplx *d=(cmplx*)malloc((N1*N2+1)*sizeof(cmplx));
     memcpy(d,D,(N1*N2+1)*sizeof(cmplx));
     long *ip=(long*)malloc(N*sizeof(long));
     zgesv_(&N, &M, d+1, &N, ip, R.D+1, &N, &info);
     free(d); free(ip);
     if (info) 
     {
	  Error_Flag(Error_Mat);
	  return false;
     }
     return true;
}

// OPTIMAL ROUTINE!!!
bool MatrixC::Solve(MatrixC &R) const
{
#ifdef DEBUG
     if (N1!=R.N1) Error("Incompatible dimensions in solver.");
#endif
     long N=N1; // number of linear equations
     long M=R.N2; // number of RHS's
     long info=0;
     cmplx *d=(cmplx*)malloc((N1*N2+1)*sizeof(cmplx));
     memcpy(d,D,(N1*N2+1)*sizeof(cmplx));
     long *ip=(long*)malloc(N*sizeof(long));
     zgesv_(&N, &M, d+1, &N, ip, R.D+1, &N, &info);
     free(ip); free(d);
     if (info) printf("Singular matrix! Result of Solve is worthless!\n");
     if (info) 
     {
	  Error_Flag(Error_Mat);
	  return false;
     }
     return true;
}

cmplx MatrixC::Det() const
{
#ifdef DEBUG
     if (!N1) return 0.0; // Non-existent matrix, zero det
     if (N1!=N2) Error("Trying to get det of non-square matrix.");
#endif
     int *I=(int*)malloc(N1*sizeof(int));
     MatrixC B(*this);
     bool error=B.LU_Decomp(I);
     if (error)
     {
	  Error_Flag(Error_Mat);
	  return 0.0;
     }
     cmplx det=1.0;
     for (long i=1;i<=N1;i++)
     	  if (i!=I[i-1]) det*=-1.0; // Now it has the right sign
     for (long i=1;i<=N1;i++)
	  det*=B(i,i);
     free(I);
     return det;
}

// so as not to make checks for "squarity", uses minimal
cmplx MatrixC::Trace() const
{
     cmplx result=0.0;
     long N=::Min(N1,N2);
     for (long i=1;i<=N;i++)
	  result+=Elem(i,i);
     return result;
}

// Full diagonalization: all eigenvalues, all eigenvectors
// NOT OPTIMAL: matrix is preserved 
// CAUTION: Hermitian matrix is assumed!!!!!!!!!!!!
void MatrixC::Diagonalize(MatrixC &B, Vector &E) const
{
#ifdef DEBUG
     if (N1!=N2) Error("Trying to diagonalize non-square matrix.");
#endif
     char jobz='V', range='A', uplo='U';
     long N=N1;
     cmplx *A=(cmplx*)malloc(N*N*sizeof(cmplx));
     memcpy(A,D+1,N*N*sizeof(cmplx));
     long lda=N;
     double vl, vu;
     long il, iu;
     char askmachine='S';
     double abstol=2.0*dlamch_(&askmachine); // This is the optimal absolute
     // tolerance to give, dlamch('S') is the minimum number whose inverse
     // does not overflow
     long m; // total number of eigenvalues found
     E.Create(N); double *W=E.D+1;
     B.Create(N); cmplx *Z=B.D+1;
     long ldz=N;
     cmplx *work=(cmplx*)malloc(8*N*sizeof(cmplx));
     long lwork=8*N;
     double *rwork=(double*)malloc(7*N*sizeof(double));
     long *iwork=(long*)malloc(5*N*sizeof(long));
     long *ifail=(long*)malloc(N*sizeof(long));
     long info=0;
     zheevx_(&jobz,&range,&uplo,&N,A,&lda,&vl,&vu,&il,&iu,&abstol,
	     &m,W,Z,&ldz,work,&lwork,rwork,iwork,ifail,&info);
     if (info)
	  Error_Flag(Error_Mat);
     free(iwork);
     free(rwork);
     free(ifail);
     free(work);
     free(A);
}

void MatrixC::Spectrum(Vector &E) const
{
#ifdef DEBUG
     if (N1!=N2) Error("Trying to diagonalize non-square matrix.");
#endif
     char jobz='N', range='A', uplo='U';
     long N=N1;
     cmplx *A=(cmplx*)malloc(N*N*sizeof(cmplx));
     memcpy(A,D+1,N*N*sizeof(cmplx));
     long lda=N;
     double vl, vu;
     long il, iu;
     char askmachine='S';
     double abstol=2.0*dlamch_(&askmachine); // This is the optimal abstol
     long m; // total number of eigenvalues found
     E.Create(N); double *W=E.D+1;
     cmplx *Z=NULL;
     long ldz=N;
     cmplx *work=(cmplx*)malloc(8*N*sizeof(cmplx));
     long lwork=8*N;
     long *iwork=(long*)malloc(5*N*sizeof(long));
     double *rwork=(double*)malloc(7*N*sizeof(double));     
     long *ifail=(long*)malloc(N*sizeof(long));
     long info=0;
     zheevx_(&jobz,&range,&uplo,&N,A,&lda,&vl,&vu,&il,&iu,&abstol,
	     &m,W,Z,&ldz,work,&lwork,rwork,iwork,ifail,&info);
     free(A); free(work); free(iwork); free(rwork); free(ifail);
     if (info) Error_Flag(Error_Mat);
}

void MatrixC::Tridiagonalize(MatrixC &B, VectorC &Diag, VectorC &S) const
{
     // Reduce to tri-diagonal form, basis is in "strange form"
     char uplo='U';
     long N=N1;
     cmplx *A=(cmplx*)malloc((N*N+1)*sizeof(cmplx));
     memcpy(A+1,D+1,N*N*sizeof(cmplx));
     long lda=N;
     Diag.Create(N); cmplx* d=Diag.D+1;
     S.Create(N-1); cmplx* e=S.D+1;
     cmplx* tau=(cmplx*)malloc(N*sizeof(cmplx));
     long lwork=N*N;
     cmplx* work=(cmplx*)malloc(lwork*sizeof(cmplx));
     long info=0;
     zhetrd_(&uplo, &N, A+1, &lda, d, e,
	     tau, work, &lwork, &info);
     if (info) Error_Flag(Error_Mat);
     // Compute the basis in "normal form" for tri-diagonal reduction
     zungtr_(&uplo, &N, A+1, &lda, tau, work, &lwork, &info);
     if (info) Error_Flag(Error_Mat);
     B.Create(N);
     B.Load(A,N,N);
     free(work); free(tau);
}

void MatrixC::NH_Diagonalize(MatrixC &B, VectorC &E) const
{
#ifdef DEBUG
     if (N1!=N2) Error("Can't NH_Diagonalize non-square matrix.");
#endif
     
     char jobl='V';
     char jobr='N';
     long N=N1;
     cmplx *A=(cmplx*)malloc(N*N*sizeof(cmplx));
     memcpy(A,D+1,N*N*sizeof(cmplx));
     long lda=N;
     E.Create(N); 
     B.Create(N); 
     cmplx* w=E.D+1;
     cmplx* vl=B.D+1;
     long ldvl=N;
     cmplx* vr=(cmplx*)NULL;
     long ldvr=1;
     long lwork=6*N;
     cmplx* work=(cmplx*)malloc(lwork*sizeof(cmplx));
     double *rwork=(double*)malloc(2*N*sizeof(double));
     long info=0;
 
     zgeev_(&jobl,&jobr,&N,A,&lda,
	    w, vl, &ldvl, vr, &ldvr, 
	    work, &lwork, rwork, &info);
     free(A); free(work); free(rwork);
     if (info) Error("Error diagonalizing non-symmetric matrix.");
}

void MatrixC::SVD(MatrixC &BU, MatrixC &BVt, Vector &SV) const
{
     char jobu='S', jobvt='S';
     long M=N1, N=N2;
     cmplx *A=(cmplx*)malloc(M*N*sizeof(cmplx));
     memcpy(A,D+1,M*N*sizeof(cmplx));
     long lda=M;
     SV.Create(min(M,N)); double *S=SV.D+1;
     BU.Create(M,M); cmplx *U=BU.D+1;
     BVt.Create(N,N); cmplx *Vt=BVt.D+1;
     BU.Zero(); BVt.Zero();
     long ldu=M; long ldvt=N;
     long lwork=10*max(M,N);
     cmplx *work=(cmplx*)malloc(lwork*sizeof(cmplx));
     double *rwork=(double*)malloc(lwork*sizeof(double));
     long info=0;
     zgesvd_(&jobu,&jobvt,&M,&N,A,&lda,S,U,&ldu,Vt,&ldvt,work,&lwork, 
	     rwork,&info);
     if (info) Error_Flag(Error_Mat);
     free(A);
     free(work);
     free(rwork);
}

void MatrixC::Write(int prec) const
{
     cmplx x;
     char form[40];
     if (!prec)
	  sprintf(form,"(%%12.8g,%%12.8g) ");
     else
	  sprintf(form,"(%%%d.%dg,%%%d.%dg) ",prec+6,prec,prec+6,prec);
     for (long i=1;i<=N1;i++)
     {
         for (long j=1;j<=N2;j++)
	 {
	      x=Elem(i,j);
	      printf(form,real(x),imag(x));
	 }
         printf("\n");
     }
     printf("\n");
}

bool MatrixC::Save_Binary(FILE *fich) const
{
     int nwrite=fwrite(&N1,sizeof(long),1,fich);
     if (nwrite!=1) 
     { 
	  Error_Flag(Error_IO); return false; 
     }
     nwrite=fwrite(&N2,sizeof(long),1,fich);
     if (nwrite!=1) 
     { 
	  Error_Flag(Error_IO); return false; 
     }
     if (!(N1*N2)) return true;
     nwrite=fwrite(D+1,sizeof(cmplx),N1*N2,fich);
     if (nwrite!=N1*N2) 
     { 
	  Error_Flag(Error_IO); return false;
     }
     if (ferror(fich)) return false;
     return true;
}

bool MatrixC::Save_Binary(const char *s) const
{
     FILE *fich=fopen(s,"wb");
     if (!fich) return false;
     bool status=Save_Binary(fich);
     fclose(fich);
     return status;
}

bool MatrixC::Load_Binary(FILE *fich)
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
	  D=(cmplx*)NULL; 
	  return true;
     }
     Create(N1,N2);
     nread=fread(D+1,sizeof(cmplx),N1*N2,fich);
     if ((nread!=N1*N2) || (ferror(fich)))
     { 
	  Error_Flag(Error_IO); 
	  return false; 
     }
     return true;
}

bool MatrixC::Load_Binary(const char *s) 
{
     FILE *fich=fopen(s,"rb");
     if (!fich) return false;
     bool status=Load_Binary(fich);
     fclose(fich);
     return status;
}

bool MatrixC::Save(const char *name) const
{
     FILE *fich=fopen(name,"wt");
     if (!fich) return false;
     Save(fich);
     fclose(fich);
     return true;
}

bool MatrixC::Load(const char *name)
{
     FILE *fich=fopen(name,"rt");
     if (!fich) return false;
     bool status=Load(fich);
     fclose(fich);
     return status;
}

bool MatrixC::Save(FILE *fich) const
{
     if (fprintf(fich,"# %ld %ld\n",N1,N2)<0)
     {
	  Error_Flag(Error_IO); return false; 
     }
     for (long i=1;i<=N1;i++)
     {
	  for (long j=1;j<=N2;j++)
	       fprintf(fich,"%16.12g %16.12g   ",
		       real(Elem(i,j)),imag(Elem(i,j)));
	  fprintf(fich,"\n");
     }
     if (ferror(fich)) 
     {
	  Error_Flag(Error_IO); return false;
     }
     return true;
}

bool MatrixC::Load(FILE *fich)
{
     double x, y; long n1, n2;
     if (!fscanf(fich,"# %ld %ld\n",&n1,&n2)) 
     {
	  Error_Flag(Error_IO); return false;
     }
     Create(n1,n2);
     for (long i=1;i<=n1;i++)
	  for (long j=1;j<=n2;j++)
	  {
	       if (feof(fich) || !fscanf(fich,"%lg %lg",&x,&y)) 
	       {
		    Error_Flag(Error_IO); return false;
	       }
	       Elem(i,j)=x+M_I*y;
	  }
     if (ferror(fich)) 
     {
	  Error_Flag(Error_IO);
	  return false;
     }
     return true;
}

cmplx& MatrixC::operator()(long i, long j)
{
#ifdef DEBUG
     if (i<1 || i>N1 || j<1 || j>N2) 
	  Error("Error accessing matrix element.");
#endif
     return D[N1*(j-1)+i];
}

cmplx MatrixC::operator() (long i, long j) const
{
#ifdef DEBUG
     if (i<1 || i>N1 || j<1 || j>N2) 
	  Error("Error accessing matrix element.");
#endif
     return D[N1*(j-1)+i];
}

MatrixC& MatrixC::operator=(const MatrixC& M)
{
     if (!(M.N1*M.N2)) 
     { 
	  Destroy(); 
	  return(*this); 
     }
     if (this==&M) return (*this);     
     Copy(*this,M);
     return (*this);
}

MatrixC& MatrixC::operator=(const Matrix& M)
{
     if (!(M.N1*M.N2)) 
     { 
	  Destroy(); 
	  return(*this); 
     }
     Copy(*this,M);
     return (*this);
}

void MatrixC::operator+=(const MatrixC& M)
{
     if (!N1) { Create(M.N1,M.N2); Zero(); }
#ifdef DEBUG
     if (N1!=M.N1 || N2!=M.N2) 
	  Error("Adding matrices with different dims.");
#endif
     cmplx alpha=1.0;
     long ix=1, n=N1*N2;
     zaxpy_(&n,&alpha,M.D+1,&ix,D+1,&ix);
}

void MatrixC::operator+=(const Matrix& M)
{
     if (!N1) { Create(M.N1,M.N2); Zero(); }
#ifdef DEBUG
     if (N1!=M.N1 || N2!=M.N2) 
	  Error("Adding matrices with different dims.");
#endif
     long N=N1*N2;
     for (long I=1;I<=N;I++)
	  D[I]+=(cmplx)M.D[I];
}

void MatrixC::operator-=(const MatrixC& M)
{
     if (!N1) { Create(M.N1,M.N2); Zero(); }
#ifdef DEBUG
     if (N1!=M.N1 || N2!=M.N2) 
	  Error("Substracting matrices with different dims.");
#endif
     cmplx alpha=-1.0;
     long ix=1, n=N1*N2;
     zaxpy_(&n,&alpha,M.D+1,&ix,D+1,&ix);
}

void MatrixC::operator-=(const Matrix& M)
{
     if (!N1) { Create(M.N1,M.N2); Zero(); }
#ifdef DEBUG
     if (N1!=M.N1 || N2!=M.N2) 
	  Error("Substracting matrices with different dims.");
#endif
     long N=N1*N2;
     for (long I=1;I<=N;I++)
	  D[I]-=(cmplx)M.D[I];
}

void MatrixC::operator*=(cmplx x)  
{
#ifdef DEBUG
    if (!N1) return;
#endif
    long ix=1, n=N1*N2;
    cmplx cosa=x;
    zscal_(&n,&cosa,D+1,&ix);
}

void MatrixC::operator*=(double x)
{
     (*this)*=(cmplx)x;
}

void MatrixC::operator*=(const MatrixC& M)
{
#ifdef DEBUG
     if (N2!=M.N1)
	  Error("Incompatible dimensions in *=");
     if (!N1 || !N2 || !M.N1 || !M.N2)
	  Error("Zero matrix in *=\n");
#endif
     cmplx alpha=1.0, beta=0.0;
     char ca='N', cb='N';
     long n1=N1, n2=M.N2, k=M.N1;
     cmplx *D2=(cmplx*)malloc((n1*n2+1)*sizeof(cmplx));
     zgemm_(&ca,&cb,&n1,&n2,&k,&alpha,D+1,&n1,M.D+1,&k,&beta,D2+1,&n1);
     Load(D2,n1,n2);     
}

void MatrixC::operator*=(const Matrix& M)
{
     (*this)*=Cmplx(M);
}

void MatrixC::operator&=(const VectorC &L)
{
     Append_Row(L);
}

void MatrixC::operator&=(const MatrixC &T)
{
     Append_Row(T);
}

void MatrixC::operator|=(const VectorC &L)
{
     Append_Col(L);
}

void MatrixC::operator|=(const MatrixC &T)
{
     Append_Col(T);
}

/////////////////////////////////////////////////////
// MatrixC: external functions
/////////////////////////////////////////////////////

void Copy(MatrixC& B, const MatrixC& A) // B <- A raw and strict copy. 
{
     B.Destroy();
     if (!A.N1) return;
     B.Create(A.N1,A.N2);
     memcpy(B.D+1,A.D+1,A.N1*A.N2*sizeof(cmplx));
}

void Copy(MatrixC &B, const Matrix &A)
{
     B.Destroy();
     if (!A.N1) return;
     B.Create(A.N1,A.N2);
     for (long i=1;i<=A.N1;i++)
	  for (long j=1;j<=A.N2;j++)
	       B(i,j)=A(i,j); // sorry, no other way...
}

MatrixC Zero_C(long N1, long N2)
{
     MatrixC Z(N1,N2);
     return Z;
}

MatrixC Unit_C(long N1, long N2)
{
     MatrixC U(N1,N2);
     U.Unit();
     return U;
}

MatrixC Diag(const VectorC &E)
{
     long N=E.N;
     MatrixC R(N);
     for (long i=1;i<=N;i++)
	  R(i,i)=E(i);
     return R;
}

MatrixC Diag_C(const Vector &E)
{
     long N=E.N;
     MatrixC R(N);
     for (long i=1;i<=N;i++)
	  R(i,i)=E(i);
     return R;
}

MatrixC operator-(const MatrixC &M)
{
     MatrixC R(M);
     R*=(cmplx)-1.0;
     return R;
}

MatrixC operator+(const MatrixC &A, const MatrixC &B)
{
     MatrixC R(A);
     R+=B;
     return R;
}

MatrixC operator+(const MatrixC &A, const Matrix &B)
{
     MatrixC R(A);
     R+=B;
     return R;
}

MatrixC operator-(const MatrixC &A, const MatrixC &B)
{
     MatrixC R(A);
     R-=B;
     return R;
}

MatrixC operator*(cmplx K, const MatrixC &A)
{
     MatrixC R(A);
     R*=K;
     return R;
}

MatrixC operator*(double K, const MatrixC &A)
{
     return ((cmplx)K)*A;
}

VectorC operator*(const MatrixC &A, const VectorC &V)
{
     VectorC R;
     Multiply(R,A,V);
     return R;
}

VectorC operator*(const MatrixC &A, const Vector &V)
{
     VectorC R;
     Multiply(R,A,Cmplx(V));
     return R;
}

MatrixC operator*(const MatrixC &A, const MatrixC &B)
{
     MatrixC R;
     Multiply(R,A,B);
     return R;
}

MatrixC operator*(const MatrixC &A, const Matrix &B)
{
     MatrixC R;
     Multiply(R,A,Cmplx(B));
     return R;
}

MatrixC operator|(const VectorC &A, const VectorC &B)
{
     MatrixC R;
     R.Append_Col(A);
     R.Append_Col(B);
     return R;
}

MatrixC operator|(const MatrixC &M, const VectorC &B)
{
     MatrixC R(M);
     R.Append_Col(B);
     return R;
}

MatrixC operator|(const VectorC &B, const MatrixC &M)
{
     MatrixC R(M);
     R.Insert_Col(B,1);
     return R;
}

MatrixC operator|(const MatrixC &M1, const MatrixC &M2)
{
     MatrixC R(M1);
     R.Append_Col(M2);
     return R;
}

MatrixC operator&(const MatrixC &M, const VectorC &B)
{
     MatrixC R(M);
     R.Append_Row(B);
     return R;
}

MatrixC operator&(const VectorC &B, const MatrixC &M)
{
     MatrixC R(M);
     R.Append_Row(B);
     return R;
}

MatrixC operator&(const MatrixC &M1, const MatrixC &M2)
{
     MatrixC R(M1);
     R.Append_Row(M2);
     return R;
}

void Multiply(VectorC &R, const MatrixC &M, const VectorC &V)
{
#ifdef DEBUG
     if (M.N2!=V.N) Error("Wrong dimensions in Multiply\n");
#endif
     cmplx alpha=1.0, beta=0.0;
     long incx=1, n1=M.N1, n2=M.N2;
     char c='N';
     R.Create(n1);
     zgemv_(&c, &n1, &n2, &alpha, M.D+1, &n1, V.D+1, &incx, 
	    &beta, R.D+1, &incx);
}

void Multiply(MatrixC &R, const MatrixC &M1, const MatrixC &M2)
{
#ifdef DEBUG
     if (M1.N2!=M2.N1) 
	  Error("Incompatible dimensions in Multiply\n");
#endif
     cmplx alpha=1.0, beta=0.0;
     if (M1.N1==0 || M1.N2==0 || M2.N1==0 || M2.N2==0)
     {
	  R.Destroy();
	  return;
     }
     char ca='N', cb='N';
     long n1=M1.N1, n2=M2.N2, k=M2.N1;
     cmplx *D2=(cmplx*)malloc((n1*n2+1)*sizeof(cmplx));
     if (n1<0 || n2<0) Error("Nonexistent matrix B!\n");
     zgemm_(&ca,&cb,&n1,&n2,&k,&alpha,M1.D+1,&n1,M2.D+1,&k,&beta,D2+1,&n1);
     R.Load(D2,n1,n2);     
}

// The most general MatrixC-MatrixC product routine
// R <- alpha*M1*M2 + beta*R
// M1 is transposed if T1=true
// M2 is transposed if T2=true
void Multiply_Add(MatrixC &R, const MatrixC &M1, const MatrixC &M2,
		  cmplx alpha, cmplx beta, bool T1, bool T2)
{
     long n1=(T1 ? M1.N2 : M1.N1), 
	  n2=(T2 ? M2.N1 : M2.N2), 
	  k=(T1 ? M1.N1 : M1.N2);
#ifdef DEBUG
     long kk=(T2 ? M2.N2 : M2.N1);
     if (k!=kk) Error("Wrong internal dim in Multiply_Add");
     if (R.N1!=n1 || R.N2!=n2) Error("Wrong external dim in Multiply_Add");
#endif
     if (M1.N1==0 || M1.N2==0 || M2.N1==0 || M2.N2==0)
	  return;
     char ca=(T1==true ? 'C' : 'N');
     char cb=(T2==true ? 'C' : 'N');
     long N1=M1.N1;
     long N2=M2.N1;
     zgemm_(&ca,&cb,&n1,&n2,&k,&alpha,M1.D+1,&N1,M2.D+1,&N2,&beta,R.D+1,&n1);
}

MatrixC Tens_Prod(const MatrixC &A, const MatrixC &B)
{
     MatrixC R;
     Tens_Prod(R,A,B);
     return R;
}

void Tens_Prod(MatrixC &R, const MatrixC &A, const MatrixC &B)
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

MatrixC Tens_Prod_Unit(const MatrixC &A, long m, Side s)
{
     MatrixC R;
     Tens_Prod_Unit(R,A,m,s);
     return R;
}

void Tens_Prod_Unit(MatrixC &R, const MatrixC& A, long m, Side s)
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

MatrixC Tens_Prod_Diag(const MatrixC &A, const VectorC &D, Side s)
{
     MatrixC R;
     Tens_Prod_Diag(R,A,D,s);
     return R;
}

void Tens_Prod_Diag(MatrixC &R, const MatrixC &A, const VectorC &D, Side s)
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

void Ket_Bra(MatrixC &R, const VectorC &V, const VectorC &W)
{
     long N1=V.N, N2=W.N;
     R.Create(N1,N2);
     for (long i=1;i<=N1;i++)
	  for (long j=1;j<=N2;j++)
	       R(i,j)=V(i)*conj(W(j));
}

MatrixC Ket_Bra(const VectorC &V, const VectorC &W)
{
     long N1=V.N, N2=W.N;
     MatrixC R(N1,N2);
     Ket_Bra(R,V,W);
     return R;
}

MatrixC Projector(const VectorC &V)
{
     return Ket_Bra(V,V);
}

// Very dirty trick: read a Matrix with the real and imaginary
// parts, Sort it, read it back as MatrixC
void Sort(Vector &V, MatrixC &M) 
{
     Matrix A;
     A.D=(double*)M.D;
     A.N1=2*M.N1;
     A.N2=M.N2;
     A.Insert_Row(V,1);
     A.Sort_Cols(V);
     V=A.Row(1);
     A.Remove_Row(1);
     M.D=(cmplx*)A.D;
     A.N1=A.N2=0;
}

// Unitary basis change
MatrixC Change_Basis(const MatrixC &M, const MatrixC &B)
{
     long ncolB=B.N2; // number of columns in B
     MatrixC R(ncolB);
     for (long i=1;i<=ncolB;i++)
          for (long j=1;j<=ncolB;j++)
               R(i,j)=M.Elem(B,i,j);
     return R;
}

MatrixC T(const MatrixC &M)
{
     MatrixC R(M);
     R.T();
     return R;
}

MatrixC Herm(const MatrixC &M)
{
     MatrixC R(M);
     R.Herm();
     return R;
}

MatrixC Part(const MatrixC &M, long i0,long j0,long i1,long j1)
{
     MatrixC R;
     Part(R,M,i0,j0,i1,j1);
     return R;
}

void Part(MatrixC &R, const MatrixC &M, long n10, long n20, long n1f, long n2f)
{
     long m1=n1f-n10+1;
     long m2=n2f-n20+1;
     R.Create(m1,m2);
     R.Zero();
     for (long i=1;i<=m2;i++)
          memcpy(R.D+m1*(i-1)+1,M.D+M.N1*(n20+i-2)+n10,m1*sizeof(cmplx));
}

MatrixC Part(const MatrixC &A, const List &L1, const List &L2)
{
     MatrixC B(L1.N,L2.N);
     for (long i=1;i<=L1.N;i++)
	  for (long j=1;j<=L2.N;j++)
	       B(i,j)=A(L1(i),L2(j));
     return B;
}

MatrixC Invert(const MatrixC &M)
{
#ifdef DEBUG
     if (M.N1!=M.N2) Error("Can't invert non-square matrix!");
#endif     
     MatrixC R(M.N1);
     R.Unit();
     M.Solve(R);
     return R;
}

VectorC Solve(const MatrixC &M, const VectorC &b)
{
     VectorC B(b);
     M.Solve(B);
     return B;
}

cmplx Det(const MatrixC &M)
{
     return M.Det();
}

cmplx Trace(const MatrixC &M)
{
     return M.Trace();
}

double Norm(const MatrixC &M)
{
     double norma=0.0;
     for (long i=1;i<=M.N1;i++)
	  for( long j=1;j<=M.N2;j++)
	       norma+=norm(M(i,j));
     return sqrt(norma);
}

Matrix Re(const MatrixC &M)
{
     Matrix A(M.N1,M.N2);
     for (long i=1;i<=M.N1;i++)
	  for (long j=1;j<=M.N2;j++)
	       A(i,j)=real(M(i,j));
     return A;
}

Matrix Im(const MatrixC &M)
{
     Matrix A(M.N1,M.N2);
     for (long i=1;i<=M.N1;i++)
	  for (long j=1;j<=M.N2;j++)
	       A(i,j)=imag(M(i,j));
     return A;
}

Matrix Abs(const MatrixC &M)
{
     Matrix A(M.N1,M.N2);
     for (long i=1;i<=M.N1;i++)
	  for (long j=1;j<=M.N2;j++)
	       A(i,j)=abs(M(i,j));
     return A;
}

Matrix Sqr(const MatrixC &M)
{
     Matrix A(M.N1,M.N2);
     for (long i=1;i<=M.N1;i++)
	  for (long j=1;j<=M.N2;j++)
	       A(i,j)=norm(M(i,j));
     return A;
}

MatrixC Conj(const MatrixC &M)
{
     MatrixC A(M.N1,M.N2);
     for (long i=1;i<=M.N1;i++)
	  for (long j=1;j<=M.N2;j++)
	       A(i,j)=conj(M(i,j));
     return A;
}

MatrixC Cmplx(const Matrix &M)
{
     MatrixC A(M.N1,M.N2);
     for (long i=1;i<=M.N1;i++)

	  for (long j=1;j<=M.N2;j++)
	       A(i,j)=(cmplx)M(i,j);
     return A;
}

MatrixC Cmplx(const Matrix &R, const Matrix &I)
{
     MatrixC A(R.N1,R.N2);
     for (long i=1;i<=R.N1;i++)
	  for (long j=1;j<=R.N2;j++)
	       A(i,j)=(cmplx)R(i,j)+M_I*(cmplx)I(i,j);
     return A;
}
