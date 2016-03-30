////////////////////////////////////////////////////////
// Hvb++
// Copyleft: Javier Rodríguez Laguna
// Common routines for Hvb
// 060404-080725-110930-150125-150510

#include"Common.h"

cmplx M_I(0.0,1.0);

cmplx operator+(long i,cmplx z)
{
     return (cmplx)i+z;
}

cmplx operator+(cmplx z,long i)
{
     return z+(cmplx)i;
}

cmplx operator-(long i,cmplx z)
{
     return (cmplx)i-z;
}

cmplx operator-(cmplx z,long i)
{
     return z-(cmplx)i;
}

double Sqr(double a)
{
     return a*a;
}

long Max(long a, long b)
{
     return (a>b ? a:b);
}

double Max(double a, double b)
{
     return (a>b ? a:b);
}

long Min(long a, long b)
{
     return (a<b ? a:b);
}

double Min(double a, double b)
{
     return (a<b ? a:b);
}

void Swap(int &a, int &b)
{
     int c; c=a; a=b; b=c;
}

void Swap(long &a, long &b)
{
     long c; c=a; a=b; b=c;
}

void Swap(double &a, double &b)
{
     double c; c=a; a=b; b=c;
}

void Swap(cmplx &a, cmplx &b)
{
     cmplx c; c=a; a=b; b=c;
}

// return (-1)^k (minus one power k)
long Mop(long k)
{
     return (k%2 ? -1 : 1);
}

double Sign(double x)
{
     if (x==0.0) return 0.0;
     if (x<0.0) return -1.0;
     else return 1.0;
}

// returns x with the sign of y
double Sign(double x, double y)
{
     return ((y)>=0.0 ? fabs(x) : -fabs(x));
}

// factorial, x!
double Fact(double x)
{
     return exp(lgamma(x+1));
}

long Fact(long n)
{
     if (!n) return 1;
     long F=1;
     for (long i=1;i<=n;i++)
          F*=i; 
     return F;
}

// Combinatorial numbers, the pedestrian way
long Comb(const long n, const long m)
{
     if (!m || n==m) return 1;
     long vari=1, perm=1;
     for (long i=n;i>=n-m+1;i--)
          vari*=i; // variations Vn,m
     for (long i=m;i>1;i--)
          perm*=i; // permutations Pm
     return (vari/perm); // comb = variations / permutations
}

// Combinatorial numbers, using Gamma functions
double Comb(double x, double y)
{
     if (x==0.0 && y==0.0) return 1.0;
     double z=(x>0.0 ? x : -x+y-1.0);
     double comb= exp(lgamma(z+1.0)-lgamma(y+1.0)-lgamma(z-y+1.0));
     if (!isfinite(comb)) printf("comb(%g %g) is not finite!\n",x,y);
     if (x<0.0) comb*=Mop((long)y);
     return comb;
}

////////////////////////////////////////////////////////
// Error handling
////////////////////////////////////////////////////////

void Error(const char *s) 
{
     printf("%s\n",s);
     fflush(stdout);
     exit(1);
}

void Error_Flag(Error_Type e)
{
     Error_Flags.Append( (long)e );
}

void Error_Clean()
{
     Error_Flags.Destroy();
}

Error_Type Error_Read()
{
     long Ne=Error_Flags.N;
     if (!Ne) return No_Error;
     Error_Type e=(Error_Type)Error_Flags(Ne);
     Error_Flags.Part(1,Ne-1);
     return e;
}

List Error_Flags;

////////////////////////////////////////////////////////
// CPU time, Time in seconds since starting machine
double Clock()
{
     struct timespec ts;
     clock_gettime(CLOCK_MONOTONIC, &ts);
     return (double)ts.tv_sec+ts.tv_nsec/1e9;
}

void Delay(double Dt)
{
     double start=Clock();
     do {} while(Clock()-start < Dt);
}

////////////////////////////////////////////////////////
// RANDOM NUMBERS
// We'll use Mersenne Twister 19937 
#define MTRNG_N 624
#define MTRNG_M 397
#define MTRNG_RAND_MAX 0xffffffffUL
static const unsigned long MTRNG_UPPER_MASK = 0x80000000UL;
static const unsigned long MTRNG_LOWER_MASK = 0x7fffffffUL;

Rand_State_Type *Rand_State;
bool Rand_Started=false;

void Rand_Open(unsigned long int s)
{
     Rand_State=(Rand_State_Type*)malloc(sizeof(Rand_State_Type));
     if (s == 0)
	  s = 4357;   /* the default seed is 4357 */
     Rand_State->mt[0]= s & 0xffffffffUL;
     int i;
     for (i=1; i<MTRNG_N; i++)
     {
	 /* See Knuth's "Art of Computer Programming" Vol. 2, 3rd
	    Ed. p.106 for multiplier. */
	  Rand_State->mt[i] =
	       (1812433253UL * 
		(Rand_State->mt[i-1] ^ (Rand_State->mt[i-1] >> 30)) + i);
	  Rand_State->mt[i] &= 0xffffffffUL;
     }
     Rand_State->mti = i;
#ifdef DEBUG
     Rand_Started=true;
#endif
}

void Rand_Close()
{
     free(Rand_State);
#ifdef DEBUG
     Rand_Started=false;
#endif
}

unsigned long Rand_Full()
{
#ifdef DEBUG
     if (!Rand_Started) Error("Using Rand without Rand_Open\n");
#endif
     unsigned long k;
     unsigned long *const mt = Rand_State->mt;
#define MTRNG_MAGIC(y) (((y)&0x1) ? 0x9908b0dfUL : 0)
     if (Rand_State->mti >= MTRNG_N)
     {   // generate N words at one time 
	  int kk;
	  for (kk=0; kk<MTRNG_N-MTRNG_M; kk++)
	  {
	       unsigned long y = (mt[kk] & MTRNG_UPPER_MASK) | 
		    (mt[kk + 1] & MTRNG_LOWER_MASK);
	       mt[kk] = mt[kk + MTRNG_M] ^ (y >> 1) ^ MTRNG_MAGIC(y);
	  }
	  for (; kk<MTRNG_N-1; kk++)
	  {
	       unsigned long y = (mt[kk] & MTRNG_UPPER_MASK) | 
		    (mt[kk + 1] & MTRNG_LOWER_MASK);
	       mt[kk] = mt[kk + (MTRNG_M - MTRNG_N)] ^ (y >> 1) ^ 
		    MTRNG_MAGIC(y);
	  }
	  unsigned long y = (mt[MTRNG_N - 1] & MTRNG_UPPER_MASK) | 
	       (mt[0] & MTRNG_LOWER_MASK);
	  mt[MTRNG_N - 1] = mt[MTRNG_M - 1] ^ (y >> 1) ^ MTRNG_MAGIC(y);
	  Rand_State->mti = 0;
     }
     // Tempering
     k = mt[Rand_State->mti];
     k ^= (k >> 11);
     k ^= (k << 7) & 0x9d2c5680UL;
     k ^= (k << 15) & 0xefc60000UL;
     k ^= (k >> 18);
     Rand_State->mti++;
     return k;
}

double Rand()
{
     return Rand_Full() / 4294967296.0 ;
}

double Rand(double a, double b)
{
     double x=Rand();
     return (b-a)*x+a;
}

// Get a uniform integer deviate between a and b,
// a and b included!
long Rand_I(long a, long b)
{
     long range=b-a+1;
     return (Rand_Full() % range)+a;
}

double Rand_Gaussian(double mean, double dev)
{
     static double extra; // an extra value produced by the method
     static bool extra_ready=false;
     if (extra_ready)
     {
	  extra_ready=false;
	  return mean + dev*extra;
     }
     double u, v, s;
     do
     {
	  u=Rand(-1.0,1.0);
	  v=Rand(-1.0,1.0);
	  s=u*u+v*v;
     }while(s>=1.0 || s==0.0);
     double a=sqrt(-2.0*log(s)/s);
     extra=v*a;
     extra_ready=true;
     return mean + dev*(u*a);
}

////////////////////////////////////////////////////////
// Binary numbers as long

long Pow_2(long num)
{
     return(1<<num);
}

long Log_2(long num)
{
     return Max_Bit(num);
}

bool Bit(long X, long pos)
{
     return (X & Pow_2(pos));
}

long Max_Bit(long X)
{ 
     long i=-1;
     if (!X) return (0);
     do
     {
	  i++;
	  X=X>>1;
     }while(X);
     return i;
}

long Flip_Bit(long X, long pos) 
{
     long power=Pow_2(pos);
     if (X & power) return(X-power);
     else return(X+power);
}

long Set_Bit(long X, long pos, long value)
{
     if (Bit(X,pos)==value) return X;
     else return Flip_Bit(X,pos);
}

long Swap_Bits(long X, long pos1, long pos2)
{
     bool temp=Bit(X,pos1);
     X=Set_Bit(X,pos1,Bit(X,pos2));
     X=Set_Bit(X,pos2,temp);
     return X;
}

long Reverse_Bits(long X, long maxbit)
{
     if (!maxbit) maxbit=Max_Bit(X);
     int x=0;
     for (int i=0;i<=maxbit;i++)
	  if (Bit(X,i)) 
	       x=Set_Bit(x,maxbit-i,1);
     return x;
}

// insert bit "b" in position "j" of number "X"
// e.g.: insert(ABCD,2,0) -> AB0CD
long Insert_Bit(long X, long pos, long value)
{
     // Let X=ABCD (15), pos=2
     long I=1<<pos;   // I=100   (4)
     long I2=I-1;     // I2=11   (3)
     long low=X & I2; // low=CD  (3)
     long x=X-low;    // x=AB00  (12)
     x*=2;            // x=AB000 (24)
     x+=low;          // AB0CD   (27)
     if (value) x+=I; // AB1CD   (31)
     return x;
}

long Count_Ones(long X) // returns the number of ones in the expansion of X
{
     long number=0;
     long C=X;
     do
     {
	  if (C%2) number++;
	  C>>=1;
     }while(C);
     return number;
}

long Next_In_Sector(long X)
// returns the smallest number which is bigger than X and has the same
// number of ones in its expansion. 
{
     long i=0, n_ones=0, Y=X, R=X;
     bool bit=Y%2, nextbit;
     do
     {
	  Y>>=1;
	  nextbit=Y%2;
	  if (bit && !nextbit) break; // if no 01 pattern is found yet
	  if (bit) // if current bit is switched
	  {
	       n_ones++;
	       R-=Pow_2(i);
	  }
	  i++;
	  bit=nextbit;
     }while(Y); 
     X=Swap_Bits(R,i+1,i);
     X+=Pow_2(n_ones)-1;
     return X;	  
}

char* Bin_2_String(long X, long maxbit)
{
     if (!maxbit) maxbit=Max_Bit(X);
     char *s=new char[maxbit+1];
     for (long i=maxbit;i>=0;i--)
	  s[maxbit-i]=(Bit(X,i) ? '1' : '0');
     s[maxbit+1]='\0';
     return s;
}


////////////////////////////////////////////////////////
// Class List
////////////////////////////////////////////////////////

List::List()
{
     N=0; // non-created
     D=NULL;
}

List::List(long n)
{
     N=0;
     Create(n);
     Zero();
}

List::List(long* P, long n)
{
     N=0;
     Create(n);
     memcpy(D+1,P,n*sizeof(long));
}

List::List(const List& L)
{
     N=0;
     Create(L.N);
     if (L.N>0) memcpy(D,L.D,(N+1)*sizeof(long));
     else D=NULL;
}

List::~List()
{
     Destroy();
}

void List::Start() // n=0
{
     D=(long*)NULL;
     N=0;
}

void List::Create(long n)
{
     if (N!=0) Destroy();
     D=(long*)malloc((n+1)*sizeof(long));
     N=n;
}

void List::Load(long *d, long n)
{
     Destroy();
     N=n;
     D=d;
}

void List::Load_Copy(long *d, long n)
{
     long *d2=(long*)malloc((n+1)*sizeof(long));
     memcpy(d2,d,(n+1)*sizeof(long));
     Load(d2,n);
}

void List::Transfer(List & L)
{
     Destroy();
     D=L.D;
     N=L.N;
     L.Start();
}

void List::Destroy()
{
     if (N) free(D);
     D=NULL;
     N=0;
}

void List::Zero()
{
     if (N==0) return;
     memset(D+1,0,N*sizeof(long));
}

void List::Set(long S) 
{
     for (long k=1;k<=N;k++)
	  D[k]=S;
}

void List::Set_Part(const List &L, long i) // put List at site i, overwrite
{
#ifdef DEBUG
     if (i+L.N-1>N)
	  Error("Inserting a vector which is too big.");
#endif
     memcpy(D+i,L.D+1,L.N*sizeof(long));
}

void List::Remove(long i, long j)
{
     long n0=i, nf=(j?j:i);
     long n=nf-n0+1;
     if (n0==1 && nf==N) { Destroy(); return; }
     memmove(D+n0,D+nf+1,(N-nf)*sizeof(long));
     N-=n;
     D=(long*)realloc(D,(N+1)*sizeof(long));
}

long List::Substract(long x) // removes all appearances of "x"
{
     if (!D) return 0;
     List Indices(N);
     long n=0;
     for (long i=1;i<=N;i++)
	  if (D[i]==x) 
	  {
	       Indices(i)=1;
	       n++;
	  }
     if (n==N) 
     {
	  Destroy();
	  return n;
     }
     long *D2=(long*)malloc((N-n+1)*sizeof(long));
     long index=1;
     for (long i=1;i<=N;i++)
	  if (!Indices(i))
	  {
	       D2[index]=D[i];
	       index++;
	  }
     free(D);
     D=D2;
     N-=n;
     return n;     
}

long List::Substract(const List &L) 
{
     if (!D) return 0;
     List Indices(N);
     long n=0;
     for (long i=1;i<=N;i++)
	  if (L.Find(D[i])) 
	  {
	       Indices(i)=1;
	       n++;
	  }
     if (n==N) 
     {
	  Destroy();
	  return n;
     }
     long *D2=(long*)malloc((N-n+1)*sizeof(long));
     long index=1;
     for (long i=1;i<=N;i++)
	  if (!Indices(i))
	  {
	       D2[index]=D[i];
	       index++;
	  }
     free(D);
     D=D2;
     N-=n;
     return n;     
}

void List::Append(long n)
{ 
     D=(long*)realloc(D,(N+2)*sizeof(long));
     D[N+1]=n;
     N++;
}

void List::Append(const List &L)
{
     D=(long*)realloc(D,(N+L.N+1)*sizeof(long));
     memcpy(D+N+1,L.D+1,L.N*sizeof(long));
     N+=L.N;
}

void List::Sort(long tipo) 
// if positive: ascending order; if negative, descending
{
     if (tipo>0)
	  qsort(D+1,N,sizeof(long),Sort_Long_12);
     else qsort(D+1,N,sizeof(long),Sort_Long_21);
}

void List::Uniquify()
{
     long i,j,x,counter=0;
     long *D2, *D3;
     D2=(long*)malloc((N+1)*sizeof(long));
     for (i=1;i<=N;i++)
     {
	  D2[i]=1;
	  for (j=1;j<i;j++)
	       if (D[i]==D[j]) D2[i]=0;
	  if (D2[i]==1) counter++;
     } // now D2[i] contains "1" at element "i" only if "i" is original
     D3=(long*)malloc((counter+1)*sizeof(long));
     x=1;
     for (i=1;i<=N;i++)
	  if (D2[i]==1) 
	  {
	       D3[x]=D[i];
	       x++;
	  }
     free(D);
     free(D2);
     D=D3;
     N=counter;
}

void List::Insert(long x, long p) // Insert "x" at position "p"
{
     D=(long*)realloc(D,(N+2)*sizeof(long));
     long *F=(long*)malloc((N-p+1)*sizeof(long));
     memcpy(F,D+p,(N-p+1)*sizeof(long));
     D[p]=x;
     memcpy(D+p+1,F,(N-p+1)*sizeof(long));
     N++;
     free(F);
}

void List::Insert(const List &K, long p) // Insert List K at position p
{
     // Save the data beyond "p"
     long *F=(long*)malloc((N-p+1)*sizeof(long));
     memcpy(F,D+p,(N-p+1)*sizeof(long));
     // Now, extend D and copy K on D
     long n=K.N;
     D=(long*)realloc(D,(N+n+1)*sizeof(long));
     memcpy(D+p,K.D+1,n*sizeof(long));
     // Now, copy again the data beyond "p"
     memcpy(D+p+n,F,(N-p+1)*sizeof(long));
     N+=n;
     free(F);
}

// Reverse an interval in a list, from i to j, both included
void List::Reverse(long i, long j)
{
     long n1=(i==0 ? 1 : i);
     long n2=(j==0 ? N : j);
     long m=n2-n1+1; // number of items to reverse
     for (long k=1;k<=m/2;k++)
	  Swap( n1+k-1, n2+1-k);
}
 
void List::Swap(long i, long j)
{
     if (i>N || j>N) Error("Swap error!\n");
     long tmp=D[i];
     D[i]=D[j];
     D[j]=tmp;
}

void List::Part(long n0, long n1)
{
     N=n1-n0+1;
     if (N<1)
     {
	  N=0;
	  free(D); D=NULL;
	  return;
     }
     long* D2=(long*)malloc((N+1)*sizeof(long));
     memcpy(D2+1,D+n0,N*sizeof(long));
     free(D);
     D=D2;
}

void List::Rand(long a, long b)
{
     for (long i=1;i<=N;i++)
	  D[i]=Rand_I(a,b);
}

long List::Find(long n, long k) const // finds first "n", starting from pos k
{
     long i=k-1;
     if (N==0) return 0;
     do
     {
	 i++; 
     }while(i<N && D[i]!=n);
     if (i<=N && D[i]==n) return i;
     else return 0;
}

long List::Find_Non(long x, long k) const
{
     long i=k-1;
     do
	  i++;
     while(i<N && D[i]==x);
     if (i<=N && D[i]!=x) return i;
     else return 0;
}

List List::Find_All(long k) const
{
     List L;
     for (long i=1;i<=N;i++)
	  if (D[i]==k) L.Append(i);
     return L;
}

long List::Count(long x) const // counts the number of appearances of x
{
     long c=0;
     for (long i=1;i<=N;i++)
	  if (D[i]==x) c++;
     return c;
}

long List::Min() const
{
     long nmin=D[1];
     for (long k=2;k<=N;k++)
	  if (D[k]<nmin) nmin=D[k];
     return nmin;
}

long List::Max() const
{
     long nmax=D[1];
     for (long k=2;k<=N;k++)
	  if (D[k]>nmax) nmax=D[k];
     return nmax;
}

long List::Min(long &pmin) const
{
     long nmin=D[1];
     pmin=1;
     for (long p=2;p<=N;p++)
 	  if (D[p]<nmin) { nmin=D[p]; pmin=p; }
     return nmin;
}

long List::Max(long &pmax) const
{
      long nmax=D[1];
      pmax=1;
      for (long p=2;p<=N;p++)
 	  if (D[p]>nmax) { nmax=D[p]; pmax=p; }
      return nmax;
}

long List::Sum() const
{
     return Sum(1,N);
}

long List::Sum(long i0, long i1) const
{
     long r=0;
     for (long k=i0;k<=i1;k++)
	  r+=D[k];
     return r;
}

long List::Prod() const
{
     return Prod(1,N);
}

long List::Prod(long i0, long i1) const
{
     long r=1;
     for (long k=i0;k<=i1;k++)
	  r*=D[k];
     return r;
}

void List::Write() const
{
     if (!N) printf("Empty list\n");
     long l=10;
     for (long i=1;i<=N;i++)
     {
	  printf("%5ld ",D[i]);
	  if (i%l==0) printf("\n");
     }
     if (N%l!=0) printf("\n");
}

bool List::Save_Binary(FILE *arch) const
{
     int nitems=fwrite(&N,sizeof(long),1,arch);
     if (nitems!=1) return false;
     if (!N) return true;
     nitems=fwrite(D+1,sizeof(long),N,arch);
     if (nitems!=N) return false;
     return true;
}

bool List::Save_Binary(const char *name) const
{
     FILE* fich=fopen(name,"wb");
     bool status=Save_Binary(fich);
     fclose(fich);
     return status;
}

bool List::Load_Binary(FILE *arch)
{
     int nitems=fread(&N,sizeof(long),1,arch);
     if (nitems!=1) return false;
     if (!N) return true;
     Create(N);
     nitems=fread(D+1,sizeof(long),N,arch);
     if (nitems!=N) return false;
     return true;
}

bool List::Load_Binary(const char *name)
{
     FILE* fich=fopen(name,"rb");
     bool status=Load_Binary(fich);
     fclose(fich);
     return status;
}

bool List::Load(const char *s)
{
     FILE *fich=fopen(s,"rt");
     bool status=Load(fich);
     fclose(fich);
     return status;
}

bool List::Load(FILE *fich)
{
     if (!fich) return false;
     long n, x;
     if (!fscanf(fich,"# %ld\n",&n)) return false;
     Create(n);
     Zero();
     for (long i=1;i<=n;i++)
     {
	  if (!fscanf(fich,"%ld\n",&x)) return false;
	  D[i]=x;
     }
     return true;
}

bool List::Save(const char *name) const
{
     FILE *fich=fopen(name,"wt");
     bool status=Save(fich);
     fclose(fich);
     return status;
}

bool List::Save(FILE *fich) const
{
     if (!fich) return false;
     if (fprintf(fich,"# %ld\n",N)<0) return false;
     for (long i=1;i<=N;i++)
	  if (fprintf(fich,"%ld\n",D[i])<0) return false;
     return true;
}

long List::operator() (long n1) const
{
    if (n1<0 || n1>N) 
	 Error("Error getting list component.\n");
    return D[n1];
}

long& List::operator() (long n1)
{
    if (n1<0 || n1>N) 
	 Error("Error putting list component.\n");
    return D[n1];
}

List& List::operator=(const List& L0) 
{
     if (this==&L0) return *this;
     Copy(*this,L0);
     return *this;
}

void List::operator+=(const List &L)
{
     if (!N) { Create(L.N); Zero(); }
#ifdef DEBUG
      if (N!=L.N) Error("Incompatible sizes in List +=");
#endif
      for (long i=1;i<=N;i++)
	   D[i]+=L(i);
}

void List::operator-=(const List &L)
{
     if (!N) { Create(L.N); Zero(); }
#ifdef DEBUG
      if (N!=L.N) Error("Incompatible sizes in List +=");
#endif
      for (long i=1;i<=N;i++)
	   D[i]-=L(i);
}

void List::operator*=(long p)
{
     if (!N) return;
     for (long i=1;i<=N;i++)
	  D[i]*=p;
}

void List::operator/=(long p)
{
     if (!N) return;
     for (long i=1;i<=N;i++)
	  D[i]/=p;
}

void List::operator&=(long p)
{
     Append(p);
}

void List::operator&=(const List &L)
{
     Append(L);
}

void List::operator+=(const long p)
{
     if (!N) return;
     for (long i=1;i<=N;i++)
	  D[i]+=p;
}

void List::operator-=(const long p)
{
     if (!N) return;
     for (long i=1;i<=N;i++)
	  D[i]-=p;
}

List operator+(const List &L1, const List &L2) 
{
     List R(L1);
     R+=L2;
     return R;
}

List operator-(const List &L1, const List &L2)
{
     List R(L1);
     R-=L2;
     return R;
}

List operator+(const List &L, const long p) 
{
     List R(L);
     R+=p;
     return R;
}

List operator-(const List &L, const long p)
{
     List R(L);
     R-=p;
     return R;
}

List operator+(const long p, const List &L) 
{
     List R(L);
     R+=p;
     return R;
}

List operator-(const long p, const List &L)
{
     List R(L);
     R*=-1;
     R+=p;
     return R;
}

List operator-(const List &L)
{
     List R(L);
     R*=-1;
     return L;
}

List operator*(const List &L, const long p)
{
     List R(L);
     R*=p;
     return R;
}

List operator/(const List &L, const long p)
{
     List R(L);
     R/=p;
     return R;
}

List operator*(const long p, const List &L)
{
     List R(L);
     R*=p;
     return R;
}

List operator&(const List &L1, const List &L2)
{
     List L(L1);
     L.Append(L2);
     return L;
}

List operator&(long p, const List &L2)
{
     List L(1); L(1)=p; 
     L.Append(L2);
     return L;
}

List operator&(const List &L2, long p)
{
     List L(L2);
     L.Append(p);
     return L;
}

void Copy(List& L1, const List &L2)
{
     L1.Destroy();
     L1.N=L2.N;
     if (L2.D!=NULL)
     {
	  L1.D=(long*)malloc((L1.N+1)*sizeof(long));
	  memcpy(L1.D,L2.D,(L1.N+1)*sizeof(long));
     }
}

// Returns true if lists are completely equal
bool Is_Equal(const List &L1, const List &L2)
{
     if (L1.N!=L2.N) return false;
     for (long i=1;i<=L1.N;i++)
	  if (L1(i)!=L2(i)) return false;
     return true;
}

// returns true if L1 is contained in L2
// more precisely, if all elements of L1 are also in L2
bool Is_Subset(const List &L1, const List &L2)
{
     for (int i=1;i<=L1.N;i++)
	  if (!L2.Find(L1(i))) return false;
     return true;
}

List Intersect(const List &L1, const List &L2)
{
     List L;
     if (L1.N*L2.N==0) return L;
     for (long i=1;i<=L1.N;i++)
	  if (L2.Find(L1(i))>0) L.Append(L1(i));
     return L;
}

// The resulting list is as large as L1
// For each element of L1, it looks up its value in L2
List Combine(const List &L1, const List &L2)
{
     List L(L1);
     for (long i=1;i<=L1.N;i++)
	  L(i)=L2(L1(i));
     return L;
}

List List_Range(long i1, long i2)
{
     List L(labs(i2-i1)+1);
     if (i2>i1)
	  for (long i=i1;i<=i2;i++) L(i-i1+1)=i;
     else
	  for (long i=i1;i>=i2;i--) L(i1-i+1)=i;
     return L;
}

List Reverse(const List &L, long i, long j)
{
     List L2(L);
     L2.Reverse(i,j);
     return L2;
}

List Part(const List &L, long i1, long i2)
{
     List L2(L);
     L2.Part(i1,i2);
     return L2;
}

List Substract(const List &L1, const List &L2) // OPTIMIZE
{
     List L(L1);
     L.Substract(L2);
     return L;
}

// Convert long integer "j" into the list of bits "1" in its binary expansion
List List_Ones(long j)
{
     List L;
     for (long k=0;k<=Max_Bit(j);k++)
          if (Bit(j,k)) L.Append(k);
     return L;
}

List  Bin_2_List(long X, int maxbit) // conversion of a binary number to a list
{
     List L;
     long Y=X;
     while(Y)
     {
	  if (Y%2) L.Append(1);
	  else L.Append(0);
	  Y/=2;
     }
     // now if size of list is nbits (or nbits is 0) done, otherwise, complete
     while (L.N<maxbit+1)
	    L.Append(0);
     return L;
}

long List_2_Bin(const List &L)
{
     long X=0;
     for (long i=1;i<=L.N;i++)
	  X+=L(i)*Pow_2(i-1);
     return X;
}

// Generic number system function
// Takes integer and list Ls giving the number system, returns repr
List Int_2_List(long ix, const List &Ls)
{
     long i=ix;
     List X(Ls.N);
     for (long k=Ls.N;k>=1;k--)
     {
	  X(k)=i % Ls(k);
	  i/=Ls(k);
     }
     return X;
}

// Takes number representation and number system, returns int
long List_2_Int(const List &X, const List &Ls)
{
     long N=Ls.N;
     long ix=X(N);
     if (Ls.N<=1) return ix;
     long prod=1;
     for (long k=Ls.N-1;k>=1;k--)
     {
	  prod*=Ls(k+1);
	  ix+=X(k)*prod;
     }
     return ix;
}

int Sort_Long_12 (const void *a, const void *b)
{
    const long *da = (const long *) a;
    const long *db = (const long *) b;
    return (*da > *db) - (*da < *db);
}

int Sort_Long_21(const void* x, const void* y) // descending order
{
    const long *xx=(const long *) x;
    const long *yy=(const long *) y;
    return (*yy > *xx) - (*xx > *yy);
}

int Sort_Double_12(const void* x, const void* y) // ascending order
{
     const double *dx=(const double *)x;
     const double *dy=(const double *)y;
     return (*dx > *dy) - (*dx < *dy);
}

int Sort_Double_21(const void* x, const void* y) // descending order
{
     const double *dx=(const double *)x;
     const double *dy=(const double *)y;
     return (*dx < *dy) - (*dx > *dy);
}

List Random_Permutation(long N) // random perm of {1...L}
{
     List L(2*N);
     for (long i=1;i<=N;i++)
     {
	  L(2*i-1)=Rand_Full();
	  L(2*i)=(double)i;
     }
     qsort(L.D+1,N,2*sizeof(long),Sort_Long_12);
     List P(N);
     for (long i=1;i<=N;i++)
	  P(i)=L(2*i);
     return P;
}

// Algorithm to get the next permutation in lexicographical order
// find the longest non-increasing tail, let k=last element of head
// find the smallest element in tail which is largest than k
// swap those elements, revert resulting tail
// the return value is \pm 1 depending on the sign of the relative permutation
// return is zero if we reached the last permutation
long Next_Permutation(List &L)
{
     if (L.N<=1) return 0;
     long k=L.N; // the index separating head and tail
     while(L(k-1)>=L(k) && k>1) k--; // now k marks the beginning of tail
     if (k==1) return false; // same permutation, meaning that we're done
     long ip=L(k-1); // pivot element
     // Now, we find first element in tail which is larger than ip
     long k2=L.N;
     while(L(k2)<=ip && k2>k) k2--;
     L.Swap(k-1,k2);
     L.Reverse(k,L.N);
     return Mop(1+(L.N-k+1)/2);
}


/////////////////////////////////////////////////////////////////////
// Class Table, array of long int
////////////////////////////////////////////////////////////////////

Table::Table()
{
     Start();
}

Table::Table(long n1, long n2) : N1(n1), N2(n2)
{
     Start();
     Create(n1,n2);
     Zero();
}

Table::Table(const Table & M)
{
     Start();
     Copy(*this,M);
}

Table::~Table()
{
     Destroy();
}

void Table::Start()
{
     N1=N2=0;
     D=(long*)NULL;
}

void Table::Create(long n1, long n2)
{
     Destroy();
     N1=n1; N2=n2;
     if (!N2) N2=N1;
     if (!N1)
     {
          D=NULL;
          return;
     }
     D=(long*)malloc((N1*N2+1)*sizeof(long));
#ifdef DEBUG
     if (!D) Error("Error allocating Table\n");
#endif
}

void Table::Load(long *d, long n1, long n2)
{
     Destroy();
     N1=n1;
     N2=n2; if (!N2) N2=N1;
     D=d;
}

void Table::Load_Copy(long *d1, long n1, long n2)
{
     if (!n2) n2=n1;
     long *d2=(long*)malloc((n1*n2+1)*sizeof(long));
     memcpy(d2,d1,(n1*n2+1)*sizeof(long));
     Load(d2,n1,n2);
}

void Table::Transfer(Table &T)
{
     Destroy();
     D=T.D;
     N1=T.N1;
     N2=T.N2;
     T.Start();
}

void Table::Destroy()
{
     if (N1) free(D);
     D=NULL;
     N1=N2=0;
}

void Table::Resize(long n1, long n2)
{
     if (!n2) n2=n1;
     if (n1==N1 && n2==N2) return; 
     long nr1=::Min(N1,n1), nr2=::Min(N2,n2);
     long oldN1=N1, oldN2=N2;
     
     long *D2=(long*)malloc((N1*N2+1)*sizeof(long));
     memcpy(D2+1,D+1,N1*N2*sizeof(long));
     
     Create(n1,n2);
     if (n1>oldN1 || n2>oldN2)
	  Zero(); // In case the new matrix is bigger
     for (long i=1;i<=nr2;i++)
	  memcpy(D+N1*(i-1)+1,D2+oldN1*(i-1)+1,nr1*sizeof(long));
     free(D2);

}

void Table::Zero()
{
     memset(D+1,0,N1*N2*sizeof(long));
}

void Table::Set(long x)
{
     for (long i=1;i<=N1;i++)
	  for (long j=1;j<=N2;j++)
	       Elem(i,j)=x;
}

long& Table::Elem(long i,long j)
{
#ifdef DEBUG
     if (i<1 || i>N1 || j<1 || j>N2) 
	  Error("Error accessing Table element.");
#endif
     return D[N1*(j-1)+i];
}

void Table::Set_Col(const List &L, long n)
{
#ifdef DEBUG
     if ((L.N!=N1) || (n>N2)) Error("Set_Col is impossible.");
#endif
     memcpy(D+(n-1)*N1+1,L.D+1,L.N*sizeof(long));
   
}

void Table::Set_Row(const List &L, long n)
{
     for (long i=1;i<=N2;i++)
	  Elem(n,i)=L(i);
}

void Table::Set_Diag(const List &L)
{
#ifdef DEBUG
     if (L.N>::Min(N1,N2)) Error("Error in Set_Diag.");
#endif
     for (long i=1;i<=L.N;i++)
	  Elem(i,i)=L(i);
}

void Table::Append_Col(const List &L)
{
     if (N1)
     	  Resize(N1,N2+1);
     else
     	  Create(L.N,1);
     Set_Col(L,N2);
}

void Table::Append_Col(const Table &T)
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

void Table::Append_Row(const List &L)
{
     if (N2)
     	  Resize(N1+1,N2);
     else
     	  Create(1,L.N);
     Set_Row(L,N1);
}

void Table::Append_Row(const Table &T)
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

void Table::Insert_Col(const List &L, long p)
{
#ifdef DEBUG
     if (L.N!=N1) Error("Insert_Col with wrong dimensions\n");
#endif
     if (!N1)
     {
	  Create(L.N,1);
	  Set_Col(L,1);
     }
     Resize(N1,N2+1);
     for (long i=N2;i>=p+1;i--)
	  memmove(D+1+(i-1)*N1,D+1+(i-2)*N1,N1*sizeof(long));
     Set_Col(L,p);
}

void Table::Insert_Row(const List &L, long p)
{
#ifdef DEBUG
     if (L.N!=N2) Error("Insert_Row with wrong dimensions\n");
#endif
     if (!N1)
     {
	  Create(1,L.N);
	  Set_Row(L,1);
     }
     Resize(N1+1,N2);
     for (long i=N1;i>=p+1;i--)
	  for (long j=1;j<=N2;j++)
	       Elem(i,j)=Elem(i-1,j);
     Set_Row(L,p);
}

void Table::Remove_Col(long p)
{
#ifdef DEBUG
     if (p>N2) Error("Remove non-existent col\n");
#endif
     if (N2==1) { Destroy(); return; }
     for (long i=p+1;i<=N2;i++)
	  memcpy(D+(i-2)*N1+1,D+(i-1)*N1+1,N1*sizeof(long));
     Resize(N1,N2-1);
}

void Table::Remove_Row(long p)
{
#ifdef DEBUG
     if (p>N1) Error("Remove non-existent row\n");
#endif
     if (N1==1) { Destroy(); return; }
     Table R(N1-1,N2);
     for (long j=1;j<=N2;j++)
     {
	  for (long i=1;i<=p-1;i++)
	       R(i,j)=Elem(i,j);
	  for (long i=p+1;i<=N1;i++)
	       R(i-1,j)=Elem(i,j);
     }
     Transfer(R);
}

void Table::Swap_Cols(long k1,long k2)
{
     double acum;
     for(long i=1;i<=N1;i++)
     {
	  acum=Elem(i,k1);
	  Elem(i,k1)=Elem(i,k2);
	  Elem(i,k2)=acum;
     }
}

void Table::Swap_Rows(long k1,long k2)
{
    long acum;
    for(long i=1;i<=N2;i++)
    {
	 acum=Elem(k1,i);
	 Elem(k1,i)=Elem(k2,i);
	 Elem(k2,i)=acum;
    }
}

void Table::Permute_Cols(const List &P)
{
     Table A(N1,N2);
     for (long i=1;i<=N2;i++)
	  A.Set_Col(Col(P(i)),i);
     Transfer(A);
}

void Table::Permute_Rows(const List &P)
{
     Table A(N1,N2);
     for (long i=1;i<=N1;i++)
	  A.Set_Row(Row(P(i)),i);
     Transfer(A);
}

void Table::Sort_Cols()
{
     qsort(D+1,N2,N1*sizeof(long),Sort_Long_12);
}

void Table::T()
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

void Table::Part(long n10,long n20,long n1f,long n2f)
{
     Table R;
     ::Part(R,*this,n10,n20,n1f,n2f);
     Transfer(R);
}

void Table::Set_Part(const Table &T, long i, long j)
{
#ifdef DEBUG
     if (i+T.N1-1>N1) Error("Incompatible dimensions in Set_Part.");
     if (j+T.N2-1>N2) Error("Incompatible dimensions in Set_Part.");
#endif
     long m=T.N2;
     for (long k=1;k<=m;k++)
	  memcpy(D+N1*(j+k-2)+i,T.D+T.N1*(k-1)+1,T.N1*sizeof(long));
}

void Table::Rand(long a, long b)
{
     for (long i=1;i<=N1*N2;i++)
	  D[i]=Rand_I(a,b);
}

bool Table::Is_Zero() const
{
     for (long i=1;i<=N1*N2;i++)
	  if (D[i]) return false;
     return true;
}

long Table::Elem(long i, long j) const
{
#ifdef DEBUG
     if (i<1 || i>N1 || j<1 || j>N2) 
	  Error("Error accessing matrix element.");
#endif
     return D[N1*(j-1)+i];
}

List Table::Col(long n) const
{
     List R;
     Col(R,n);
     return R;
}

List Table::Row(long n) const
{
     List R;
     Row(R,n);
     return R;
}

List Table::Diag() const
{
     long M=::Min(N1,N2);
     List V(M);
     for (long i=1;i<=M;i++)
	  V(i)=Elem(i,i);
     return V;
}

void Table::Col(List &R, long n) const
{
     long *d=(long*)malloc((N1+1)*sizeof(long));
     memcpy(d+1,D+1+(n-1)*N1,N1*sizeof(long));
     R.Load(d,N1);
}

void Table::Row(List &R, long n) const
{
     R.Create(N2);
     R.Zero();
     for (long i=1;i<=N2;i++)
	  R(i)=Elem(n,i);
}

long Table::Max() const
{
     long X=LONG_MIN;
     for (long i=1;i<=N1;i++)
	  for (long j=1;j<=N2;j++)
	  {
	       long x=Elem(i,j);
	       if (x>X) X=x;
	  }
     return X;
}

long Table::Max(long &im, long &jm) const
{
     long X=LONG_MIN;
     for (long i=1;i<=N1;i++)
	  for (long j=1;j<=N2;j++)
	  {
	       long x=Elem(i,j);
	       if (x>X) 
	       {
		    X=x;
		    im=i; jm=j;
	       }
	  }
     return X;
}

long Table::Min() const
{
     long X=LONG_MAX;
     for (long i=1;i<=N1;i++)
	  for (long j=1;j<=N2;j++)
	  {
	       long x=Elem(i,j);
	       if (x<X) X=x;
	  }
     return X;
}

long Table::Min(long &im, long &jm) const
{
     long X=LONG_MAX;
     for (long i=1;i<=N1;i++)
	  for (long j=1;j<=N2;j++)
	  {
	       long x=Elem(i,j);
	       if (x<X) 
	       {
		    X=x;
		    im=i; jm=j;
	       }
	  }
     return X;
}

void Table::Write() const
{
     for (long i=1;i<=N1;i++)
     {
         for (long j=1;j<=N2;j++)
	 {
	      long x=Elem(i,j);
	      printf("%5ld ",x);
	 }
         printf("\n");
     }
     printf("\n");
}

bool Table::Save_Binary(const char *name) const
{
     FILE *fich=fopen(name,"wb");
     if (!fich) return false;
     bool status=Save_Binary(fich);
     fclose(fich);
     return status;
}

bool Table::Load_Binary(const char *name)
{
     FILE *fich=fopen(name,"rb");
     if (!fich) return false;
     bool status=Load_Binary(fich);
     fclose(fich);
     return status;
}

bool Table::Save_Binary(FILE *fich) const
{
     int nwrite=fwrite(&N1,sizeof(long),1,fich);
     if (nwrite!=1) { Error_Flag(Error_IO); return false; }
     nwrite=fwrite(&N2,sizeof(long),1,fich);
     if (nwrite!=1) { Error_Flag(Error_IO); return false; }
     if (!(N1*N2)) return true;
     nwrite=fwrite(D+1,sizeof(long),N1*N2,fich);
     if (nwrite!=N1*N2) { Error_Flag(Error_IO); return false;}
     if (ferror(fich)) return false;
     return true;
}

bool Table::Load_Binary(FILE *fich)
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
	  D=(long*)NULL; 
	  return true;
     }
     Create(N1,N2);
     Zero();
     nread=fread(D+1,sizeof(long),N1*N2,fich);
     if ((nread!=N1*N2) || (ferror(fich)))
     { 
	  Error_Flag(Error_IO); 
	  return false; 
     }
     return true;
}

bool Table::Save(const char *name) const
{
     FILE *fich=fopen(name,"wt");
     if (!fich) return false;
     Save(fich);
     fclose(fich);
     return true;

}

bool Table::Load(const char *name)
{
     FILE *fich=fopen(name,"rt");
     if (!fich) return false;
     bool status=Load(fich);
     fclose(fich);
     return status;
}

bool Table::Save(FILE *fich) const
{
     if (fprintf(fich,"# %ld %ld\n",N1,N2)<0)
     {
	  Error_Flag(Error_IO); return false; 
     }
     for (long i=1;i<=N1;i++)
     {
	  for (long j=1;j<=N2;j++)
	       fprintf(fich,"%ld ",Elem(i,j));
	  fprintf(fich,"\n");
     }
     if (ferror(fich)) 
     {
	  Error_Flag(Error_IO); return false;
     }
     return true;
}

bool Table::Load(FILE *fich)
{
     long n1, n2;
     if (fscanf(fich,"# %ld %ld\n",&n1,&n2)<2)
     {
	  Error_Flag(Error_IO); return false;
     }
     Destroy();
     Create(n1,n2);
     for (long i=1;i<=n1;i++)
	  for (long j=1;j<=n2;j++)
	  {
	       long x;
	       if (fscanf(fich,"%ld ",&x)<1) 
	       {
		    Error_Flag(Error_IO); return false;
	       }
	       Elem(i,j)=x;
	  }
     return true;
}

long& Table::operator()(long i, long j)
{
#ifdef DEBUG
     if (i<1 || i>N1 || j<1 || j>N2) 
	  Error("Error accessing table element.");
#endif
     return D[N1*(j-1)+i];
}

long Table::operator() (long i, long j) const
{
#ifdef DEBUG
     if (i<1 || i>N1 || j<1 || j>N2) 
	  Error("Error accessing table element.");
#endif
     return D[N1*(j-1)+i];
    
}

Table& Table::operator=(const Table& T)
{
     if (this==&T) return *this;
     Copy(*this,T);
     return (*this);
}

void Table::operator+=(const Table& T)
{
     if (!N1) { Create(T.N1,T.N2); Zero(); }
#ifdef DEBUG
     if (N1!=T.N1 || N2!=T.N2) 
	  Error("Adding matrices with different dims.");
#endif
     for (long i=1;i<=T.N1;i++)
	  for (long j=1;j<=T.N2;j++)
	       Elem(i,j)+=T(i,j);
}

void Table::operator-=(const Table& T)
{
     if (!N1) { Create(T.N1,T.N2); Zero(); }
#ifdef DEBUG
     if (N1!=T.N1 || N2!=T.N2) 
	  Error("Adding matrices with different dims.");
#endif
     for (long i=1;i<=T.N1;i++)
	  for (long j=1;j<=T.N2;j++)
	       Elem(i,j)-=T(i,j);
}

void Table::operator+=(long K)
{
#ifdef DEBUG
     if (!N1) return;
#endif
     for (long i=1;i<=N1;i++)
	  for (long j=1;j<=N2;j++)
	       Elem(i,j)+=K;
}

void Table::operator-=(long K)
{
#ifdef DEBUG
     if (!N1) return;
#endif
     for (long i=1;i<=N1;i++)
	  for (long j=1;j<=N2;j++)
	       Elem(i,j)-=K;
}

void Table::operator*=(long K)
{
#ifdef DEBUG
     if (!N1) return;
#endif
     for (long i=1;i<=N1;i++)
	  for (long j=1;j<=N2;j++)
	       Elem(i,j)*=K;
}
  
void Table::operator/=(long K)
{
#ifdef DEBUG
     if (!N1) return;
#endif
     for (long i=1;i<=N1;i++)
	  for (long j=1;j<=N2;j++)
	       Elem(i,j)/=K;
}

void Table::operator&=(const List &L)
{
     Append_Row(L);
}

void Table::operator&=(const Table &T)
{
     Append_Row(T);
}

void Table::operator|=(const List &L)
{
     Append_Col(L);
}

void Table::operator|=(const Table &T)
{
     Append_Col(T);
}

void Copy(Table& B, const Table& A) // B <- A raw and strict copy. 
{
     if (!A.N1) return;
     B.Create(A.N1,A.N2);
     memcpy(B.D+1,A.D+1,A.N1*A.N2*sizeof(long));
}

void Write(const Table &T)
{
     T.Write();
}

void Part(Table &R, const Table &M,long n10,long n20,long n1f,long n2f)
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
		 m1*sizeof(long));
}

Table operator-(const Table& T)
{
     Table R(T);
     R*=-1;
     return R;
}

Table operator+(const Table& A, const Table& B)
{
#ifdef DEBUG
     if (A.N1!=B.N1 || A.N2!=B.N2) 
	  Error("Adding tables with different dimensions.");
#endif
     Table R(A.N1,A.N2);
     for (long i=1;i<=A.N1;i++)
	  for (long j=1;j<=A.N2;j++)
	       R(i,j)=A(i,j)+B(i,j);
     return R;
}

Table operator-(const Table& A, const Table& B)
{
#ifdef DEBUG
     if (A.N1!=B.N1 || A.N2!=B.N2) 
	  Error("Subtracting tables with different dimensions.");
#endif
     Table R(A.N1,A.N2);
     for (long i=1;i<=A.N1;i++)
	  for (long j=1;j<=A.N2;j++)
	       R(i,j)=A(i,j)-B(i,j);
     return R;

}

Table operator*(long K, const Table &A)
{
     Table T(A);
     T*=K;
     return T;
}

Table operator*(const Table &A, long K)
{
     Table T(A);
     T*=K;
     return T;
}

Table operator/(const Table &A, long K)
{
     Table T(A);
     T/=K;
     return T;
}

Table operator|(const List &L1, const List &L2)
{
     Table R;
     R.Append_Col(L1);
     R.Append_Col(L2);
     return R;
}

Table operator|(const Table &T, const List &L)
{
     Table R(T);
     R.Append_Col(L);
     return R;
}

Table operator|(const List &L, const Table &T)
{
     Table R(T);
     R.Insert_Col(L,1);
     return R;
}

Table operator|(const Table &T1, const Table &T2)
{
     Table R(T1);
     R.Append_Col(T2);
     return R;
}

Table operator&(const Table &T, const List &L)
{
     Table R(T);
     R.Append_Row(L);
     return R;
}

Table operator&(const List &L, const Table &T)
{
     Table R(T);
     R.Append_Row(L);
     return R;
}

Table operator&(const Table &T1, const Table &T2)
{
     Table R(T1);
     R.Append_Row(T2);
     return R;
}
