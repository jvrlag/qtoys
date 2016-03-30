////////////////////////////////////////////////////////
// hvb++ 1.0
// Copyleft: Javier Rodríguez Laguna
// 080725-150125
// Common routines library
// JaviRL, 060404
#ifndef COMMON_HEADER
#define COMMON_HEADER
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/times.h>
#include <unistd.h>
#include <complex>
#include <iostream>
#include <ctype.h>
#include <limits.h>
#include <errno.h>

////////////////////////////////////////////////////////
// General definitions
////////////////////////////////////////////////////////

using namespace std;
#define cmplx complex<double>
extern cmplx M_I;
cmplx operator+(long,cmplx);
cmplx operator+(cmplx,long);
cmplx operator-(long,cmplx);
cmplx operator-(cmplx,long);

enum Side {Left, Right};

double Sqr(double a);
long   Max(long a, long b);
double Max(double a, double b);
long   Min(long a, long b);
double Min(double a, double b);
void   Swap(int &a, int &b);
void   Swap(long &a, long &b);
void   Swap(double &a, double &b);
void   Swap(cmplx &a, cmplx &b);
long   Mop(long n); // minus one power n, (-1)^n
double Sign(double x);
double Sign(double x, double y);
// combinatorial functions
double Fact(double x);
long   Fact(long n);
double Comb(double x, double y);
long   Comb(long n, long m);

////////////////////////////////////////////////////////
// Error handling
////////////////////////////////////////////////////////

class List;

enum   Error_Type { No_Error=0, Error_IO, Error_Mem, Error_Mat, Error_Dev };
void   Error(const char *) __attribute__ ((noinline));
void   Error_Flag(Error_Type);
Error_Type Error_Read();
void   Error_Clean();
extern List Error_Flags;

////////////////////////////////////////////////////////
// Time measurement
////////////////////////////////////////////////////////

double Clock();
void   Delay(double Dt);

////////////////////////////////////////////////////////
// Random number generator
////////////////////////////////////////////////////////
#define MTRNG_N 624
typedef struct
{
     unsigned long mt[MTRNG_N];
     int mti;
}Rand_State_Type;
extern Rand_State_Type *Rand_State;
extern bool Rand_Started;
void   Rand_Open(unsigned long s);
unsigned long Rand_Full();
double Rand();
double Rand(double a, double b);
long   Rand_I(long a, long b);
double Rand_Gaussian(double mu, double sigma);
void   Rand_Close();


////////////////////////////////////////////////////////
// Binary numbers
////////////////////////////////////////////////////////

long  Pow_2(long num);
long  Log_2(long num);
bool  Bit(long X, long pos);
long  Max_Bit(long X);
long  Flip_Bit(long X, long pos); 
long  Put_Bit(long X, long pos, long value);
long  Swap_Bits(long X, long pos1, long pos2);
long  Reverse_Bits(long X, long maxbit=0);
long  Insert_Bit(long X, long pos, long value);
long  Count_Ones(long X); 
long  Next_In_Sector(long X);
char* Bin_2_String(long X, long maxbit=0);

//////////////////////////////////////////////////////////////////////

class List
{
public:
     long N; // size of list
     long *D;

     List();     // empty initializer
     List(long); // initialize a list with a given size
     List(long*,long); // initialize a list copying a number of items from long*
     List(const List&); // initialize a list from another
     ~List(); // destructor
     void Start();      // void allocation
     void Create(long); // allocates
     void Load(long *, long);
     void Load_Copy(long *, long);
     void Transfer(List &);
     void Destroy();    // deallocates
     
     void Zero();        // puts all elements to zero
     void Set(long);     // initialize the list to a given value
     void Set_Part(const List &L, long p); // put List at site p, overwriting
     void Remove(long i,long j=0);  // remove from position i to j
     long Substract(long);    // eliminates all occurrences of (long)
     long Substract(const List &); // remove all elements from List
     void Append(long);        // appends an element at the end
     void Append(const List &);      // appends a full list at the end
     void Sort(long);                // sorts the list
     void Uniquify();                // removes all repeated values
     void Insert(long, long);        // insert a number at a given position
     void Insert(const List&, long); // insert a list at a given position
     void Reverse(long i=0, long j=0); // reverse the order of (part) of a list
     void Swap(long,long); // swaps the values at two positions
     void Part(long,long); // returns a part of the list
     void Rand(long,long);

     long Find(long,long k=1) const; // finds first appearance from pos k
     long Find_Non(long,long k=1) const;
     List Find_All(long) const;
     long Count(long x) const; // counts the number of appearances of x
     long Min() const;
     long Max() const;
     long Min(long &imin) const;
     long Max(long &imax) const;
     long Sum() const;
     long Sum(long,long) const; // sum of elements between two given
     long Prod() const; // product of all elements
     long Prod(long,long) const; // takes the product of all elements

     void Write() const; // writes List to stdout
     bool Save_Binary(FILE *) const;
     bool Save_Binary(const char *name) const;
     bool Load_Binary(FILE *);
     bool Load_Binary(const char *name);
     bool Load(const char *);
     bool Load(FILE *);
     bool Save(const char *) const;
     bool Save(FILE *) const;

     long operator() (long) const;
     long& operator() (long);
     List& operator=(const List&);
     void operator+=(const List&);
     void operator-=(const List&);
     void operator*=(long);
     void operator/=(long);
     void operator+=(long);
     void operator-=(long);

     void operator&=(const List &);
     void operator&=(long);
};

List operator+(const List &L1, const List &L2); 
List operator-(const List &L1, const List &L2);
List operator+(const List &L, long p); 
List operator-(const List &L, long p);
List operator+(long p, const List &L); 
List operator-(long p, const List &L);
List operator-(const List &L);
List operator*(const List &L, long p);
List operator/(const List &L, long p);
List operator*(long p, const List &L);

List operator&(const List &L1, const List &L2);
List operator&(long, const List &L);
List operator&(const List &L, long);

void Copy(List &L1, const List &L2);
bool Is_Equal(const List &L1, const List &L2);
bool Is_Subset(const List &L1, const List &L2);
List Intersect(const List &L1, const List &L2);
List Combine(const List &, const List &);
List List_Range(long, long);
List Reverse(const List &L, long i=0, long j=0);
List Part(const List &L, long i1, long i2);
List Substract(const List &L1, const List &L2);

List List_Ones(long X); // returns the list of positions of the 1's in X
List Bin_2_List(long X, int maxbit=0); // conversion of a binary to a list
long List_2_Bin(const List &L);
List Int_2_List(long ix, const List &Ls);
long List_2_Int(const List &X, const List &Ls);

int Sort_Double_12(const void* x, const void* y); // ascending order
int Sort_Double_21(const void* x, const void* y); // descending order
int Sort_Long_12(const void* x, const void* y); // ascending order
int Sort_Long_21(const void* x, const void* y); // descending order

List Random_Permutation(long N); // random perm of {1...L}
long Next_Permutation(List &L); // returns relative sign


/////////////////////////////////////////////////////////////////////
// Class Table, 2D array of long int
/////////////////////////////////////////////////////////////////////

class Table
{
public:
     long *D;
     long N1, N2;

     Table();
     Table(long n, long m=0);
     Table(const Table &);
     ~Table();
     void Start();
     void Create(long n, long m=0);
     void Load(long *, long, long=0);
     void Load_Copy(long *, long, long=0);
     void Transfer(Table &);
     void Destroy();
     void Resize(long n, long m=0);

     void Zero();
     void Set(long x);
     long &Elem(long,long);
     void Set_Col(const List &, long);     
     void Set_Row(const List &, long);
     void Set_Diag(const List &);
     void Append_Col(const List &);
     void Append_Col(const Table &);
     void Append_Row(const List &);
     void Append_Row(const Table &);
     void Insert_Col(const List &, long);
     void Insert_Row(const List &, long);
     void Remove_Col(long);
     void Remove_Row(long);
     void Swap_Cols(long,long);
     void Swap_Rows(long,long);
     void Permute_Cols(const List &L);
     void Permute_Rows(const List &L);
     void Sort_Cols();
     void T();
     void Part(long,long,long,long);
     void Set_Part(const Table &, long, long);
     void Rand(long,long);

     bool Is_Zero() const;
     long Elem(long, long) const;
     List Col(long) const;
     List Row(long) const;
     List Diag() const;
     void Col(List &, long) const;
     void Row(List &, long) const;
     long Max() const;
     long Max(long &, long &) const;
     long Min() const;
     long Min(long &, long &) const;

     void Write() const;
     bool Save_Binary(const char *name) const;
     bool Load_Binary(const char *name);
     bool Save_Binary(FILE *fich) const;
     bool Load_Binary(FILE *fich);
     bool Save(const char *name) const;
     bool Load(const char *name);
     bool Save(FILE *fich) const;
     bool Load(FILE *fich);

     long& operator()(long, long);
     long operator() (long, long) const;
     Table& operator=(const Table&);
     void operator+=(const Table&);
     void operator-=(const Table&);
     void operator+=(long);
     void operator-=(long);
     void operator*=(long);  
     void operator/=(long);
     void operator&=(const List &);
     void operator&=(const Table &);
     void operator|=(const List &);
     void operator|=(const Table &);
};

void Copy(Table& B, const Table& A);
void Write(const Table &T);
void Part(Table &R, const Table &M,long,long,long,long);

Table operator-(const Table &);
Table operator+(const Table &, const Table &);
Table operator-(const Table &, const Table &);
Table operator*(long K, const Table &A);
Table operator*(const Table &A, long K);
Table operator/(const Table &A, long K);

Table operator|(const List &, const List &);
Table operator|(const Table &, const List &);
Table operator|(const List &, const Table &);
Table operator|(const Table &, const Table &);
Table operator&(const Table &, const List &);
Table operator&(const List &, const Table &);
Table operator&(const Table &, const Table &);


#endif













