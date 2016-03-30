// Parsing routines
// 111130-150510

#include"Text.h"

Text::Text()
{
     N=0; D=(char*)NULL;
}

Text::Text(const char *z)
{
     Start();
     Create(z);
}

Text::Text(const Text &S)
{
     Start();
     Create(S.D);
}

void Text::Start()
{
     N=0; D=(char*)NULL;
}

Text::~Text()
{
     Destroy();
}

void Text::Destroy()
{
     if (N) free(D);
     N=0;
}

void Text::Create(const char *z)
{
     Destroy();
     if (!z) return;
     N=strlen(z);
     D=(char*)malloc(N+1);
     strncpy(D,z,N);
     D[N]=0;
}

void Text::Load(char *z) // just take the pointer
{
     D=z;
     N=strlen(z);
}

// returns number of read characters, if 0 then EOF, if -1, then error
long Text::Get_Line(FILE *fich)
{
     Destroy();
     char *z=NULL;
     ssize_t Nr; size_t N=0;
     Nr=getline(&z,&N,fich);
     if (Nr==-1)  // smth went wrong, EOF or other thing
     { 
	  free(z); 
	  if (feof(fich)) return -2;  // plain EOF, no big deal
	  Error_Flag(Error_IO); // something worse
	  return -1; 
     }
     if (Nr) Load(z); 
     return Nr;
}
 
void Text::Write() const
{
     if (!D) { printf("Null Text\n"); return; }
     printf("%s\n",D);
}

// Append a string, until you get \0
void Text::Append(const char *s)
{
     if (!s) return;
     long n=strlen(s);
     if (!D) D=(char*)malloc(n+1);
     else D=(char*)realloc(D,N+n+1);
     strncpy(D+N,s,n);
     N+=n;
     D[N]=0;
}

// Append n chars from a string
void Text::Append(const char *s, long n)
{
      if (!s) return;
      D=(char*)realloc(D,N+n+1);
      strncpy(D+N,s,n);
      N+=n;
      D[N]=0;
}

// Append another string
void Text::Append(const Text &S)
{
     Append(S.D);
}

// Append the string resulting from printing "i" with a given format
void Text::Append_F(const char *formato, long i)
{
     char *s=(char*)malloc(100);
     sprintf(s,formato,i);
     Append(s);
     free(s);
}

// Append the string resulting from printing "i" with a given format
void Text::Append_F(const char *formato, double x)
{
     char *s=(char*)malloc(100);
     sprintf(s,formato,x);
     Append(s);
     free(s);
}

// Return S from chars i0 to i1, both included (first char is 0!!!)
Text Text::Part(long i0, long i1) const
{
     Text S2;
     if (i0>i1) return S2; 
     long n=i1-i0+1;
     S2.Append(D+i0,n);
     return S2;
}

Text& Text::operator=(const Text &S)
{
     //if (this==&S) return *this;
     Copy(*this,S);
     return *this;
}

Text& Text::operator=(const char *z)
{
     Create(z);
     return *this;
}

ostream& operator<< (ostream &out, Text &S)
{
     out << S.D;
     return out;
}


void Copy(Text &S, const Text &S1)
{
     S.Destroy();
     S.N=S1.N;
     if (S1.D!=NULL)
     {
	  S.D=(char*)malloc(S1.N+1);
	  strncpy(S.D,S1.D,S1.N);
	  S.D[S.N]=0;
     }
}

void Copy(Text &S, const char *z)
{
     S.Destroy();
     if (!z) return;
     long n=strlen(z);
     S.D=(char*)malloc(n+1);
     strcpy(S.D,z);
     S.N=n;
}

void Copy(Text &S, const char *z, long n)
{
     S.Destroy();
     if (!z) return;
     S.D=(char*)malloc(n+1);
     strncpy(S.D,z,n);
     S.D[n]=0;
     S.N=n;
}

// remove blanks from beginning and end of char, free the pointer
void Text::Strip_Blanks()
{
     if (N==0) return;
     long i0=0;
     while( (D[i0]==' ' || D[i0]==10) && i0<N)
	  i0++;  // now D[i0]!=' '
     long i1=strlen(D)-1;
     while( (D[i1]==' ' || D[i1]==10) && i1>i0)
	  i1--;
     long n=i1-i0+1;
     if (n<=0) 
     {
	  Destroy();
	  return;
     }
     char *z=(char*)malloc(n+1);
     strncpy(z,D+i0,n);
     z[n]=0;
     Destroy();
     D=z;
     N=n;
}

void Text::Append(const char q)
{
     D=(char*)realloc(D,N+2);
     D[N]=q;
     D[N+1]=0;
     N++;
}

bool operator==(const Text &S1, const Text &S2)
{
     if (S1.N!=S2.N) return false;
     return !strcmp(S1.D,S2.D);
}

bool operator!=(const Text &S1, const Text &S2)
{
     return !(S1==S2);
}

bool Text::Is_Here(const char *z)
{
     if (!N) return false;
     return(memcmp(D,z,strlen(z))==0);
}

bool Text::Is_Here(const Text &S)
{
     return Is_Here(S.D);
}

bool Text::Is_There(const char *z)
{
     if (!N) return false;
     if (strstr(D,z)!=NULL) return true;
     else return false;
}

bool Text::Is_There(const Text &S)
{
     return Is_There(S.D);
}

// return the substring from the n-th appearance of q0 
// until the next appearance of q1
// if q0==0, it means "from the beginning"; if q1==0, "to the end"
// if no appropriate q0 is found, returns empty Text
// if no appropriate q1 is found, returns from q0 to end of the Text
Text Text::Token(char q0, char q1, long n) const
{
     long i0, i1; // limiting indices
     if (!q0) i0=-1;
     else
     {
	  i0=Find_Nth(q0,n);
	  if (i0==-1) { Text S; return S; } // return empty string  
     }
     Text S=Part(i0+1,N);
     if (!q1) i1=S.N;
     else
	  i1=S.Find_Nth(q1,1);
     if (i1==-1) i1=S.N;
     return S.Part(0,i1-1);
}

long Text::Get_Int(char q0, char q1, long n) const
{

     Text S=Token(q0,q1,n);
     return S.To_Int();
}

double Text::Get_Real(char q0, char q1, long n) const
{
     Text S=Token(q0,q1,n);
     return S.To_Real();
}

// count how many appearances of the given character
long Text::Count(char q) const
{
     long n=0;
     for (long i=0;i<=N;i++) 
	  if (D[i]==q) n++;
     return n;
}

// find the n-th appearance of the character q in the string
// -1 means not found
long Text::Find_Nth(char q, long n) const
{
     long vez=0, i;
     for (i=0;i<N;i++) 
     {
	  if (D[i]==q) vez++;
	  if (vez==n) break;
     }
     if (i==N) return -1;
     return i;
}

long Text::To_Int() const
{
     if (!D)
     {
	  Error_Flag(Error_IO);
	  return 0;
     }
     errno=0;
     double x=strtod(D,NULL); 
     if (errno) Error_Flag(Error_IO);
     return x;     
}

double Text::To_Real() const
{
     if (!D) 
     {
	  Error_Flag(Error_IO);
	  return 0.0;
     }
     errno=0;
     double x=strtod(D,NULL); 
     if (errno) Error_Flag(Error_IO);
     return x;     
}

void Text::To_LowerCase()
{
     for (long i=0;i<N;i++)
	  D[i]=tolower(D[i]);
}

void Text::To_UpperCase()
{
     for (long i=0;i<N;i++)
	  D[i]=toupper(D[i]);
}

Text Text::Get_Field(long n) const
{
     long count=0;
     bool spacing_old=true;
     long i,j;
     for (i=0;i<N;i++)
     {
          bool spacing=Is_Space(D[i]);
          if (spacing_old && !spacing) count++;
          if (count==n) break;
          spacing_old=spacing;
     }
     if (count!=n)  // error!
     {
          Error_Flag(Error_IO);
          Text Zp;
          return Zp;
     }
     for (j=i;j<N;j++)
     {
          bool spacing=Is_Space(D[i]);
          if (spacing) break;
     }
     return Part(i,j);
}

double Text::Get_Real(long n) const
{
     Text Z=Get_Field(n);
     return Z.To_Real();
}

long Text::Get_Int(long n) const
{
     Text Z=Get_Field(n);
     return Z.To_Int();
}

long Text::Count_Fields() const
{
     long count=0;
     bool spacing_old=true;
     for (long i=0;i<N;i++)
     {
          bool spacing=Is_Space(D[i]);
          if (spacing_old && !spacing) count++;
          spacing_old=spacing;
     }
     return count;
}

Vector Text::To_Vector() const
{
     List S; // starting chars for fields
     bool spacing_old=true;
     for (long i=0;i<N;i++)
     {
          bool spacing=Is_Space(D[i]);
          if (spacing_old && !spacing)
	       S.Append(i);
          spacing_old=spacing;
     }
     long Nfields=S.N;
     S.Append(N);
     Vector V(Nfields);
     for (long k=1;k<=Nfields;k++)
     {
	  double x=Part(S(k),S(k+1)-1).To_Real();
	  V(k)=x;
     }
     return V;
}

List Text::To_List() const
{
     List S; // starting chars for fields
     bool spacing_old=true;
     for (long i=0;i<N;i++)
     {
          bool spacing=Is_Space(D[i]);
          if (spacing_old && !spacing)
	       S.Append(i);
          spacing_old=spacing;
     }
     long Nfields=S.N;
     S.Append(N);
     List R(Nfields);
     for (long k=1;k<=Nfields;k++)
     {
	  long x=Part(S(k),S(k+1)-1).To_Int();
	  R(k)=x;
     }
     return R;
}

bool Is_Space(const char q)
{
     return (q==' ' || q=='\t' || q=='\n');
}
