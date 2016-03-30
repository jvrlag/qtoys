// Parsing routines
// 111130
#ifndef TEXT_H
#define TEXT_H
#include"Common.h"
#include"Matrix.h"

class Text
{
public:
     long N; // number of chars, not including the final 0
     char *D;
     Text();
     Text(const char*);
     Text(const Text&);
     ~Text();
     void Start();
     void Create(const char*);
     void Load(char*);
     void Destroy();

     long Get_Line(FILE*); // read a line from a file, ret -1 is failure
     void Write() const;

     void Append(const char*);
     void Append(const char*, long);
     void Append(const Text&);
     void Append(const char);

     void Append_F(const char*, long i);
     void Append_F(const char*, double f);

     Text& operator=(const Text&);
     Text& operator=(const char*);

     void Strip_Blanks();

     bool Is_Here(const char *);
     bool Is_Here(const Text &);
     bool Is_There(const char *);
     bool Is_There(const Text &);

     Text Token(char q0, char q1, long n) const;
     Text Part(long i0, long i1) const;
     long   Get_Int(char q0, char q1, long n) const;
     double Get_Real(char q0, char q1, long n) const;

     long   Count(char q) const;
     long   Find_Nth(char q, long n) const;
     long   To_Int() const;
     double To_Real() const;
     
     void To_LowerCase(); 
     void To_UpperCase();

     Text   Get_Field(long n) const;
     double Get_Real(long n) const;
     long   Get_Int(long n) const;
     long   Count_Fields() const;
     Vector To_Vector() const;
     List   To_List() const;
};

void Copy(Text &S, const Text &S1);
void Copy(Text &S, const char *z);
void Copy(Text &S, const char *z, long n);

ostream& operator<< (ostream &out, Text &S);
bool operator==(const Text &S1, const Text &S2);

bool Is_Space(const char q);

#endif
