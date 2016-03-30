////////////////////////////////////////////////////////
// hvb++ 1.0
// Copyleft: Javier Rodríguez Laguna
// Simple driver to create PostScript files from C
// 0004-0208-0505-0607-0807-1501-1505
#ifndef EASYPS_HEADER
#define EASYPS_HEADER
#include"Matrix.h"

extern FILE *PS_Current_File;
FILE *PS_Open(const char *name, int x0, int y0, int x1, int y1);

void PS_Close(FILE *psfile);

void PS_Translate(FILE *psfile, double x, double y);
void PS_Translate(double x, double y);

void PS_Line(FILE *psfile, double x0, double y0, double x1, double y1);
void PS_Line(double x0, double y0, double x1, double y1);

void PS_Circle(FILE *psfile, double cx0, double cy0, double R);
void PS_Circle(double cx0, double cy0, double R);
     
void PS_Fill_Circle(FILE *psfile, double cx0, double cy0, double R);
void PS_Fill_Circle(double cx0, double cy0, double R);

void PS_Rectangle(FILE *psfile, double x0, double y0, double wx, double wy);
void PS_Rectangle(double x0, double y0, double wx, double wy);

void PS_Fill_Rectangle(FILE *psfile, double x0, double y0, double wx, double wy);
void PS_Fill_Rectangle(double x0, double y0, double wx, double wy);

void PS_Polygon(FILE *psfile, const Vector &X, const Vector &Y);
void PS_Polygon(const Vector &X, const Vector &Y);

void PS_Fill_Polygon(FILE *psfile, const Vector &X, const Vector &Y);
void PS_Fill_Polygon(const Vector &X, const Vector &Y);

void PS_Normal_Line(FILE *psfile);
void PS_Normal_Line();

void PS_Dashed_Line(FILE *psfile);
void PS_Dashed_Line();

void PS_Set_Line_Width(FILE *psfile, double f);
void PS_Set_Line_Width(double f);

void PS_Color(FILE *psfile, double red, double green, double blue);
void PS_Color(double red, double green, double blue);

extern const char *PS_Default_Font;

// Remember: angle in degrees!!!!
void PS_Arc(FILE *psfile, double cx0, double cy0, double R,
	    double a0, double a1);
void PS_Arc(double cx0, double cy0, double R,double a0, double a1);

/* void PS_Curve(FILE *psfile, double x1, double y1, double x2, double y2, */
/* 	      double x3, double y3, double x4, double y4); */
/* void PS_Curve(double x1, double y1, double x2, double y2, */
/* 	      double x3, double y3, double x4, double y4); */

void PS_Prepare_Font(FILE *psfile, const char *fontname, int size);
void PS_Prepare_Font(const char *fontname, int size);

void PS_Text(FILE *psfile, double x0, double y0, const char *S);
void PS_Text(double x0, double y0, const char *S);

void PS_Box_Text(FILE *psfile, double x0, double y0, const char *S, double b);
void PS_Box_Text(double x0, double y0, const char *S, double b);

#endif


