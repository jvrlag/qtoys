// 050511, first version, used Imlib11
// 100914, second version, uses Imlib2
// 141024, exporting/importing RGB matrices
// Works with a "current image"
#ifndef EASYIM_HEADER
#define EASYIM_HEADER

#include"EX.h"
#include"Matrix.h"
#include"Imlib2.h"

typedef Imlib_Image EI_Image;

void EI_Start();
EI_Image EI_Load(const char *name); // also sets current image
int EI_Get_Width();
int EI_Get_Height();
void EI_Render(int x, int y);
void EI_Render_Scaled(int x, int y, int w, int h);
void EI_Free();
EI_Image EI_Capture(int x, int y, int w, int h); // also sets current image
void EI_Save(const char *name);

// Save the current image into matrices R,G,B and (alpha channel) A
void EI_To_Matrices(Matrix &R, Matrix &G, Matrix &B, Matrix &A);
// Put matrices (R,G,B) and alpha-channel A as current image
void EI_From_Matrices(const Matrix &R, const Matrix &G, const Matrix &B,
		      const Matrix &A);

#endif
