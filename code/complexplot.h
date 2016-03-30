// Complex plot,
// Routines in order to show complex functions, with a color code
// 110912-160327

#include"MatrixC.h"
#include"EX.h"

// returns a different "top" color, periodic-continuous in alpha
// alpha in [0,1] !!!
Vector Get_Color_From_Phase(double alpha);

long Find_Color_Index(cmplx z, long Nint, long Nphase, double Rmax);

List Build_Palette_With_Phases(long Nc_int, long Nc_phase, bool white);

void Complex_Plot(const MatrixC &V, long Nc_int, long Nc_phase, 
		  long xsize, bool white, double saturation);
