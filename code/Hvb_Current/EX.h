////////////////////////////////////////////////////////
// hvb++
// Copyleft: Javier Rodríguez Laguna
// EasyX library, part of HVB
// JaviRL, 2001-2006-2015
#ifndef EASYX_HEADER
#define EASYX_HEADER
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/keysym.h>
#include "Common.h"


typedef struct
{
     Display *display;
     int screen;
     Window root;
     Visual *visual;
     int depth;
     GC gc;
     XFontStruct *font;
     Colormap colormap;
     long black, white;
} EX_Info_Type;

typedef struct
{
     bool use_buffer;
     Window window;
     Pixmap buffer;
     long width, height;
}EX_Window;

extern EX_Info_Type EX_Info;
extern EX_Window* EX_CW; // EX Current Window

void EX_Start();
EX_Window* EX_Start(int x, int y, int width, int height);
EX_Window* EX_Create_Window(int x, int y, int width, int height);
void EX_Set_Name(const char *name);
void EX_Destroy(EX_Window*);
void EX_Close();

void EX_Enable_Buffer();
void EX_Disable_Buffer();

void EX_Get_Window_Size(long &width, long &height);

void EX_Pixel(int x, int y);
void EX_Flush();
void EX_Line(int x0, int y0, int x, int y);
// angles in degrees!
void EX_Arc(int x0, int y0, int r, double a0, double a1);
void EX_Fill_Arc(int x0, int y0, int r, double a0, double a1);
void EX_Circle(int x, int y, int R);
void EX_Fill_Circle(int x, int y, int R);
void EX_Rectangle(int x0, int y0, int w, int h);
void EX_Fill_Rectangle(int x0, int y0, int w, int h);
void EX_Polygon(const List &X, const List &Y);
void EX_Fill_Polygon(const List &X, const List &Y);
void EX_Clear();

void EX_Set_Line_Width(int lw);

void EX_Get_Event(XEvent *);
KeySym EX_Key_2_Keysym (XEvent *event);
char EX_Key_Pressed();
char EX_Read_Key();

// Struct used to return the position where pointer button was pressed
typedef struct
{
     int x;
     int y;
     int button;
} EX_Pointer_State;

int EX_Pointer();
int EX_Pointer_Pressed();
EX_Pointer_State EX_Read_Pointer();

extern const char EX_Default_Font[];
void EX_Prepare_Font(const char *font_id);
void EX_Draw_String(int x, int y, const char *cadena);
long EX_Text_Width(const char *cadena, long n);
long EX_Text_Ascent();
long EX_Text_Descent();

XImage* EX_Get_Image(int x0, int y0, int xsize, int ysize);
void EX_Put_Image(XImage* imagen, int x0, int y0, int xsize, int ysize);

void EX_Color(const char *colorname);
void EX_Color(double R, double G, double B);
void EX_Set_Color(long colornum);
long EX_Alloc_Named_Color(const char *colorname);
long EX_Alloc_RGB_Color(double red, double green, double blue);
List EX_Palette(double R1, double G1, double B1,
		double R2, double G2, double B2, long Ncolors);

void EX_Mode(int mode);
void EX_Mode_Or();
void EX_Mode_Copy(); 
	 
#endif
