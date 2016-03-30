////////////////////////////////////////////////////////
// hvb++ 1.0
// Copyleft: Javier Rodríguez Laguna
// 080725-140328-150127
// Routines for easy connection and usage of X-Windows server. 
// JaviRL Jan 1999. 
// All you have to do is: include "easyx.h" in your header, and use the 
// appropriate compilation order adding -L/usr/X11R6/lib -lX11.  
// large modifications: April 16, 2000; Feb 26, 2004.
// last modification: March 28, 2014
#ifndef EASYX
#define EASYX

#include"EX.h"

EX_Info_Type EX_Info;
EX_Window* EX_CW; // EX Current Window
const char EX_Default_Font[]="-misc-fixed-medium-*-*-*-*-*-*-*-c-100-iso8859-*";

void EX_Start()
{
     EX_Info.display=XOpenDisplay((char*)NULL);
     if (EX_Info.display==(Display*)NULL)
	  Error("I couldn't connect to Xserver");
     EX_Info.screen=DefaultScreen(EX_Info.display);
     EX_Info.root=RootWindow(EX_Info.display,EX_Info.screen);

     EX_Info.visual=DefaultVisual(EX_Info.display,EX_Info.screen);
     EX_Info.depth=DefaultDepth(EX_Info.display,EX_Info.screen);
     EX_Info.colormap=DefaultColormap(EX_Info.display,EX_Info.screen);
     EX_Info.black=BlackPixel(EX_Info.display,EX_Info.screen);
     EX_Info.white=WhitePixel(EX_Info.display,EX_Info.screen);
}

EX_Window* EX_Create_Window(int x, int y, int width, int height)
{
     unsigned long event_mask; 
     event_mask=StructureNotifyMask;

     unsigned long attr_mask=CWEventMask | CWBackPixel | 
	  CWBorderPixel | CWOverrideRedirect;
     
     XSetWindowAttributes attributes;
     attributes.event_mask=event_mask;
     attributes.border_pixel=WhitePixel(EX_Info.display,EX_Info.screen);
     attributes.background_pixel=BlackPixel(EX_Info.display,EX_Info.screen);
     attributes.override_redirect=false;
          
     Window window=XCreateWindow(EX_Info.display, EX_Info.root, 
				 x, y, width, height,
				 2, CopyFromParent, InputOutput,
				 EX_Info.visual, attr_mask, &attributes);
     
     //printf("created window: %d\n",window);
     // Size hints
     XSizeHints *size_hints = XAllocSizeHints();
     size_hints->x = x; 
     size_hints->y = y;
     size_hints->height = height; 
     size_hints->width = width;
     size_hints->min_height = height;
     size_hints->min_width = width;
     size_hints->flags = USPosition | USSize | PMinSize;
     size_hints->base_width = width;
     size_hints->base_height = height;
     size_hints->flags |= PBaseSize;
     XSetWMNormalHints(EX_Info.display, window, size_hints);
     XFree((char*)size_hints);

     // Class hints (name)
     XClassHint class_hints;
     class_hints.res_class=(char*)NULL;
     class_hints.res_name=(char*)NULL;
     XSetClassHint (EX_Info.display, window, &class_hints);
 
     // Window-manager hints
     XWMHints wm_hints;
     wm_hints.flags = InputHint | StateHint;
     wm_hints.initial_state = NormalState;
     wm_hints.input = true;
     XSetWMHints (EX_Info.display, window, &wm_hints);

     XMapRaised(EX_Info.display,window);

     event_mask|=ExposureMask | PointerMotionMask | KeyPressMask | 
     	  StructureNotifyMask | ButtonPressMask | ButtonMotionMask | 
     	  ButtonReleaseMask;

     XSelectInput(EX_Info.display,window,event_mask);
     
     XEvent evt;
     do
     {
      	  XNextEvent( EX_Info.display , &evt );   // calls XFlush
     }while(evt.type != Expose);
     usleep(1);

     XSync(EX_Info.display, false);

     EX_Window* W=(EX_Window*)malloc(sizeof(EX_Window));
     W->window=window;
     W->use_buffer=false;
     W->width=width;
     W->height=height;	  
     EX_CW=W;
     return W;
}

void EX_Set_Name(const char *name)
{
     XStoreName(EX_Info.display, EX_CW->window, name);
     XMapRaised(EX_Info.display, EX_CW->window);
     XFlush(EX_Info.display);
}

// returns the index of the window
EX_Window* EX_Start(int x, int y, int width, int height, const char *name)
{
     EX_Start();
     EX_Window* W = EX_Create_Window(x, y, width, height); 
     XGCValues xgcvalues;
     xgcvalues.foreground=WhitePixel(EX_Info.display,EX_Info.screen);
     xgcvalues.background=BlackPixel(EX_Info.display,EX_Info.screen);
     EX_Info.gc=XCreateGC(EX_Info.display, W->window, 
			 (GCForeground | GCBackground), &xgcvalues);
     EX_CW=W;
//     EX_Prepare_Font(EX_Default_Font);
     return W;
}

// For backwards-compatibility reasons
EX_Window* EX_Start(int x, int y, int width, int height)
{
     return EX_Start(x,y,width,height,"");
}

void EX_Enable_Buffer()
{
     if (EX_CW->use_buffer) return;
     XWindowAttributes wa;
     XGetWindowAttributes(EX_Info.display, EX_CW->window, &wa);
     EX_CW->buffer = XCreatePixmap(EX_Info.display, EX_Info.root,
				  wa.width, wa.height, wa.depth);
     EX_CW->use_buffer=true;
     EX_Clear();
}

void EX_Disable_Buffer()
{
     XFreePixmap(EX_Info.display, EX_CW->buffer);
     EX_CW->use_buffer=false;
}

// returns the width of the WINDOW we're using
void EX_Get_Window_Size(long &width, long &height)
{
     Window rootw;
     int x, y; unsigned int w, h, border, depth;
     XGetGeometry(EX_Info.display,EX_CW->window,&rootw,&x,&y,
		  &w,&h,&border,&depth);
     width=w;
     height=h;
}

void EX_Pixel(int x, int y)
{
     if (EX_CW->use_buffer)
	  XDrawPoint(EX_Info.display,EX_CW->buffer,EX_Info.gc,x,y);
     else
	  XDrawPoint(EX_Info.display,EX_CW->window,EX_Info.gc,x,y);
}

void EX_Flush()
{
     if (EX_CW->use_buffer)
     {
	  XCopyArea(EX_Info.display, EX_CW->buffer, 
		    EX_CW->window, EX_Info.gc,
	   	    0, 0, EX_CW->width, EX_CW->height, 0, 0);
	  // Now, capture the Expose event that will be created
	  XEvent evento;
	  XCheckWindowEvent(EX_Info.display,EX_CW->window,
	  		    ExposureMask,&evento);

     }
     XFlush(EX_Info.display);
}

void EX_Close()
{
//     XFreeFont(EX_Info.display,EX_Info.font);
     XFree(EX_Info.gc);
     EX_Destroy(EX_CW);
     XSetCloseDownMode(EX_Info.display,DestroyAll);
     XCloseDisplay(EX_Info.display);
}

void EX_Destroy(EX_Window *EW)
{
     if (EW->use_buffer)
	  XFreePixmap(EX_Info.display,EW->buffer);
     XDestroyWindow(EX_Info.display,EW->window);
     XFree(EW);
}

void EX_Line(int x0, int y0, int x, int y)
{
     if (EX_CW->use_buffer)
	  XDrawLine(EX_Info.display,EX_CW->buffer,EX_Info.gc,x0,y0,x,y);
     else
	  XDrawLine(EX_Info.display,EX_CW->window,EX_Info.gc,x0,y0,x,y);
}

void EX_Circle(int x, int y, int R)
{
     if (EX_CW->use_buffer)
	  XDrawArc(EX_Info.display,EX_CW->buffer,EX_Info.gc,
		   x-R,y-R,2*R,2*R,0,64*360);
     else
	  XDrawArc(EX_Info.display,EX_CW->window,EX_Info.gc,
		   x-R,y-R,2*R,2*R,0,64*360);
}

void EX_Fill_Circle(int x, int y, int R)
{
     if (EX_CW->use_buffer)
	  XFillArc(EX_Info.display,EX_CW->buffer,EX_Info.gc,
		   x-R,y-R,2*R,2*R,0,64*360);
     else
	  XFillArc(EX_Info.display,EX_CW->window,EX_Info.gc,
		   x-R,y-R,2*R,2*R,0,64*360);
}

// angles in degrees!
void EX_Arc(int x0, int y0, int r, double a0, double a1)
{
     if (EX_CW->use_buffer)
	  XDrawArc(EX_Info.display,EX_CW->buffer,EX_Info.gc,
		   x0-r,y0-r,2*r,2*r,(int)round(a0*64),(int)round(a1*64));
     else
	  XDrawArc(EX_Info.display,EX_CW->window,EX_Info.gc,
		   x0-r,y0-r,2*r,2*r,(int)round(a0*64),(int)round(a1*64));

}

void EX_Fill_Arc(int x0, int y0, int r, double a0, double a1)
{
     if (EX_CW->use_buffer)
	  XFillArc(EX_Info.display,EX_CW->buffer,EX_Info.gc,
		   x0-r,y0-r,2*r,2*r,(int)round(a0*64),(int)round(a1*64));
     else
	  XFillArc(EX_Info.display,EX_CW->window,EX_Info.gc,
		   x0-r,y0-r,2*r,2*r,(int)round(a0*64),(int)round(a1*64));
}

void EX_Fill_Rectangle(int x0, int y0, int x1, int y1)
{
     if (EX_CW->use_buffer)
	  XFillRectangle(EX_Info.display,EX_CW->buffer,EX_Info.gc,
			 x0,y0,x1,y1);
     else
	  XFillRectangle(EX_Info.display,EX_CW->window,EX_Info.gc,
			 x0,y0,x1,y1);
}

void EX_Rectangle(int x0, int y0, int x1, int y1)
{
     if (EX_CW->use_buffer)
	  XDrawRectangle(EX_Info.display,EX_CW->buffer,EX_Info.gc,
			 x0,y0,x1,y1);
     else
	  XDrawRectangle(EX_Info.display,EX_CW->window,EX_Info.gc,
			 x0,y0,x1,y1);
}

void EX_Clear()
{
     if (!EX_CW->use_buffer)
	  XClearWindow(EX_Info.display,EX_CW->window);
     else
     {
	  XWindowAttributes wa;
	  XGetWindowAttributes(EX_Info.display, EX_CW->window, &wa);
	  XSetForeground(EX_Info.display, EX_Info.gc, EX_Info.black);
	  XFillRectangle(EX_Info.display, EX_CW->buffer, EX_Info.gc, 
			 0, 0, wa.width, wa.height);
	  XSetForeground(EX_Info.display, EX_Info.gc, EX_Info.white);
     }
}

void EX_Set_Line_Width(int lw)
{
     XSetLineAttributes(EX_Info.display,EX_Info.gc,lw,
			LineSolid,CapNotLast,JoinMiter);
}

// Maybe the interface is not as efficient as it should...
void EX_Fill_Polygon(const List &X, const List &Y)
{
     XPoint *P=(XPoint*)malloc(X.N*sizeof(XPoint));
     for (long i=0;i<X.N;i++)
     {
	  P[i].x=X(i+1);
	  P[i].y=Y(i+1);
     }
     if (EX_CW->use_buffer)
	  XFillPolygon(EX_Info.display,EX_CW->buffer,EX_Info.gc,
		       P,X.N,Convex,CoordModeOrigin);
     else
	  XFillPolygon(EX_Info.display,EX_CW->window,EX_Info.gc,
		        P,X.N,Convex,CoordModeOrigin);
     free(P);
}

// Maybe the interface is not as efficient as it should...
void EX_Polygon(const List &X, const List &Y)
{
     if (!X.N || !Y.N) return;
     XPoint *P=(XPoint*)malloc(X.N*sizeof(XPoint));
     for (long i=0;i<X.N;i++)
     {
	  P[i].x=X(i+1);
	  P[i].y=Y(i+1);
     }
     if (EX_CW->use_buffer)
	  XDrawLines(EX_Info.display,EX_CW->buffer,EX_Info.gc,
		     P,X.N,CoordModeOrigin);
     else
	  XDrawLines(EX_Info.display,EX_CW->window,EX_Info.gc,
		     P,X.N,CoordModeOrigin);
     free(P);
}


void EX_Get_Event(XEvent *evento)
{
     XNextEvent(EX_Info.display,evento);
}

long EX_Alloc_Named_Color(const char *colorname)
{
     XColor hardwarecolor, exactcolor;
     long color=0;
     int status;
     
     status=XAllocNamedColor(EX_Info.display,EX_Info.colormap,
			     colorname, &hardwarecolor,&exactcolor);
     if (status!=0) color=hardwarecolor.pixel;
     else printf("Error allocating color %s\n",colorname);
     return(color);
}

long EX_Alloc_RGB_Color(double red, double green, double blue)
{
     XColor search_color;
     long color=0;
     int status;
     
     search_color.red=(int)(65535*red);
     search_color.green=(int)(65535*green);
     search_color.blue=(int)(65535*blue);
     status=XAllocColor(EX_Info.display,EX_Info.colormap,&search_color);
     if (status==0) printf("Error allocating RGB color\n");
     else color=search_color.pixel;
     return(color);
}

void EX_Set_Color(long colornum)
{
     XSetForeground(EX_Info.display,EX_Info.gc,colornum);
}

void EX_Color(const char *colorname)
{
     long color=EX_Alloc_Named_Color(colorname);
     EX_Set_Color(color);
}

void EX_Color(double R, double G, double B)
{
     long color=EX_Alloc_RGB_Color(R,G,B);
     EX_Set_Color(color);
}

char EX_Key_Pressed() // if a key is pressed, return it, otherwise return 0
{
     // XEvent evento;
     // bool pressed;
     // pressed=XCheckWindowEvent(EX_Info.display,EX_CW->window,
     // 			       KeyPressMask,&evento);
     // XSync(EX_Info.display,true); // discard events
     // if (!pressed) return 0;
     // KeySym tecla=EX_Key_2_Keysym(&evento);
     // return (char)tecla;
     XEvent evento;
     int pressed;
     pressed=XCheckWindowEvent(EX_Info.display,EX_CW->window, 
			       KeyPressMask, &evento);
     if (!pressed) return 0;
     KeySym tecla=EX_Key_2_Keysym(&evento);
     return (char)tecla;

}    

char EX_Read_Key() // Wait until a key is pressed, return it
{
     XEvent evento;
     KeySym tecla;
     
     do
     {
	  XNextEvent(EX_Info.display,&evento);
     }while(evento.type!=KeyPress);
     
     tecla=EX_Key_2_Keysym(&evento);
     return (char)tecla;
}

int EX_Pointer() // if the pointer has moved or clicked, returns 1
// Should be followed by a EX_Read_Pointer
{
     XEvent evento;
     int pressed;
     pressed=XCheckWindowEvent(EX_Info.display,EX_CW->window, 
			       ButtonPressMask | 
			       ButtonReleaseMask | ButtonMotionMask, &evento);
     if (pressed) XPutBackEvent(EX_Info.display,&evento);
     return pressed;
}

int EX_Pointer_Pressed() // if the pointer has been clicked, returns 1
// Should be followed by a EX_Read_Pointer
{
     XEvent evento;
     int pressed=XCheckWindowEvent(EX_Info.display,EX_CW->window,
				   ButtonPressMask,&evento);
     if (pressed) XPutBackEvent(EX_Info.display, &evento);
     return pressed;
}
     	
EX_Pointer_State EX_Read_Pointer()
{
     XEvent evento;
     EX_Pointer_State C;
     int button=0;
     
     do
     {	
	  XNextEvent(EX_Info.display,&evento);
     }while(evento.type!=MotionNotify && evento.type!=ButtonPress && 
	    evento.type!=ButtonRelease && evento.xany.window!=EX_CW->window);
     if (evento.type==MotionNotify)
     {
	  C.x=evento.xmotion.x;
	  C.y=evento.xmotion.y;
	  C.button=button;
     }
     if (evento.type==ButtonPress)
     {
	  C.x=evento.xbutton.x;
	  C.y=evento.xbutton.y;
	  C.button=evento.xbutton.button;
	  button=C.button;
     }
     if (evento.type==ButtonRelease)
     {
	  C.x=evento.xbutton.x;
	  C.y=evento.xbutton.y;
	  button=0;
	  C.button=0;
     }
     XFlush(EX_Info.display);
     return C;
}

void EX_Prepare_Font(const char *font_id)
{
     EX_Info.font=XLoadQueryFont(EX_Info.display,font_id);     
     if (EX_Info.font == (XFontStruct*)NULL)
     {
	  /* if this fails, go to nofont, which is the default */
	  fprintf(stderr,"Error loading font [%s]\n",font_id);
	  EX_Info.font=XLoadQueryFont(EX_Info.display,"fixed");
	  if (EX_Info.font ==(XFontStruct *)NULL)
	  {
	       fprintf(stderr,"Error loading font fixed\n");
	       XCloseDisplay(EX_Info.display);
	       exit(1);
	  }
     }
     XSetFont(EX_Info.display,EX_Info.gc,EX_Info.font->fid);
}

void EX_Draw_String(int x, int y, const char *cadena)
{
     if (EX_CW->use_buffer)
	  XDrawString(EX_Info.display,EX_CW->buffer,EX_Info.gc,
		      x,y,cadena,strlen(cadena));
     else
	  XDrawString(EX_Info.display,EX_CW->window,EX_Info.gc,
		      x,y,cadena,strlen(cadena));
}

long EX_Text_Ascent()
{
     return EX_Info.font->ascent;
}

long EX_Text_Descent()
{
     return EX_Info.font->descent;
}

long EX_Text_Width(const char *cadena, long n)
{
     long nr=n;
     if (!n) nr=strlen(cadena);
     return XTextWidth(EX_Info.font, cadena, nr);
}

XImage* EX_Get_Image(int x0, int y0, int xsize, int ysize)
{
     if (EX_CW->use_buffer)
	  return XGetImage(EX_Info.display,EX_CW->buffer,
			   x0,y0,xsize,ysize,AllPlanes,XYPixmap);
     else
	  return XGetImage(EX_Info.display,EX_CW->window,
			   x0,y0,xsize,ysize,AllPlanes,XYPixmap);
}

void EX_Put_Image(XImage* imagen, int x0, int y0, int xsize, int ysize)
{
     if (EX_CW->use_buffer)
	  XPutImage(EX_Info.display,EX_CW->buffer,EX_Info.gc,
		    imagen,0,0,x0,y0,xsize,ysize);
     else
	  XPutImage(EX_Info.display,EX_CW->window,EX_Info.gc,
		    imagen,0,0,x0,y0,xsize,ysize);
}

List EX_Palette(double R1, double G1, double B1,
		double R2, double G2, double B2, long Ncolors)
{
     List Pal(Ncolors);
     for (long i=1;i<=Ncolors;i++)
     {
	 double s=(double)i/(double)Ncolors;
	 double rojo=R1+s*(R2-R1);
	 double verde=G1+s*(G2-G1);
	 double azul=B1+s*(B2-B1);
	 Pal(i)=EX_Alloc_RGB_Color(rojo,verde,azul);
     }
     return Pal;
}     

KeySym EX_Key_2_Keysym (XEvent *event)
{
     XComposeStatus compose;
     KeySym keysym;
     XKeyEvent *keyevent;
     char cad[20];
     
     keyevent=(XKeyEvent*) event;
     XLookupString(keyevent,cad,19,&keysym,&compose);
     return keysym;
}

// Set the mode for the graphics context, i.e.: how to "mix" new bits coming
// from the painting orders with the old ones, already in the screen
// 3 -> COPY, usual one
// 7 -> OR, 
void EX_Mode(int mode)
{
     XGCValues xgcvalues;
     xgcvalues.function=mode;
     XChangeGC(EX_Info.display,EX_Info.gc,GCFunction,&xgcvalues);
//     XFlushGC(EX_Info.display,EX_Info.gc);
}

void EX_Mode_Or()
{
     EX_Mode(7);
}

void EX_Mode_Copy()
{
     EX_Mode(3);
}
	 
#endif
