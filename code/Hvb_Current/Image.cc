// 100914, easyim
#include"Image.h"

void EI_Start()
{
     imlib_set_cache_size(2048 * 1024);
     imlib_context_set_display(EX_Info.display);
     imlib_context_set_visual(EX_Info.visual);
     imlib_context_set_colormap(EX_Info.colormap);
}

EI_Image EI_Load(const char *name) // also sets current image!
{
     EI_Image image;
     image = imlib_load_image(name);
     imlib_context_set_image(image);
     return image;
}

int EI_Get_Width()
{
     return imlib_image_get_width();
}

int EI_Get_Height()
{
     return imlib_image_get_height();
}

void EI_Render(int x, int y)
{
     if (EX_CW->use_buffer)
	  imlib_context_set_drawable(EX_CW->buffer); 
     else
	  imlib_context_set_drawable(EX_CW->window); 
     imlib_render_image_on_drawable(x,y);
}

void EI_Render_Scaled(int x, int y, int w, int h)
{
     if (EX_CW->use_buffer)
	  imlib_context_set_drawable(EX_CW->buffer); 
     else
	  imlib_context_set_drawable(EX_CW->window); 
     EI_Image nueva=imlib_create_image(w,h);
     EI_Image antigua=imlib_context_get_image();
     int oldw=EI_Get_Width();
     int oldh=EI_Get_Height();
     imlib_context_set_image(nueva);
     imlib_blend_image_onto_image(antigua,0,0,0,oldw,oldh,0,0,w,h);
     imlib_render_image_on_drawable(x,y);
     imlib_free_image();
     imlib_context_set_image(antigua);
}

void EI_Free()
{
    imlib_free_image();
}

EI_Image EI_Capture(int x, int y, int w, int h) // also sets current image
{
     if (EX_CW->use_buffer)
	  imlib_context_set_drawable(EX_CW->buffer); 
     else
	  imlib_context_set_drawable(EX_CW->window); 
     EI_Image image=imlib_create_image_from_drawable(0,0,0,w,h,1);
     imlib_context_set_image(image);
     return image;
}

void EI_Save(const char *name)
{
     imlib_save_image(name);
}

// Current image to matrices
void EI_To_Matrices(Matrix &R, Matrix &G, Matrix &B, Matrix &A)
{
     long width=EI_Get_Width();
     long height=EI_Get_Height();
     unsigned char* data=(unsigned  char*)imlib_image_get_data();
     long size=width*height;
     A.Create(width,height);
     R.Create(width,height); 
     G.Create(width,height); 
     B.Create(width,height);
     for (long i=0;i<size;i++)
     {
	  long y=i/width+1;
	  long x=i%width+1;
	  B(x,y)=(double)data[4*i]/256.0;
	  G(x,y)=(double)data[4*i+1]/256.0;
	  R(x,y)=(double)data[4*i+2]/256.0;
	  A(x,y)=(double)data[4*i+3]/256.0;
	  if (B(x,y)<0.0) printf("Caution! %d -> %g\n",data[4*i],B(x,y));
     }
}

// tranform x \in [0,1] into char in [0,255]
char Double_To_Char(double x)
{
     long i=(long)round(x*256.0);
     if (i<0) i=0;
     if (i>255) i=255;
     return (char)i;
}

// Matrices to current image
void EI_From_Matrices(const Matrix &R, const Matrix &G, const Matrix &B,
		      const Matrix &A)
{
     long width=R.N1;
     long height=R.N2;
     long size=width*height;
     char* data=(char*)malloc(size*4);
     for (long i=0;i<size;i++)
     {
	  long y=i/width+1;
	  long x=i%width+1;
	  data[4*i]=Double_To_Char(B(x,y));
	  data[4*i+1]=Double_To_Char(G(x,y));
	  data[4*i+2]=Double_To_Char(R(x,y));
	  data[4*i+3]=Double_To_Char(A(x,y));
     }
     Imlib_Image im=imlib_create_image_using_data(width,height,(DATA32*)data);
     imlib_context_set_image(im);     
}

