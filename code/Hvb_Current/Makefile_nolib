# Hvb generic makefile

hvb=$(HOME)/lib/Hvb_Current/

FLAGS=-pedantic -Wall -g -I$(hvb) -O3 -DDEBUG
LIBS=-llapack -lblas -lm -L/usr/X11R6/lib -lX11 -lImlib2 

HVBS= $(hvb)Matrix.o $(hvb)Text.o $(hvb)Common.o $(hvb)MatrixC.o \
$(hvb)EX.o $(hvb)Image.o $(hvb)Postscript.o $(hvb)Optimize.o \
$(hvb)Graph.o $(hvb)Manybody.o $(hvb)Calculus.o

OBJECTS = 

objects: $(HVBS) $(OBJECTS)

%.o: %.cc
	g++ -c $< $(FLAGS) -o $@

x%: x%.o objects
	g++ $(FLAGS) -o $@ $< $(OBJECTS) $(HVBS) $(LIBS)

clean:
	ls x* | grep -v "\." | xargs rm -f
	rm *.o 

mrproper: 	
	rm *.o $(hvb)/*.o

