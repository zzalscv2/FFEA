gsl=-I "/home/DS/pyrar/gsl/usr/include" -L "/home/DS/pyrar/gsl/usr/lib/" -lgsl -lgslcblas
glut =-I "./" -L . -lglut -lGLU -lGL
opt = -lm -lpthread -finline-functions -O2 -Wall -Winline
gtk = -I "./gtk-2.0/" -lgtk -lgdk

all: wd make_cubic_mesh

wd: walrus_display.c 
	cc -g walrus_display.c $(opt) $(glut) -o wd

make_cubic_mesh: make_cubic_mesh.c
	cc make_cubic_mesh.c -o make_cubic_mesh

gtk_hello_world: gtk_hello_world.c
	cc gtk_hello_world.c $(gtk) -o gtk_hello_world

clean:
	rm -f *.o
