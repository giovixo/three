#
# Makefile for OS X
#

# Compiler
CPP = g++ -Wall -Wno-deprecated-declarations

all: three

three: obj/main.o obj/particle.o obj/particle_system.o obj/clock.o 
	$(CPP) -o three obj/main.o obj/particle.o obj/particle_system.o obj/clock.o -framework OpenGL -framework GLUT

three++: obj/main_gl.o obj/particle.o obj/particle_system.o obj/clock.o 
	$(CPP) -o three++ obj/main_gl.o obj/particle.o obj/particle_system.o obj/clock.o -framework OpenGL -framework GLUT

obj/main.o: code/main.cpp
	$(CPP) -c code/main.cpp -o obj/main.o

obj/main_gl.o: code/main_gl.cpp
	$(CPP) -c code/main_gl.cpp -o obj/main_gl.o

obj/particle.o: code/particle.cpp
	$(CPP) -c code/particle.cpp -o obj/particle.o

obj/particle_system.o: code/particle_system.cpp
	$(CPP) -c code/particle_system.cpp -o obj/particle_system.o

obj/clock.o: code/clock.cpp
	$(CPP) -c code/clock.cpp -o obj/clock.o
	
clean: 
	rm -f obj/*o
