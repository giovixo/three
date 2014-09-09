

#include <iostream>
#include <ctime>
using namespace std;

#include "particle_system.h"

#define PROGRAM_TITLE           "Three bodies simulation"
#define WINDOW_SIZE             696
#define FRAME_RATE_SAMPLES      50

// Function headers
void Display();
void Idle();

// Break time function
void breakTime( int seconds);

// Global variable (three particles in phase space)
ParticleSystem PSystem;


int main(int argc, char *argv[]) {
//
// Three bodies simulation
//

	double dt = 5.;
	PSystem.set_time_step(dt);
	PSystem.printInfo();

	// Start GLUT and negotiate a window with the operating system
	glutInit(&argc, argv);
	// Select RGB color scheme and use double buffering (see Wikipedia)
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
	// Set initial window size
	glutInitWindowSize(WINDOW_SIZE, WINDOW_SIZE);
	// Create GLUT window
	glutCreateWindow(PROGRAM_TITLE);
	
	// We will be working in the [-1.0, 1.0]x[-1.0, 1.0] square
	gluOrtho2D(-10.0, 10.0, -10.0, 10.0);
	
	// Tell GLUT witch functions to call to do the drawing, keyboard
	// management, etc
    glutDisplayFunc(Display);
    glutIdleFunc(Idle);
    
    // Pass control to GLUT
	glutMainLoop();

	return 0;
}


void Display() {       
        // Clear the screen
        glClear(GL_COLOR_BUFFER_BIT);

        // Make legal area white
        float x_max, x_min, y_max, y_min;

        x_max =  10.; 
        x_min = -10.;
        y_max =  10.;
        y_min = -10.;

        glColor3f(1.0, 1.0, 1.0);
        glBegin(GL_QUADS);
                glVertex2f(x_max, y_max);
                glVertex2f(x_min, y_max);
                glVertex2f(x_min, y_min);
                glVertex2f(x_max, y_min);
        glEnd();

        // Draw Particle System 
        // glColor3f(1.0, 0.0, 0.);
        PSystem.Draw();

        // We finished drawing. Display what's on this buffer now.
        glutSwapBuffers();
}

void Idle() {
	// PSystem.print(); // For test
	PSystem.step_rk();	
	// Call the display function 
	glutPostRedisplay();
}

void breakTime( int seconds)
{
    clock_t temp;
    temp = clock () + seconds * CLOCKS_PER_SEC ;
    while (clock() < temp) {}
}


