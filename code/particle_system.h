
#include <cmath>
#include <vector>
#include <GLUT/glut.h>
#include "particle.h"

#include <fstream>

// Macros
        #define FOR_EACH_PARTICLE       for (unsigned int k=0; k<particles.size(); k++) {
        #define END_FOR_EACH_PARTICLE   }

class ParticleSystem {
public: 
	// Constructor
	ParticleSystem();
	
	// Set the time step
    void set_time_step(double time_step){dt = time_step;}
			
	// Step forward by dt (Runge-Kutta method) 
	void step_rk();
	
	// Step forward by dt (Euler method)
	void step_euler();

	// Print some infos on the particles system 
	void printInfo(); 	
	// Print the position and velocity at a given time
	void print();
	// Print the header before print() [OBSOLETE]
	void printHeader();
	
	// Open the output file
	void open_output();
	// Print on file the position and velocity at a given time
	void printf();
	// Close the output file
	void close_output();
	
	// Draw the particles with GLUT library
	void Draw();

private:
	// Define a vector of points
	std::vector< Particle > particles;
	
	// Position of the center of mass (obsolete)
	double cm[3];
	
	// Position and velocity of the center of mass 
	double x_cm, y_cm, z_cm;
	double vx_cm, vy_cm, vz_cm;
	
	// The energy
	double energy;
	
	double time, dt; // Time and time step 
	
	// Initial condition (positions and velocities)
	void initialCondition();
	
	// Estimate d(xv)/dt at given time
	void eval_dxvdt();
	
	// Estimate the force between particle i and j
	void eval_force(double force[3], int pi, int pj);
	
	// Estimate the energy
    void eval_energy();
    void eval_pot(double *pen, int pi, int pj);
    
    // Evaluate the center of mass position
    void eval_cm();
    
    // Pointer to the output file
    ofstream file_w; 
};

