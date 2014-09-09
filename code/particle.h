
#include <cmath>

class Particle {
public: 
	// Constructors
	Particle();
	Particle(double xv[6]);
	
// The index of a given particle 
	int n;
		
	// Set the time step
	void set_time_step(double time_step){dt = time_step;}
	
	// Evaluate the energy
	void eval_energy();
	
	// Update the state with the Euler method
	void step();
	
	// Update the state with the RK method
	void step_rk();

	// Print the state (position and velocity)
	void print(int format);
	
	// Parameters: time step dt
	double dt;
	
	// Particle mass, state, acceleration, and energy at a given time
	// xv_i is the state at the start point of a given time interval; 
	// this value is used by Runge Kutta method
	double mass;
	double xv_i[6], xv[6]; // position (0:2) and velocity (3:5)
	double energy; // energy
	double time; // Time
	
	// dx1,2,3,4 are the xv steps used by the unge Kutta method
	double dx1[6], dx2[6], dx3[6], dx4[6];
	
	// time derivative of vector [x(3),v(3)]
	double dxvdt[6]; 
	
	// Evaluate the time derivative of vector xv=[x(3),v(3)]
	void eval_dxvdt();
	
};

