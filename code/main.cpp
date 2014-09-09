

#include <iostream>
using namespace std;

#include "particle_system.h"
#include "clock.h"

int main() {
//
// Three bodies simulation 
//
	// Time duration of this run
	Clock run_time;

	// 	Define three particles in phase space and print position and velocity
	ParticleSystem* particles; 
	particles = new ParticleSystem();
	double dt = 0.01;
	particles->set_time_step(dt);
	
	// Open the output file
	particles->open_output();
	
	// Print some basic infos
	particles->printInfo();
	
	// Print the header for particle list and the initial state
	// particles->printHeader();
	particles->printf();
	
	// Run

	run_time.start();
	for (int i=0; i < 400; i++) {
		particles->step_rk();
		particles->printf(); // Print 
	}
	run_time.stop();
	run_time.print();
	
	// Close the output file
	particles->close_output();
	
	delete particles;

	return 0;
}
