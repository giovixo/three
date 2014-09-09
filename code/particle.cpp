


#include <iostream>
using namespace std;

#include "particle.h"

#include <cmath>

Particle::Particle() {

//  Initial state
//  ==============	

	xv[0] = 1.;
	xv[1] = 0.;
	xv[2] = 0.;
	
	xv[3] = 0.;
	xv[4] = 0.;
	xv[5] = 0.;
	
	
	// Time step
	time = 0;
	dt = 0.1;
}

Particle::Particle(double xv0[6]) {

//  Initial state
//  ==============
	for (int i=0; i<6; i++)
		xv[i] = xv0[i];

	// Time step
	time = 0;
	dt = 0.1;
}

void Particle::print(int format) {
	if (format == 0) {
		cout << time << " " << xv[0] << " " << xv[1] << " " << xv[2];
		cout << " " << xv[3] << " " << xv[4] << " " << xv[5];
		cout << " " << energy << endl;
	} else {
		cout << "Time: " << time;
		cout << "   Position: " << xv[0] << " " << xv[1] << " " << xv[2];
		cout << "   Velocity: " << xv[3] << " " << xv[4] << " " << xv[5];
		cout << "   Energy: " << energy << endl;
	}
}

void Particle::eval_energy() {
	energy = 0.5*(xv[3]*xv[3] + xv[4]*xv[4] + xv[5]*xv[5]) 
	         - 1./sqrt(xv[0]*xv[0] + xv[1]*xv[1] + xv[2]*xv[2]);
}

void Particle::eval_dxvdt() {

	// Evaluate the acceleration (dxvdt[3:5])
	// ======================================
	//
	// Change this funtion to set another force field

	dxvdt[0] = xv[3]; // dx/dt = v = xv[3:5]
	dxvdt[1] = xv[4];
	dxvdt[2] = xv[5];
	
	
	// -> First estimate the versor r directed to the center
	double r[3];
	double r2;
	r2 = xv[0]*xv[0] + xv[1]*xv[1] + xv[2]*xv[2];
	for (int i=0; i<3; i++)
		r[i] = xv[i]/sqrt(r2);
	// -> Then estimate the force (acceleration)
	dxvdt[3] = -r[0]/r2;
	dxvdt[4] = -r[1]/r2;
	dxvdt[5] = -r[2]/r2;
}


void Particle::step() {

	eval_dxvdt();

	for (int i=0; i<6; i++) 
		xv[i] += dxvdt[i]*dt;
	
	eval_energy();
	time += dt;
}



void Particle::step_rk() {
	
	for (int i=0; i<6; i++) { 
		xv_i[i] = xv[i];
	}
	
	eval_dxvdt();
	for (int i=0; i<6; i++) {
		dx1[i] = dt*dxvdt[i];
		xv[i] = xv_i[i] + 0.5*dx1[i];
	}

	eval_dxvdt();
	for (int i=0; i<6; i++) {
		dx2[i] = dt*dxvdt[i];
		xv[i] = xv_i[i] + 0.5*dx2[i];
	}
	
	eval_dxvdt();
	for (int i=0; i<6; i++) {
		dx3[i] = dt*dxvdt[i];
		xv[i]  = xv_i[i] + dx3[i];
	}
	
	eval_dxvdt();
	for (int i=0; i<6; i++) {
		dx4[i] = dt*dxvdt[i];
	}
	
	for (int i=0; i<6; i++) {
		xv[i] = xv_i[i] + (dx1[i] + 2*dx2[i] + 2*dx3[i] + dx4[i])/6.0;
	}
	
	eval_energy();
	time += dt;
}



