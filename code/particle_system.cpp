


#include <iostream>
#include <fstream>
#include <cassert>

using namespace std;

#include <cmath>
#include "particle_system.h"

// K2 = Gaussian gravitational constant squared
#define K2 2.959122082855911e-4

ParticleSystem::ParticleSystem() {

//  Initial state
//  ==============	
	particles.clear();
	particles.resize(3); // Three particles
	//cout << "Number of particles: " <<  particles.size() << endl; Test only
	
	unsigned int i;
	for (i=0; i<particles.size(); i++) particles[i].n = i; // assign the particle index
	
	initialCondition(); // Assign the initial conditions in phase space
	time = 0; // Set the initial time
	dt = 0.04;       // Set the time step
}

void ParticleSystem::printInfo() {
	cout << endl << " =========== INITIAL STATE (in the reference of center of mass) ===============" <<  endl << endl;

	cout << " No. particles: " << particles.size() << endl;
	cout << " Total mass: " << particles[0].mass + particles[1].mass + particles[2].mass << endl; 
	
	for (unsigned int i=0; i<particles.size(); i++) {
		cout << endl;
		cout << " Particle index: " << particles[i].n << endl;
		
		cout << " Mass: " << particles[i].mass << endl;
		
		cout << " Position: ";
		for (int j=0; j<3; j++)  cout << particles[i].xv[j] << " ";
		cout << endl;
		
		cout << " Velocity: ";
		for (int j=3; j<6; j++)  cout << particles[i].xv[j] << " ";
		cout << endl;			
	}
		
	cout << endl << " ==============================================================================" <<  endl << endl;
}

void ParticleSystem::initialCondition() {
	
	Particle *pcurr;
	double m, x[3], v[3];
	
	// Read the initial condition (mass, position, velocity) from file init.dat 
	ifstream read_file("init.dat");
	assert(read_file.is_open()); // Check if file exists
	
	string dummyLine;	
	getline(read_file, dummyLine);

	FOR_EACH_PARTICLE
		pcurr = &(particles[k]);
		read_file >> m >> x[0] >> x[1] >> x[2] >> v[0] >> v[1] >> v[2];
		pcurr->mass = m;
		pcurr->xv[0] = x[0];
		pcurr->xv[1] = x[1];
		pcurr->xv[2] = x[2]; 
		pcurr->xv[3] = v[0];
		pcurr->xv[4] = v[1];
		pcurr->xv[5] = v[2];
	END_FOR_EACH_PARTICLE
	
	// Center of mass: position and velocity
	x_cm = 0.; y_cm = 0.; z_cm = 0.;
	vx_cm = 0.; vy_cm = 0.; vz_cm = 0.;
	double total_mass( particles[0].mass + particles[1].mass + particles[2].mass );
	FOR_EACH_PARTICLE
		x_cm = x_cm +  ( particles[k].mass * particles[k].xv[0] ) / total_mass;
		y_cm = y_cm +  ( particles[k].mass * particles[k].xv[1] ) / total_mass;
		z_cm = z_cm +  ( particles[k].mass * particles[k].xv[2] ) / total_mass;
		vx_cm = vx_cm +  ( particles[k].mass * particles[k].xv[3] ) / total_mass;
		vy_cm = vy_cm +  ( particles[k].mass * particles[k].xv[4] ) / total_mass;
		vz_cm = vz_cm +  ( particles[k].mass * particles[k].xv[5] ) / total_mass;
	END_FOR_EACH_PARTICLE

	// Convert coordinates and velocity in the center of mass reference
	FOR_EACH_PARTICLE
		pcurr = &(particles[k]);
		pcurr->xv[0] = pcurr->xv[0] - x_cm; 
		pcurr->xv[1] = pcurr->xv[1] - y_cm;
		pcurr->xv[2] = pcurr->xv[2] - z_cm;
		pcurr->xv[3] = pcurr->xv[3] - vx_cm;
		pcurr->xv[4] = pcurr->xv[4] - vy_cm;
		pcurr->xv[5] = pcurr->xv[5] - vz_cm;
	END_FOR_EACH_PARTICLE
	
	read_file.close();
	
}

void ParticleSystem::printHeader() {
	cout << "time   p_index  mass   x1   x2   x3   v1   v2   v3 energy x_cm y_cm z_cm" << endl;
}

void ParticleSystem::print() {
	
	Particle *pcurr;
	
	eval_energy();
	eval_cm();
	
	FOR_EACH_PARTICLE
		pcurr = &(particles[k]);
		cout << time << " " << pcurr->n << " " << pcurr->mass << " "; 
		for (int i=0; i<6; i++) cout << pcurr->xv[i] << " ";
		cout << energy << " ";
		for (int i=0; i<3; i++) cout << cm[i] << " ";
		cout << endl;
	END_FOR_EACH_PARTICLE
}

/*
void ParticleSystem::eval_dxvdt() {
//
//   Evaluate dxv/dt (the force) on each particles
// ================================================
//-> <Force on particles 0 by particle 1
	double force_01[3], r[3];
	double dist2=0; // Square distance
	for (int i=0; i<3; i++) dist2 += (particles[0].xv[i] - particles[1].xv[i])*
	(particles[0].xv[i] - particles[1].xv[i]);
	for (int i=0; i<3; i++) {
		r[i] = (particles[0].xv[i] - particles[1].xv[i])/sqrt(dist2); // versor r01
		force_01[i] = -r[i]*( (particles[0].mass*particles[1].mass)/dist2 );
	}
	
	Particle *pcurr;
	pcurr = &(particles[0]);
	for (int i=0; i<3; i++) {
		pcurr->dxvdt[i] = pcurr->xv[i+3];
		pcurr->dxvdt[i+3] = force_01[i]/pcurr->mass;
	}
	pcurr = &(particles[1]);
	for (int i=0; i<3; i++) {
		pcurr->dxvdt[i] = pcurr->xv[i+3];
		pcurr->dxvdt[i+3] = -force_01[i]/pcurr->mass;
	}
}
*/

void ParticleSystem::eval_dxvdt() {
//
//   Evaluate dxv/dt (the force) on each particles
// ================================================

	//-> First evaluate the forces f01, f12, f20 between the three particles
	//   Fij is the force on particle i due tuo particle j ... so fji = -fij
	double f01[3], f12[3], f20[3];
	eval_force(f01, 0, 1);
	eval_force(f12, 1, 2);
	eval_force(f20, 2, 0);
	
	//-> Evaluate dxv/dt for each particle
	Particle *pcurr;
	pcurr = &(particles[0]);
	for (int i=0; i<3; i++) {
		pcurr->dxvdt[i] = pcurr->xv[i+3];
		pcurr->dxvdt[i+3] = ( f01[i] -f20[i] )/pcurr->mass;
	}
	pcurr = &(particles[1]);
	for (int i=0; i<3; i++) {
		pcurr->dxvdt[i] = pcurr->xv[i+3];
		pcurr->dxvdt[i+3] = ( -f01[i] + f12[i] )/pcurr->mass;
	}
	pcurr = &(particles[2]);
	for (int i=0; i<3; i++) {
		pcurr->dxvdt[i] = pcurr->xv[i+3];
		pcurr->dxvdt[i+3] = ( -f12[i] + f20[i] )/pcurr->mass;
	}

}

void ParticleSystem::eval_force(double force[3], int pi, int pj) {

	//-> evaluate the distance between particles i and j
	double dist2=0; // Square distance
	double r[3]; // versor from i to j
	for (int l=0; l<3; l++) {
		dist2 += (particles[pi].xv[l] - particles[pj].xv[l])*(particles[pi].xv[l] - particles[pj].xv[l]);
	}
	//-> evaluate the distance between particles i and j
	for (int l=0; l<3; l++) {
		r[l] = (particles[pj].xv[l] - particles[pi].xv[l])/sqrt(dist2); // versor rij
		force[l] = K2 * r[l]*( (particles[pi].mass*particles[pj].mass)/dist2 );
	}
	
}


void ParticleSystem::step_rk() {

	Particle *pcurr;
	
	// Set xv_i to the initial value 	
	FOR_EACH_PARTICLE
		pcurr = &(particles[k]);
		for (int i=0; i<6; i++) pcurr->xv_i[i] = pcurr->xv[i];
	END_FOR_EACH_PARTICLE
	
	// First step: dx1
	eval_dxvdt();
	FOR_EACH_PARTICLE
		pcurr = &(particles[k]);
		for (int i=0; i<6; i++) {
        	pcurr->dx1[i] = dt*pcurr->dxvdt[i];
        	pcurr->xv[i]  = pcurr->xv_i[i] + 0.5*pcurr->dx1[i];
        }
	END_FOR_EACH_PARTICLE
	
	// Second step: dx2
	eval_dxvdt();
	FOR_EACH_PARTICLE
		pcurr = &(particles[k]);
		for (int i=0; i<6; i++) {
        	pcurr->dx2[i] = dt*pcurr->dxvdt[i];
        	pcurr->xv[i]  = pcurr->xv_i[i] + 0.5*pcurr->dx2[i];
        }
	END_FOR_EACH_PARTICLE
	
	// Third step: dx3
	eval_dxvdt();
	FOR_EACH_PARTICLE
		pcurr = &(particles[k]);
		for (int i=0; i<6; i++) {
        	pcurr->dx3[i] = dt*pcurr->dxvdt[i];
        	pcurr->xv[i]  = pcurr->xv_i[i] + pcurr->dx3[i];
        }
	END_FOR_EACH_PARTICLE
	
	// Fourth step: dx4
	eval_dxvdt();
	FOR_EACH_PARTICLE
		pcurr = &(particles[k]);
		for (int i=0; i<6; i++) {
        	pcurr->dx4[i] = dt*pcurr->dxvdt[i];
        }
	END_FOR_EACH_PARTICLE
	
	// Finally update the state
	FOR_EACH_PARTICLE
		pcurr = &(particles[k]);
		for (int i=0; i<6; i++) {
        	pcurr->xv[i] = pcurr->xv_i[i] + 
        		(pcurr->dx1[i] + 2*pcurr->dx2[i] + 2*pcurr->dx3[i] + pcurr->dx4[i])/6.0;
        }
	END_FOR_EACH_PARTICLE
	
	// Update the time
	time += dt;
	//print();
}

void ParticleSystem::step_euler()  {
	
	Particle *pcurr;
	
	eval_dxvdt();
	FOR_EACH_PARTICLE
		pcurr = &(particles[k]);
		for (int i=0; i<6; i++) {
       		pcurr->xv[i] += pcurr->dxvdt[i]*dt;
      	}
	END_FOR_EACH_PARTICLE
		
	// Update the time
	time += dt;
	print();		
}

void ParticleSystem::eval_cm() {
	double inv_mtot;
	inv_mtot = 1./( particles[0].mass + particles[1].mass + particles[2].mass);
	for (int i=0; i<3; i++) {
		cm[i] = inv_mtot*( (particles[0].mass)*(particles[0].xv[i]) + 
			(particles[1].mass)*(particles[1].xv[i]) + (particles[2].mass)*(particles[2].xv[i])); 
	}
}

void ParticleSystem::eval_energy() {

	Particle *pcurr;

	//-> kinetic energy
	double kin=0.; 
	FOR_EACH_PARTICLE
		pcurr = &(particles[k]);
		for (int i=3; i<6; i++) {
			kin += 0.5*pcurr->mass*pcurr->xv[i]*pcurr->xv[i];
		}
	END_FOR_EACH_PARTICLE
	
	//-> Potential energy
	double pot=0., couple=0.;
	
	eval_pot(&couple, 0, 1);
	pot += couple;
	
	eval_pot(&couple, 0, 2);
	pot += couple;
	
	eval_pot(&couple, 1, 2);
	pot += couple;

	
	energy = kin + pot;

}

void ParticleSystem::eval_pot(double *pen, int pi, int pj) {
	//
	// Evaluate the potential energy between particles pi and pj
	//
	
	double dist2=0.;
	for (int i=0; i<3; i++) {
		dist2 += (particles[pi].xv[i]-particles[pj].xv[i])*(particles[pi].xv[i]-particles[pj].xv[i]);
	}
	*pen = -(particles[pi].mass)*(particles[pj].mass)/sqrt(dist2);
}

void ParticleSystem::open_output() {
	file_w.open("output.dat");
	file_w << "time   p_index  mass   x1   x2   x3   v1   v2   v3 energy x_cm y_cm z_cm" << endl;
}

void ParticleSystem::printf() {
	Particle *pcurr;
	
	eval_energy();
	eval_cm();
	
	FOR_EACH_PARTICLE
		pcurr = &(particles[k]);
		file_w << time << " " << pcurr->n << " " << pcurr->mass << " "; 
		for (int i=0; i<6; i++) file_w << pcurr->xv[i] << " ";
		file_w << energy << " ";
		for (int i=0; i<3; i++) file_w << cm[i] << " ";
		file_w << endl;
	END_FOR_EACH_PARTICLE	
}

void ParticleSystem::close_output() {
	file_w.close();
}

void ParticleSystem::Draw()
{
        float size = 0.126; // 0.042; 
        Particle* pcurr;

        FOR_EACH_PARTICLE
                pcurr = &(particles[k]);
                if (k == 0) {
                	glColor3f(1.0, 0.0, 0.0);
                } else if (k == 1) {
                	glColor3f(0.0, 0.0, 1.0);
                }
                else {
                	glColor3f(0.0, 1.0, 0.0);
                }
                
                // Draw circles
                glBegin(GL_POLYGON);
                        glVertex2f(pcurr->xv[0] + size       , pcurr->xv[1]);
                        glVertex2f(pcurr->xv[0] + 0.707*size , pcurr->xv[1] + 0.707*size);
                        glVertex2f(pcurr->xv[0]              , pcurr->xv[1] + size);
                        glVertex2f(pcurr->xv[0] - 0.707*size , pcurr->xv[1] + 0.707*size);
                        glVertex2f(pcurr->xv[0] - size       , pcurr->xv[1]);
                        glVertex2f(pcurr->xv[0] - 0.707*size , pcurr->xv[1] - 0.707*size);
                        glVertex2f(pcurr->xv[0]              , pcurr->xv[1] - size);
                        glVertex2f(pcurr->xv[0] + 0.707*size , pcurr->xv[1] - 0.707*size);
                glEnd();
        END_FOR_EACH_PARTICLE
}


