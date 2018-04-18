#pragma once
#include "stdafx.h"
#include "constant.h"
#include "Particle.h"

/*
Species - data structure to store vector and basic info for particular particle species
*/

class Species
{
public:

	/* variable list */
	double mass;                          // particle's mass
	double charge;                        // particle's charge
	double specific_charge;               // particle's specific charge
	vector<Particle<double>> part;        // particle vector to store particle information

	/* function list */
	Species(double m, double c) :mass(m), charge(c), specific_charge(c / m) {};  // constructor function
 	void initialMaxwell(int num, double temp, double xbound1, double xbound2);   // initialize Maxwell distributed particles
	void saveParticle(string directory);                                         // save particle information in given directory
	void addParticle(double x, double v);                                        // add a particle with x posiiton, v velocity
	void remParticle(int index);                                                 // remove a particle with given index number 
	double kinetic();                                                            // return the total kinetic energy for this species
};	
