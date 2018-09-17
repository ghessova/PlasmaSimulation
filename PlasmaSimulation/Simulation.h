

#ifndef PLASMASIMULATION_SIMULATION_H
#define PLASMASIMULATION_SIMULATION_H


#include <vector>
#include "Particle.h"



class Simulation {

	// properties:
	//std::vector<Particle *> electrons;
	//std::vector<Particle *> ions;
	//int gridSize;           // size of grid site
	//int** gridCharges;

private:
	std::vector<Particle *> Simulation::initialize(bool isSource, double velocityRange);
public:
	Simulation(); // constructor
	void simulate();

};
void printParticles(std::vector<Particle *> v, const char* fileName);
//void printPotential(double** phi, const char* fileName);


#endif //PLASMASIMULATION_SIMULATION_H
