//
// Created on 23.06.2018.
//

#include <iostream>
#include <fstream>
#include "Simulation.h"
#include <random>
#include <time.h>
#include <math.h>
#include <stdlib.h>



// simulation parameters and constants
const int grid = 100;
const long int macroPartsCount = 10000;  // number of macro particles
const int steps = 100;					// number of simulation steps

const double gridWidth = 1.5e-3;		// grid width (m) // or 1e-2
const double gridHeight = 1.5e-3;		// grid height (m)
const double particlesDensity = 1e12;   // electrons per meter^2 // or 10^8                              
const double particlesCount = particlesDensity * gridWidth * gridHeight; // computed number of particles
const double macroPartSize = particlesCount / macroPartsCount; // number of particles (ions / electrons) in macro particle (real/used particles) .. or 1e8

const double q = 1.60217662e-19;                        // charge of electron
const double macroQ = macroPartSize * q;
const double maxVelocity = 1;
const double epsilon = 8.854e-12;	// permitivity
const double elMass = 9.1e-31;		// electron mass (kg)
const double ionMass = 1.673e-27;   // ion mass (kg)

const double PI = 3.141592653589793;
double const dt = 1e-5;

//int crashCount;
int plateCrashCount;	// number of particles that crashed on the plates
int gapCount;			// number of particles in the gap


const double gap = 0.1 * gridWidth;				// gap width: -----------    ----------
const double plateHeight = 0.3 * gridHeight;	// height of the plates surroundind the gap

double elPotential[grid + 1][grid + 1];     // phi - electron potential field (nodes)
double iontPotential[grid + 1][grid + 1];   // rho - ion potential field

double cellWidth = gridWidth / grid;
double cellHeight = gridHeight / grid;

// Constructor														
Simulation::Simulation() {
	//crashCount = 0;
}


// initialize specified 2D field to 0
void init2DField(double field[grid + 1][grid + 1]) {
	for (int x = 0; x < grid + 1; x++) {	
		for (int y = 0; y < grid + 1; y++) {		
			field[x][y] = 0;
		}
	}

}


/**
* Generates macroparticles in the beginning of the simulation.
* This method is called twice - once for electrons and then for ions.
*/
std::vector<Particle *> Simulation::initialize() {
	std::vector<Particle *> particles; // vector of macro-particles (ions / electrons)
	double rm = double(RAND_MAX); 
	time_t t;
	srand((unsigned)time(&t));

	for (int i = 0; i < macroPartsCount; i++) {
		 // coordinates
		 double x = rand() / rm; // number between 0 and 1
		 double y = rand() / rm;
		 double *coordinates = (double *)malloc(2 * sizeof(double));
		 coordinates[0] = x * gridWidth;
		 coordinates[1] = plateHeight + y * (gridHeight - plateHeight); // the particles are not generated below the plate level


		 // velocity vector - auxiliary variables
		 double u1 = (rand() / rm);
		 double u2 = (rand() / rm);
		 double u3 = (rand() / rm);
		 double u4 = (rand() / rm);

		 if (u1 == 0) {
			 u1++;
		 }
		 if (u2 == 0) {
			 u2++;
		 }
		 if (u3 == 0) {
			 u3++;
		 }
		 if (u4 == 0) {
			 u4++;
		 }

		 // velocity vector 
		 double *velocity = (double *)malloc(3 * sizeof(double));

		 velocity[0] = sqrt(-2 * log(u1)) * cos(2 * PI * u2) * maxVelocity; // x
		 velocity[1] = sqrt(-2 * log(u1)) * sin(2 * PI * u2) * maxVelocity; // y
		 velocity[2] = sqrt(-2 * log(u3)) * cos(2 * PI * u4) * maxVelocity; // z

		 //double *elForce = (double *)malloc(3 * sizeof(double));
		 Particle *particle = createParticle(coordinates, velocity);
		 particles.push_back(particle);
		 //std::cout << particle->coords[0] << " " << particle->coords[1] << std::endl;

	}
	return particles;	
}

bool intersectsThePlate(double x, double y) {
	return y <= plateHeight && (x <= gridWidth / 2 - gap / 2 || x >= gridWidth / 2 + gap / 2);
}

void countCoordsAndVelocity(std::vector<Particle *> particles) {
	for (int i = 0; i < macroPartsCount; i++) { // iteration through particles
		Particle *particle = particles.at(i);
		double x = particle->coords[0];
		double y = particle->coords[1];

		double vx = particle->velocity[0];
		double vy = particle->velocity[1];

		// without forces
		x = x + dt * vx;
		y = y + dt * vy;

		// position is out of range -> the particle is reflected
		if (x < 0) {
			x = -x;
			vx = -vx;
		}
		if (x > gridWidth) {
			x = 2 * gridWidth - x;
			vx = -vx;
		}
		if (y > gridHeight) {
			y = 2 * gridHeight - y;
			vy = -vy;
		}
		if (y < plateHeight) {	// special case - below the plate level
			if (x > gridWidth / 2 - gap / 2 && x < gridWidth / 2 + gap / 2) { // particle is in the gap																			  
				gapCount++;
				if (y < 0) { // particle is in the bottom
					x = rand() / RAND_MAX * gridWidth;
					y = y + gridHeight - plateHeight; // new particle is generated on the opposite side
				}

			}
			else { // particle is in on the plates
				y = y + gridHeight - plateHeight; // new particle is generated on the opposite side
				plateCrashCount++;
			}
		}

		particle->coords[0] = x;
		particle->coords[1] = y;

		particle->velocity[0] = vx;
		particle->velocity[1] = vy;
	}
}

void printParticles(std::vector<Particle *> v, const char* fileName) {
	//typedef std::vector<Particle *>::iterator it_type;	
	std::ofstream outputFile;
	outputFile.open(fileName);
	for (int i = 0; i < v.size(); i++) {
		Particle *particle = v.at(i);
		outputFile << particle->coords[0] << " " << particle->coords[1] << std::endl;
	}
	outputFile.close();
}

void printMatrix (double phi[grid +1][grid+1], const char* fileName) {
	//typedef std::vector<Particle *>::iterator it_type;	
	std::ofstream outputFile;
	outputFile.open(fileName);
	//for (int j = grid; j >= 0; j--) { // inner nodes
	for (int j = 0; j < grid + 1; j++) { // inner nodes

		for (int i = 0; i < grid + 1; i++) {
			outputFile << phi[i][j] << " ";

		}
		outputFile << std::endl;
	}
	outputFile.close();
}



double getCellIndex(double position, double range, int cellsCount) {
	return position * cellsCount / range;
}

/*

	This method counts electric charge in the nodes of the grid using the cloud-in-cell algorithm.
	The values are later used for getting the electric potential.

	see https://pdfs.semanticscholar.org/5010/d47d9fcc9539cc315a54400cae2ce17eb1e2.pdf

	rho - chargeField - output parameter
	Q - charge
*/
void countCharge(std::vector<Particle *> particles, double rho[grid + 1][grid + 1], double Q) {
	
	for (int i = 0; i < macroPartsCount; i++) { // iteration through particles
		Particle *particle = particles.at(i);
		double x = particle->coords[0];
		double y =  particle->coords[1];

		// indexes of node in the grid (node matrix) // todo rename

		int cellX = getCellIndex(x, gridWidth, grid);
		int cellY = getCellIndex(y, gridHeight, grid);
		double dx = x - cellX  * cellWidth;		// x distance from the left bottom corner of the cell
		double dy = y - cellY  * cellHeight;		// y distance from the left bottom corner of the cell

		rho[cellX][cellY] += Q * (cellWidth - dx) * (cellHeight - dy) / (cellWidth * cellHeight);
		rho[cellX + 1][cellY] += Q * dx * (cellHeight - dy) / (cellWidth * cellHeight);
		rho[cellX][cellY + 1] += Q * (cellWidth - dx) *  dy / (cellWidth * cellHeight);
		rho[cellX + 1][cellY + 1] += Q * dx * dy / (cellWidth * cellHeight);
	}

}

/*
	iters - number of iterations
*/
double getOmega(int iters) {
	return 2 / (1 + sin(PI * iters / (iters + 1)));
}

/*
	rho - matrix of charges - input parameter
	phi - matrix of potentials - output parameter
	iters - number of iterations of SOR
	i - simulation step
*/
void countPotential(double rho[grid + 1][grid + 1], double phi[grid + 1][grid + 1], int iters, int i) {
	
	double omega = getOmega(iters);
	for (int k = 0; k < iters; k++) {
		for (int i = 0; i < grid + 1; i++) { // outer nodes
			phi[i][0] = phi[i][1];
			phi[i][grid] = phi[i][grid - 1];
			phi[0][i] = phi[1][i];
			phi[grid][i] = phi[grid - 1][i];
		}

		// SOR - Successive over-relaxation - iterative method
		// see https://en.wikipedia.org/wiki/Successive_over-relaxation
		for (int i = 1; i < grid; i++) { // inner nodes
			for (int j = 1; j < grid; j++) {
				double x = i * cellWidth;		// absolute coordinates of the nodes
				double y = j * cellHeight;
				if (intersectsThePlate(x, y)) { // on the plates, the potential is zero
					phi[i][j] = 0;
				}
				else {
					phi[i][j] = (1 - omega) * phi[i][j] + (omega / 4) * (phi[i][j + 1] + phi[i + 1][j] + phi[i - 1][j] + phi[i][j - 1] + rho[i][j] * cellWidth * cellHeight / epsilon);
				}
			}
		}
		//getchar();
	}

}

/*
phi - matrix of potentials - input parameter
e - electric field - output parameter 
*/
void countElectricField(double e[grid + 1][grid + 1][2], double phi[grid + 1][grid + 1]) {
	for (int i = 1; i < grid; i++) {
		for (int j = 1; j < grid; j++) {
				e[i][j][0] = (phi[i][j - 1] - phi[i][j + 1]) / (2 * cellWidth);			// ex
			e[i][j][1] = (phi[i - 1][j] - phi[i + 1][j]) / (2 * cellHeight);		// ey
		}
	}
	// border values
	for (int i = 0; i <= grid; i++) {
		// ex
		e[0][i][0] = e[1][i][0];
		e[grid][i][0] = e[grid - 1][i][0];
		e[i][0][0] = e[i][1][0];
		e[i][grid][0] = e[i][grid - 1][0];
		// ey
		e[0][i][1] = e[1][i][1];
		e[grid][i][1] = e[grid - 1][i][1];
		e[i][0][1] = e[i][1][1];
		e[i][grid][1] = e[i][grid - 1][1];

	}
}

double billinearInterpolation(double x, double y, double x1, double x2, double y1, double y2, double fQ11, double fQ12, double fQ21, double fQ22) {
	double fxy1 = (x2 - x) * fQ11 / (x2 - x1) + (x - x1) * fQ21 / (x2 - x1);
	double fxy2 = (x2 - x) * fQ12 / (x2 - x1) + (x - x1) * fQ22 / (x2 - x1);
	return (y2 - y) * fxy1 / (y2 - y1) + (y - y1) * fxy2 / (y2 - y1);
}

void interpolateElectricField(double e[grid + 1][grid + 1][2], std::vector<Particle *> particles, double *exArray, double *eyArray) {
	double ex, ey;	

	for (int i = 0; i < macroPartsCount; i++) { // iteration through particles
		Particle *particle = particles.at(i);
		double x = particle->coords[0];
		double y = particle->coords[1];

		// indexes of node in the grid (node matrix) 
		int cellX = getCellIndex(x, gridWidth, grid);
		int cellY = getCellIndex(y, gridHeight, grid);

		exArray[i] = billinearInterpolation(x, y, cellX * cellWidth, (cellX + 1) * cellWidth, cellY * cellHeight, (cellY + 1) * cellHeight, e[cellX][cellY][0], e[cellX][cellY + 1][0], e[cellX + 1][cellY][0], e[cellX + 1][cellY + 1][0]);
		eyArray[i] = billinearInterpolation(x, y, cellX * cellWidth, (cellX + 1) * cellWidth, cellY * cellHeight, (cellY + 1) * cellHeight, e[cellX][cellY][1], e[cellX][cellY + 1][1], e[cellX + 1][cellY][1], e[cellX + 1][cellY + 1][1]);

		//particle->ex = ex;
		//particle->ey = ey;	
		//particle->ez = 0; // ez
	}

	//getchar();
}







void Simulation::simulate() {
	std::vector<Particle *> electrons = initialize();
	std::vector<Particle *> ions = initialize();

	//printParticles(electrons, 0);
	double rho[grid + 1][grid + 1];	// charge matrix
	double phi[grid + 1][grid + 1]; // potential matrix
	double elField[grid + 1][grid + 1][2]; // matrix of electric field vectors
	double exArray[macroPartsCount];
	double eyArray[macroPartsCount];

	for (int t = 0; t < steps; t = t++ /*+ dt*/) { // time iteration

		gapCount = 0;
		// 1. coordinates and velocity
		countCoordsAndVelocity(ions);
		

		// 2. charge in grid nodes (both ions and electrons contribute)
		init2DField(rho);
		countCharge(ions, rho,  macroQ);
		countCharge(electrons, rho, -macroQ);

		// 3. potential in grid nodes (is obtained from charge)
		init2DField(phi);
		countPotential(rho, phi, 1000, t);
		countElectricField(elField, phi);
		interpolateElectricField(elField, electrons, exArray, eyArray);
		//interpolateElectricField(elField, ions);

		if (t == 0) {
			printMatrix(phi, "potential0.txt");
			printMatrix(rho, "charge0.txt");
			printParticles(electrons, "electrons0.txt");
		}
		countCoordsAndVelocity(electrons);
		if (t % 10 == 0) {
			std::cout << "Step " << t << " completed." << std::endl;
		}

		//std::cout << particles.at(0)->coords[0] << " " << particles.at(0)->coords[1] << std::endl;

	}
	printParticles(electrons, "electrons.txt");
	printParticles(ions, "ions.txt");
	printMatrix(rho, "charge.txt");
	printMatrix(phi, "potential.txt");
	


	std::cout << plateCrashCount << std::endl;
	std::cout << gapCount << std::endl;


}






