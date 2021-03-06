#include <iostream>
#include <fstream>
#include "Simulation.h"
#include <random>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <map>



// simulation parameters and constants
const int grid = 100;
const long int macroPartsCount = 1e5;  // number of macro particles
const int steps = 100;					// number of simulation steps

const double gridWidth = 1e-2;		// grid width (m) // or 1e-2
const double gridHeight = 1e-2;		// grid height (m)
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



// magnetic field
const double theta = PI / 4;	// angle
const double B0 = 1.15;			// size

//int crashCount;
int plateCrashCount;	// number of particles that crashed on the plates
int gapCount;			// number of particles in the gap


const double gap = 0.2 * gridWidth;				// gap width: -----------    ----------
const double plateHeight = 0.23 * gridHeight;	// height of the plates surroundind the gap

double elPotential[grid + 1][grid + 1];     // phi - electron potential field (nodes)
double iontPotential[grid + 1][grid + 1];   // rho - ion potential field

double cellWidth = gridWidth / grid;
double cellHeight = gridHeight / grid;

const double T_e = 1e7; // electron temperature (Kelvins)
const double T_i = 1e7; // ion temperature (Kelvins)

const double  k_B = 1.38064852e-23;                        // m2 kg s-2 K-1 (Boltzman)
double cMult = macroPartsCount / macroPartSize;
const double  v_Te = sqrt((k_B * T_e) / (cMult * elMass));          // m s-1
const double  v_Ti = sqrt((k_B * T_i) / (cMult * ionMass));

//double const dt = 1e-12;
// Half of maximum allowed time not to croos double cell (conditions are allright)
double const dt = 5e-12;

// Constructor														
Simulation::Simulation() {
}


// initialize specified 2D field to 0
void init2DField(double field[grid + 1][grid + 1]) {
	for (int x = 0; x < grid + 1; x++) {	
		for (int y = 0; y < grid + 1; y++) {		
			field[x][y] = 0;
		}
	}

}

std::vector<double> operator+(std::vector<double> a, std::vector<double> b) {
	std::vector<double> res(a.size());
	for (unsigned int i = 0; i < a.size(); i++) {
		res[i] = a[i] + b[i];
	}
	return res;
}

std::vector<double> operator*(double a, std::vector<double> b) {
	std::vector<double> res(b.size());
	for (unsigned int i = 0; i < b.size(); i++) {
		res[i] = a*b[i];
	}
	return res;
}


/*
Utility function for filling 3D vector.
*/
std::vector<double> fillVector(double a, double b, double c) {
	std::vector<double> res(3);
	res[0] = a;
	res[1] = b;
	res[2] = c;
	return res;
}

std::vector<double> crossProduct(std::vector<double> a, std::vector<double> b) {
	std::vector<double> result(3);
	int j, k;
	for (int i = 0; i < 3; i++) {
		j = (i + 1) % 3;
		k = (i + 2) % 3;
		result[i] = a[j] * b[k] - a[k] * b[j];
	}
	return result;
}
double getVectorSize(std::vector<double> vec) {
	double size = 0;
	for (int i = 0; i < vec.size(); i++) {
		size += vec.at(i)*vec.at(i);
	}
	return sqrt(size);
}



/**
* Generates macroparticles in the beginning of the simulation.
* This method is called twice - once for electrons and then for ions.
*/
std::vector<Particle *> Simulation::initialize(bool isSource, double velocityRange) {
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
		 if (isSource) {
			 coordinates[1] = y * gridHeight;
		 }
		 else {
			 coordinates[1] = plateHeight + y * (gridHeight - plateHeight); // the particles are not generated below the plate level
		 }


		 // velocity vector 
		 double *velocity = (double *)malloc(3 * sizeof(double));

		 //do {
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

			 velocity[0] = sqrt(-2 * log(u1)) * cos(2 * PI * u2) * (velocityRange/sqrt(3)); // x
			 velocity[1] = sqrt(-2 * log(u1)) * sin(2 * PI * u2) * (velocityRange/sqrt(3)); // y
			 velocity[2] = sqrt(-2 * log(u3)) * cos(2 * PI * u4) * (velocityRange/sqrt(3)); // z

		 //} while (getVectorSize(fillVector(velocity[0], velocity[1], velocity[2])) > velocityRange);

		 Particle *particle = createParticle(coordinates, velocity);
		 particles.push_back(particle);

	}
	return particles;	
}


bool intersectsThePlate(double x, double y) {
	return y <= plateHeight && (x <= gridWidth / 2 - gap / 2 || x >= gridWidth / 2 + gap / 2);
}

/*
	Evaluates coordinates and velocities of the particles using Boris scheme.
*/
void boris(std::vector<Particle *> *particles, std::vector<double> *ex, std::vector<double> *ey, bool electrons) {
	std::vector<double> B(3);  // magnetic field
	B = fillVector(B0 * sin(theta), -B0 *cos(theta), 0);

	std::vector<double> ef(3), vm(3), s(3), T(3), E, coords, v;

	double m;   // mass
	int charge; // +-1

	if (electrons) {
		m = elMass;
		charge = -1;
	}
	else {
		m = ionMass;
		charge = 1;
	}


	for (int i = 0; i < particles->size(); i++) { // iteration through particles
		Particle *particle = particles->at(i);
		double x = particle->coords[0];
		double y = particle->coords[1];

		double vx = particle->velocity[0];
		double vy = particle->velocity[1];
		double vz = particle->velocity[2];

		coords = fillVector(x, y, 0);
		v = fillVector(vx, vy, vz);
		E = fillVector(ex->at(i), ey->at(i), 0);


		for (int j = 0; j < 3; j++) {
			ef[j] = charge * q / m*dt / 2 * E[j];
			T[j] = charge * q*B[j] * dt / m / 2;
		}
		double T_vel = getVectorSize(T);

		for (int j = 0; j < 3; j++) {
			s[j] = 2 * T[j] / (1 + T_vel * T_vel);
		}
		vm = v + ef;
		v = vm + crossProduct(vm + crossProduct(vm, T), s) + ef;
		coords = coords + dt * v;

		particle->coords[0] = coords.at(0);
		particle->coords[1] = coords.at(1);

		particle->velocity[0] = v.at(0);
		particle->velocity[1] = v.at(1);
		particle->velocity[2] = v.at(2);

	}
}


/*
   Evaluates coordinates and velocities of the source particles.
   No forces affect these particles, only the particles that reach the area borders are reflected.
   No plates included.
   Returns vector of reflected particles.
*/
std::vector<Particle *> source(std::vector<Particle *> *particles, bool areElectrons) {
	std::vector<Particle *> reflectedParticles;

	double dtx = areElectrons ? dt*2 : dt;

	for (int i = 0; i < particles->size(); i++) { // iteration through particles
		Particle *particle = particles->at(i);
		double x = particle->coords[0];
		double y = particle->coords[1];

		double vx = particle->velocity[0];
		double vy = particle->velocity[1];
		double vz = particle->velocity[2];


		// without forces
		x = x + dtx * vx;
		y = y + dtx * vy;

		bool isReflected = false;

		// position is out of range -> the particle is reflected
		if (x < 0) {
			x = -x;
			vx = -vx;
			isReflected = true;			
		}
		if (x > gridWidth) {
			x = 2 * gridWidth - x;
			vx = -vx;
			isReflected = true;
		}
		if (y > gridHeight) {
			y = 2 * gridHeight - y;
			vy = -vy;
			isReflected = true;
		}
		if (y < 0) {	// special case 
			y = gridHeight + y;
			isReflected = true;
		}

		particle->coords[0] = x;
		particle->coords[1] = y;

		particle->velocity[0] = vx;
		particle->velocity[1] = vy;
		particle->velocity[2] = vz;

		if (isReflected) {
			Particle *copy = createParticleCopy(particle);
			reflectedParticles.push_back(copy);
		}		

	}
	//std::cout << reflectedParticles.size() << std::endl;
	return reflectedParticles;

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

void printVector(std::vector<double> *numbers, const char *fileName) {
	std::ofstream outputFile;
	outputFile.open(fileName);
	for (int i = 0; i < numbers->size(); i++) {
		outputFile << numbers->at(i) << std::endl;
	}
	outputFile.close();
}


void printFrequencies(std::map<double, int> *map, const char* fileName) {
	std::ofstream outputFile;
	outputFile.open(fileName);
	typedef std::map<double, int>::iterator it_type;

	for (it_type iterator = map->begin(); iterator != map->end(); iterator++) {
		outputFile << iterator->first << " " << iterator->second << std::endl;;
	}
	outputFile.close();
}



double getCellIndex(double position, double range, int cellsCount) {
	if (position == range) { // the particle is on the border
		return cellsCount - 1;
	}
	return position * cellsCount / range;
}

/*

	This method counts electric charge in the nodes of the grid using the cloud-in-cell algorithm.
	The values are later used for getting the electric potential.

	see https://pdfs.semanticscholar.org/5010/d47d9fcc9539cc315a54400cae2ce17eb1e2.pdf

	rho - chargeField - output parameter
	Q - charge
*/
void countCharge(std::vector<Particle *> *particles, double rho[grid + 1][grid + 1], double Q) {
	
	for (int i = 0; i < particles->size(); i++) { // iteration through particles
		Particle *particle = particles->at(i);
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



double getIntervalSize(double min, double max, int intervals) {
	return (max - min) / intervals;
}

double getMapKey(double velocity, double intervalSize) {
	int intervalIndex = velocity / intervalSize;
	return intervalIndex * intervalSize;
}

/*
	Adds velocity to map of their frequencies.
	Map key is lower boundary of velocity interval.
*/


void addToMap(std::map<double, int> *velocities, double velocity, double intervalSize) {
	// first determine map key
	double key = getMapKey(velocity, intervalSize);
	// insert
	if (velocities->count(key) == 0) {
		velocities->insert(std::pair<double, int>(key, 1));
	}
	else {
		double freq = velocities->at(key);
		freq++;
		velocities->erase(key);
		velocities->insert(std::pair<double, int>(key, freq));
	}

}

void getVelocities(std::vector<Particle *> *particles, std::vector<double> *velocities, double *min, double *max) {
	*min = v_Te * 10;
	*max = 0;
	double velocity;
	for (int i = 0; i < particles->size(); i++) {
		Particle *particle = particles->at(i);
			std::vector<double> vel = fillVector(particle->velocity[0], particle->velocity[1], particle->velocity[2]);
			velocity = getVectorSize(vel);
			if (velocity < *min) {
				*min = velocity;
			}
			if (velocity > *max) {
				*max = velocity;
			}
			velocities->push_back(velocity);
	}
}

void createFreqMap(std::map<double, int> *map, std::vector<double> *velocities, double min, double max, int intervals) {
	double intervalSize = getIntervalSize(min, max, intervals);
	for (int i = 0; i < velocities->size(); i++) {
		addToMap(map, velocities->at(i), intervalSize);
	}
}

void velocityCheck(std::vector<Particle *> *particles, const char* outputFile) {
	int intervals = 100;
	std::vector<double> velocities;
	double min = 0;
	double max = 0;
	getVelocities(particles, &velocities, &min, &max);
	std::map<double, int> map;
	createFreqMap(&map, &velocities, min, max, intervals);
	printFrequencies(&map, outputFile);
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

double bilinearInterpolation(double x, double y, double x1, double x2, double y1, double y2, double fQ11, double fQ12, double fQ21, double fQ22) {
	double fxy1 = (x2 - x) * fQ11 / (x2 - x1) + (x - x1) * fQ21 / (x2 - x1);
	double fxy2 = (x2 - x) * fQ12 / (x2 - x1) + (x - x1) * fQ22 / (x2 - x1);
	return (y2 - y) * fxy1 / (y2 - y1) + (y - y1) * fxy2 / (y2 - y1);
}

/*
 e - electric field in grid nodes
 particles - ions or electrons
 ex - 
*/
void interpolateElectricField(double e[grid + 1][grid + 1][2], std::vector<Particle *> *particles, std::vector<double> *ex, std::vector<double> *ey) {
	for (int i = 0; i < particles->size(); i++) { // iteration through particles
		Particle *particle = particles->at(i);
		double x = particle->coords[0];
		double y = particle->coords[1];

		// indexes of node in the grid (node matrix) 
		int cellX = getCellIndex(x, gridWidth, grid);
		int cellY = getCellIndex(y, gridHeight, grid);

		ex->push_back(bilinearInterpolation(x, y, cellX * cellWidth, (cellX + 1) * cellWidth, cellY * cellHeight, (cellY + 1) * cellHeight, e[cellX][cellY][0], e[cellX][cellY + 1][0], e[cellX + 1][cellY][0], e[cellX + 1][cellY + 1][0]));
		ey->push_back(bilinearInterpolation(x, y, cellX * cellWidth, (cellX + 1) * cellWidth, cellY * cellHeight, (cellY + 1) * cellHeight, e[cellX][cellY][1], e[cellX][cellY + 1][1], e[cellX + 1][cellY][1], e[cellX + 1][cellY + 1][1]));

		//particle->ex = ex;
		//particle->ey = ey;	
		//particle->ez = 0; // ez
	}

	//getchar();
}


std::vector<Particle *> deterministicInitialize() {
	std::vector<Particle *> particles; // vector of macro-particles (ions / electrons)

	for (int i = 0; i < macroPartsCount; i++) {
		double *coordinates = (double *)malloc(2 * sizeof(double));
		coordinates[0] = gridWidth / 2;
		coordinates[1] = gridHeight / 2;

		double *velocity = (double *)malloc(3 * sizeof(double));
		velocity[0] = 0;
		velocity[1] = 0;
		velocity[2] = 0;


		Particle *particle = createParticle(coordinates, velocity);
		particles.push_back(particle);
	}
	return particles;
}

void initElField(double elField[grid + 1][grid + 1][2]) {
	for (int i = 0; i < grid + 1; i++) {
		for (int j = 0; j < grid + 1; j++) {
			elField[i][j][0] = 0;
			elField[i][j][1] = 0;
		}
	}
}

void init1DArray(double eArray[macroPartsCount]) {
	for (int i = 0; i < grid + 1; i++) {
		eArray[i] = 0;
	}
}

void initVector(std::vector<double> *vec, int count) {
	for (int i = 0; i < count; i++) {
		vec->push_back(0);
	}
}

bool isOutOfTheBox(double x, double y, bool isSource) {
	if (x < 0 || x > gridWidth || y < 0 || y > gridHeight) {
		return true;
	}
	if (!isSource && intersectsThePlate(x, y)) {
		plateCrashCount++;
		return true;
	}
	return false;	
}

void checkParticles(std::vector<Particle *> *particles, std::vector<Particle* > *out, bool erase) {
	out->clear();
	int count = 0;
	for (int i = particles->size() - 1; i >= 0 ; i--) {
		Particle *particle = particles->at(i);
		double x = particle->coords[0];
		double y = particle->coords[1];
		if (isOutOfTheBox(x,y,!erase)) {
			if (erase) {
				particles->erase(particles->begin() + i);
			}	
			//free(particle);
			count++ ;
			out->push_back(particle);
		}
	}
	count = count;
}


/*
	Adds reflected source particles to the vector of particles in the simulated area.
*/
void addParticles(std::vector<Particle *> *particles, std::vector<Particle *> *reflected) {
	for (int i = 0; i < reflected->size(); i++) {
		if (particles->size() >= macroPartsCount) {
			break;
		}
		particles->push_back(reflected->at(i));
	}
}


void Simulation::simulate() {

	// particles in external source
	std::vector<Particle *> sourceElectrons = initialize(true, v_Te);
	std::vector<Particle *> sourceIons = initialize(true, v_Ti);

	// particles in simulated area
	std::vector<Particle *> electrons = initialize(false, v_Te);
	std::vector<Particle *> ions = initialize(false, v_Ti);


	double rho[grid + 1][grid + 1];	// charge matrix
	double phi[grid + 1][grid + 1]; // potential matrix
	double elField[grid + 1][grid + 1][2]; // matrix of electric field vectors
	initElField(elField);
	std::vector<double> ex;
	std::vector<double> ey;

	std::vector<Particle *> outIonsSource;
	std::vector<Particle *> outElectronsSource;
	std::vector<Particle *> outIonsSimulation;
	std::vector<Particle *> outElectronsSimulation;
	std::vector<double> potentialVector;


	for (int t = 0; t < steps; t = t++ /*+ dt*/) { // time iteration
		plateCrashCount = 0;
		gapCount = 0;


		outIonsSource.clear();
		outElectronsSource.clear();
		outIonsSource = source(&sourceIons, false);
		outElectronsSource = source(&sourceElectrons, true);

		// add reflected particles from source to simulation
		addParticles(&ions, &outIonsSource);
		addParticles(&electrons, &outElectronsSource);



		// 2. charge in grid nodes (both ions and electrons contribute)
		init2DField(rho);

		countCharge(&ions, rho, macroQ);
		countCharge(&electrons, rho, -macroQ);

		// 3. potential in grid nodes (is obtained from charge)
		init2DField(phi);
		countPotential(rho, phi, 400, t);
		potentialVector.push_back(phi[10][25]);
		countElectricField(elField, phi);

		ex.clear();
		ey.clear();
		interpolateElectricField(elField, &electrons, &ex, &ey);
		interpolateElectricField(elField, &ions, &ex, &ey);
		boris(&electrons, &ex, &ey, true);
		boris(&ions, &ex, &ey, false);

		// check which particles are out of the box
		checkParticles(&ions, &outIonsSimulation, true);
		checkParticles(&electrons, &outElectronsSimulation, true);



		if (t == 0) {
			printMatrix(phi, "potential0.txt");
			printMatrix(rho, "charge0.txt");
			printParticles(electrons, "electrons0.txt");
		}


		if (t % 5 == 0) {
			std::cout << "Step " << t << " completed." << std::endl;
			std::cout << "Simulation Electrons out " << outElectronsSimulation.size() << std::endl;
			std::cout << "Source Electrons out " << outElectronsSource.size() << std::endl;
			std::cout << "CrashPlates " << plateCrashCount << std::endl;
		}

	}
	printParticles(electrons, "electrons.txt");
	printParticles(ions, "ions.txt");
	printParticles(sourceElectrons, "electronsSource.txt");
	printParticles(sourceIons, "ionsSource.txt");
	
	printMatrix(rho, "charge.txt");
	printMatrix(phi, "potential.txt");
	
	printVector(&potentialVector, "potentialZone.txt");

	velocityCheck(&sourceElectrons, "checksource.txt");
	velocityCheck(&electrons, "elVelocities.txt");
	velocityCheck (&ions, "ionVelocities.txt");

	velocityCheck(&outIonsSimulation, "outIonsSimulation.txt");
	velocityCheck(&outElectronsSimulation, "outElectronsSimulation.txt");
	velocityCheck(&outIonsSource, "outIonsSource.txt");
	velocityCheck(&outElectronsSource, "outElectronsSource.txt");

	//std::cout << plateCrashCount << std::endl;
	//std::cout << gapCount << std::endl;


}






