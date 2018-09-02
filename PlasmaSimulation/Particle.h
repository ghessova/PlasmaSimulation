#ifndef PLASMASIMULATION_PARTICLE_H
#define PLASMASIMULATION_PARTICLE_H

typedef struct Particle {
	double *coords;      // coordinates (2 elements)
	double *velocity;     // velocity (3 elements)
	//double *elForce;
	/*double ex;
	double ey;
	double ez;*/
};

Particle *createParticle(double coords[2], double velocity[3]);

Particle *createParticleCopy(Particle *particle);

void freeParticle(Particle *particle);


#endif //PLASMASIMULATION_PARTICLE_H#pragma once
