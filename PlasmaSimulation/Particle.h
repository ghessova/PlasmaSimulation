#ifndef PLASMASIMULATION_PARTICLE_H
#define PLASMASIMULATION_PARTICLE_H

typedef struct Particle {
	double *coords;      // coordinates (2 elements)
	double *velocity;     // velocity (3 elements)
};

Particle *createParticle(double coords[2], double velocity[3]);
#endif //PLASMASIMULATION_PARTICLE_H#pragma once
