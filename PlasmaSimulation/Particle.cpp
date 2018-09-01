//
// Created on 23.06.2018.
//

#include <stdlib.h>
#include <malloc.h>
#include "Particle.h"

Particle *createParticle(double coords[2], double velocity[3]) {
	Particle *particle = (Particle *)malloc(sizeof(Particle *));
	particle->velocity = velocity;
	particle->coords = coords;
	return particle;
}

Particle *createParticleCopy(Particle *particle) {
	Particle *copy = (Particle *)malloc(sizeof(Particle *));
	copy->coords = (double *)malloc(2 * sizeof(double));
	copy->coords[0] = particle->coords[0];
	copy->coords[1] = particle->coords[1];

	copy->velocity = (double *)malloc(3 * sizeof(double));
	copy->velocity[0] = 0;
	copy->velocity[1] = 0;
	copy->velocity[2] = 0;
	return copy;
}



