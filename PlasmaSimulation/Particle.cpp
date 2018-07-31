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


