// Main.cpp : Defines the entry point for the console application.
//
#include "Particle.h"
#include "Simulation.h"

#include <iostream>
#include <fstream>
#include <random>
#include <time.h>
#include <math.h>
#include <stdlib.h>


int main()
{
	/*std::vector<int> numbers;
	for (int i = 0; i < 10; i++) {
		numbers.push_back(i);
	}
	for (int i = 0; i < numbers.size(); i++) {
		std::cout << numbers[i] << std::endl;
	}
	for (int i = 0; i < numbers.size(); i++) {
		if (numbers[i] < 1 == 0) {
			numbers.erase(numbers.begin() + i);
		}
	}
	for (int i = 0; i < numbers.size(); i++) {
		std::cout << numbers[i] << std::endl;
	}*/
	
	std::cout << "Simulation started." << std::endl;
	Simulation simulation;
	clock_t start = clock();
	simulation.simulate();
	clock_t end = clock();
	float seconds = (float)(end - start) / CLOCKS_PER_SEC;

	int min = (int)seconds / 60;
	int sec = (int)seconds % 60;

	std::cout << "Simulation completed in " << min <<  " min " << sec << " s." <<std::endl;

	getchar();
	return 0;
}
