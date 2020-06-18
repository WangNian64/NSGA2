#pragma once
#include "GAPara.h"
#include "Individual.h"
#ifndef POPULATION_H
#define POPULATION_H
using namespace std;

class Population
{
public:
	vector<Individual> individualSet;
	GAPara gaPara_pop;

	Population(GAPara gaPara);
	void initialize();
	void clear();
	Population copy();
	Population copy_all();
	Population combination(Population);

	void evaluation(int fitnessCount, ProblemParas proParas);
};

#endif // !POPULATION_H