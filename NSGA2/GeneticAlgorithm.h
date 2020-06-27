#pragma once
#include "Population.h"
#include "ProblemParas.h"
#include "BestPathInfo.h"
#ifndef GENETIC_ALGORITHM_H
#define GENETIC_ALGORITHM_H

typedef vector<Individual*> Front;

using namespace std;
class GeneticAlgorithm
{
public:
	int curIterNum;
	GAPara gaPara;
	vector<BestPathInfo> bestPathInfoList;		//最优路径信息

	void wheelSelection(Population&);
	void crossover(Population&);
	void mutation(Population&);
	void SBX(int, Population&);
	void PLM(int index, double* lowBounds, double* upBounds, Population& p);
	vector<Front> fastNonDominatedSort(Population*);
	void crowdingDistanceAssignment(Front);

	GeneticAlgorithm(GAPara gaPara);
	void setLowBounds();
	void setUpperBounds();

	Population InitialPopulation();
	Population GetChildPopulation(Population parentPop);
	//Population UpdatePopulation();
	void normalize();
	void associate();
	void niching();
};

#endif // !GENETIC_ALGORITHM_H
