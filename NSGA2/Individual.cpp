#include "Individual.h"
#include "Tools.h"
#include "ProblemParas.h"
#include "FitnessFunction.h"
//个体集合
Individual::Individual() {}

//初始化一个个体（geneSize是染色体基因的长度）
Individual::Individual(int geneSize, int fitnessCount, double* lower_bounds, double* upper_bounds)
{
	this->geneSize = geneSize;
	this->fitnessCount = fitnessCount;
	genes = new double[geneSize];
	for (int i = 0; i < geneSize; i++)
	{
		double tempGene = GetDoubleRand() * (upper_bounds[i] - lower_bounds[i]) + lower_bounds[i];
		genes[i] = tempGene;//构造初始种群
	}
}

//计算个体的所有适应度（这里的适应度函数很简单，后面要改）
void Individual::evaluation(int fitnessCount, ProblemParas proParas)
{
	Individual individual = *this;
	objectiveSet = FitnessFunction(individual, proParas);
}

//判断q是否支配当前染色体
bool Individual::dominate(Individual q)
{
	for (int i = 0; i < fitnessCount; i++)
	{
		if (objectiveSet[i] > q.objectiveSet[i])
		{
			return false;
		}
	}
	return true;
}