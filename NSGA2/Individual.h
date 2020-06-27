#pragma once
#include <random>
#include "BestPathInfo.h"
#include "ProblemParas.h"
using namespace std;

class Individual
{
public:
	using Comparator = bool(Individual::*)(Individual, Individual);

	double* genes;//存储染色体的数组（一个染色体是一堆基因组成的）
	int geneSize;//基因数组大小
	double* objectiveSet;//目标函数的值的集合
	int fitnessCount;//适应度数目
	vector<Individual*> dominatedSet;//支配其他解的解集
	int dominatedCount;//被支配解的数目
	int rank;//Pareto等级
	double distance;//解的距离（和拥挤度相反，距离越大拥挤度越小）

	Individual();
	void evaluation(int curIterNum, int maxIterNum, vector<BestPathInfo>& bestPathInfoList, ProblemParas& proParas);

	Individual(int geneSize, int fitnessCount, double* lower_bounds, double* upper_bounds, ProblemParas& problemParas);
	Individual(int s, int geneSize, int fitnessCount, vector<double> lower_bounds, vector<double> upper_bounds, ProblemParas& problemParas);
	bool dominate(Individual);
};
