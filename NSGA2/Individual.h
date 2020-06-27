#pragma once
#include <random>
#include "BestPathInfo.h"
#include "ProblemParas.h"
using namespace std;

class Individual
{
public:
	using Comparator = bool(Individual::*)(Individual, Individual);

	double* genes;//�洢Ⱦɫ������飨һ��Ⱦɫ����һ�ѻ�����ɵģ�
	int geneSize;//���������С
	double* objectiveSet;//Ŀ�꺯����ֵ�ļ���
	int fitnessCount;//��Ӧ����Ŀ
	vector<Individual*> dominatedSet;//֧��������Ľ⼯
	int dominatedCount;//��֧������Ŀ
	int rank;//Pareto�ȼ�
	double distance;//��ľ��루��ӵ�����෴������Խ��ӵ����ԽС��

	Individual();
	void evaluation(int curIterNum, int maxIterNum, vector<BestPathInfo>& bestPathInfoList, ProblemParas& proParas);

	Individual(int geneSize, int fitnessCount, double* lower_bounds, double* upper_bounds, ProblemParas& problemParas);
	Individual(int s, int geneSize, int fitnessCount, vector<double> lower_bounds, vector<double> upper_bounds, ProblemParas& problemParas);
	bool dominate(Individual);
};
