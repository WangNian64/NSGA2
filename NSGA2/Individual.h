#pragma once
#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H
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
	vector<Individual*> dominatedSet;//��֧��⼯
	int dominatedCount;//��֧������Ŀ
	int rank;//Pareto�ȼ�
	double distance;//��ľ��루��ӵ�����෴��

	Individual();
	Individual(int geneSize, int fitnessCount, double* lower_bounds, double* upper_bounds);
	void evaluation(int fitnessCount, ProblemParas proParas);
	bool dominate(Individual);
};

#endif // !INDIVIDUAL_H
