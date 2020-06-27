#include "Population.h"
#include "BestPathInfo.h"
#include "ProblemParas.h"
#include "Individual.h"

Population::Population(GAPara gaPara)
{
	this->gaPara_pop = gaPara;
}

//��ʼ��һ����Ⱥ
void Population::initialize()
{
	for (int i = 0; i < gaPara_pop.populationSize; i++)
	{
		individualSet.push_back(Individual(gaPara_pop.geneSize, gaPara_pop.fitnessCount, gaPara_pop.lower_bound_, gaPara_pop.upper_bound_, gaPara_pop.problemParas));
	}
}

void Population::clear()
{
	individualSet.clear();
}

Population Population::copy()
{
	return Population(gaPara_pop);
}

//����һ����ǰ��Ⱥ�Ŀ���
Population Population::copy_all()
{
	Population _tmpPop(gaPara_pop);
	_tmpPop.individualSet = individualSet;
	return _tmpPop;
}

Population Population::combination(Population& q)
{
	GAPara gaParaCopy = gaPara_pop.Copy();
	gaParaCopy.populationSize = gaPara_pop.populationSize * 2;
	gaParaCopy.geneSize = gaPara_pop.geneSize;
	//double newCrossoverProb = gaPara_pop.crossoverProb;
	//double newMutationProb = gaPara_pop.mutationProb;
	Population _tmp(gaParaCopy);
	_tmp.individualSet = individualSet;

	for (auto ind : q.individualSet)//�ں���׷���µ�Ԫ��
	{
		_tmp.individualSet.push_back(ind);
	}
	return _tmp;
}

void Population::evaluation(int curIterNum, int maxIterNum, vector<BestPathInfo>& bestPathInfoList, ProblemParas& proParas)
{
	for (int i = 0; i < individualSet.size(); i++)
	{
		individualSet[i].evaluation(curIterNum, maxIterNum, bestPathInfoList, proParas);
	}
}

