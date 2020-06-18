#include "Individual.h"
#include "Tools.h"
#include "ProblemParas.h"
#include "FitnessFunction.h"
//���弯��
Individual::Individual() {}

//��ʼ��һ�����壨geneSize��Ⱦɫ�����ĳ��ȣ�
Individual::Individual(int geneSize, int fitnessCount, double* lower_bounds, double* upper_bounds)
{
	this->geneSize = geneSize;
	this->fitnessCount = fitnessCount;
	genes = new double[geneSize];
	for (int i = 0; i < geneSize; i++)
	{
		double tempGene = GetDoubleRand() * (upper_bounds[i] - lower_bounds[i]) + lower_bounds[i];
		genes[i] = tempGene;//�����ʼ��Ⱥ
	}
}

//��������������Ӧ�ȣ��������Ӧ�Ⱥ����ܼ򵥣�����Ҫ�ģ�
void Individual::evaluation(int fitnessCount, ProblemParas proParas)
{
	Individual individual = *this;
	objectiveSet = FitnessFunction(individual, proParas);
}

//�ж�q�Ƿ�֧�䵱ǰȾɫ��
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