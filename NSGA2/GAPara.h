#pragma once
#include "ProblemParas.h"
struct GAPara
{
	int geneSize;						//������Ŀ
	int fitnessCount;					//��Ӧ����Ŀ
	int populationSize;					//��Ⱥ��С
	int maxIterNum;						//����������

	double mutationProb;				//�������
	double crossoverProb;				//�������

	double* lower_bound_;				//�����������
	double* upper_bound_;				//�����������
	ProblemParas problemParas;			//��Ӧ���豸���ֲ���

	GAPara Copy()
	{
		GAPara copyOne(geneSize);
		copyOne.fitnessCount = fitnessCount;
		copyOne.populationSize = populationSize;
		copyOne.maxIterNum = maxIterNum;
		copyOne.mutationProb = mutationProb;
		copyOne.crossoverProb = crossoverProb;

		for (int i = 0; i < geneSize; i++) {
			copyOne.lower_bound_[i] = lower_bound_[i];
		}
		for (int i = 0; i < geneSize; i++) {
			copyOne.upper_bound_[i] = upper_bound_[i];
		}
		copyOne.problemParas = problemParas;
		return copyOne;
	}
	GAPara() {}
	GAPara(int geneSize)
	{
		this->geneSize = geneSize;
		lower_bound_ = new double[geneSize];
		upper_bound_ = new double[geneSize];
	}
	// �����������ͷŶ��ڴ�
	~GAPara()
	{
		//if (lower_bound_) { delete[] lower_bound_; }
		//if (upper_bound_) { delete[] upper_bound_; }
	}
	void setLowBounds(double lowBoundX, double lowBoundY, int lowBoundDirect)
	{
		for (int i = 0; i < geneSize; i++) {
			if (i % 3 == 0)
				lower_bound_[i] = lowBoundX + problemParas.deviceParaList[i / 3].size.x * 0.5 + problemParas.deviceParaList[i / 3].spaceLength;
			else if (i % 3 == 1)
				lower_bound_[i] = lowBoundY + problemParas.deviceParaList[i / 3].size.y * 0.5 + problemParas.deviceParaList[i / 3].spaceLength;
			else
				lower_bound_[i] = lowBoundDirect;
		}
	}
	void setUpBounds(double upBoundX, double upBoundY, int upBoundDirect)
	{
		for (int i = 0; i < geneSize; i++) {
			if (i % 3 == 0)
				upper_bound_[i] = upBoundX - problemParas.deviceParaList[i / 3].size.x * 0.5 - problemParas.deviceParaList[i / 3].spaceLength;
			else if (i % 3 == 1)
				upper_bound_[i] = upBoundY - problemParas.deviceParaList[i / 3].size.y * 0.5 - problemParas.deviceParaList[i / 3].spaceLength;
			else
				upper_bound_[i] = upBoundDirect;
		}
	}
};

