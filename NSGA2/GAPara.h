#pragma once
#include "ProblemParas.h"
struct GAPara
{
	int geneSize;						//������Ŀ
	int fitnessCount;					//��Ӧ����Ŀ
	int populationSize;					//��Ⱥ��С
	int maxIterNum;						//����������

	double mutationProb;
	double crossoverProb;

	double* lowerBounds;				//�����������
	double* upBounds;					//�����������
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
			copyOne.lowerBounds[i] = lowerBounds[i];
		}
		for (int i = 0; i < geneSize; i++) {
			copyOne.upBounds[i] = upBounds[i];
		}
		copyOne.problemParas = problemParas;
		return copyOne;
	}
	GAPara() {}
	GAPara(int geneSize)
	{
		this->geneSize = geneSize;
		lowerBounds = new double[geneSize];
		upBounds = new double[geneSize];
	}
	// �����������ͷŶ��ڴ�
	~GAPara()
	{
		//if (lowerBounds) { delete[]lowerBounds; }
		//if (upBounds) { delete[]upBounds; }
	}
	void setLowBounds(double lowBoundX, double lowBoundY)
	{
		for (int i = 0; i < geneSize; i++) {
			if (i % 2 == 0)
				lowerBounds[i] = lowBoundX + problemParas.DeviceSizeArray[i / 2].x * 0.5;
			else
				lowerBounds[i] = lowBoundY + problemParas.DeviceSizeArray[i / 2].y * 0.5;
		}
	}
	void setUpBounds(double upBoundX, double upBoundY)
	{
		for (int i = 0; i < geneSize; i++) {
			if (i % 2 == 0)
				upBounds[i] = upBoundX - problemParas.DeviceSizeArray[i / 2].x * 0.5;
			else
				upBounds[i] = upBoundY - problemParas.DeviceSizeArray[i / 2].y * 0.5;
		}
	}
};

