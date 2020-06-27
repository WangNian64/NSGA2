#pragma once
#include "ProblemParas.h"
struct GAPara
{
	int geneSize;						//基因数目
	int fitnessCount;					//适应度数目
	int populationSize;					//种群大小
	int maxIterNum;						//最大迭代次数

	double mutationProb;				//变异概率
	double crossoverProb;				//交叉概率

	double* lower_bound_;				//基因随机下限
	double* upper_bound_;				//基因随机上限
	ProblemParas problemParas;			//对应的设备布局参数

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
	// 析构函数：释放堆内存
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

