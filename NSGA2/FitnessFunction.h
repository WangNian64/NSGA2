#pragma once
#include <math.h>
#include <cmath>
#define PI 3.1415926
#define DOUBLE_MAX 1.7976931348623158e+308
#define DOUBLE_MIN 2.2250738585072014e-308
#define MAX_FITNESS 100000000.0
double* FitnessFunction(Individual& individual, ProblemParas proParas);
double CalcuTotalArea(Individual& individual, ProblemParas proParas);
//默认的适应度计算函数，可以替换
double* FitnessFunction(Individual& individual, ProblemParas proParas)
{
	double* fitness = new double[individual.fitnessCount];
	fitness[0] = fitness[1] = 0.0;
	double deviceDist;

	for (int i = 0; i < individual.geneSize; i += 2) {
		for (int j = 0; j < individual.geneSize; j += 2) {
			deviceDist = abs(individual.genes[i] - individual.genes[j])
				+ abs(individual.genes[i + 1] - individual.genes[j + 1]);
			//if (deviceDist < proParas.MinDistArray[i / 2][j / 2])
			//{
			//	fitness[0] = fitness[1] = MAX_FITNESS;
			//	return fitness;
			//}
			//添加成本作为适应度值（考虑最小距离）
			fitness[0] += deviceDist * proParas.CostParaArray[i / 2][j / 2].UnitCost * proParas.CostParaArray[i / 2][j / 2].MatFlow;
		}
	}
	//如果重叠，加入惩罚
	for (int i = 0; i < individual.geneSize; i += 2) {
		double firstLowX = individual.genes[i] - 0.5 * proParas.DeviceSizeArray[i / 2].x;
		double firstUpX = individual.genes[i] + 0.5 * proParas.DeviceSizeArray[i / 2].x;
		double firstLowY = individual.genes[i + 1] - 0.5 * proParas.DeviceSizeArray[i / 2].y;
		double firstUpY = individual.genes[i + 1] + 0.5 * proParas.DeviceSizeArray[i / 2].y;
		for (int j = i + 2; j < individual.geneSize; j += 2) {
			double secondLowX = individual.genes[j] - 0.5 * proParas.DeviceSizeArray[j / 2].x;
			double secondUpX = individual.genes[j] + 0.5 * proParas.DeviceSizeArray[j / 2].x;
			double secondLowY = individual.genes[j + 1] - 0.5 * proParas.DeviceSizeArray[j / 2].y;
			double secondUpY = individual.genes[j + 1] + 0.5 * proParas.DeviceSizeArray[j / 2].y;
			if (IsRangeOverlap(firstLowX, firstUpX, secondLowX, secondUpX) && IsRangeOverlap(firstLowY, firstUpY, secondLowY, secondUpY)) {
				fitness[0] = fitness[1] = MAX_FITNESS;
				return fitness;
			}
		}
	}
	//面积也作为适应度值
	fitness[1] = CalcuTotalArea(individual, proParas);
	//fitness[0] = 1.0 / fitness[0];
	return fitness;
}
//计算占地面积
double CalcuTotalArea(Individual& individual, ProblemParas proParas) {
	double area = 0;
	double min_X, min_Y, max_X, max_Y;
	int min_X_index, min_Y_index, max_X_index, max_Y_index;
	min_X = min_Y = DOUBLE_MAX;
	max_X = max_Y = DOUBLE_MIN;
	min_X_index = min_Y_index = max_X_index = max_Y_index = 0;
	for (int i = 0; i < individual.geneSize; i += 2) {
		if (individual.genes[i] - proParas.DeviceSizeArray[i / 2].x * 0.5 < min_X) {
			min_X = individual.genes[i] - proParas.DeviceSizeArray[i / 2].x * 0.5;
			min_X_index = i / 2;
		}
		if (individual.genes[i + 1] - proParas.DeviceSizeArray[i / 2].y * 0.5 < min_Y) {
			min_Y = individual.genes[i + 1] - proParas.DeviceSizeArray[i / 2].y * 0.5;
			min_Y_index = i / 2;
		}
		if (individual.genes[i] + proParas.DeviceSizeArray[i / 2].x * 0.5 > max_X) {
			max_X = individual.genes[i] + proParas.DeviceSizeArray[i / 2].x * 0.5;
			max_X_index = i / 2;
		}
		if (individual.genes[i + 1] - proParas.DeviceSizeArray[i / 2].y * 0.5 > max_Y) {
			max_Y = individual.genes[i + 1] - proParas.DeviceSizeArray[i / 2].y * 0.5;
			max_Y_index = i / 2;
		}
	}
	//计算总面积
	area = (max_X - min_X) * (max_Y - min_Y);
	return area;
}
