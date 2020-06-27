#pragma once
#include "DevicePara.h"
#include <set>
//存储最优粒子的输送线信息
struct BestPathInfo
{
	double curBestFitnessVal;//当前目标的最优值，用于判断是否更新
	vector<double> devicePosList;//设备位置列表
	vector<InoutPoint> inoutPoints;//出入口集合
	set<StraightConveyorInfo> strConveyorList;//直线输送机信息列表
	set<Vector2Int> curveConveyorList;//转弯输送机信息列表
	BestPathInfo() {
		curBestFitnessVal = INT_MAX;
	}
	BestPathInfo(int geneSize) {
		curBestFitnessVal = INT_MAX;
		devicePosList = vector<double>(geneSize);
	}
	void Clear() {
		devicePosList.clear();
		inoutPoints.clear();
		strConveyorList.clear();
		curveConveyorList.clear();
	}
};