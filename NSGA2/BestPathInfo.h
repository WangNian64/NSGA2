#pragma once
#include "DevicePara.h"
#include <set>
//�洢�������ӵ���������Ϣ
struct BestPathInfo
{
	double curBestFitnessVal;//��ǰĿ�������ֵ�������ж��Ƿ����
	vector<double> devicePosList;//�豸λ���б�
	vector<InoutPoint> inoutPoints;//����ڼ���
	set<StraightConveyorInfo> strConveyorList;//ֱ�����ͻ���Ϣ�б�
	set<Vector2Int> curveConveyorList;//ת�����ͻ���Ϣ�б�
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