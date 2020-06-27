#include "Individual.h"
#include "Tools.h"
#include "ProblemParas.h"
#include "FitnessFunction.h"
//个体集合
Individual::Individual() {}

//初始化一个个体（geneSize是染色体基因的长度）
Individual::Individual(int geneSize, int fitnessCount, double* lower_bounds, double* upper_bounds, ProblemParas& problemParas)
{
	this->geneSize = geneSize;
	this->fitnessCount = fitnessCount;
	this->objectiveSet = new double[fitnessCount];
	for (int i = 0; i < fitnessCount; i++) {
		this->objectiveSet[i] = INT_MAX;
	}
	genes = new double[geneSize];

	//复制一次上下边界
	double* lowerBounds = new double[geneSize];
	double* upperBounds = new double[geneSize];
	for (int i = 0; i < geneSize; ++i) {
		lowerBounds[i] = lower_bounds[i];
		upperBounds[i] = upper_bounds[i];
	}
	//直接随机的，不够
	//for (int i = 0; i < geneSize; i++)
	//{
	//	double tempGene = GetDoubleRand() * (upper_bounds[i] - lower_bounds[i]) + lower_bounds[i];
	//	genes[i] = tempGene;//构造初始种群
	//}

	#pragma region 初始化一个染色体
	//先随机朝向，然后根据朝向调整设备尺寸范围
	for (int i = 2; i < geneSize; i += 3)
	{
		genes[i] = GetDoubleRand() * (upper_bounds[i] - lower_bounds[i]) + lower_bounds[i];
	}
	//根据朝向修改设备上下界范围&设备坐标
	vector<Size> deviceSizeCopy;
	for (int i = 0; i < problemParas.DeviceSum; i++)
	{
		deviceSizeCopy.push_back(Size(problemParas.deviceParaList[i].size.x, problemParas.deviceParaList[i].size.y));
	}
	for (int i = 2; i < geneSize; i += 3)
	{
		//double转int，转换为Direction，然后根据朝向重新计算设备尺寸和出入口
		//Rotate90或者Rotate270,修改上下限
		DeviceDirect curDirect = (DeviceDirect)(int)genes[i];
		if (curDirect == DeviceDirect::Rotate90 || curDirect == DeviceDirect::Rotate270)
		{
			//x和y
			lowerBounds[i - 2] = 0 + problemParas.deviceParaList[i / 3].size.y * 0.5 + problemParas.deviceParaList[i / 3].spaceLength;
			lowerBounds[i - 1] = 0 + problemParas.deviceParaList[i / 3].size.x * 0.5 + problemParas.deviceParaList[i / 3].spaceLength;

			upperBounds[i - 2] = problemParas.workShopLength - problemParas.deviceParaList[i / 3].size.y * 0.5 - problemParas.deviceParaList[i / 3].spaceLength;
			upperBounds[i - 1] = problemParas.workShopWidth - problemParas.deviceParaList[i / 3].size.x * 0.5 - problemParas.deviceParaList[i / 3].spaceLength;

			//size的x和y需要互换
			swap(deviceSizeCopy[i / 3].x, deviceSizeCopy[i / 3].y);
		}
		else
		{
			//x和y
			lowerBounds[i - 2] = 0 + problemParas.deviceParaList[i / 3].size.x * 0.5 + problemParas.deviceParaList[i / 3].spaceLength;
			lowerBounds[i - 1] = 0 + problemParas.deviceParaList[i / 3].size.y * 0.5 + problemParas.deviceParaList[i / 3].spaceLength;

			upperBounds[i - 2] = problemParas.workShopLength - problemParas.deviceParaList[i / 3].size.x * 0.5 - problemParas.deviceParaList[i / 3].spaceLength;
			upperBounds[i - 1] = problemParas.workShopWidth - problemParas.deviceParaList[i / 3].size.y * 0.5 - problemParas.deviceParaList[i / 3].spaceLength;

		}
	}

	#pragma region 考虑非重叠约束，这里分块产生随机点
	//(每隔1米产生一个随机点，只要找到一个随机点满足非重叠约束，就采用）朝向默认为0
	//新的随机：随机设备的摆放顺序
	vector<int> unmakeDeviceIndexVec;
	vector<int> madeDeviceIndexVec;
	for (int i = 0; i < problemParas.DeviceSum; i++)
	{
		unmakeDeviceIndexVec.push_back(i);
	}
	default_random_engine e;
	while (unmakeDeviceIndexVec.size() > 0)
	{
		//微秒级精度的随机数种子
		e.seed(GetRamdonSeed());
		uniform_int_distribution<unsigned> u(0, unmakeDeviceIndexVec.size() - 1);
		int randomVecIndex = u(e);
		int randomDeviceIndex = unmakeDeviceIndexVec[randomVecIndex];//得到设备的index
		int j = randomDeviceIndex * 3;

		double Xstart = lowerBounds[j];
		double Ystart = lowerBounds[j + 1];

		double tempPositionX = 0;
		double tempPositionY = 0;

		bool findParticle = false;
		while (Ystart <= upperBounds[j + 1] - 1 && findParticle == false) {//X和Y要在范围内
			Xstart = lowerBounds[j];
			while (Xstart <= upperBounds[j] - 1 && findParticle == false) {
				tempPositionX = GetDoubleRand() * 1.0 + Xstart;//得到Xstart到Xstart+1之间的一个随机数
				tempPositionY = GetDoubleRand() * 1.0 + Ystart;//得到Ystart到Ystart+1之间的一个随机数
				double halfX = deviceSizeCopy[j / 3].x * 0.5 + problemParas.deviceParaList[j / 3].spaceLength;
				double halfY = deviceSizeCopy[j / 3].y * 0.5 + problemParas.deviceParaList[j / 3].spaceLength;
				double tempLowX = tempPositionX - halfX;
				double tempUpX = tempPositionX + halfX;
				double tempLowY = tempPositionY - halfY;
				double tempUpY = tempPositionY + halfY;

				bool IsCross = false;
				//检查当前设备是否与其他重叠
				for (int k = 0; k < madeDeviceIndexVec.size(); k++)
				{
					int curDeviceIndex = madeDeviceIndexVec[k];
					int curDimIndex = curDeviceIndex * 3;

					double halfX1 = deviceSizeCopy[curDeviceIndex].x * 0.5 + problemParas.deviceParaList[curDeviceIndex].spaceLength;
					double halfY1 = deviceSizeCopy[curDeviceIndex].y * 0.5 + problemParas.deviceParaList[curDeviceIndex].spaceLength;

					double curLowX = genes[curDimIndex] - halfX1;
					double curUpX = genes[curDimIndex] + halfX1;
					double curLowY = genes[curDimIndex + 1] - halfY1;
					double curUpY = genes[curDimIndex + 1] + halfY1;
					//若发生重叠，退出
					if (IsRangeOverlap(tempLowX, tempUpX, curLowX, curUpX) && IsRangeOverlap(tempLowY, tempUpY, curLowY, curUpY)) {
						IsCross = true;
						break;
					}
				}
				//全部不重叠，给粒子赋值
				if (IsCross == false) {
					findParticle = true;
					genes[j] = tempPositionX;
					genes[j + 1] = tempPositionY;

					//更新vec
					madeDeviceIndexVec.push_back(randomDeviceIndex);
					unmakeDeviceIndexVec.erase(unmakeDeviceIndexVec.begin() + randomVecIndex);

				}
				Xstart++;
				if (Xstart >= upperBounds[j] - 1) {
					Ystart++;
				}
			}
		}
	}
	#pragma endregion

	#pragma endregion
}

//初始化一个个体(完全随机)
//这个也要注意上下界范围
Individual::Individual(int s, int geneSize, int fitnessCount, vector<double> lower_bounds, vector<double> upper_bounds, ProblemParas& problemParas)
{
	this->geneSize = geneSize;
	this->fitnessCount = fitnessCount;
	this->objectiveSet = new double[fitnessCount];
	for (int i = 0; i < fitnessCount; i++) {
		this->objectiveSet[i] = INT_MAX;
	}
	genes = new double[geneSize];
	
	#pragma region 初始化一个染色体
	//先随机朝向，然后根据朝向调整设备尺寸范围
	for (int i = 2; i < geneSize; i += 3)
	{
		genes[i] = GetDoubleRand() * (upper_bounds[i] - lower_bounds[i]) + lower_bounds[i];
	}
	//根据朝向修改设备上下界范围&设备坐标
	for (int i = 2; i < geneSize; i += 3)
	{
		//double转int，转换为Direction，然后根据朝向重新计算设备尺寸和出入口
		//Rotate90或者Rotate270,修改上下限
		DeviceDirect curDirect = (DeviceDirect)(int)genes[i];
		if (curDirect == DeviceDirect::Rotate90 || curDirect == DeviceDirect::Rotate270)
		{
			//x和y
			lower_bounds[i - 2] = 0 + problemParas.deviceParaList[i / 3].size.y * 0.5 + problemParas.deviceParaList[i / 3].spaceLength;
			lower_bounds[i - 1] = 0 + problemParas.deviceParaList[i / 3].size.x * 0.5 + problemParas.deviceParaList[i / 3].spaceLength;

			upper_bounds[i - 2] = problemParas.workShopLength - problemParas.deviceParaList[i / 3].size.y * 0.5 - problemParas.deviceParaList[i / 3].spaceLength;
			upper_bounds[i - 1] = problemParas.workShopWidth - problemParas.deviceParaList[i / 3].size.x * 0.5 - problemParas.deviceParaList[i / 3].spaceLength;

		}
		else
		{
			//x和y
			lower_bounds[i - 2] = 0 + problemParas.deviceParaList[i / 3].size.x * 0.5 + problemParas.deviceParaList[i / 3].spaceLength;
			lower_bounds[i - 1] = 0 + problemParas.deviceParaList[i / 3].size.y * 0.5 + problemParas.deviceParaList[i / 3].spaceLength;

			upper_bounds[i - 2] = problemParas.workShopLength - problemParas.deviceParaList[i / 3].size.x * 0.5 - problemParas.deviceParaList[i / 3].spaceLength;
			upper_bounds[i - 1] = problemParas.workShopWidth - problemParas.deviceParaList[i / 3].size.y * 0.5 - problemParas.deviceParaList[i / 3].spaceLength;

		}
	}

	#pragma region 设备坐标直接随机
	for (int i = 0; i < geneSize; i++)
	{
		if (i % 3 != 2) {
			double tempGene = GetDoubleRand() * (upper_bounds[i] - lower_bounds[i]) + lower_bounds[i];
			genes[i] = tempGene;//构造初始种群
		}
	}
	#pragma endregion

	#pragma endregion
}
//计算个体的所有适应度（这里的适应度函数很简单，后面要改）
void Individual::evaluation(int curIterNum, int maxIterNum, vector<BestPathInfo>& bestPathInfoList, ProblemParas& proParas)
{
	FitnessFunction(curIterNum, maxIterNum, bestPathInfoList, proParas, *this);
}

//判断q是否支配当前染色体
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