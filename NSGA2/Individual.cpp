#include "Individual.h"
#include "Tools.h"
#include "ProblemParas.h"
#include "FitnessFunction.h"
//���弯��
Individual::Individual() {}

//��ʼ��һ�����壨geneSize��Ⱦɫ�����ĳ��ȣ�
Individual::Individual(int geneSize, int fitnessCount, double* lower_bounds, double* upper_bounds, ProblemParas& problemParas)
{
	this->geneSize = geneSize;
	this->fitnessCount = fitnessCount;
	this->objectiveSet = new double[fitnessCount];
	for (int i = 0; i < fitnessCount; i++) {
		this->objectiveSet[i] = INT_MAX;
	}
	genes = new double[geneSize];

	//����һ�����±߽�
	double* lowerBounds = new double[geneSize];
	double* upperBounds = new double[geneSize];
	for (int i = 0; i < geneSize; ++i) {
		lowerBounds[i] = lower_bounds[i];
		upperBounds[i] = upper_bounds[i];
	}
	//ֱ������ģ�����
	//for (int i = 0; i < geneSize; i++)
	//{
	//	double tempGene = GetDoubleRand() * (upper_bounds[i] - lower_bounds[i]) + lower_bounds[i];
	//	genes[i] = tempGene;//�����ʼ��Ⱥ
	//}

	#pragma region ��ʼ��һ��Ⱦɫ��
	//���������Ȼ����ݳ�������豸�ߴ緶Χ
	for (int i = 2; i < geneSize; i += 3)
	{
		genes[i] = GetDoubleRand() * (upper_bounds[i] - lower_bounds[i]) + lower_bounds[i];
	}
	//���ݳ����޸��豸���½緶Χ&�豸����
	vector<Size> deviceSizeCopy;
	for (int i = 0; i < problemParas.DeviceSum; i++)
	{
		deviceSizeCopy.push_back(Size(problemParas.deviceParaList[i].size.x, problemParas.deviceParaList[i].size.y));
	}
	for (int i = 2; i < geneSize; i += 3)
	{
		//doubleתint��ת��ΪDirection��Ȼ����ݳ������¼����豸�ߴ�ͳ����
		//Rotate90����Rotate270,�޸�������
		DeviceDirect curDirect = (DeviceDirect)(int)genes[i];
		if (curDirect == DeviceDirect::Rotate90 || curDirect == DeviceDirect::Rotate270)
		{
			//x��y
			lowerBounds[i - 2] = 0 + problemParas.deviceParaList[i / 3].size.y * 0.5 + problemParas.deviceParaList[i / 3].spaceLength;
			lowerBounds[i - 1] = 0 + problemParas.deviceParaList[i / 3].size.x * 0.5 + problemParas.deviceParaList[i / 3].spaceLength;

			upperBounds[i - 2] = problemParas.workShopLength - problemParas.deviceParaList[i / 3].size.y * 0.5 - problemParas.deviceParaList[i / 3].spaceLength;
			upperBounds[i - 1] = problemParas.workShopWidth - problemParas.deviceParaList[i / 3].size.x * 0.5 - problemParas.deviceParaList[i / 3].spaceLength;

			//size��x��y��Ҫ����
			swap(deviceSizeCopy[i / 3].x, deviceSizeCopy[i / 3].y);
		}
		else
		{
			//x��y
			lowerBounds[i - 2] = 0 + problemParas.deviceParaList[i / 3].size.x * 0.5 + problemParas.deviceParaList[i / 3].spaceLength;
			lowerBounds[i - 1] = 0 + problemParas.deviceParaList[i / 3].size.y * 0.5 + problemParas.deviceParaList[i / 3].spaceLength;

			upperBounds[i - 2] = problemParas.workShopLength - problemParas.deviceParaList[i / 3].size.x * 0.5 - problemParas.deviceParaList[i / 3].spaceLength;
			upperBounds[i - 1] = problemParas.workShopWidth - problemParas.deviceParaList[i / 3].size.y * 0.5 - problemParas.deviceParaList[i / 3].spaceLength;

		}
	}

	#pragma region ���Ƿ��ص�Լ��������ֿ���������
	//(ÿ��1�ײ���һ������㣬ֻҪ�ҵ�һ�������������ص�Լ�����Ͳ��ã�����Ĭ��Ϊ0
	//�µ����������豸�İڷ�˳��
	vector<int> unmakeDeviceIndexVec;
	vector<int> madeDeviceIndexVec;
	for (int i = 0; i < problemParas.DeviceSum; i++)
	{
		unmakeDeviceIndexVec.push_back(i);
	}
	default_random_engine e;
	while (unmakeDeviceIndexVec.size() > 0)
	{
		//΢�뼶���ȵ����������
		e.seed(GetRamdonSeed());
		uniform_int_distribution<unsigned> u(0, unmakeDeviceIndexVec.size() - 1);
		int randomVecIndex = u(e);
		int randomDeviceIndex = unmakeDeviceIndexVec[randomVecIndex];//�õ��豸��index
		int j = randomDeviceIndex * 3;

		double Xstart = lowerBounds[j];
		double Ystart = lowerBounds[j + 1];

		double tempPositionX = 0;
		double tempPositionY = 0;

		bool findParticle = false;
		while (Ystart <= upperBounds[j + 1] - 1 && findParticle == false) {//X��YҪ�ڷ�Χ��
			Xstart = lowerBounds[j];
			while (Xstart <= upperBounds[j] - 1 && findParticle == false) {
				tempPositionX = GetDoubleRand() * 1.0 + Xstart;//�õ�Xstart��Xstart+1֮���һ�������
				tempPositionY = GetDoubleRand() * 1.0 + Ystart;//�õ�Ystart��Ystart+1֮���һ�������
				double halfX = deviceSizeCopy[j / 3].x * 0.5 + problemParas.deviceParaList[j / 3].spaceLength;
				double halfY = deviceSizeCopy[j / 3].y * 0.5 + problemParas.deviceParaList[j / 3].spaceLength;
				double tempLowX = tempPositionX - halfX;
				double tempUpX = tempPositionX + halfX;
				double tempLowY = tempPositionY - halfY;
				double tempUpY = tempPositionY + halfY;

				bool IsCross = false;
				//��鵱ǰ�豸�Ƿ��������ص�
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
					//�������ص����˳�
					if (IsRangeOverlap(tempLowX, tempUpX, curLowX, curUpX) && IsRangeOverlap(tempLowY, tempUpY, curLowY, curUpY)) {
						IsCross = true;
						break;
					}
				}
				//ȫ�����ص��������Ӹ�ֵ
				if (IsCross == false) {
					findParticle = true;
					genes[j] = tempPositionX;
					genes[j + 1] = tempPositionY;

					//����vec
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

//��ʼ��һ������(��ȫ���)
//���ҲҪע�����½緶Χ
Individual::Individual(int s, int geneSize, int fitnessCount, vector<double> lower_bounds, vector<double> upper_bounds, ProblemParas& problemParas)
{
	this->geneSize = geneSize;
	this->fitnessCount = fitnessCount;
	this->objectiveSet = new double[fitnessCount];
	for (int i = 0; i < fitnessCount; i++) {
		this->objectiveSet[i] = INT_MAX;
	}
	genes = new double[geneSize];
	
	#pragma region ��ʼ��һ��Ⱦɫ��
	//���������Ȼ����ݳ�������豸�ߴ緶Χ
	for (int i = 2; i < geneSize; i += 3)
	{
		genes[i] = GetDoubleRand() * (upper_bounds[i] - lower_bounds[i]) + lower_bounds[i];
	}
	//���ݳ����޸��豸���½緶Χ&�豸����
	for (int i = 2; i < geneSize; i += 3)
	{
		//doubleתint��ת��ΪDirection��Ȼ����ݳ������¼����豸�ߴ�ͳ����
		//Rotate90����Rotate270,�޸�������
		DeviceDirect curDirect = (DeviceDirect)(int)genes[i];
		if (curDirect == DeviceDirect::Rotate90 || curDirect == DeviceDirect::Rotate270)
		{
			//x��y
			lower_bounds[i - 2] = 0 + problemParas.deviceParaList[i / 3].size.y * 0.5 + problemParas.deviceParaList[i / 3].spaceLength;
			lower_bounds[i - 1] = 0 + problemParas.deviceParaList[i / 3].size.x * 0.5 + problemParas.deviceParaList[i / 3].spaceLength;

			upper_bounds[i - 2] = problemParas.workShopLength - problemParas.deviceParaList[i / 3].size.y * 0.5 - problemParas.deviceParaList[i / 3].spaceLength;
			upper_bounds[i - 1] = problemParas.workShopWidth - problemParas.deviceParaList[i / 3].size.x * 0.5 - problemParas.deviceParaList[i / 3].spaceLength;

		}
		else
		{
			//x��y
			lower_bounds[i - 2] = 0 + problemParas.deviceParaList[i / 3].size.x * 0.5 + problemParas.deviceParaList[i / 3].spaceLength;
			lower_bounds[i - 1] = 0 + problemParas.deviceParaList[i / 3].size.y * 0.5 + problemParas.deviceParaList[i / 3].spaceLength;

			upper_bounds[i - 2] = problemParas.workShopLength - problemParas.deviceParaList[i / 3].size.x * 0.5 - problemParas.deviceParaList[i / 3].spaceLength;
			upper_bounds[i - 1] = problemParas.workShopWidth - problemParas.deviceParaList[i / 3].size.y * 0.5 - problemParas.deviceParaList[i / 3].spaceLength;

		}
	}

	#pragma region �豸����ֱ�����
	for (int i = 0; i < geneSize; i++)
	{
		if (i % 3 != 2) {
			double tempGene = GetDoubleRand() * (upper_bounds[i] - lower_bounds[i]) + lower_bounds[i];
			genes[i] = tempGene;//�����ʼ��Ⱥ
		}
	}
	#pragma endregion

	#pragma endregion
}
//��������������Ӧ�ȣ��������Ӧ�Ⱥ����ܼ򵥣�����Ҫ�ģ�
void Individual::evaluation(int curIterNum, int maxIterNum, vector<BestPathInfo>& bestPathInfoList, ProblemParas& proParas)
{
	FitnessFunction(curIterNum, maxIterNum, bestPathInfoList, proParas, *this);
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