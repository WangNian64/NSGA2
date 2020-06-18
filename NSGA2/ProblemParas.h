#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include "Tools.h"
using namespace std;
//�豸�ߴ�
struct Size {
	double x;
	double y;
};
//�ɱ��������
struct CostPara {
	double MatFlow;		//��������
	double UnitCost;	//��λ�������ĳɱ�
};
struct ProblemParas
{
	int DeviceSum;					//�豸����
	Size* DeviceSizeArray;			//�豸�ߴ�����

	Size WorkshopSize;				//����ߴ�

	CostPara** CostParaArray;		//�ɱ������������

	double** MinDistArray;			//�豸��С��������
	ProblemParas() {}

	ProblemParas(int deviceNum)
	{
		DeviceSum = deviceNum;

		DeviceSizeArray = new Size[deviceNum];

		CostParaArray = new CostPara * [DeviceSum];
		for (int i = 0; i < DeviceSum; i++) {
			CostParaArray[i] = new CostPara[DeviceSum];
		}

		MinDistArray = new double* [DeviceSum];
		for (int i = 0; i < DeviceSum; i++) {
			MinDistArray[i] = new double[DeviceSum];
		}




		ifstream fileIn("../../Para.txt");
		string line;
		if (fileIn) // �и��ļ�
		{
			//����ߴ�
			getline(fileIn, line);//��һ��
			getline(fileIn, line);
			vector<string> shopSizeStr = split(line, ",");
			WorkshopSize.x = atof(shopSizeStr[0].c_str());
			WorkshopSize.y = atof(shopSizeStr[1].c_str());


			//�豸�ߴ�
			getline(fileIn, line);//��һ��
			getline(fileIn, line);
			vector<string> strSplit = split(line, " ");
			for (int i = 0; i < DeviceSum; i++)
			{
				vector<string> deviceSizeStr = split(strSplit[i], ",");
				DeviceSizeArray[i].x = atof(deviceSizeStr[0].c_str());
				DeviceSizeArray[i].y = atof(deviceSizeStr[1].c_str());
			}


			//��λ����ɱ�����
			getline(fileIn, line);//��һ��
			for (int i = 0; i < DeviceSum; i++) {
				getline(fileIn, line);
				vector<string> strSplit = split(line, ",");
				for (int j = 0; j < DeviceSum; j++)
				{
					CostParaArray[i][j].UnitCost = atof(strSplit[j].c_str());
				}
			}

			//������������
			getline(fileIn, line);//��һ��
			for (int i = 0; i < DeviceSum; i++) {
				getline(fileIn, line);
				vector<string> strSplit = split(line, ",");
				for (int j = 0; j < DeviceSum; j++)
				{
					CostParaArray[i][j].MatFlow = atof(strSplit[j].c_str());
				}
			}

			//�豸��С��������
			getline(fileIn, line);//��һ��
			for (int i = 0; i < DeviceSum; i++) {
				getline(fileIn, line);
				vector<string> strSplit = split(line, ",");
				for (int j = 0; j < DeviceSum; j++)
				{
					MinDistArray[i][j] = atof(strSplit[j].c_str());
					if (i != j) {
						MinDistArray[i][j] += 2;
					}
				}
				cout << endl;
			}
		}
		else // û�и��ļ�
		{
			cout << "no such file" << endl;
		}
	}
};
