// NSGA2.cpp : 此文件包含“main”函数。在那里开始执行程序结束
#include <iomanip>
#include <ctime>
#include "Population.h"
#include "Individual.h"
#include "GeneticAlgorithm.h"
#include "ProblemParas.h"

typedef vector<Individual*> Front;

bool descending(Individual*, Individual*);

int main()
{
	int problemSum = 1;//问题的数目
	int testSum = 1;//每个实验跑的实验次数
	for (int curProblem = 0; curProblem < problemSum; curProblem++)//跑多个问题
	{
		cout << "跑第" + to_string(curProblem + 1) + "个问题" << endl;
		#pragma region 设置nsga2参数
		ifstream inputFile;
		inputFile.open("../../InputParas/InputPara" + to_string(curProblem + 1) + ".txt");

		//Para.txt参数读取
		ProblemParas proParas(inputFile);

		//设置遗传算法参数
		int geneSize = proParas.DeviceSum * 3;
		GAPara gaPara(geneSize);
		gaPara.fitnessCount = 2;
		gaPara.populationSize = 100;
		gaPara.maxIterNum = 400;
		gaPara.crossoverProb = 0.9;
		gaPara.mutationProb = 0.1;
		gaPara.problemParas = proParas;
		gaPara.setLowBounds(0, 0, DeviceDirect::Default);
		gaPara.setUpBounds(proParas.workShopLength, proParas.workShopWidth, DeviceDirect::Rotate270 + 1);

		#pragma endregion

		#pragma region 调用算法，并输出结果
		GeneticAlgorithm ga(gaPara);
		ga.curIterNum = 0;
		std::srand(GetRamdonSeed());

		string curProblemFolderName = "Problem" + to_string(curProblem + 1);
		for (int curTest = 0; curTest < testSum; curTest++) {//每个问题跑多次
			clock_t startTime, endTime;//记录调用时间

			startTime = clock();//计时开始
			#pragma region 初始化
			Population _tmpParent = ga.InitialPopulation();
			#pragma endregion

			#pragma region 迭代
			//目标1的值放在archiveList1中，目标2的值放在archiveList2中
			//第n次实验放到文件里去
			//文件夹的名字叫Testn
			ofstream OutFile;
			ofstream OutFile1;
			string curTestFolderName = "Test" + to_string(curTest + 1);
			OutFile.open("../../Results/" + curProblemFolderName + "/" + curTestFolderName + "/archiveList1.txt");
			OutFile1.open("../../Results/" + curProblemFolderName + "/" + curTestFolderName + "/archiveList2.txt");
			//开始遗传算法的迭代
			for (int i = 0; i < gaPara.maxIterNum; i++)
			{
				cout << (i + 1) << endl;
				Population _tmpChild = ga.GetChildPopulation(_tmpParent);//得到子种群
				Population _tmpComb = _tmpParent.combination(_tmpChild);//子代与父代混合得到总的种群
				vector<Front> front = ga.fastNonDominatedSort(&_tmpComb);//通过快速非支配排序得到新的Pareto前沿的数组
				_tmpParent.clear();//释放父类的空间

				//拥挤度排序
				int j = 0;
				//front[0]最好，后面越来越差
				while (_tmpParent.individualSet.size() + front[j].size() <= gaPara.populationSize)
				{
					if (front[j].size() == 0) { cin.get(); break; }
					//计算当前种群所有pareto前沿中所有染色体的拥挤度，并按照拥挤度由低到高排序
					//拥挤度越低，被选中的概率越高，反之越低，这里使用的是类似于轮盘赌的方法
					ga.crowdingDistanceAssignment(front[j]);
					//将个体加入到父代的种群中，按照rank递增的顺序
					for (Individual* ind : front[j])
					{
						_tmpParent.individualSet.push_back(*ind);
					}
					j++;
				}
				++ga.curIterNum;//迭代数目+1
				//对所有的Pareto前沿进行排序，排序规则：rank低的在前，或者rank相等但是解的距离大
				sort(front[j].begin(), front[j].end(), descending);

				int size = _tmpParent.individualSet.size();//size是更新后的父代种群的个体数目
				for (int k = 0; k < (gaPara.populationSize - size); k++)//再次选择一些个体加入到父种群中，直到数目达到上限
				{
					_tmpParent.individualSet.push_back(*front[j][k]);
				}
				//存储每次迭代的最优个体（分为目标1和目标2）
				string f1line = to_string(ga.bestPathInfoList[0].curBestFitnessVal) + "\n";
				OutFile << f1line;
				string f2line = to_string(ga.bestPathInfoList[1].curBestFitnessVal) + "\n";
				OutFile1 << f2line;

			}
			OutFile.close();
			OutFile1.close();
			#pragma endregion


			endTime = clock();
			cout << "迭代" << gaPara.maxIterNum << "次的最终用时:" << static_cast<double>(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;

			#pragma region 保存设备尺寸&最终布局结果&连线点的坐标
			OutFile.open("../../Results/" + curProblemFolderName + "/" + curTestFolderName + "/FinalResult.txt");
			#pragma region 记录最终布局结果	
			int resultIndex = 0;
			int minHandleCost = INT_MAX;
			int minConveyValue = INT_MAX;
			//优先选物料运输成本最低的
			for (int i = 0; i < _tmpParent.individualSet.size(); i++)
			{
				//cout << _tmpParent.individualSet[i].objectiveSet[0] << endl;
				if (_tmpParent.individualSet[i].objectiveSet[0] < minHandleCost)
				{
					minHandleCost = _tmpParent.individualSet[i].objectiveSet[0];
					resultIndex = i;
				}
			}
			//优先选输送机成本低的
			//for (int i = 0; i < _tmpParent.individualSet.size(); i++)
			//{
			//	if (_tmpParent.individualSet[i].objectiveSet[1] < minHandleCost)
			//	{
			//		minConveyValue = _tmpParent.individualSet[i].objectiveSet[1];
			//		resultIndex = i;
			//	}
			//}

			//for (int i = 0; i < geneSize; i += 3)
			//{
			//	OutFile /*<< fixed << setprecision(1)*/ << _tmpParent.individualSet[resultIndex].genes[i];
			//	OutFile << ",";
			//	OutFile /*<< fixed << setprecision(1)*/ << _tmpParent.individualSet[resultIndex].genes[i + 1];
			//	OutFile << "\n";
			//}
			for (int i = 0; i < geneSize; i += 3)
			{
				OutFile /*<< fixed << setprecision(1)*/ << ga.bestPathInfoList[0].devicePosList[i];
				OutFile << ",";
				OutFile /*<< fixed << setprecision(1)*/ << ga.bestPathInfoList[0].devicePosList[i + 1];
				OutFile << "\n";
			}
			#pragma endregion

			#pragma region 记录设备尺寸
			for (int i = 2; i < geneSize; i += 3)
			{
				//DeviceDirect direct = (DeviceDirect)(int)_tmpParent.individualSet[resultIndex].genes[i];
				DeviceDirect direct = (DeviceDirect)(int)ga.bestPathInfoList[0].devicePosList[i];
				string line = "";
				if (direct == DeviceDirect::Rotate90 || direct == DeviceDirect::Rotate270)
				{
					line = to_string(ga.gaPara.problemParas.deviceParaList[i / 3].size.y) + "," +
						to_string(ga.gaPara.problemParas.deviceParaList[i / 3].size.x);
				}
				else {
					line = to_string(ga.gaPara.problemParas.deviceParaList[i / 3].size.x) + "," +
						to_string(ga.gaPara.problemParas.deviceParaList[i / 3].size.y);
				}
				OutFile << line + "\n";
			}
			#pragma endregion

			int fitnessIndex = 0;
			#pragma region 记录出入口坐标（旋转之后的，不带设备坐标）
			vector<InoutPoint> ioPoints = ga.bestPathInfoList[fitnessIndex].inoutPoints;
			OutFile << to_string(ioPoints.size()) + "\n";//出入口数目
			for (int i = 0; i < ioPoints.size(); i++)
			{
				//string line = "";
				if (ioPoints[i].pointDirect == PointDirect::Up || ioPoints[i].pointDirect == PointDirect::Down)
				{
					//line += "Vertical ";
					OutFile << "Vertical ";
				}
				else
				{
					//line += "Horizon ";
					OutFile << "Horizon ";
				}
				//line += to_string(ioPoints[i].pointAxis.x) + " " + to_string(ioPoints[i].pointAxis.y) + " \n";
				//OutFile << line;
				OutFile /*<< fixed << setprecision(1)*/ << ioPoints[i].pointAxis.x;
				OutFile << " ";
				OutFile /*<< fixed << setprecision(1)*/ << ioPoints[i].pointAxis.y;
				OutFile << "\n";
			}
			#pragma endregion

			#pragma region 记录出入口路径
			#pragma endregion

			#pragma region 记录直线输送机和转弯输送机参数
			cout << ga.bestPathInfoList[fitnessIndex].curBestFitnessVal << endl;
			set<StraightConveyorInfo> strInfoList = ga.bestPathInfoList[fitnessIndex].strConveyorList;
			set<Vector2Int> curveInfoList = ga.bestPathInfoList[fitnessIndex].curveConveyorList;
			OutFile << strInfoList.size() << "\n";
			for (StraightConveyorInfo sci : strInfoList)
			{
				OutFile << to_string(sci.startPos.x) << "," << to_string(sci.startPos.y)
					<< ";" << to_string(sci.endPos.x) << "," << to_string(sci.endPos.y)
					<< ";" << to_string(sci.startHnum) << ";" << to_string(sci.startVnum)
					<< ";" << to_string(sci.endHnum) << ";" << to_string(sci.endVnum)
					<< "\n";
			}
			OutFile << curveInfoList.size() << "\n";
			for (Vector2Int v : curveInfoList)
			{
				OutFile << to_string(v.x) << "," << to_string(v.y) << "\n";
			}


			//set<SegPath> segPathSet = _tmpParent.individualSet[resultIndex].segPathSet;
			//OutFile << segPathSet.size() << "\n";
			//for (SegPath sp : segPathSet)
			//{
			//	OutFile << to_string(sp.p1.x) << "," << to_string(sp.p1.y)
			//			<< ";" << to_string(sp.p2.x) << "," << to_string(sp.p2.y) << "\n";
			//}

			//map<Vector2, PointInfo> pathPointInfoMap = _tmpParent.individualSet[resultIndex].pathPointInfoMap;
			//OutFile << pathPointInfoMap.size() << "\n";
			//for (auto it = pathPointInfoMap.begin(); it != pathPointInfoMap.end(); it++) 
			//{
			//	OutFile << to_string(it->first.x) << "," << to_string(it->first.y)
			//		<< ";" << to_string(it->second.horiDirNum) << ";" << to_string(it->second.vertDirNum) << "\n";
			//}
			#pragma endregion

			OutFile.close();
			#pragma endregion
		}
		
		
		#pragma endregion
	}

	system("pause");
}

//不同的Pareto集合，按照Rank值进行排序
bool descending(Individual* a, Individual* b)
{
	if ((a->rank < b->rank) || ((a->rank == b->rank) && (a->distance > b->distance)))
	{
		return true;
	}
	else
	{
		return false;
	}
}