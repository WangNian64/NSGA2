// NSGA2.cpp : 此文件包含“main”函数。在那里开始执行程序结束
//

#include "Population.h"
#include "Individual.h"
#include "GeneticAlgorithm.h"
#include "ProblemParas.h"

#define POPULATION_SIZE 100
#define GENE_SIZE 1//染色体的基因数目
#define CROSSOVER_PROB 0.9
#define MUTATION_PROB (1 / GENE_SIZE)

typedef vector<Individual*> Front;

bool descending(Individual*, Individual*);

////保存结果到txt
//void SaveLayoutResults(PSOOptimizer psooptimizer, double* result) {
//	ofstream OutFile;
//	//存每一次的适应度值
//	OutFile.open("/多目标算法集合/IterateResult.txt");
//	for (int i = 0; i < psooptimizer.max_iter_num_; i++)
//	{
//		string line = to_string(1.0 / result[i]) + "\n";
//		OutFile << line;
//	}
//	OutFile.close();
//	//存最后的布局结果
//	OutFile.open("/多目标算法集合/LayoutResult.txt");
//	for (int i = 0; i < psooptimizer.dim_; i += 2) {
//		string line = to_string(psooptimizer.all_best_position_[i]) + "," + to_string(psooptimizer.all_best_position_[i + 1]) + "\n";
//		OutFile << line;
//	}
//	OutFile.close();
//}
int main()
{

	//Para.txt参数读取
	int deviceNum = 6;
	ProblemParas proParas(deviceNum);

	//设置遗传算法参数
	int geneSize = deviceNum * 2;
	GAPara gaPara(geneSize);
	gaPara.fitnessCount = 2;
	gaPara.populationSize = 100;
	gaPara.maxIterNum = 100;

	gaPara.crossoverProb = CROSSOVER_PROB;
	gaPara.mutationProb = MUTATION_PROB;
	gaPara.problemParas = proParas;
	gaPara.setLowBounds(0, 0);
	gaPara.setUpBounds(proParas.WorkshopSize.x, proParas.WorkshopSize.y);

	//初始化算法
	GeneticAlgorithm ga(gaPara);

	//初始化种群
	Population _tmpParent = ga.InitialPopulation();
	//for (int i = 0; i < ga.gaPara.populationSize; i++) 
	//{
	//	cout << _tmpParent.individualSet[i].genes[0] << ",";
	//}

	double* object1Results = new double[gaPara.maxIterNum];
	double* object2Results = new double[gaPara.maxIterNum];
	//开始遗传算法的迭代
	for (int i = 0; i < gaPara.maxIterNum; i++)
	{
		//得到子种群
		Population _tmpChild = ga.GetChildPopulation(_tmpParent);

		//子代与父代混合得到总的种群
		Population _tmpComb = _tmpParent.combination(_tmpChild);

		//通过快速非支配排序得到新的Pareto前沿的数组
		vector<Front> front = ga.fastNonDominatedSort(&_tmpComb);
		_tmpParent.clear();//释放空间

		//拥挤度排序
		int j = 0;
		//front[0]最好，后面越来越差
		while (_tmpParent.individualSet.size() + front[j].size() <= POPULATION_SIZE)//
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
		//对所有的Pareto前沿进行排序
		//排序规则：rank低的在前，或者rank相等但是解的距离大
		sort(front[j].begin(), front[j].end(), descending);

		int size = _tmpParent.individualSet.size();//size是更新后的父代种群的个体数目
		for (int k = 0; k < (POPULATION_SIZE - size); k++)//再次选择一些个体加入到父种群中，直到数目达到上限
		{
			_tmpParent.individualSet.push_back(*front[j][k]);
		}
		//for (int i = 0; i < _tmpParent.individualSet.size(); i++) 
		//{
			//for (int j = 0; j < _tmpParent.geneSize; j++)
			//{
			//	cout << _tmpParent.individualSet[0].genes[j] << ",";
			//}
			//cout << endl;
		//}
		//cout << endl;
		//cout << "Parent size: " << _tmpParent.individualSet.size() << endl;
		//cout << "Complete generation [" << i << "]" << endl;

		//找出当前种群中适应度1最小的
		int smallIndex = 0;
		double obj1Small = _tmpParent.individualSet[0].objectiveSet[0];
		for (int i = 1; i < gaPara.populationSize; i++)
		{
			if (_tmpParent.individualSet[i].objectiveSet[0] < obj1Small)
			{
				smallIndex = i;
				obj1Small = _tmpParent.individualSet[i].objectiveSet[0];
			}
		}
		object1Results[i] = _tmpParent.individualSet[smallIndex].objectiveSet[0];
		object2Results[i] = _tmpParent.individualSet[smallIndex].objectiveSet[1];

	}


#pragma region 输出结果到文件中
	//输出每一代的迭代结果
	ofstream os("IterateResult.txt");
	if (os.is_open())
	{
		for (int i = 0; i < gaPara.maxIterNum; i++)
		{
			os << object1Results[i] << "," << object2Results[i] << endl;
		}
	}
	os.close();
	//输出最后的布局结果
	os.open("LayoutResult.txt");
	if (os.is_open())
	{
		for (int i = 0; i < _tmpParent.individualSet.size(); i++)
		{
			if (i == 0)
			{
				for (int j = 0; j < _tmpParent.individualSet[i].geneSize; j += 2)
				{
					os << _tmpParent.individualSet[i].genes[j] << "," << _tmpParent.individualSet[i].genes[j + 1] << endl;
				}
				//os << i.objectiveSet[0] << "," << i.objectiveSet[1] << endl;
			}
		}
	}
	os.close();
#pragma endregion
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