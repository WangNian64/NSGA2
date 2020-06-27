#include "GeneticAlgorithm.h"
#include <random>
GeneticAlgorithm::GeneticAlgorithm(GAPara gaPara)
{
	this->gaPara = gaPara;
	this->bestPathInfoList = vector<BestPathInfo>(gaPara.fitnessCount);
	for (int i = 0; i < gaPara.fitnessCount; ++i) {
		bestPathInfoList[i] = BestPathInfo(gaPara.geneSize);//初始化
	}
}

//初始化种群
Population GeneticAlgorithm::InitialPopulation()
{
	Population initPop(gaPara);
	initPop.initialize();
	initPop.evaluation(curIterNum, gaPara.maxIterNum, bestPathInfoList, gaPara.problemParas);//计算初始种群的适应度数组
	return initPop;
}

//根据父代得到交叉变异之后的子种群
Population GeneticAlgorithm::GetChildPopulation(Population parentPop)
{
	Population _tmpChild = parentPop.copy_all();//拷贝当前种群 作为子代的初始
	SBX(20, _tmpChild);//模拟二进制交叉
	PLM(20, gaPara.lower_bound_, gaPara.upper_bound_, _tmpChild);//多项式变异
	//计算子代的适应度
	_tmpChild.evaluation(curIterNum, gaPara.maxIterNum, bestPathInfoList, gaPara.problemParas);
	return _tmpChild;
}


//轮盘赌选择
void GeneticAlgorithm::wheelSelection(Population& p)
{

}
//交叉操作
void GeneticAlgorithm::crossover(Population& p)
{
	for (int i = 0; i < p.gaPara_pop.populationSize / 2; i++)
	{
		std::random_device rd;
		std::mt19937 rng(rd());
		std::uniform_real_distribution<> dist_real(0.0, 1.0);
		std::uniform_int_distribution<std::mt19937::result_type> dist_int(0, p.gaPara_pop.geneSize - 2);

		if (dist_real(rng) <= p.gaPara_pop.crossoverProb)
		{
			int pos[3];
			for (int j = 0; j < 3; j++)
			{
				pos[j] = dist_int(rng);
			}
			std::sort(pos, pos + 3);

			for (int j = pos[0] + 1; j < pos[1]; j++)
			{
				int aux = p.individualSet[i].genes[j];
				p.individualSet[i].genes[j] = p.individualSet[i + p.gaPara_pop.populationSize / 2].genes[j];
				p.individualSet[i + p.gaPara_pop.populationSize / 2].genes[j] = aux;
			}
			for (int j = pos[2] + 1; j < p.gaPara_pop.geneSize; j++)
			{
				int aux = p.individualSet[i].genes[j];
				p.individualSet[i].genes[j] = p.individualSet[i + p.gaPara_pop.populationSize / 2].genes[j];
				p.individualSet[i + p.gaPara_pop.populationSize / 2].genes[j] = aux;
			}
		}
	}
}

//变异操作
void GeneticAlgorithm::mutation(Population& p)
{
	for (int i = 0; i < p.gaPara_pop.populationSize; i++)
	{
		for (int j = 0; j < p.gaPara_pop.geneSize; j++)
		{
			std::random_device rd;
			std::mt19937 rng(rd());
			std::uniform_real_distribution<> dist(0, 1);

			if (dist(rng) < p.gaPara_pop.fitnessCount) {
				if (p.individualSet[i].genes[j] == 0) {
					p.individualSet[i].genes[j] = 1;
				}
				else {
					p.individualSet[i].genes[j] = 0;
				}
			}
		}
	}
}

//模拟二进制交叉
void GeneticAlgorithm::SBX(int index, Population& p)
{
	vector<Individual> _tmpSet;
	std::random_device rd;
	std::mt19937 rng(rd());
	std::uniform_real_distribution<> dist_real(0.0, 1.0);

	//不能改动p的gaPara
	vector<double> lowerBound_Ind1(p.gaPara_pop.geneSize, 0);
	vector<double> lowerBound_Ind2(p.gaPara_pop.geneSize, 0);
	vector<double> upperBound_Ind1(p.gaPara_pop.geneSize, 0);
	vector<double> upperBound_Ind2(p.gaPara_pop.geneSize, 0);
	for (int i = 0; i < p.gaPara_pop.geneSize; ++i) {
		lowerBound_Ind1[i] = p.gaPara_pop.lower_bound_[i];
		lowerBound_Ind2[i] = p.gaPara_pop.lower_bound_[i];
	}
	for (int i = 0; i < p.gaPara_pop.geneSize; ++i) {
		upperBound_Ind1[i] = p.gaPara_pop.upper_bound_[i];
		upperBound_Ind2[i] = p.gaPara_pop.upper_bound_[i];
	}

	for (int i = 0; i < p.gaPara_pop.populationSize / 2; i++)//一次处理两个个体
	{
		Individual& ind1 = p.individualSet[i];//待交叉的第一个
		Individual& ind2 = p.individualSet[i + (p.gaPara_pop.populationSize / 2)];//待交叉的第二个

		//需要产生新的个体吗？
		Individual c1(0, p.gaPara_pop.geneSize, p.gaPara_pop.fitnessCount, lowerBound_Ind1, upperBound_Ind1, p.gaPara_pop.problemParas);
		Individual c2(0, p.gaPara_pop.geneSize, p.gaPara_pop.fitnessCount, lowerBound_Ind2, upperBound_Ind2, p.gaPara_pop.problemParas);

		//cout << "ind1:" << endl;
		//for (int i = 0; i < p.gaPara_pop.geneSize; i++) {
		//	cout << ind1.genes[i] << ", ";
		//}
		//cout << endl;
		//cout << "ind2:" << endl;
		//for (int i = 0; i < p.gaPara_pop.geneSize; i++) {
		//	cout << ind2.genes[i] << ", ";
		//}
		//cout << endl;
		//cout << endl;

		//for (int j = 0; j < p.gaPara_pop.geneSize; j++)
		//{
		//	//随机数小于交叉概率，说明应该交叉
		//	if (dist_real(rng) <= p.gaPara_pop.crossoverProb)
		//	{
		//		double u = dist_real(rng);
		//		double x_bar = (p1.genes[j] + p2.genes[j]) / 2;
		//		double b_bar;
		//		if (u <= 0.5)
		//		{
		//			b_bar = pow(2.0 * u, 1.0 / ((double)index + 1.0));
		//		}
		//		else
		//		{
		//			b_bar = pow(1.0 / 2.0 * (1.0 - u), 1.0 / ((double)index + 1.0));
		//		}
		//		p1.genes[j] = x_bar - b_bar * abs(p2.genes[j] - p1.genes[j]) / 2;
		//		p2.genes[j] = x_bar + b_bar * abs(p2.genes[j] - p1.genes[j]) / 2;
		//	}
		//}

		#pragma region 方法1
		//先修改设备朝向
		for (int j = 2; j < p.gaPara_pop.geneSize; j += 3) {
			//随机数小于交叉概率，说明应该交叉
			if (dist_real(rng) <= p.gaPara_pop.crossoverProb)
			{
				double u = dist_real(rng);
				double x_bar = (ind1.genes[j] + ind2.genes[j]) / 2;
				double b_bar;
				if (u <= 0.5)
				{
					b_bar = pow(2.0 * u, 1.0 / ((double)index + 1.0));
				}
				else
				{
					b_bar = pow(1.0 / 2.0 * (1.0 - u), 1.0 / ((double)index + 1.0));
				}
				c1.genes[j] = x_bar - b_bar * abs(ind1.genes[j] - ind2.genes[j]) / 2;
				c2.genes[j] = x_bar + b_bar * abs(ind1.genes[j] - ind2.genes[j]) / 2;

				//处理朝向越界
				if (c1.genes[j] >= p.gaPara_pop.upper_bound_[j]) {//注意，=也不行
					double thre = GetDoubleRand(99);
					if (ind1.genes[j] >= p.gaPara_pop.upper_bound_[j] - 1)
					{
						c1.genes[j] = GetDoubleRand() *
							(p.gaPara_pop.upper_bound_[j] - p.gaPara_pop.lower_bound_[j]) + p.gaPara_pop.lower_bound_[j];
					}
					else if (thre < 0.5)
					{
						c1.genes[j] = p.gaPara_pop.upper_bound_[j]
							- (p.gaPara_pop.upper_bound_[j] - ind1.genes[j]) * GetDoubleRand();
					}
					else
					{
						c1.genes[j] = p.gaPara_pop.upper_bound_[j] - 0.5;
					}
				}
				if (c1.genes[j] < p.gaPara_pop.lower_bound_[j]) {
					double thre = GetDoubleRand(99);
					if (ind1.genes[j] == p.gaPara_pop.lower_bound_[j])
					{
						c1.genes[j] = GetDoubleRand() *
							(p.gaPara_pop.upper_bound_[j] - p.gaPara_pop.lower_bound_[j]) + p.gaPara_pop.lower_bound_[j];
					}
					else if (thre < 0.5)
					{
						c1.genes[j] = p.gaPara_pop.lower_bound_[j]
							+ (ind1.genes[j] - p.gaPara_pop.lower_bound_[j]) * GetDoubleRand();
					}
					else
					{
						c1.genes[j] = p.gaPara_pop.lower_bound_[j];
					}
				}

				if (c2.genes[j] >= p.gaPara_pop.upper_bound_[j]) {
					double thre = GetDoubleRand(99);
					if (ind2.genes[j] >= p.gaPara_pop.upper_bound_[j] - 1) 
					{
						ind2.genes[j] = GetDoubleRand() *
							(p.gaPara_pop.upper_bound_[j] - p.gaPara_pop.lower_bound_[j]) + p.gaPara_pop.lower_bound_[j];
					}
					else if (thre < 0.5)
					{
						ind2.genes[j] = p.gaPara_pop.upper_bound_[j]
							- (p.gaPara_pop.upper_bound_[j] - ind1.genes[j]) * GetDoubleRand();
					}
					else
					{
						ind2.genes[j] = p.gaPara_pop.upper_bound_[j] - 0.5;
					}
				}
				if (c2.genes[j] < p.gaPara_pop.lower_bound_[j]) {
					double thre = GetDoubleRand(99);
					if (ind1.genes[j] == p.gaPara_pop.lower_bound_[j])
					{
						c2.genes[j] = GetDoubleRand() *
							(p.gaPara_pop.upper_bound_[j] - p.gaPara_pop.lower_bound_[j]) + p.gaPara_pop.lower_bound_[j];
					}
					else if (thre < 0.5)
					{
						c2.genes[j] = p.gaPara_pop.lower_bound_[j]
							+ (ind2.genes[j] - p.gaPara_pop.lower_bound_[j]) * GetDoubleRand();
					}
					else
					{
						c2.genes[j] = p.gaPara_pop.lower_bound_[j];
					}
				}
			}
		}

		//根据朝向修改设备上下界范围
		for (int j = 2; j < p.gaPara_pop.geneSize; j += 3) {
			DeviceDirect curDirect = (DeviceDirect)(int)ind1.genes[j];
			if (curDirect == DeviceDirect::Rotate90 || curDirect == DeviceDirect::Rotate270)
			{
				//x和y
				lowerBound_Ind1[j - 2] = 0 + p.gaPara_pop.problemParas.deviceParaList[j / 3].size.y * 0.5 + p.gaPara_pop.problemParas.deviceParaList[j / 3].spaceLength;
				lowerBound_Ind1[j - 1] = 0 + p.gaPara_pop.problemParas.deviceParaList[j / 3].size.x * 0.5 + p.gaPara_pop.problemParas.deviceParaList[j / 3].spaceLength;

				upperBound_Ind1[j - 2] = p.gaPara_pop.problemParas.workShopLength - p.gaPara_pop.problemParas.deviceParaList[j / 3].size.y * 0.5 - p.gaPara_pop.problemParas.deviceParaList[j / 3].spaceLength;
				upperBound_Ind1[j - 1] = p.gaPara_pop.problemParas.workShopWidth - p.gaPara_pop.problemParas.deviceParaList[j / 3].size.x * 0.5 - p.gaPara_pop.problemParas.deviceParaList[j / 3].spaceLength;
			}
			else
			{
				//x和y
				lowerBound_Ind1[j - 2] = 0 + p.gaPara_pop.problemParas.deviceParaList[j / 3].size.x * 0.5 + p.gaPara_pop.problemParas.deviceParaList[j / 3].spaceLength;
				lowerBound_Ind1[j - 1] = 0 + p.gaPara_pop.problemParas.deviceParaList[j / 3].size.y * 0.5 + p.gaPara_pop.problemParas.deviceParaList[j / 3].spaceLength;

				upperBound_Ind1[j - 2] = p.gaPara_pop.problemParas.workShopLength - p.gaPara_pop.problemParas.deviceParaList[j / 3].size.x * 0.5 - p.gaPara_pop.problemParas.deviceParaList[j / 3].spaceLength;
				upperBound_Ind1[j - 1] = p.gaPara_pop.problemParas.workShopWidth - p.gaPara_pop.problemParas.deviceParaList[j / 3].size.y * 0.5 - p.gaPara_pop.problemParas.deviceParaList[j / 3].spaceLength;
			}
		}
		for (int j = 2; j < p.gaPara_pop.geneSize; j += 3) {
			DeviceDirect curDirect = (DeviceDirect)(int)ind2.genes[j];
			if (curDirect == DeviceDirect::Rotate90 || curDirect == DeviceDirect::Rotate270)
			{
				//x和y
				lowerBound_Ind2[j - 2] = 0 + p.gaPara_pop.problemParas.deviceParaList[j / 3].size.y * 0.5 + p.gaPara_pop.problemParas.deviceParaList[j / 3].spaceLength;
				lowerBound_Ind2[j - 1] = 0 + p.gaPara_pop.problemParas.deviceParaList[j / 3].size.x * 0.5 + p.gaPara_pop.problemParas.deviceParaList[j / 3].spaceLength;

				upperBound_Ind2[j - 2] = p.gaPara_pop.problemParas.workShopLength - p.gaPara_pop.problemParas.deviceParaList[j / 3].size.y * 0.5 - p.gaPara_pop.problemParas.deviceParaList[j / 3].spaceLength;
				upperBound_Ind2[j - 1] = p.gaPara_pop.problemParas.workShopWidth - p.gaPara_pop.problemParas.deviceParaList[j / 3].size.x * 0.5 - p.gaPara_pop.problemParas.deviceParaList[j / 3].spaceLength;

			}
			else
			{
				//x和y
				lowerBound_Ind2[j - 2] = 0 + p.gaPara_pop.problemParas.deviceParaList[j / 3].size.x * 0.5 + p.gaPara_pop.problemParas.deviceParaList[j / 3].spaceLength;
				lowerBound_Ind2[j - 1] = 0 + p.gaPara_pop.problemParas.deviceParaList[j / 3].size.y * 0.5 + p.gaPara_pop.problemParas.deviceParaList[j / 3].spaceLength;

				upperBound_Ind2[j - 2] = p.gaPara_pop.problemParas.workShopLength - p.gaPara_pop.problemParas.deviceParaList[j / 3].size.x * 0.5 - p.gaPara_pop.problemParas.deviceParaList[j / 3].spaceLength;
				upperBound_Ind2[j - 1] = p.gaPara_pop.problemParas.workShopWidth - p.gaPara_pop.problemParas.deviceParaList[j / 3].size.y * 0.5 - p.gaPara_pop.problemParas.deviceParaList[j / 3].spaceLength;

			}
		}

		//处理坐标值
		for (int j = 0; j < p.gaPara_pop.geneSize; j++)
		{
			if (j % 3 != 2)
			{
				if (dist_real(rng) <= p.gaPara_pop.crossoverProb)
				{
					double u = dist_real(rng);
					double x_bar = (ind1.genes[j] + ind2.genes[j]) / 2;
					double b_bar;
					if (u <= 0.5)
					{
						b_bar = pow(2.0 * u, 1.0 / ((double)index + 1.0));
					}
					else
					{
						b_bar = pow(1.0 / 2.0 * (1.0 - u), 1.0 / ((double)index + 1.0));
					}
					c1.genes[j] = x_bar - b_bar * abs(ind1.genes[j] - ind2.genes[j]) / 2;
					c2.genes[j] = x_bar + b_bar * abs(ind1.genes[j] - ind2.genes[j]) / 2;
				}

				//防止越界
				if (c1.genes[j] > upperBound_Ind1[j]) {
					double thre = GetDoubleRand(99);
					if (ind1.genes[j] >= upperBound_Ind1[j])
					{
						c1.genes[j] = GetDoubleRand() *
							(upperBound_Ind1[j] - lowerBound_Ind1[j]) + lowerBound_Ind1[j];
					}
					else if (thre < 0.5)
					{
						c1.genes[j] = upperBound_Ind1[j]
							- abs(upperBound_Ind1[j] - ind1.genes[j]) * GetDoubleRand();
					}
					else
					{
						c1.genes[j] = upperBound_Ind1[j];
					}
				}
				if (c1.genes[j] < lowerBound_Ind1[j]) {
					double thre = GetDoubleRand(99);
					if (ind1.genes[j] <= lowerBound_Ind1[j])
					{
						c1.genes[j] = GetDoubleRand() *
							(upperBound_Ind1[j] - lowerBound_Ind1[j]) + lowerBound_Ind1[j];
					}
					else if (thre < 0.5)
					{
						c1.genes[j] = lowerBound_Ind1[j]
							+ abs(ind1.genes[j] - lowerBound_Ind1[j]) * GetDoubleRand();
					}
					else
					{
						c1.genes[j] = lowerBound_Ind1[j];
					}
				}


				if (c2.genes[j] > upperBound_Ind2[j]) {
					double thre = GetDoubleRand(99);
					if (ind2.genes[j] >= upperBound_Ind2[j])
					{
						c2.genes[j] = GetDoubleRand() *
							(upperBound_Ind2[j] - lowerBound_Ind2[j]) + lowerBound_Ind2[j];
					}
					else if (thre < 0.5)
					{
						c2.genes[j] = upperBound_Ind2[j]
							- abs(upperBound_Ind2[j] - ind2.genes[j]) * GetDoubleRand();
					}
					else
					{
						c2.genes[j] = upperBound_Ind2[j];
					}
				}
				if (c2.genes[j] < lowerBound_Ind2[j]) {
					double thre = GetDoubleRand(99);
					if (ind2.genes[j] <= lowerBound_Ind2[j])
					{
						c2.genes[j] = GetDoubleRand() *
							(upperBound_Ind2[j] - lowerBound_Ind2[j]) + lowerBound_Ind2[j];
					}
					else if (thre < 0.5)
					{
						c2.genes[j] = lowerBound_Ind2[j]
							+ abs(ind2.genes[j] - lowerBound_Ind2[j]) * GetDoubleRand();
					}
					else
					{
						c2.genes[j] = lowerBound_Ind2[j];
					}
				}
			}
		}

		#pragma endregion


		//#pragma region 方法2
		////先交叉
		//for (int j = 0; j < p.gaPara_pop.geneSize; j++)
		//{
		//	//随机数小于交叉概率，说明应该交叉
		//	if (dist_real(rng) <= p.gaPara_pop.crossoverProb)
		//	{
		//		double u = dist_real(rng);
		//		double x_bar = (ind1.genes[j] + ind2.genes[j]) / 2;
		//		double b_bar;
		//		if (u <= 0.5)
		//		{
		//			b_bar = pow(2.0 * u, 1.0 / ((double)index + 1.0));
		//		}
		//		else
		//		{
		//			b_bar = pow(1.0 / 2.0 * (1.0 - u), 1.0 / ((double)index + 1.0));
		//		}
		//		c1.genes[j] = x_bar - b_bar * abs(ind1.genes[j] - ind2.genes[j]) / 2;
		//		c2.genes[j] = x_bar + b_bar * abs(ind1.genes[j] - ind2.genes[j]) / 2;
		//	}
		//}

		//#pragma region 处理P1
		////先处理朝向越界
		//for (int j = 2; j < p.gaPara_pop.geneSize; j += 3) {
		//	if (c1.genes[j] > p.gaPara_pop.upper_bound_[j]) {
		//		double thre = GetDoubleRand(99);
		//		if (ind1.genes[j] == p.gaPara_pop.upper_bound_[j])
		//		{
		//			c1.genes[j] = GetDoubleRand() *
		//				(p.gaPara_pop.upper_bound_[j] - p.gaPara_pop.lower_bound_[j]) + p.gaPara_pop.lower_bound_[j];
		//		}
		//		else if (thre < 0.5)
		//		{
		//			c1.genes[j] = p.gaPara_pop.upper_bound_[j]
		//				- (p.gaPara_pop.upper_bound_[j] - ind1.genes[j]) * GetDoubleRand();
		//		}
		//		else
		//		{
		//			c1.genes[j] = p.gaPara_pop.upper_bound_[j];
		//		}
		//	}
		//	if (c1.genes[j] < p.gaPara_pop.lower_bound_[j]) {
		//		double thre = GetDoubleRand(99);
		//		if (ind1.genes[j] == p.gaPara_pop.lower_bound_[j])
		//		{
		//			c1.genes[j] = GetDoubleRand() *
		//				(p.gaPara_pop.upper_bound_[j] - p.gaPara_pop.lower_bound_[j]) + p.gaPara_pop.lower_bound_[j];
		//		}
		//		else if (thre < 0.5)
		//		{
		//			c1.genes[j] = p.gaPara_pop.lower_bound_[j]
		//				+ (ind1.genes[j] - p.gaPara_pop.lower_bound_[j]) * GetDoubleRand();
		//		}
		//		else
		//		{
		//			c1.genes[j] = p.gaPara_pop.lower_bound_[j];
		//		}
		//	}
		//}
		////根据朝向修改设备上下界范围
		//for (int j = 2; j < p.gaPara_pop.geneSize; j += 3) {
		//	DeviceDirect curDirect = (DeviceDirect)(int)c1.genes[j];
		//	if (curDirect == DeviceDirect::Rotate90 || curDirect == DeviceDirect::Rotate270)
		//	{
		//		//x和y
		//		p.gaPara_pop.lower_bound_[j - 2] = 0 + p.gaPara_pop.problemParas.deviceParaList[j / 3].size.y * 0.5 + p.gaPara_pop.problemParas.deviceParaList[j / 3].spaceLength;
		//		p.gaPara_pop.lower_bound_[j - 1] = 0 + p.gaPara_pop.problemParas.deviceParaList[j / 3].size.x * 0.5 + p.gaPara_pop.problemParas.deviceParaList[j / 3].spaceLength;

		//		p.gaPara_pop.upper_bound_[j - 2] = p.gaPara_pop.problemParas.workShopLength - p.gaPara_pop.problemParas.deviceParaList[j / 3].size.y * 0.5 - p.gaPara_pop.problemParas.deviceParaList[j / 3].spaceLength;
		//		p.gaPara_pop.upper_bound_[j - 1] = p.gaPara_pop.problemParas.workShopWidth - p.gaPara_pop.problemParas.deviceParaList[j / 3].size.x * 0.5 - p.gaPara_pop.problemParas.deviceParaList[j / 3].spaceLength;

		//	}
		//	else
		//	{
		//		//x和y
		//		p.gaPara_pop.lower_bound_[j - 2] = 0 + p.gaPara_pop.problemParas.deviceParaList[j / 3].size.x * 0.5 + p.gaPara_pop.problemParas.deviceParaList[j / 3].spaceLength;
		//		p.gaPara_pop.lower_bound_[j - 1] = 0 + p.gaPara_pop.problemParas.deviceParaList[j / 3].size.y * 0.5 + p.gaPara_pop.problemParas.deviceParaList[j / 3].spaceLength;

		//		p.gaPara_pop.upper_bound_[j - 2] = p.gaPara_pop.problemParas.workShopLength - p.gaPara_pop.problemParas.deviceParaList[j / 3].size.x * 0.5 - p.gaPara_pop.problemParas.deviceParaList[j / 3].spaceLength;
		//		p.gaPara_pop.upper_bound_[j - 1] = p.gaPara_pop.problemParas.workShopWidth - p.gaPara_pop.problemParas.deviceParaList[j / 3].size.y * 0.5 - p.gaPara_pop.problemParas.deviceParaList[j / 3].spaceLength;

		//	}
		//}
		////处理坐标值
		//for (int j = 0; j < p.gaPara_pop.geneSize; j++)
		//{
		//	if (j % 3 != 2)
		//	{
		//		if (c1.genes[j] > p.gaPara_pop.upper_bound_[j]) {
		//			double thre = GetDoubleRand(99);
		//			if (ind1.genes[j] >= p.gaPara_pop.upper_bound_[j])
		//			{
		//				c1.genes[j] = GetDoubleRand() *
		//					(p.gaPara_pop.upper_bound_[j] - p.gaPara_pop.lower_bound_[j]) + p.gaPara_pop.lower_bound_[j];
		//			}
		//			else if (thre < 0.5)
		//			{
		//				c1.genes[j] = p.gaPara_pop.upper_bound_[j]
		//					- abs(p.gaPara_pop.upper_bound_[j] - ind1.genes[j]) * GetDoubleRand();
		//			}
		//			else
		//			{
		//				c1.genes[j] = p.gaPara_pop.upper_bound_[j];
		//			}
		//		}
		//		if (c1.genes[j] < p.gaPara_pop.lower_bound_[j]) {
		//			double thre = GetDoubleRand(99);
		//			if (ind1.genes[j] <= p.gaPara_pop.lower_bound_[j])
		//			{
		//				c1.genes[j] = GetDoubleRand() *
		//					(p.gaPara_pop.upper_bound_[j] - p.gaPara_pop.lower_bound_[j]) + p.gaPara_pop.lower_bound_[j];
		//			}
		//			else if (thre < 0.5)
		//			{
		//				c1.genes[j] = p.gaPara_pop.lower_bound_[j]
		//					+ abs(ind1.genes[j] - p.gaPara_pop.lower_bound_[j]) * GetDoubleRand();
		//			}
		//			else
		//			{
		//				c1.genes[j] = p.gaPara_pop.lower_bound_[j];
		//			}
		//		}
		//	}
		//}
		//#pragma endregion

		//#pragma region 处理P2
		////先处理朝向越界
		//for (int j = 2; j < p.gaPara_pop.geneSize; j += 3) {
		//	if (c2.genes[j] > p.gaPara_pop.upper_bound_[j]) {
		//		double thre = GetDoubleRand(99);
		//		if (ind2.genes[j] == p.gaPara_pop.upper_bound_[j])
		//		{
		//			c2.genes[j] = GetDoubleRand() *
		//				(p.gaPara_pop.upper_bound_[j] - p.gaPara_pop.lower_bound_[j]) + p.gaPara_pop.lower_bound_[j];
		//		}
		//		else if (thre < 0.5)
		//		{
		//			c2.genes[j] = p.gaPara_pop.upper_bound_[j]
		//				- (p.gaPara_pop.upper_bound_[j] - ind1.genes[j]) * GetDoubleRand();
		//		}
		//		else
		//		{
		//			c2.genes[j] = p.gaPara_pop.upper_bound_[j];
		//		}
		//	}
		//	if (c2.genes[j] < p.gaPara_pop.lower_bound_[j]) {
		//		double thre = GetDoubleRand(99);
		//		if (ind1.genes[j] == p.gaPara_pop.lower_bound_[j])
		//		{
		//			c2.genes[j] = GetDoubleRand() *
		//				(p.gaPara_pop.upper_bound_[j] - p.gaPara_pop.lower_bound_[j]) + p.gaPara_pop.lower_bound_[j];
		//		}
		//		else if (thre < 0.5)
		//		{
		//			c2.genes[j] = p.gaPara_pop.lower_bound_[j]
		//				+ (ind2.genes[j] - p.gaPara_pop.lower_bound_[j]) * GetDoubleRand();
		//		}
		//		else
		//		{
		//			c2.genes[j] = p.gaPara_pop.lower_bound_[j];
		//		}
		//	}
		//}
		////根据朝向修改设备上下界范围
		//for (int j = 2; j < p.gaPara_pop.geneSize; j += 3) {
		//	DeviceDirect curDirect = (DeviceDirect)(int)ind2.genes[j];
		//	if (curDirect == DeviceDirect::Rotate90 || curDirect == DeviceDirect::Rotate270)
		//	{
		//		//x和y
		//		p.gaPara_pop.lower_bound_[j - 2] = 0 + p.gaPara_pop.problemParas.deviceParaList[j / 3].size.y * 0.5 + p.gaPara_pop.problemParas.deviceParaList[j / 3].spaceLength;
		//		p.gaPara_pop.lower_bound_[j - 1] = 0 + p.gaPara_pop.problemParas.deviceParaList[j / 3].size.x * 0.5 + p.gaPara_pop.problemParas.deviceParaList[j / 3].spaceLength;

		//		p.gaPara_pop.upper_bound_[j - 2] = p.gaPara_pop.problemParas.workShopLength - p.gaPara_pop.problemParas.deviceParaList[j / 3].size.y * 0.5 - p.gaPara_pop.problemParas.deviceParaList[j / 3].spaceLength;
		//		p.gaPara_pop.upper_bound_[j - 1] = p.gaPara_pop.problemParas.workShopWidth - p.gaPara_pop.problemParas.deviceParaList[j / 3].size.x * 0.5 - p.gaPara_pop.problemParas.deviceParaList[j / 3].spaceLength;

		//	}
		//	else
		//	{
		//		//x和y
		//		p.gaPara_pop.lower_bound_[j - 2] = 0 + p.gaPara_pop.problemParas.deviceParaList[j / 3].size.x * 0.5 + p.gaPara_pop.problemParas.deviceParaList[j / 3].spaceLength;
		//		p.gaPara_pop.lower_bound_[j - 1] = 0 + p.gaPara_pop.problemParas.deviceParaList[j / 3].size.y * 0.5 + p.gaPara_pop.problemParas.deviceParaList[j / 3].spaceLength;

		//		p.gaPara_pop.upper_bound_[j - 2] = p.gaPara_pop.problemParas.workShopLength - p.gaPara_pop.problemParas.deviceParaList[j / 3].size.x * 0.5 - p.gaPara_pop.problemParas.deviceParaList[j / 3].spaceLength;
		//		p.gaPara_pop.upper_bound_[j - 1] = p.gaPara_pop.problemParas.workShopWidth - p.gaPara_pop.problemParas.deviceParaList[j / 3].size.y * 0.5 - p.gaPara_pop.problemParas.deviceParaList[j / 3].spaceLength;

		//	}
		//}
		////处理坐标值
		//for (int j = 0; j < p.gaPara_pop.geneSize; j++)
		//{
		//	if (j % 3 != 2)
		//	{
		//		if (c2.genes[j] > p.gaPara_pop.upper_bound_[j]) {
		//			double thre = GetDoubleRand(99);
		//			if (ind1.genes[j] >= p.gaPara_pop.upper_bound_[j])
		//			{
		//				c2.genes[j] = GetDoubleRand() *
		//					(p.gaPara_pop.upper_bound_[j] - p.gaPara_pop.lower_bound_[j]) + p.gaPara_pop.lower_bound_[j];
		//			}
		//			else if (thre < 0.5)
		//			{
		//				c2.genes[j] = p.gaPara_pop.upper_bound_[j]
		//					- abs(p.gaPara_pop.upper_bound_[j] - ind1.genes[j]) * GetDoubleRand();
		//			}
		//			else
		//			{
		//				c2.genes[j] = p.gaPara_pop.upper_bound_[j];
		//			}
		//		}
		//		if (ind2.genes[j] < p.gaPara_pop.lower_bound_[j]) {
		//			double thre = GetDoubleRand(99);
		//			if (ind2.genes[j] <= p.gaPara_pop.lower_bound_[j])
		//			{
		//				c2.genes[j] = GetDoubleRand() *
		//					(p.gaPara_pop.upper_bound_[j] - p.gaPara_pop.lower_bound_[j]) + p.gaPara_pop.lower_bound_[j];
		//			}
		//			else if (thre < 0.5)
		//			{
		//				c2.genes[j] = p.gaPara_pop.lower_bound_[j]
		//					+ abs(ind2.genes[j] - p.gaPara_pop.lower_bound_[j]) * GetDoubleRand();
		//			}
		//			else
		//			{
		//				c2.genes[j] = p.gaPara_pop.lower_bound_[j];
		//			}
		//		}
		//	}
		//}
		//#pragma endregion
		//#pragma endregion







		//for (int j = 0; j < p.gaPara_pop.geneSize; j++)
		//{
		//	//随机数小于交叉概率，说明应该交叉
		//	if (dist_real(rng) <= p.gaPara_pop.crossoverProb)
		//	{
		//		double u = dist_real(rng);
		//		double x_bar = (p1.genes[j] + p2.genes[j]) / 2;
		//		double b_bar;
		//		if (u <= 0.5)
		//		{
		//			b_bar = pow(2.0 * u, 1.0 / ((double)index + 1.0));
		//		}
		//		else
		//		{
		//			b_bar = pow(1.0 / 2.0 * (1.0 - u), 1.0 / ((double)index + 1.0));
		//		}
		//		p1.genes[j] = x_bar - b_bar * abs(p2.genes[j] - p1.genes[j]) / 2;
		//		p2.genes[j] = x_bar + b_bar * abs(p2.genes[j] - p1.genes[j]) / 2;
		//	}
		//}



		//cout << "ind1:" << endl;
		//for (int i = 0; i < p.gaPara_pop.geneSize; i++) {
		//	cout << ind1.genes[i] << ", ";
		//}
		//cout << endl;
		//cout << "ind2:" << endl;
		//for (int i = 0; i < p.gaPara_pop.geneSize; i++) {
		//	cout << ind2.genes[i] << ", ";
		//}
		//cout << endl;
		//cout << endl;


		_tmpSet.push_back(c1);
		_tmpSet.push_back(c2);
	}
	p.individualSet = _tmpSet;//更新种群
}

//变异操作
void GeneticAlgorithm::PLM(int index, double* lowBounds, double* upBounds, Population& p)
{
	std::random_device rd;
	std::mt19937 rng(rd());
	std::uniform_real_distribution<> dist_real(0, 1);

	double last_gene = 0.0;
	for (int i = 0; i < p.individualSet.size(); i++)
	{
		Individual& ind = p.individualSet[i];

		//先考虑朝向
		for (int j = 2; j < p.gaPara_pop.geneSize; j += 3) {
			if (dist_real(rng) <= p.gaPara_pop.mutationProb) {
				last_gene = ind.genes[j];
				double u = dist_real(rng);
				if (u <= 0.5)
				{
					double delta_l = pow(2.0 * u, 1.0 / (1.0 + (double)index)) - 1.0;
					ind.genes[j] = ind.genes[j] + delta_l * (ind.genes[j] - lowBounds[j]);
				}
				else
				{
					double delta_r = 1.0 - pow(2.0 * (1.0 - u), 1.0 / (1.0 + (double)index));
					ind.genes[j] = ind.genes[j] + delta_r * (upBounds[j] - ind.genes[j]);
				}

				//在这之后直接修改上下界范围
				DeviceDirect curDirect = (DeviceDirect)(int)ind.genes[j];
				if (curDirect == DeviceDirect::Rotate90 || curDirect == DeviceDirect::Rotate270)
				{
					//x和y
					p.gaPara_pop.lower_bound_[j - 2] = 0 + p.gaPara_pop.problemParas.deviceParaList[j / 3].size.y * 0.5 + p.gaPara_pop.problemParas.deviceParaList[j / 3].spaceLength;
					p.gaPara_pop.lower_bound_[j - 1] = 0 + p.gaPara_pop.problemParas.deviceParaList[j / 3].size.x * 0.5 + p.gaPara_pop.problemParas.deviceParaList[j / 3].spaceLength;

					p.gaPara_pop.upper_bound_[j - 2] = p.gaPara_pop.problemParas.workShopLength - p.gaPara_pop.problemParas.deviceParaList[j / 3].size.y * 0.5 - p.gaPara_pop.problemParas.deviceParaList[j / 3].spaceLength;
					p.gaPara_pop.upper_bound_[j - 1] = p.gaPara_pop.problemParas.workShopWidth - p.gaPara_pop.problemParas.deviceParaList[j / 3].size.x * 0.5 - p.gaPara_pop.problemParas.deviceParaList[j / 3].spaceLength;

				}
				else
				{
					//x和y
					p.gaPara_pop.lower_bound_[j - 2] = 0 + p.gaPara_pop.problemParas.deviceParaList[j / 3].size.x * 0.5 + p.gaPara_pop.problemParas.deviceParaList[j / 3].spaceLength;
					p.gaPara_pop.lower_bound_[j - 1] = 0 + p.gaPara_pop.problemParas.deviceParaList[j / 3].size.y * 0.5 + p.gaPara_pop.problemParas.deviceParaList[j / 3].spaceLength;

					p.gaPara_pop.upper_bound_[j - 2] = p.gaPara_pop.problemParas.workShopLength - p.gaPara_pop.problemParas.deviceParaList[j / 3].size.x * 0.5 - p.gaPara_pop.problemParas.deviceParaList[j / 3].spaceLength;
					p.gaPara_pop.upper_bound_[j - 1] = p.gaPara_pop.problemParas.workShopWidth - p.gaPara_pop.problemParas.deviceParaList[j / 3].size.y * 0.5 - p.gaPara_pop.problemParas.deviceParaList[j / 3].spaceLength;

				}
			}
		}
		//处理坐标值
		for (int j = 0; j < p.gaPara_pop.geneSize; j++)
		{
			if (j % 3 != 2)
			{
				//先变异
				if (dist_real(rng) <= p.gaPara_pop.mutationProb) {
					double u = dist_real(rng);
					if (u <= 0.5)
					{
						double delta_l = pow(2.0 * u, 1.0 / (1.0 + (double)index)) - 1.0;
						ind.genes[j] = ind.genes[j] + delta_l * (ind.genes[j] - lowBounds[j]);
					}
					else
					{
						double delta_r = 1.0 - pow(2.0 * (1.0 - u), 1.0 / (1.0 + (double)index));
						ind.genes[j] = ind.genes[j] + delta_r * (upBounds[j] - ind.genes[j]);
					}
				}


				last_gene = ind.genes[j];
				if (ind.genes[j] > p.gaPara_pop.upper_bound_[j]) {
					double thre = GetDoubleRand(99);
					if (last_gene >= p.gaPara_pop.upper_bound_[j])
					{
						ind.genes[j] = GetDoubleRand() *
							(p.gaPara_pop.upper_bound_[j] - p.gaPara_pop.lower_bound_[j]) + p.gaPara_pop.lower_bound_[j];
					}
					else if (thre < 0.5)
					{
						ind.genes[j] = p.gaPara_pop.upper_bound_[j]
							- abs(p.gaPara_pop.upper_bound_[j] - last_gene) * GetDoubleRand();
					}
					else
					{
						ind.genes[j] = p.gaPara_pop.upper_bound_[j];
					}
				}
				if (ind.genes[j] < p.gaPara_pop.lower_bound_[j]) {
					double thre = GetDoubleRand(99);
					if (last_gene <= p.gaPara_pop.lower_bound_[j])
					{
						ind.genes[j] = GetDoubleRand() *
							(p.gaPara_pop.upper_bound_[j] - p.gaPara_pop.lower_bound_[j]) + p.gaPara_pop.lower_bound_[j];
					}
					else if (thre < 0.5)
					{
						ind.genes[j] = p.gaPara_pop.lower_bound_[j]
							+ abs(last_gene - p.gaPara_pop.lower_bound_[j]) * GetDoubleRand();
					}
					else
					{
						ind.genes[j] = p.gaPara_pop.lower_bound_[j];
					}
				}
			}
		}

		//for (int j = 0; j < p.gaPara_pop.geneSize; j++)
		//{
		//	if (dist_real(rng) <= p.gaPara_pop.mutationProb)
		//	{
		//		double u = dist_real(rng);
		//		if (u <= 0.5)
		//		{
		//			double delta_l = pow(2.0 * u, 1.0 / (1.0 + (double)index)) - 1.0;
		//			ind.genes[j] = ind.genes[j] + delta_l * (ind.genes[j] - lowBounds[j]);
		//		}
		//		else
		//		{
		//			double delta_r = 1.0 - pow(2.0 * (1.0 - u), 1.0 / (1.0 + (double)index));
		//			ind.genes[j] = ind.genes[j] + delta_r * (upBounds[j] - ind.genes[j]);
		//		}
		//	}
		//}
	}
}

//快速非支配排序算法，得到r的所有级别的pareto解集
vector<Front> GeneticAlgorithm::fastNonDominatedSort(Population* r)
{
	vector<Front> front;
	front.resize(1);

	//找到Rank0的Pareto解集
	for (int i = 0; i < r->gaPara_pop.populationSize; i++)
	{
		//每一个个体都有一个被支配解集
		Individual* p = &r->individualSet[i];
		p->dominatedSet.clear();
		p->dominatedCount = 0;

		for (int j = 0; j < r->gaPara_pop.populationSize; j++)
		{
			if (i == j) { continue; }
			Individual* q = &r->individualSet[j];
			if (p->dominate(*q))//p支配q
			{
				p->dominatedSet.push_back(q);
			}
			else if (q->dominate(*p))//q支配p，p的被支配数目+1
			{
				p->dominatedCount++;
			}
		}
		//没有能支配p的，说明p是最优的那个集合
		if (p->dominatedCount == 0)
		{
			p->rank = 0;
			front[0].push_back(p);
		}
	}

	int i = 0;
	while (front[i].size() != 0)
	{
		vector<Individual*> _tmpSet;
		_tmpSet.resize(0);
		for (int j = 0; j < front[i].size(); j++)
		{
			Individual* p = front[i][j];
			for (int k = 0; k < p->dominatedSet.size(); k++)
			{
				Individual* q = p->dominatedSet[k];
				q->dominatedCount--;
				if (q->dominatedCount == 0)
				{
					q->rank = i + 1;
					_tmpSet.push_back(q);
				}
			}
		}
		i++;
		front.push_back(_tmpSet);
	}

	return front;
}

bool comparator_0(Individual* a, Individual* b)
{
	return a->objectiveSet[0] < b->objectiveSet[0];
}

bool comparator_1(Individual* a, Individual* b)
{
	return a->objectiveSet[1] < b->objectiveSet[1];
}

//计算pareto前沿的拥挤距离
void GeneticAlgorithm::crowdingDistanceAssignment(Front f)
{
	int objectiveSize = 2;
	//Individual f_max(gaPara.geneSize, gaPara.fitnessCount, gaPara.lower_bound_, gaPara.upper_bound_, gaPara.problemParas);
	//Individual f_min(gaPara.geneSize, gaPara.fitnessCount, gaPara.lower_bound_, gaPara.upper_bound_, gaPara.problemParas);
	//f_max.evaluation(curIterNum, gaPara.maxIterNum, bestPathInfoList, gaPara.problemParas);
	//f_min.evaluation(curIterNum, gaPara.maxIterNum, bestPathInfoList, gaPara.problemParas);

	vector<double> f_max(objectiveSize, 0);
	vector<double> f_min(objectiveSize, MAX_FITNESS);
	for (int i = 0; i < f.size(); ++i) {
		for (int j = 0; j < objectiveSize; ++j) {
			f_max[j] = max(f_max[j], f[i]->objectiveSet[j]);
			f_min[j] = min(f_min[j], f[i]->objectiveSet[j]);
		}
	}
	int l = f.size();
	for (int i = 0; i < l; i++)
	{
		f[i]->distance = 0;
	}
	for (int m = 0; m < objectiveSize; m++)
	{
		if (m == 0)
		{
			sort(f.begin(), f.end(), comparator_0);
		}
		else if (m == 1)
		{
			sort(f.begin(), f.end(), comparator_1);
		}
		f[0]->distance = numeric_limits<double>::infinity();
		f[l - 1]->distance = numeric_limits<double>::infinity();
		for (int j = 1; j < l - 1; j++)
		{
			if (f[j + 1]->objectiveSet[m] == f[j - 1]->objectiveSet[m]) {
				f[j]->distance += 0;
			}
			else {
				f[j]->distance += (f[j + 1]->objectiveSet[m] - f[j - 1]->objectiveSet[m]) / (f_max[m] - f_min[m]);
				//f[j]->distance += (f[j + 1]->objectiveSet[m] - f[j - 1]->objectiveSet[m]) / (f_max.objectiveSet[m] - f_min.objectiveSet[m]);
			}
		}
	}
}