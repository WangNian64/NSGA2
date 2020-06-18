#include "GeneticAlgorithm.h"
#include <random>
GeneticAlgorithm::GeneticAlgorithm(GAPara gaPara)
{
	this->gaPara = gaPara;
}

//初始化种群
Population GeneticAlgorithm::InitialPopulation()
{
	Population initPop(gaPara);
	initPop.initialize();
	initPop.evaluation(gaPara.fitnessCount, gaPara.problemParas);//计算初始种群的适应度数组
	return initPop;
}

//根据父代得到交叉变异之后的子种群
Population GeneticAlgorithm::GetChildPopulation(Population parentPop)
{
	Population _tmpChild = parentPop.copy_all();//拷贝当前种群 作为子代的初始
	//交叉、变异操作，得到新的子代
	SBX(20, _tmpChild);//模拟二进制交叉
	PLM(20, gaPara.lowerBounds, gaPara.upBounds, _tmpChild);//多项式变异
	//计算子代的适应度
	_tmpChild.evaluation(gaPara.fitnessCount, gaPara.problemParas);
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

	for (int i = 0; i < p.gaPara_pop.populationSize / 2; i++)
	{
		Individual p1 = p.individualSet[i];//待交叉的第一个
		Individual p2 = p.individualSet[i + (p.gaPara_pop.populationSize / 2)];//待交叉的第二个

		Individual c1(p.gaPara_pop.geneSize, p.gaPara_pop.fitnessCount, gaPara.lowerBounds, gaPara.upBounds);
		Individual c2(p.gaPara_pop.geneSize, p.gaPara_pop.fitnessCount, gaPara.lowerBounds, gaPara.upBounds);

		for (int j = 0; j < p.gaPara_pop.geneSize; j++)
		{
			if (dist_real(rng) <= p.gaPara_pop.crossoverProb)//随机数小于交叉概率，说明应该交叉
			{
				double u = dist_real(rng);
				double x_bar = (p1.genes[j] + p2.genes[j]) / 2;
				double b_bar;
				if (u <= 0.5)
				{
					b_bar = pow(2.0 * u, 1.0 / ((double)index + 1.0));
				}
				else
				{
					b_bar = pow(1.0 / 2.0 * (1.0 - u), 1.0 / ((double)index + 1.0));
				}
				c1.genes[j] = x_bar - b_bar * abs(p2.genes[j] - p1.genes[j]) / 2;
				c2.genes[j] = x_bar + b_bar * abs(p2.genes[j] - p1.genes[j]) / 2;
			}
		}
		//交叉后放回种群
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

	for (int i = 0; i < p.individualSet.size(); i++)
	{
		Individual& ind = p.individualSet[i];
		for (int j = 0; j < p.gaPara_pop.geneSize; j++)
		{
			if (dist_real(rng) <= p.gaPara_pop.mutationProb)
			{
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
		}
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
	return a->objectiveSet[0] > b->objectiveSet[0];
}

bool comparator_1(Individual* a, Individual* b)
{
	return a->objectiveSet[1] > b->objectiveSet[1];
}

//计算pareto前沿的拥挤距离
void GeneticAlgorithm::crowdingDistanceAssignment(Front f)
{
	int objectiveSize = 2;
	Individual f_max(gaPara.geneSize, gaPara.fitnessCount, gaPara.lowerBounds, gaPara.upBounds);
	Individual f_min(gaPara.geneSize, gaPara.fitnessCount, gaPara.lowerBounds, gaPara.upBounds);
	f_max.genes[0] = gaPara.upBounds[0];
	f_min.genes[0] = gaPara.lowerBounds[0];
	f_max.evaluation(gaPara.fitnessCount, gaPara.problemParas);
	f_min.evaluation(gaPara.fitnessCount, gaPara.problemParas);

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
		for (int j = 2; j < l - 1; j++)
		{
			f[j]->distance = f[j]->distance + (f[j + 1]->objectiveSet[m] - f[j - 1]->objectiveSet[m]) / (f_max.objectiveSet[m] - f_min.objectiveSet[m]);
		}
	}
}