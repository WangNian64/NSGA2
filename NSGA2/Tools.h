#pragma once
#include <vector>
#include <algorithm>
using namespace std;
//�ж����������Ƿ��ص�
static bool IsRangeOverlap(double low1, double upper1, double low2, double upper2) {
	if (max(low1, low2) <= min(upper1, upper2)) {
		return true;
	}
	else
	{
		return false;
	}
}

// ��ȡ˫���������(����һ��0-1֮���С����
static double GetDoubleRand(int N = 99)
{
	double temp = rand() % (N + 1) / (double)(N + 1);
	return temp;
}

//�ַ����ָ��
static vector<string> split(const string& str, const string& pattern)
{
	vector<string> res;
	if ("" == str)
		return res;

	string strs = str + pattern;

	size_t pos = strs.find(pattern);
	size_t size = strs.size();

	while (pos != string::npos)
	{
		string x = strs.substr(0, pos);
		res.push_back(x);//stoi(x)ת����
		strs = strs.substr(pos + 1, size);
		pos = strs.find(pattern);
	}
	return res;
}