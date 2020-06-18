#pragma once
#include <vector>
#include <algorithm>
using namespace std;
//判断两个区间是否重叠
static bool IsRangeOverlap(double low1, double upper1, double low2, double upper2) {
	if (max(low1, low2) <= min(upper1, upper2)) {
		return true;
	}
	else
	{
		return false;
	}
}

// 获取双精度随机数(返回一个0-1之间的小数）
static double GetDoubleRand(int N = 99)
{
	double temp = rand() % (N + 1) / (double)(N + 1);
	return temp;
}

//字符串分割函数
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
		res.push_back(x);//stoi(x)转整型
		strs = strs.substr(pos + 1, size);
		pos = strs.find(pattern);
	}
	return res;
}