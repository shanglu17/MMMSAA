#include <iostream>
#include <vector>
#include <set>
#include <tuple>
#include <algorithm>
using namespace std;

//类型别名
using Task=string; //任务类型
using TaskSet=vector<Task>; //任务集合类型
using Partition=vector<TaskSet>; //划分类型
using PartitionResult=tuple<double, int, Partition>; //谱系图类型

//函数声明
double findMax(const vector<vector<double>>& C); //查找最大值
double getMaxCommCost(const TaskSet& set1, const TaskSet& set2, const vector<Task>&D,const vector<vector<double>>&C); //计算最大通信代价
vector<PartitionResult>MMMSAA(const vector<Task>&D,const vector<vector<double>>&C,int M,int EN); //MMMSAA算法
void printPartitionResult(const vector<PartitionResult>& DE); //打印划分结果