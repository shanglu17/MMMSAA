#include <iostream>
#include <vector>
#include <set>
#include <tuple>
#include <algorithm>
#include <map>
using namespace std;

//类型别名
using Task=string; //任务类型
using TaskSet=vector<Task>; //任务集合类型
using Partition=vector<TaskSet>; //划分类型
struct PartitionResult
{
    double c; //通信代价
    int k; //类的数目
    Partition Kc; //类集合
};//谱系图结果类型

//函数声明
double findMax(const vector<vector<double>>& C); //查找最大值
double getMaxCommCost(const TaskSet& set1, const TaskSet& set2, const vector<Task>&D,const vector<vector<double>>&C); //计算最大通信代价
vector<PartitionResult>MMMSAA(vector<Task>&D,vector<vector<double>>&C,int M,int EN); //MMMSAA算法
vector<vector<double>>computeNewC(const Partition& Kc, const vector<Task>& D, const vector<vector<double>>& C); //计算新的通信代价矩阵
vector<Task>buildNewD(int count); //构建新的任务集合
double computeModuleC(const TaskSet& moduleA, const TaskSet& moduleB, const vector<Task>& D, const vector<vector<double>>& C); //计算模块间通信代价
int getTaskIndex(const vector<Task>& D, const Task& task); //获取任务在集合中的索引
void printPartitionResult(const vector<PartitionResult> &DE,const vector<Task>&D,const vector<vector<double>>&C); //打印划分结果