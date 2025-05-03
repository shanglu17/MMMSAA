#include "MMMSAA.h"
/*MMM单链接聚类算法
 * 算法输入：D={t1,t2,...,tn}，n个元素的集合；
 *          C=(Cij)m*n，元素间的通信代价矩阵；
 *          M:划分模块总数；
 *          EN:每个模块的任务数的最大上限
 * 算法输出：DE:谱系图（有序三元组<c,k,K>其中c是通信代价阈值，k是类的数目，K是类的集合）
 * 算法步骤：
  c=1+Max(C);
  k=n;
  Kc={{t1},{t2},...,{tn}};
  DE={<c,k,Kc>};
  repeat
   for each Ki,Kj in Kc do
     {
     larc=largest communication cost between all ti in Ki and tj in Kj;
     if larc >= c and |Ki union Kj|<=EN then Kc=Kc-{ki}-{Kj} union {Ki union Kj};
     }
     k=|Kc|;
   DE=DE union {<c,k,Kc>};
  c=c-1;
 until c=0 or k=M
if c=0 and k>M then calling again MMMSAA for K1 as the input.
*/

// 找到C中的最大值-用于计算通信代价阈值
double findMax(const vector<vector<double>> &C)
{
    double maxVal = 0; // 初始化为0
    for (const auto &row : C)
    {
        for (const auto &val : row)
        {
            if (val > maxVal)
            {
                maxVal = val;
            }
        }
    }
    return maxVal;
}

// 计算最大通信代价
double getMaxCommCost(const TaskSet &set1, const TaskSet &set2, const vector<Task> &D, const vector<vector<double>> &C)
{
    double maxCost = 0; // 初始化为0
    for (const auto &task1 : set1)
    {
        for (const auto &task2 : set2)
        {
            int i = find(D.begin(), D.end(), task1) - D.begin(); // 找到task1在D中的索引
            int j = find(D.begin(), D.end(), task2) - D.begin(); // 找到task2在D中的索引
            maxCost = max(maxCost, C[i][j]);                     // 更新最大通信代价
        }
    }
    return maxCost; // 返回最大通信代价
}

// MMMSAA算法
vector<PartitionResult> MMMSAA(vector<Task> &D, vector<vector<double>> &C, int M, int EN)
{
    vector<PartitionResult> DE; // 存储谱系图结果
    int n = D.size();           // 任务数
    double c = 1 + findMax(C);  // 初始化通信代价阈值
    int k = n;                  // 初始化类的数目
    Partition Kc(n);            // 初始化类的集合
    for (int i = 0; i < n; ++i)
    {
        Kc[i].push_back(D[i]); // 将每个任务放入单独的类中
    }
    DE.push_back({c, k, Kc}); // 将初始结果加入谱系图
    while (c > 0 && k > M)
    {                         // 当通信代价阈值大于0且类的数目大于M时继续迭代
        bool merged = false;  // 标记是否有类被合并
        Partition newKc = Kc; // 新的类集合

        for (size_t i = 0; i < Kc.size(); ++i)
        {
            for (size_t j = i + 1; j < Kc.size(); ++j)
            {
                const auto &Ki = Kc[i]; // 获取第i个类
                const auto &Kj = Kc[j]; // 获取第j个类

                double larc = getMaxCommCost(Ki, Kj, D, C); // 计算最大通信代价

                // 如果最大通信代价大于等于阈值且合并后的类大小不超过EN，则进行合并
                int sizeI = 0;
                int sizeJ = 0;
                for (const auto &task : Ki)
                {
                    sizeI += task.originalSize; // 计算第i个类的原始大小
                }
                for (const auto &task : Kj)
                {
                    sizeJ += task.originalSize; // 计算第j个类的原始大小
                }
                if (larc >= c && (sizeI + sizeJ) <= EN)
                {
                    TaskSet mergedSet = Ki;                                  // 创建新的合并类
                    mergedSet.insert(mergedSet.end(), Kj.begin(), Kj.end()); // 合并类

                    newKc.clear(); // 清空新的类集合
                    for (size_t m = 0; m < Kc.size(); ++m)
                    {
                        if (m != i && m != j)
                        { // 如果不是被合并的类，则加入新的类集合
                            newKc.push_back(Kc[m]);
                        }
                    }
                    newKc.push_back(mergedSet); // 将合并后的类加入新的类集合

                    k = newKc.size();         // 更新类的数目
                    Kc = newKc;               // 更新类集合
                    DE.push_back({c, k, Kc}); // 将新的结果加入谱系图
                    merged = true;            // 标记有类被合并
                    break;                    // 退出内层循环
                }
            }
            if (merged)
                break; // 如果有类被合并，则退出外层循环
        }
        if (!merged)
        {
            //测试用例(1)的步长
            c--; // 如果没有类被合并，则降低通信代价阈值
            // //测试用例(2)的步长
            // c-=0.5; // 如果没有类被合并，则降低通信代价阈值
        }
    }

    // 如果通信代价降到0且类的数目大于M，则调用MMMSAA算法进行重新划分
    if (c <= 0 && k > M)
    {
        auto newD = buildNewD(Kc);         // 构建新的任务集合
        auto newC = computeNewC(Kc, D, C); // 计算新的通信代价矩阵
        if (newD.size() == D.size())
        {
            bool same = true;
            for (size_t i = 0; i < D.size(); ++i)
            {
                if (newD[i].id != D[i].id || newD[i].originalSize != D[i].originalSize)
                {
                    same = false;
                    break;
                }
            }
            if (same)
            {
                cerr << "No actual reduction after merging, breaking recursion at size = " << D.size() << endl;
                return DE;
            }
        }
        auto subDE = MMMSAA(newD, newC, M, EN); // 递归调用MMMSAA算法
        // subDE中的任务名映射回原始任务组合
        map<string, TaskSet> groupMap;
        for (size_t i = 0; i < Kc.size(); ++i)
        {
            groupMap["G" + to_string(i)] = Kc[i];
        }

        // 将 subDE 中的任务重新映射回原任务集
        for (auto &res : subDE)
        {
            Partition realKc;
            for (const auto &group : res.Kc)
            {
                TaskSet expanded;
                for (const auto &task : group)
                {
                    const auto &original = groupMap[task.id]; // task 为 Gx
                    expanded.insert(expanded.end(), original.begin(), original.end());
                }
                realKc.push_back(expanded);
            }
            res.Kc = realKc;
            res.k = realKc.size();
        }
        DE.insert(DE.end(), subDE.begin(), subDE.end()); // 将新的谱系图结果加入到原谱系图中
    }
    return DE; // 返回谱系图结果
}

// 计算新的通信代价矩阵
vector<vector<double>> computeNewC(const Partition &Kc, const vector<Task> &oldD, const vector<vector<double>> &oldC)
{
    int n = Kc.size();                                    // 获取新的类的数目
    vector<vector<double>> newC(n, vector<double>(n, 0)); // 初始化新的通信代价矩阵
    map<string, int> taskIndexMap;                        // 任务ID到索引的映射
    for (int i = 0; i < oldD.size(); ++i)
    {
        taskIndexMap[oldD[i].id] = i; // 将任务ID和索引映射
    }
    for (int i = 0; i < n; ++i)
    {
        for (int j = i + 1; j < n; ++j)
        {
            double sumComm = 0; // 初始化最大通信代价
            for (const auto &taskA : Kc[i])
            {
                for (const auto &taskB : Kc[j])
                {
                    int indexA = taskIndexMap[taskA.id]; // 获取任务A的索引
                    int indexB = taskIndexMap[taskB.id]; // 获取任务B的索引
                    sumComm += oldC[indexA][indexB];     // 更新最大通信代价
                }
            }
            newC[i][j] = newC[j][i] = sumComm; // 更新新的通信代价矩阵
        }
    }
    return newC; // 返回新的通信代价矩阵
}

// 构建新的任务集合
vector<Task> buildNewD(const Partition &Kc)
{
    vector<Task> newD; // 新的任务集合
    for (size_t i = 0; i < Kc.size(); ++i)
    {
        int originalSize = 0; // 初始化原始大小
        for (const auto &task : Kc[i])
        {
            originalSize += task.originalSize; // 累加原始大小
        }
        newD.push_back({"G" + to_string(i), originalSize}); // 将合并后的任务加入新的任务集合
    }
    return newD; // 返回新的任务集合
}
// 计算模块间通信代价
double computeModuleC(const TaskSet &moduleA, const TaskSet &moduleB, const vector<Task> &D, const vector<vector<double>> &C)
{
    double totalCost = 0; // 初始化总通信代价
    for (const auto &taskA : moduleA)
    {
        int indexA = getTaskIndex(D, taskA); // 找到taskA在D中的索引
        for (const auto &taskB : moduleB)
        {
            int indexB = getTaskIndex(D, taskB); // 找到taskB在D中的索引
            totalCost += C[indexA][indexB];      // 累加通信代价
        }
    }
    return totalCost; // 返回总通信代价
}
// 获取任务在集合中的索引
int getTaskIndex(const vector<Task> &D, const Task &task)
{
    auto it = find(D.begin(), D.end(), task); // 查找任务在集合中的位置
    if (it != D.end())
    {
        return distance(D.begin(), it); // 返回索引
    }
    return -1; // 如果未找到，返回-1
}
// 打印划分结果
void printPartitionResult(const vector<PartitionResult> &DE, const vector<Task> &D, const vector<vector<double>> &C)
{
    for (const auto &[c, k, K] : DE)
    {
        cout << "c= " << c << ", k= " << k << ", K= {";
        for (const auto &set : K)
        {
            cout << "{";
            for (const auto &task : set)
            {
                cout << task.id << ", ";
            }
            cout << "}, ";
        }
        cout << "}" << endl;
    }
    const auto &[c, k, K] = DE.back(); // 获取最后一个结果
    double totalCost = 0;              // 初始化总通信代价
    for (size_t i = 0; i < K.size(); ++i)
    {
        for (size_t j = i + 1; j < K.size(); ++j)
        {
            double commCost = computeModuleC(K[i], K[j], D, C); // 计算模块间通信代价
            totalCost += commCost;                              // 累加通信代价
        }
    }
    cout << "Final Communication Cost: " << totalCost << endl; // 打印最终通信代价
}