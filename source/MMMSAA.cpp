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
                if (larc >= c && (Ki.size() + Kj.size()) <= EN)
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
            c--; // 如果没有类被合并，则降低通信代价阈值
        }
    }

    // 如果通信代价降到0且类的数目大于M，则调用MMMSAA算法进行重新划分
    while (k > M)
    {
        if (c == 0)
        {
            // 1.将当前Kc的每个模块视为新任务
            vector<Task> newD = buildNewD(Kc.size()); // 构建新的任务集合

            // 2.计算新的通信代价矩阵
            vector<vector<double>> newC = computeNewC(Kc, newD, C); // 计算新的通信代价矩阵

            // 3.聚类，模块元素数量限制在EN内
            auto updated_MMMSAA = [&](const vector<Task> &newD, const vector<vector<double>> &newC, const Partition &oldKc)
            {
                // 统计每个新模块的原始任务
                map<string, int> taskSize;
                for (size_t i = 0; i < oldKc.size(); ++i)
                {
                    taskSize["G" + to_string(i)] = oldKc[i].size(); // 统计每个新模块的原始任务数
                }
                vector<PartitionResult> subDE; // 存储新的谱系图结果
                double c2 = 1 + findMax(newC);    // 初始化新的通信代价阈值
                int k2 = oldKc.size();         // 初始化新的类的数目
                Partition Kc2(k2);             // 初始化新的类集合
                for (const auto &t : newD)
                {
                    Kc2.push_back({t}); // 将每个新任务放入单独的类中
                }
                subDE.push_back({c2, k2, Kc2}); // 将初始结果加入新的谱系图
                while (c2 > 0 && k2 > M)
                {                           // 当新的通信代价阈值大于0且新的类的数目大于M时继续迭代
                    bool merged = false;    // 标记是否有类被合并
                    Partition newKc2 = Kc2; // 新的类集合

                    for (size_t i = 0; i < Kc2.size(); ++i)
                    {
                        for (size_t j = i + 1; j < Kc2.size(); ++j)
                        {
                            const auto &Ki = Kc2[i]; // 获取第i个类
                            const auto &Kj = Kc2[j]; // 获取第j个类

                            double larc = getMaxCommCost(Ki, Kj, newD, newC); // 计算最大通信代价

                            // 统计合并后的类大小
                            int mergedSize = 0;
                            for (const auto &task : Ki)
                            {
                                mergedSize += taskSize[task]; // 统计合并后的类大小
                            }
                            for (const auto &task : Kj)
                            {
                                mergedSize += taskSize[task]; // 统计合并后的类大小
                            }
                            // 如果最大通信代价大于等于阈值且合并后的类大小不超过EN，则进行合并
                            if (larc >= c2 && mergedSize <= EN)
                            {
                                TaskSet mergedSet = Ki;                                  // 创建新的合并类
                                mergedSet.insert(mergedSet.end(), Kj.begin(), Kj.end()); // 合并类

                                newKc2.clear(); // 清空新的类集合
                                for (size_t m = 0; m < Kc2.size(); ++m)
                                {
                                    if (m != i && m != j)
                                    { // 如果不是被合并的类，则加入新的类集合
                                        newKc2.push_back(Kc2[m]);
                                    }
                                }
                                newKc2.push_back(mergedSet); // 将合并后的类加入新的类集合

                                k2 = newKc2.size();             // 更新新的类的数目
                                Kc2 = newKc2;                   // 更新新的类集合
                                subDE.push_back({c2, k2, Kc2}); // 将新的结果加入谱系图
                                merged = true;                  // 标记有类被合并
                                break;                          // 退出内层循环
                            }
                        }
                        if (merged)
                            break; // 如果有类被合并，则退出外层循环
                    }
                    if (!merged)
                    {
                        c2--; // 如果没有类被合并，则降低通信代价阈值
                    }
                }
                return subDE; // 返回新的谱系图结果
            };
            // 4.调用更新后的MMMSAA算法
            auto newDE = updated_MMMSAA(newD, newC, Kc);     // 调用更新后的MMMSAA算法
            if (newDE.back().k >= k) {
                break; // 没有进一步减少类的数量，跳出死循环
            }
            DE.insert(DE.end(), newDE.begin(), newDE.end()); // 将新的谱系图结果加入原有结果
            
            //更新Kc和k:根据newDE的最后聚类结果重新构造Kc
            Kc=newDE.back().Kc; // 更新Kc为最后的聚类结果
            k=newDE.back().k; // 更新k为最后的类的数目
            D = newD;   // 保证后续使用的任务集合与Kc对应
            C = newC;   // 保证通信矩阵和任务集合一致
        }
        if (k <= M)
            break; // 如果类的数目小于等于M，则退出循环
    }
    return DE; // 返回谱系图结果
}

// 计算新的通信代价矩阵
vector<vector<double>> computeNewC(const Partition &Kc, const vector<Task> &D, const vector<vector<double>> &C)
{
    int k = Kc.size();                                    // 获取类的数目
    vector<vector<double>> newC(k, vector<double>(k, 0)); // 初始化新的通信代价矩阵
    for (size_t i = 0; i < k; ++i)
    {
        for (size_t j = i + 1; j < k; ++j)
        {
            double maxCost = INT_MIN; // 初始化最大通信代价
            for (const auto &task1 : Kc[i])
            { // 遍历第i个类的任务
                for (const auto &task2 : Kc[j])
                {                                                              // 遍历第j个类的任务
                    int index_i = find(D.begin(), D.end(), task1) - D.begin(); // 找到task1在D中的索引
                    int index_j = find(D.begin(), D.end(), task2) - D.begin(); // 找到task2在D中的索引
                    maxCost = max(maxCost, C[index_i][index_j]);               // 更新最大通信代价
                }
            }
            newC[i][j] = maxCost; // 更新新的通信代价矩阵
            newC[j][i] = maxCost; // 对称赋值
        }
    }
    return newC; // 返回新的通信代价矩阵
}

// 构建新的任务集合
vector<Task> buildNewD(int count)
{
    vector<Task> newD; // 新的任务集合
    for (int i = 0; i < count; ++i)
    {
        newD.push_back("G" + to_string(i)); // 将每个类命名为G0,G1,...,Gn
    }
    return newD; // 返回新的任务集合
}
// 打印划分结果
void printPartitionResult(const vector<PartitionResult> &DE)
{
    for (const auto &[c, k, K] : DE)
    {
        cout << "c= " << c << ", k= " << k << ", K= {";
        for (const auto &set : K)
        {
            cout << "{";
            for (const auto &task : set)
            {
                cout << task << ", ";
            }
            cout << "}, ";
        }
        cout << "}" << endl;
    }
}