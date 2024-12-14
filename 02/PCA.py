import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# 读取数据
map_df = pd.read_csv('map.csv')
data_df = pd.read_csv('Pdata.csv', index_col=0)

# 映射组名
mapping = dict(zip(map_df['Sample_ID'], map_df['Tissue_Type']))
data_df.rename(columns=mapping, inplace=True)
print(data_df)
data_df.T.to_csv('mapping_data.csv', sep=',', header=True, encoding='utf-8', index=True)


# 3. 准备 PCA 数据
A = data_df.T.to_numpy()  # 转换为 NumPy 数组
MEAN = np.mean(A, axis=0)  # 求均值
X = np.subtract(A, MEAN)  # 去中心化
COV = np.dot(X.T, X)  # 协方差矩阵
W, V = np.linalg.eig(COV)  # 特征值和特征向量

# 提取前两主成分
e1, e2 = V.T[0], V.T[1]
z1 = np.dot(X, e1)
z2 = np.dot(X, e2)

# 降维结果
RES = np.array([z1, z2]).T  # 转置为样本 * 主成分矩阵

# 4. 可视化降维结果
# 提取组别信息
group_labels = data_df.columns.tolist()
unique_groups = list(set(group_labels))
colors = plt.cm.tab10.colors  # 使用 Matplotlib 默认调色板
group_to_color = {group: colors[i] for i, group in enumerate(unique_groups)}

# 绘制散点图
plt.figure(figsize=(8, 6))
for group in unique_groups:
    indices = [i for i, label in enumerate(group_labels) if label == group]
    plt.scatter(RES[indices, 0], RES[indices, 1],
                color=group_to_color[group], alpha=0.7, label=group)

# 图表设置
plt.title('Q1 a) PCA Visualization', fontsize=14)
plt.xlabel('Principal Component 1', fontsize=12)
plt.ylabel('Principal Component 2', fontsize=12)
plt.axhline(0, color='gray', linestyle='--', linewidth=0.8)
plt.axvline(0, color='gray', linestyle='--', linewidth=0.8)
plt.legend(title="Groups", fontsize=10)
plt.grid(True, linestyle='--', alpha=0.6)

# 显示图像
plt.show()