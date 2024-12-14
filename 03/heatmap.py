import pandas as pd
from scipy.stats import zscore

# 加载数据
df = pd.read_csv("Pdata.csv", index_col=0)
# df = df.T
# Z-score 标准化
zscore_df = df.apply(zscore, axis=1)

import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage

# 计算距离矩阵和层次聚类
linkage_matrix = linkage(zscore_df.T, method='ward')  # 对列（样本）聚类
sns.clustermap(zscore_df, cmap="vlag", method='ward', figsize=(10, 8))

plt.show()