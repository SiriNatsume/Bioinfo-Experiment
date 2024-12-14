import pandas as pd
from scipy.stats import ttest_ind, ttest_rel
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import numpy as np

# 读取数据
df = pd.read_csv("mapping_data.csv", index_col=0)  # 第一列为索引

# 确保表达数据为数值型
df = df.apply(pd.to_numeric)

# 将所有行的数据值转换为 2 的次方
df = df.applymap(lambda x: 2**x)

# 按第一行（组名）分组并计算每组的平均值
grouped_df = df.groupby(level=0).mean()
# 对平均值求 log2
grouped_df = grouped_df.applymap(lambda x: np.log2(x))

# 进行分组

Group_N = grouped_df.loc[['N'], :]
Group_T = grouped_df.loc[['T'], :]
Group_U = grouped_df.loc[['U'], :]
df_N = df[df.index == 'N']  # 提取所有标记为 N 的行
df_T = df[df.index == 'T']  # 提取所有标记为 T 的行
df_U = df[df.index == 'U']  # 提取所有标记为 U 的行

TNlog2FC = []
UNlog2FC = []
p_TN = []
p_UN = []


for protein in Group_N.columns:
    N_values = Group_N[protein].values
    T_values = Group_T[protein].values
    U_values = Group_U[protein].values
    N_df = df_N[protein]
    T_df = df_T[protein]
    U_df = df_U[protein]

    # TN t检验
    t_TN, p_tn = ttest_ind(T_df, N_df, equal_var=False)
    p_TN.append(p_tn)

    # UN t检验
    t_UN, p_un = ttest_ind(U_df, N_df, equal_var=False)
    p_UN.append(p_un)

    # TN log2FC
    log2FC_TN = float(T_values - N_values)

    # UN log2FC
    log2FC_UN = float(U_values - N_values)

    TNlog2FC.append(log2FC_TN)
    UNlog2FC.append(log2FC_UN)

# 调整 p 值
adjusted_p_TN = multipletests(p_TN, method='fdr_bh')[1]
adjusted_p_UN = multipletests(p_UN, method='fdr_bh')[1]

# 输出比较结果
res_df = pd.DataFrame({'TNlog2FC': TNlog2FC, 'UNlog2FC': UNlog2FC, 'TN_PValue': adjusted_p_TN, 'UN_PValue': adjusted_p_UN,})
res_df.index = Group_N.columns


# 添加 -log10 转换的列
res_df['-log10(TN_PValue)'] = -np.log10(res_df['TN_PValue'])
res_df['-log10(UN_PValue)'] = -np.log10(res_df['UN_PValue'])

# ---- 标注前10和后10 ---- #
# 对 T vs N 按 TNlog2FC 排序，提取前后10
top_TN_proteins = res_df.sort_values(by='TNlog2FC', ascending=False).head(10)
bottom_TN_proteins = res_df.sort_values(by='TNlog2FC', ascending=True).head(10)

# 对 U vs N 按 UNlog2FC 排序，提取前后10
top_UN_proteins = res_df.sort_values(by='UNlog2FC', ascending=False).head(10)
bottom_UN_proteins = res_df.sort_values(by='UNlog2FC', ascending=True).head(10)

# 绘制火山图
plt.figure(figsize=(16, 10))

# 绘制 T vs N 的散点
plt.scatter(res_df['TNlog2FC'], res_df['-log10(TN_PValue)'], color='grey', label='T vs N', alpha=0.7)

# 绘制 U vs N 的散点
plt.scatter(res_df['UNlog2FC'], res_df['-log10(UN_PValue)'], color='lightgrey', label='U vs N', alpha=0.7)

# 特殊上色 T vs N 前10和后10
plt.scatter(top_TN_proteins['TNlog2FC'], top_TN_proteins['-log10(TN_PValue)'],
            color='darkred', label='T vs N: Top 10 Proteins')
plt.scatter(bottom_TN_proteins['TNlog2FC'], bottom_TN_proteins['-log10(TN_PValue)'],
            color='red', label='T vs N: Bottom 10 Proteins')

# 特殊上色 U vs N 前10和后10
plt.scatter(top_UN_proteins['UNlog2FC'], top_UN_proteins['-log10(UN_PValue)'],
            color='darkblue', label='U vs N: Top 10 Proteins')
plt.scatter(bottom_UN_proteins['UNlog2FC'], bottom_UN_proteins['-log10(UN_PValue)'],
            color='blue', label='U vs N: Bottom 10 Proteins')

# 标注 T vs N 的前10和后10
for idx in top_TN_proteins.index:
    plt.text(res_df.loc[idx, 'TNlog2FC'], res_df.loc[idx, '-log10(TN_PValue)'], idx, color='darkred', fontsize=8)
for idx in bottom_TN_proteins.index:
    plt.text(res_df.loc[idx, 'TNlog2FC'], res_df.loc[idx, '-log10(TN_PValue)'], idx, color='pink', fontsize=8)

# 标注 U vs N 的前10和后10
for idx in top_UN_proteins.index:
    plt.text(res_df.loc[idx, 'UNlog2FC'], res_df.loc[idx, '-log10(UN_PValue)'], idx, color='darkblue', fontsize=8)
for idx in bottom_UN_proteins.index:
    plt.text(res_df.loc[idx, 'UNlog2FC'], res_df.loc[idx, '-log10(UN_PValue)'], idx, color='lightblue', fontsize=8)

# 样式
plt.axhline(-np.log10(0.05), color='green', linestyle='--', label='Significance Threshold (p=0.05)')
plt.xlabel('Log2 Fold Change')
plt.ylabel('-Log10(P-Value)')
plt.title('Q1 c) Volcano Plot: Top/Bottom 10 Proteins Highlighted (T vs N and U vs N)')

# 图例设置外置
plt.legend(loc='upper left', bbox_to_anchor=(1.05, 1), fontsize=10)
plt.tight_layout(rect=[0, 0, 0.8, 1])  # 调整布局，为图例留出空间
plt.grid(alpha=0.4)

# 展示图
plt.show()


# 提取 T vs N 中前 10 和后 10 蛋白的名字
top_TN_names = top_TN_proteins.index.tolist()
bottom_TN_names = bottom_TN_proteins.index.tolist()

# 提取 U vs N 中前 10 和后 10 蛋白的名字
top_UN_names = top_UN_proteins.index.tolist()
bottom_UN_names = bottom_UN_proteins.index.tolist()

# 打印结果
print(top_TN_proteins)
print("T vs N Top 10 Proteins:")
print(top_TN_names)
print("\nT vs N Bottom 10 Proteins:")
print(bottom_TN_names)

print("\nU vs N Top 10 Proteins:")
print(top_UN_names)
print("\nU vs N Bottom 10 Proteins:")
print(bottom_UN_names)


# 保存到新的 CSV 文件
# grouped_df.to_csv("averaged_expression.csv")