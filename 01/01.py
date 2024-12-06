import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
import gseapy as gp
import matplotlib.pyplot as plt

# --------以下为配置部分--------
# 输入任务名称（重复的任务名称将导致上一次任务的结果被覆盖）
task: str = ''  # task: str = 'B-C'
# 分组信息（保证只有 A、B 两组，且不可更改组名，仅更改 "sample"）
group_info = {
    'sample1': 'A',
    'sample2': 'A',
    'sample3': 'B',
    'sample4': 'B'
}
# group_info = {
#     # 'CS_1': 'A',
#     # 'CS_2': 'A',
#     # 'CSR100_1': 'B',
#     # 'CSR100_2': 'B',
#     'CSR100_1': 'A',
#     'CSR100_2': 'A',
#     'FBS_1': 'B',
#     'FBS_2': 'B'
# }
# 使用的基因集路径（也可使用 KEGG 等在线基因集）
gene_sets: str = 'data/h.all.v2024.1.Hs.symbols.gmt'

# --------以下为创建者实验课题目映射--------
# CS : A
# CSR100 : B
# FBS : C


# --------以下为非修改部分--------
# 读取存储数据的 CSV 文件
data = pd.read_csv(f'data/{task}.csv')
# 分组 DataFrame
group_df = pd.DataFrame.from_dict(group_info, orient='index', columns=['Group'])
# 除 0
filtered_data = data.loc[(data.iloc[:, 1:] > 0).any(axis=1)]  # 样本表达量从第2列开始
# 提取并清洗表达数据
frag_counts = filtered_data.iloc[:, 1:]
frag_counts = frag_counts.dropna()
frag_counts = frag_counts.groupby(frag_counts.index).mean()
frag_counts.columns = group_info.keys()
# 分离分组
Group_A = frag_counts.loc[:, group_df['Group'] == 'A']
Group_B = frag_counts.loc[:, group_df['Group'] == 'B']

# 差异表达分析 (t 检验)
p_values = []
log_fold_changes = []

for gene in frag_counts.index:
    A_values = Group_A.loc[gene].values
    B_values = Group_B.loc[gene].values

    # t 检验
    t_stat, p_val = ttest_ind(B_values, A_values, equal_var=False)
    p_values.append(p_val)

    # 计算对数倍数变化
    mean_control = np.mean(A_values)
    mean_treatment = np.mean(B_values)
    log_fc = np.log2(mean_treatment + 1) - np.log2(mean_control + 1)
    log_fold_changes.append(log_fc)

# 调整 p 值
adjusted_p_values = multipletests(p_values, method='fdr_bh')[1]

# 构建差异表达表
de_results = pd.DataFrame({
    'Gene': filtered_data['Gene_symbol'],
    'Log2FC': log_fold_changes,
    'PValue': p_values,
    'AdjustedPValue': adjusted_p_values
})

# 筛选显著基因
significant_genes = de_results[(de_results['AdjustedPValue'] < 0.05) & (abs(de_results['Log2FC']) > 0.5)]
# 排序基因用于 GSEA
gene_rank = de_results[['Gene', 'Log2FC']].set_index('Gene').sort_values('Log2FC', ascending=False)
# 清洗排序基因列表
gene_rank = gene_rank.dropna()
gene_rank = gene_rank[~gene_rank.index.duplicated(keep='first')]
# 保存排序基因列表至新文件
gene_rank.to_csv(f'{task}_gene_ranked_list.txt', sep='\t')

# 运行 GSEA
enr = gp.prerank(
    rnk=f'{task}_gene_ranked_list.txt',  # 排序文件
    gene_sets=gene_sets,  # 使用的基因集
    outdir=f'{task}_gsea_output',  # 输出目录
    permutation_num=100  # 默认置换次数
)

# terms = enr.res2d.Term
# axs = enr.plot(terms=terms[1:11], ofname=f"result.png", show_ranking=True)

# 筛选 NES 为正数的前十基因集
top10_positive_nes = enr.res2d[enr.res2d['NES'] > 0].sort_values('NES', ascending=False).head(10)
# 筛选 NES 为倒数的前十基因集
bottom10_negative_nes = enr.res2d[enr.res2d['NES'] < 0].sort_values('NES', ascending=True).head(10)

# 选择相关列
table_data = top10_positive_nes[['Term', 'NES', 'FDR q-val']].reset_index(drop=True)
table_data.columns = ['Gene Set', 'NES', 'Adjusted P-value']  # 重命名列

# 绘制正数表格
_, ax = plt.subplots(figsize=(8, 4))
ax.axis('tight')
ax.axis('off')

# 创建表格
table = ax.table(cellText=table_data.values,
                 colLabels=table_data.columns,
                 loc='center',
                 cellLoc='center',
                 colColours=['#f2f2f2'] * 3)

# 调整表格样式
table.auto_set_font_size(False)
table.set_fontsize(10)
table.auto_set_column_width(col=list(range(len(table_data.columns))))

# 保存图像
plt.savefig(f'{task}_top10_positive_nes_table.png', dpi=300, bbox_inches='tight')

# 选择相关列
table_data = bottom10_negative_nes[['Term', 'NES', 'FDR q-val']].reset_index(drop=True)
table_data.columns = ['Gene Set', 'NES', 'Adjusted P-value']  # 重命名列

# 绘制倒数表格
_, ax = plt.subplots(figsize=(8, 4))
ax.axis('tight')
ax.axis('off')

# 创建表格
table = ax.table(cellText=table_data.values,
                 colLabels=table_data.columns,
                 loc='center',
                 cellLoc='center',
                 colColours=['#f2f2f2'] * 3)

# 调整表格样式
table.auto_set_font_size(False)
table.set_fontsize(10)
table.auto_set_column_width(col=list(range(len(table_data.columns))))

# 保存图像
plt.savefig(f'{task}_bottom10_negative_nes_table.png', dpi=300, bbox_inches='tight')
