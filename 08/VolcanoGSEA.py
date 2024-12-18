import gseapy as gp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm
import os
import sys


print(f'Loading: 正在读取并检查输入文件。')
try :
    # 读取肿瘤数据
    t_df = pd.read_csv("t.csv", index_col=0)  # 第一列为索引
    t_df = t_df.groupby(t_df.index).mean()  # 合并同名项并计算平均
    t_df = t_df.T
    n_df = pd.read_csv('n.csv', index_col=0)
    n_df = n_df.groupby(n_df.index).mean()
    n_df = n_df.T
except FileNotFoundError:
    print("Error: 未检测到输入文件。\n您确认您存放的文件格式正确吗？\n您确认您存放的文件路径正确吗？\n您确认您存放的文件名称正确吗？\n您确认您存放的文件数量完整吗？")
    input(f'输入回车终止程序。')
    exit(1)


# 检测 n_df 与 t_df 的列索引是否一致
if not n_df.columns.equals(t_df.columns):
    print("Error: 输入文件的索引不一致。\n您确认您输入的文件内容正确吗？")
    input(f'输入回车终止程序。')
    exit(1)

print(f'INFO: Success\n')

# 创建输出文件夹
os.makedirs('gsea_output', exist_ok=True)
os.makedirs('output', exist_ok=True)

log2FC = []
p = []

# 遍历每个蛋白质列，计算 log2FC 和 p 值
for protein in tqdm(n_df.columns, desc="Processing T-tests: ", total=len(n_df.columns)):
    N_values = n_df[protein].values
    T_values = t_df[protein].values

    # 检查并跳过包含无效值的列
    if not (np.isfinite(N_values).all() and np.isfinite(T_values).all()):
        print(f"Invalid values detected in column: {protein}, skipping...")
        continue

    # t 检验
    _, p_value = ttest_ind(T_values, N_values, equal_var=False)
    p.append(p_value)


    # 计算 log2FC
    mean_control = np.mean(N_values)
    mean_treatment = np.mean(T_values)
    log_fc = np.log2(mean_treatment + 1) - np.log2(mean_control + 1)
    log2FC.append(log_fc)


# 构建差异表达表
de_res = pd.DataFrame({
    'Gene': n_df.columns,
    'Log2FC': log2FC,
    'PValue': p,
    # 'AdjustedPValue': adjusted_p_values
})
de_res = de_res.dropna(subset=['PValue'])
p_list = de_res['PValue'].tolist()
adj_p_value = []

# 确保 p 是一维的 NumPy 数组
p_list = np.array(p_list)
# 调整 p 值
if len(p_list) > 0:  # 确保存在有效的 p 值
    adj_p_value = multipletests(p_list, method='fdr_bh')[1]
else:
    print("No valid p-values were calculated.")

# 构建差异表达表
de_res['AdjustedPValue'] = adj_p_value
de_res.to_csv('output/de_results.csv', index=False)
# 筛选显著基因
significant_genes = de_res[(de_res['AdjustedPValue'] < 0.05) & (abs(de_res['Log2FC']) > 1.0)]


# 识别前十/后十基因
top10_upregulated = significant_genes.nlargest(10, 'Log2FC')
top10_downregulated = significant_genes.nsmallest(10, 'Log2FC')

# 打印前十/后十基因
print("\n")
print("Top 10 Upregulated Genes:")
print(top10_upregulated[['Gene', 'Log2FC']].to_string(index=False))
print("\nTop 10 Downregulated Genes:")
print(top10_downregulated[['Gene', 'Log2FC']].to_string(index=False))
input(f'请确认复制以上数据后输入回车\n')

# 绘制火山图
plt.figure(figsize=(10, 7))
plt.scatter(de_res['Log2FC'], -np.log10(de_res['AdjustedPValue']), color='gray', alpha=0.5, label='All Genes')

# 重新着色显著基因
plt.scatter(significant_genes['Log2FC'], -np.log10(significant_genes['AdjustedPValue']), alpha=0.7, label='Significant')

# 重新着色前十/后十基因
plt.scatter(top10_upregulated['Log2FC'], -np.log10(top10_upregulated['AdjustedPValue']),
            color='darkred', label='Top 10 Upregulated')
plt.scatter(top10_downregulated['Log2FC'], -np.log10(top10_downregulated['AdjustedPValue']),
            color='darkblue', label='Top 10 Downregulated')

# 标记前十/后十基因
for _, row in top10_upregulated.iterrows():
    plt.text(row['Log2FC'], -np.log10(row['AdjustedPValue']), row['Gene'], fontsize=8, ha='right', color='darkred')

for _, row in top10_downregulated.iterrows():
    plt.text(row['Log2FC'], -np.log10(row['AdjustedPValue']), row['Gene'], fontsize=8, ha='right', color='darkblue')

# 添加 FDR 置信度虚线
plt.axhline(y=-np.log10(0.05), color='green', linestyle='--', linewidth=1, label='P=0.05')

# 添加 log2FC 置信度虚线
plt.axvline(x=1.0, color='gray', linestyle='--', linewidth=1, alpha=0.7, label='|Log2FC|=1.0')
plt.axvline(x=-1.0, color='gray', linestyle='--', linewidth=1, alpha=0.7)

# 添加标题和坐标轴标题
plt.title('Volcano Plot of Differentially Expressed Genes')
plt.xlabel('Log2 Fold Change (Log2FC)')
plt.ylabel('-Log10(FDR)')

# 添加图例
plt.legend()

# 调整布局并输出
plt.tight_layout()
# 保存图像
plt.savefig(f'output/volcano.png', dpi=300, bbox_inches='tight')
# plt.show()


# 排序基因用于 GSEA
gene_rank = de_res[['Gene', 'Log2FC']].set_index('Gene').sort_values('Log2FC', ascending=False)
# 清洗排序基因列表
gene_rank = gene_rank.dropna()
gene_rank = gene_rank[~gene_rank.index.duplicated(keep='first')]
# 保存排序基因列表至新文件
gene_rank.to_csv(f'output/gene_ranked_list.txt', sep='\t')

# 检查是否在 PyInstaller 打包环境
if getattr(sys, 'frozen', False) and hasattr(sys, '_MEIPASS'):
    # 如果是打包的 .exe 文件，临时路径存放在 _MEIPASS
    base_path = sys._MEIPASS
else:
    # 如果未打包，则使用脚本当前目录
    base_path = os.path.dirname(os.path.abspath(__file__))

# 基因集文件路径
gene_sets = os.path.join(base_path, 'h.all.v2024.1.Hs.symbols.gmt')

# 调用 GSEA
enr = gp.prerank(
    rnk=f'output/gene_ranked_list.txt',  # 排序文件
    gene_sets=gene_sets,  # 使用的基因集
    outdir=f'gsea_output',  # 输出目录
    permutation_num=100  # 默认置换次数
)


# 筛选 NES 为正数的前十基因集
top10_positive_nes = enr.res2d[enr.res2d['NES'] > 0].sort_values('NES', ascending=False).head(10)
# 筛选 NES 为倒数的前十基因集
bottom10_negative_nes = enr.res2d[enr.res2d['NES'] < 0].sort_values('NES', ascending=True).head(10)

# 选择相关列
top10_table_data = top10_positive_nes[['Term', 'NES', 'FDR q-val']].reset_index(drop=True)
top10_table_data.columns = ['Gene Set', 'NES', 'Adjusted P-value']  # 重命名列
bottom10_table_data = bottom10_negative_nes[['Term', 'NES', 'FDR q-val']].reset_index(drop=True)
bottom10_table_data.columns = ['Gene Set', 'NES', 'Adjusted P-value']  # 重命名列


# 合并数据
combined_table_data = pd.concat([top10_table_data, bottom10_table_data], axis=0, ignore_index=True)

# 设置颜色：深蓝表示正向（Top 10 Positive），深红表示负向（Bottom 10 Negative）
colors = ['darkred'] * len(top10_table_data) + ['darkblue'] * len(bottom10_table_data)

# 创建绘图
plt.figure(figsize=(14, 7))
plt.bar(combined_table_data['Gene Set'], combined_table_data['NES'], color=colors)

# 添加标题和坐标轴标签
plt.title('Top 10 Positive and Bottom 10 Negative Gene Sets', fontsize=16)
plt.xlabel('Gene Sets', fontsize=14)
plt.ylabel('NES', fontsize=14)

# 设置x轴标签，旋转角度防止重叠
plt.xticks(rotation=45, ha='right', fontsize=12)

# 在柱状图上添加每个条的 FDR 值
for i, row in combined_table_data.iterrows():
    plt.text(i, -2 if row['NES'] > 0 else 0.15,  # 调整到柱状图底部
             f"FDR: {row['Adjusted P-value']:.2e}",  # 标注 NES 和 Adjusted P-value
             ha='center', fontsize=10, color='black', rotation=45)  # 旋转 45 度

# 调整图表布局
plt.tight_layout()
# 保存图像
plt.savefig(f'output/NES.png', dpi=300, bbox_inches='tight')
# 显示图像
# plt.show()

# 提示处理完成
print(f'INFO: Complete\n请检查 output 文件夹查看输出结果。\n请检查 gsea_output 文件夹查看 GSEA 输出结果。')
input(f'输入回车终止程序。')
sys.exit(0)
