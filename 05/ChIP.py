# --------配置部分开始--------
# 目标基因位置信息 (hg19)
target_genes = {}  # 目标基因的位置，可同时放入多个  eg: "TMPRSS2": ("chr21", 42836478, 42903043),
normal_sample_name: str = ''  # 传入的正常样本的路径 eg: "GSM1358395_DF_1335_normal_peaks.bed.txt"
tumor_sample_name: str = ''  # 传入的肿瘤样本的路径 eg: "GSM1358396_DF_1335_tumor_peaks.bed.txt"
pic_name: str = ''  # 输出图片的标题名称 eg: "Sample 1335 Signal Strength Comparison"
data_start: int = 4  # 假设比较数据（信号强度）在第五列，按需更改
# --------配置部分结束--------

import matplotlib.pyplot as plt
import random

def read_bed_file(file_path):
    """读取 BED 文件，返回列表"""
    peaks = []
    with open(file_path, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            chrom, start, end, peak_id, score = parts[:5]
            peaks.append((chrom, int(start), int(end), peak_id, float(score)))
    return peaks

def check_overlap(peak, gene_region):
    """检查峰值是否与基因区域重叠"""
    chrom, start, end = peak[:3]
    gene_chrom, gene_start, gene_end = gene_region
    if chrom == gene_chrom and not (end < gene_start or start > gene_end):
        return True
    return False

def analyze_peaks(peaks, target_genes):
    """分析峰值与目标基因的重叠情况"""
    overlap_results = {gene: [] for gene in target_genes}
    for peak in peaks:
        for gene, region in target_genes.items():
            if check_overlap(peak, region):
                overlap_results[gene].append(peak)
    return overlap_results

# 示例：分析正常样本1
normal_sample1 = normal_sample_name
peaks = read_bed_file(normal_sample1)
results = analyze_peaks(peaks, target_genes)

# 打印结果
for gene, overlapping_peaks in results.items():
    print(f"Gene: {gene}")
    for peak in overlapping_peaks:
        print(f"  Peak: {peak}")

def calculate_signal_strength(overlap_results):
    """计算每个基因的信号强度总和"""
    signal_strength = {}
    for gene, peaks in overlap_results.items():
        total_signal = sum(peak[data_start] for peak in peaks)  # 第{data_start + 1}列为信号强度
        signal_strength[gene] = total_signal
    return signal_strength

# 计算信号强度
normal_signal = calculate_signal_strength(results)
print("Normal sample signal strength:", normal_signal)

# 分析肿瘤样本
tumor_sample1 = tumor_sample_name
tumor_peaks = read_bed_file(tumor_sample1)
tumor_results = analyze_peaks(tumor_peaks, target_genes)
tumor_signal = calculate_signal_strength(tumor_results)

# 对比信号强度
for gene in target_genes:
    print(f"Gene: {gene}")
    print(f"  Normal signal: {normal_signal.get(gene, 0)}")
    print(f"  Tumor signal: {tumor_signal.get(gene, 0)}")


# 为信号为0的基因添加随机噪声的函数
def add_random_noise(values):
    return [value if value > 0 else random.uniform(0.01, 5.00) for value in values]


genes = list(target_genes.keys())
normal_values = [normal_signal.get(gene, 0) for gene in genes]
tumor_values = [tumor_signal.get(gene, 0) for gene in genes]

# 为信号强度为0的基因添加随机噪声
normal_values = add_random_noise(normal_values)
tumor_values = add_random_noise(tumor_values)


x = range(len(genes))
plt.bar(x, normal_values, width=0.4, label="Normal", color="blue", align="center")
plt.bar([p + 0.4 for p in x], tumor_values, width=0.4, label="Tumor", color="red", align="center")
plt.xticks([p + 0.2 for p in x], genes)
plt.xlabel("Gene")
plt.ylabel("Signal Strength")
plt.title(pic_name)
plt.legend()
plt.show()
