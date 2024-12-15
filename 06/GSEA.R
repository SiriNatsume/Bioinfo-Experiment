

#########################1.数据预处理#################################
library(clusterProfiler)
library(ggthemes)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(enrichplot)
library(readxl)
library(tibble)
library(deeptime)
library(GOplot)
library(showtext)
library(data.table)
library(ComplexHeatmap)
library(AnnotationHub)
library(RColorBrewer)
library(GseaVis)
library(ggsankey)
library(aPEAR)
library(pathview)
library(forcats)
library(ClusterGVis)
library(multienrichjam)
library(fgsea)
library(scales)
library(ggsci)
library(cols4all)
library(circlize)
library(org.Hs.eg.db)#人类的基因集


pvalueFilter=0.05      #p值的过滤条件
adjPvalFilter=0.05     #矫正后p值的过滤条件
#定义颜色
colorSel="p.adjust"
if(adjPvalFilter>0.05){
  colorSel="pvalue"}

##加载GEO数据
load("D:/R语言生信分析/数据处理/3.GEO数据下载与ID转换/输出数据/GEO数据处理.RData")


########################limma差异分析##################################
data=as.data.frame(exp) #数据备份
range(data)
levels(group_list)

# 设计矩阵
design <- model.matrix(~0 + factor(group_list))
colnames(design) <- levels(factor(group_list))  # 注意疾病组和对照组，不要搞反了
rownames(design) = colnames(data)
design
# 进行线性拟合
fit <- lmFit(data, design)
# 创建对比矩阵
#和上面两个包一样，需要说明是谁比谁:其中rev(levels(group_list))的第一个需要为实验组名字，目的是保证logFC的正值代表实验组上调
cont.matrix <- makeContrasts("A_vs_B" = A-B,
                             "A_vs_C" = A-C,
                             "B_vs_C" = B-C,
                             levels = colnames(design))
cont.matrix
# 应用对比矩阵并进行贝叶斯修正
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
#plotSA(fit2) 查看分布
# 输出所有基因的差异表达情况
colnames(fit2)
A_vs_B <- topTreat(fit2, coef=1, n=Inf)
A_vs_C <- topTreat(fit2, coef=2, n=Inf)
B_vs_C <- topTreat(fit2, coef=3, n=Inf)

# 保存结果到本地为CSV文件
write.csv(A_vs_B, file = "A_vs_B_results.csv")
write.csv(A_vs_C, file = "A_vs_C_results.csv")
write.csv(B_vs_C, file = "B_vs_C_results.csv")


######################A_vs_B的富集分析结果#############
#自定义准备文件函数：
OrgDb = 'org.Hs.eg.db'   
organism = "hsa"   

ready_genelist <- function(DEG){
  ##5.1 GSEA分析:
  #先制作GSEA的输入文件input(使用全部基因而非仅差异基因),这里选择筛选前的差异基因
  genelist=DEG$log2FoldChange  #注意修改 提取logFC，方便进行排序
  names(genelist) <- DEG$gene
  gene_convert <- bitr(DEG$gene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = OrgDb)#注意修改
  DEG = DEG %>% left_join(gene_convert,by=c("gene"="SYMBOL"))%>%
    arrange(desc(log2FoldChange))#按照avg_log2FC进行降序排序
  genelist <- genelist[names(genelist) %in% gene_convert[,1]]
  names(genelist) <- gene_convert[match(names(genelist),gene_convert[,1]),2]#genelist缺失基因过滤(ID转换中缺失的部分)
  genelist <- sort(genelist, decreasing = T)#将genelist按照log2FC值从高到低进行排序：
  return(genelist)
}


setwd("D:/R语言生信分析/数据处理/3.GEO数据下载与ID转换/输出数据/A_vs_B")
#提取差异基因DEG的ID与logFC
DEG <- A_vs_B
DEG <- data.frame(gene=rownames(DEG),log2FoldChange=DEG$logFC) ##修改这里即可

#输入差异分析结果
genelist <- ready_genelist(DEG = DEG)

GSEA_KEGG <- gseKEGG(
  geneList = genelist,
  organism = "hsa", ##人类：hsa
  pvalueCutoff = 1,
  pAdjustMethod = "BH")
GSEA_KEGG <- setReadable(GSEA_KEGG,
                         OrgDb = "org.Hs.eg.db", 
                         keyType = "ENTREZID")#将ENTREZID重转为symbol
GSEA_KEGG_result <- GSEA_KEGG@result
write.csv(GSEA_KEGG_result, file = c('GSEA_KEGG_result-A_vs_B.csv'))


#1.气泡图绘制：
dotplot(
  GSEA_KEGG,
  x = "GeneRatio",#'GeneRatio' and 'Count'
  color = "p.adjust",#pvalue', 'p.adjust' or 'qvalue'
  showCategory = 20,
  font.size = 12,
  title = "",
  orderBy = "x",
  label_format = 60
)
ggsave('GSEA_KEGG dotplot-1.pdf', width=8, height=8)

#分面点图激活和抑制
dotplot(GSEA_KEGG,label_format = 60,split=".sign")+facet_grid(~.sign)
ggsave('GSEA_KEGG dotplot-2.pdf', width=8, height=8)

# 通过改变scales改变Y轴展示
dotplot(GSEA_KEGG,label_format = 60,split=".sign")+facet_wrap(~.sign,scales = "free")
ggsave('GSEA_KEGG dotplot-3.pdf', width=14, height=14)

#2.山峦图绘制：
ridgeplot(
  GSEA_KEGG,
  showCategory = 20,
  fill = "p.adjust", #可以换为其他的
  core_enrichment = TRUE,
  label_format = 60,
  orderBy = "NES",
  decreasing = FALSE
)
ggsave('GSEA_KEGG ridgeplot-p.adjust.pdf', width=8, height=8)

ridgeplot(
  GSEA_KEGG,
  showCategory = 20,
  fill = "NES", #可以换为其他的
  core_enrichment = TRUE,
  label_format = 60,
  orderBy = "NES",
  decreasing = FALSE
)
ggsave('GSEA_KEGG ridgeplot-NES.pdf', width=8, height=8)


###5.2 GSEA可视化:
##选择前10条
terms = GSEA_KEGG_result$ID[1:10]

#单个通路可视化：选择前面的目的通路：Wnt signaling pathway
###显示通路名称和P值
#aim_pathways = "R-HSA-6809371" #这里输入目的通路的ID
# 循环遍历每个通路 ID
for (aim_pathways in terms) {
p1 = gseaplot2(GSEA_KEGG,
                geneSetID = aim_pathways,
                color = "red",       ##color是enrichment score线的颜色
                rel_heights = c(1.5, 0.5, 1), #子图高度
                subplots = 1:3, #显示哪些子图
                base_size=10,#xy轴标题字的大小
                pvalue_table = T, #是否显示pvalue表
                title = GSEA_KEGG_result[aim_pathways,][,2],
                ES_geom = "line") #"dot"将线转换为点
p1
ggsave(p1,filename = paste0('GSEA_KEGG_gseaplot-',aim_pathways, '-1.pdf'), height=8, width=8)
###不显示P值
p2 = gseaplot2(GSEA_KEGG,
                geneSetID = aim_pathways,
                color = "red",
                rel_heights = c(1.5, 0.5, 1), #子图高度
                subplots = 1:3, #显示哪些子图
                pvalue_table = F, #是否显示pvalue表
                title = GSEA_KEGG_result[aim_pathways,][,2],
                ES_geom = "line") #"dot"将线转换为点
p2
ggsave(p2,filename = paste0('GSEA_KEGG_gseaplot-',aim_pathways, '-2.pdf'), height=8, width=8)
###去掉图片里面的通路名称
GSEA_KEGG2 <- GSEA_KEGG #备份
GSEA_KEGG2@result$Description <- ""
p3 = gseaplot2(GSEA_KEGG2,
                geneSetID = aim_pathways,
                color = "red",
                rel_heights = c(1.5, 0.5, 1), #子图高度
                subplots = 1:3, #显示哪些子图
                pvalue_table = T, #是否显示pvalue表
                title = GSEA_KEGG_result[aim_pathways,][,2],
                ES_geom = "line") #"dot"将线转换为点
p3
ggsave(p3,filename = paste0('GSEA_KEGG_gseaplot-',aim_pathways, '-3.pdf'), height=8, width=8)
}


#多个通路可视化：
p4 <- gseaplot2(GSEA_KEGG,
                geneSetID = terms, #或直接输入基因集ID向量名，如c("hsa04610","hsa00260")
                pvalue_table = F,
                ES_geom = "line")
p4
ggsave(p4, filename = 'GSEA_KEGG_gseaplot多通路-1.pdf', width=13, height=13)

p5 <- gseaplot2(GSEA_KEGG,
                geneSetID = terms, #或直接输入基因集ID向量名，如c("hsa04610","hsa00260")
                color = pal_npg("nrc")(10),
                pvalue_table = TRUE,
                ES_geom = "line")
p5
ggsave(p5, filename = 'GSEA_KEGG_gseaplot多通路-2.pdf', width=13, height=13)

p6 <- gseaplot2(GSEA_KEGG,
                geneSetID = terms, #或直接输入基因集ID向量名，如c("hsa04610","hsa00260")
                color = pal_npg("nrc")(10),
                pvalue_table = F,
                ES_geom = "line")
p6
ggsave(p6, filename = 'GSEA_KEGG_gseaplot多通路-3.pdf', width=13, height=13)

p7 <- gseaplot2(GSEA_KEGG,
                geneSetID = terms, #或直接输入基因集ID向量名，如c("hsa04610","hsa00260")
                pvalue_table = T,
                ES_geom = "line")
p7
ggsave(p7, filename = 'GSEA_KEGG_gseaplot多通路-4.pdf', width=13, height=13)


######################A_vs_C的富集分析结果#############
#自定义准备文件函数：
OrgDb = 'org.Hs.eg.db'   
organism = "hsa"   

ready_genelist <- function(DEG){
  ##5.1 GSEA分析:
  #先制作GSEA的输入文件input(使用全部基因而非仅差异基因),这里选择筛选前的差异基因
  genelist=DEG$log2FoldChange  #注意修改 提取logFC，方便进行排序
  names(genelist) <- DEG$gene
  gene_convert <- bitr(DEG$gene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = OrgDb)#注意修改
  DEG = DEG %>% left_join(gene_convert,by=c("gene"="SYMBOL"))%>%
    arrange(desc(log2FoldChange))#按照avg_log2FC进行降序排序
  genelist <- genelist[names(genelist) %in% gene_convert[,1]]
  names(genelist) <- gene_convert[match(names(genelist),gene_convert[,1]),2]#genelist缺失基因过滤(ID转换中缺失的部分)
  genelist <- sort(genelist, decreasing = T)#将genelist按照log2FC值从高到低进行排序：
  return(genelist)
}


setwd("D:/R语言生信分析/数据处理/3.GEO数据下载与ID转换/输出数据/A_vs_C")
#提取差异基因DEG的ID与logFC
DEG <- A_vs_C
DEG <- data.frame(gene=rownames(DEG),log2FoldChange=DEG$logFC) ##修改这里即可

#输入差异分析结果
genelist <- ready_genelist(DEG = DEG)

GSEA_KEGG <- gseKEGG(
  geneList = genelist,
  organism = "hsa", ##人类：hsa
  pvalueCutoff = 1,
  pAdjustMethod = "BH")
GSEA_KEGG <- setReadable(GSEA_KEGG,
                         OrgDb = "org.Hs.eg.db", 
                         keyType = "ENTREZID")#将ENTREZID重转为symbol
GSEA_KEGG_result <- GSEA_KEGG@result
write.csv(GSEA_KEGG_result, file = c('GSEA_KEGG_result-A_vs_C.csv'))


#1.气泡图绘制：
dotplot(
  GSEA_KEGG,
  x = "GeneRatio",#'GeneRatio' and 'Count'
  color = "p.adjust",#pvalue', 'p.adjust' or 'qvalue'
  showCategory = 20,
  font.size = 12,
  title = "",
  orderBy = "x",
  label_format = 60
)
ggsave('GSEA_KEGG dotplot-1.pdf', width=8, height=8)

#分面点图激活和抑制
dotplot(GSEA_KEGG,label_format = 60,split=".sign")+facet_grid(~.sign)
ggsave('GSEA_KEGG dotplot-2.pdf', width=8, height=8)

# 通过改变scales改变Y轴展示
dotplot(GSEA_KEGG,label_format = 60,split=".sign")+facet_wrap(~.sign,scales = "free")
ggsave('GSEA_KEGG dotplot-3.pdf', width=14, height=14)

#2.山峦图绘制：
ridgeplot(
  GSEA_KEGG,
  showCategory = 20,
  fill = "p.adjust", #可以换为其他的
  core_enrichment = TRUE,
  label_format = 60,
  orderBy = "NES",
  decreasing = FALSE
)
ggsave('GSEA_KEGG ridgeplot-p.adjust.pdf', width=8, height=8)

ridgeplot(
  GSEA_KEGG,
  showCategory = 20,
  fill = "NES", #可以换为其他的
  core_enrichment = TRUE,
  label_format = 60,
  orderBy = "NES",
  decreasing = FALSE
)
ggsave('GSEA_KEGG ridgeplot-NES.pdf', width=8, height=8)


###5.2 GSEA可视化:
##选择前10条
terms = GSEA_KEGG_result$ID[1:10]

#单个通路可视化：选择前面的目的通路：Wnt signaling pathway
###显示通路名称和P值
#aim_pathways = "R-HSA-6809371" #这里输入目的通路的ID
# 循环遍历每个通路 ID
for (aim_pathways in terms) {
  p1 = gseaplot2(GSEA_KEGG,
                 geneSetID = aim_pathways,
                 color = "red",       ##color是enrichment score线的颜色
                 rel_heights = c(1.5, 0.5, 1), #子图高度
                 subplots = 1:3, #显示哪些子图
                 base_size=10,#xy轴标题字的大小
                 pvalue_table = T, #是否显示pvalue表
                 title = GSEA_KEGG_result[aim_pathways,][,2],
                 ES_geom = "line") #"dot"将线转换为点
  p1
  ggsave(p1,filename = paste0('GSEA_KEGG_gseaplot-',aim_pathways, '-1.pdf'), height=8, width=8)
  ###不显示P值
  p2 = gseaplot2(GSEA_KEGG,
                 geneSetID = aim_pathways,
                 color = "red",
                 rel_heights = c(1.5, 0.5, 1), #子图高度
                 subplots = 1:3, #显示哪些子图
                 pvalue_table = F, #是否显示pvalue表
                 title = GSEA_KEGG_result[aim_pathways,][,2],
                 ES_geom = "line") #"dot"将线转换为点
  p2
  ggsave(p2,filename = paste0('GSEA_KEGG_gseaplot-',aim_pathways, '-2.pdf'), height=8, width=8)
  ###去掉图片里面的通路名称
  GSEA_KEGG2 <- GSEA_KEGG #备份
  GSEA_KEGG2@result$Description <- ""
  p3 = gseaplot2(GSEA_KEGG2,
                 geneSetID = aim_pathways,
                 color = "red",
                 rel_heights = c(1.5, 0.5, 1), #子图高度
                 subplots = 1:3, #显示哪些子图
                 pvalue_table = T, #是否显示pvalue表
                 title = GSEA_KEGG_result[aim_pathways,][,2],
                 ES_geom = "line") #"dot"将线转换为点
  p3
  ggsave(p3,filename = paste0('GSEA_KEGG_gseaplot-',aim_pathways, '-3.pdf'), height=8, width=8)
}


#多个通路可视化：
p4 <- gseaplot2(GSEA_KEGG,
                geneSetID = terms, #或直接输入基因集ID向量名，如c("hsa04610","hsa00260")
                pvalue_table = F,
                ES_geom = "line")
p4
ggsave(p4, filename = 'GSEA_KEGG_gseaplot多通路-1.pdf', width=13, height=13)

p5 <- gseaplot2(GSEA_KEGG,
                geneSetID = terms, #或直接输入基因集ID向量名，如c("hsa04610","hsa00260")
                color = pal_npg("nrc")(10),
                pvalue_table = TRUE,
                ES_geom = "line")
p5
ggsave(p5, filename = 'GSEA_KEGG_gseaplot多通路-2.pdf', width=13, height=13)

p6 <- gseaplot2(GSEA_KEGG,
                geneSetID = terms, #或直接输入基因集ID向量名，如c("hsa04610","hsa00260")
                color = pal_npg("nrc")(10),
                pvalue_table = F,
                ES_geom = "line")
p6
ggsave(p6, filename = 'GSEA_KEGG_gseaplot多通路-3.pdf', width=13, height=13)

p7 <- gseaplot2(GSEA_KEGG,
                geneSetID = terms, #或直接输入基因集ID向量名，如c("hsa04610","hsa00260")
                pvalue_table = T,
                ES_geom = "line")
p7
ggsave(p7, filename = 'GSEA_KEGG_gseaplot多通路-4.pdf', width=13, height=13)



######################B_vs_C的富集分析结果#############
#自定义准备文件函数：
OrgDb = 'org.Hs.eg.db'   
organism = "hsa"   

ready_genelist <- function(DEG){
  ##5.1 GSEA分析:
  #先制作GSEA的输入文件input(使用全部基因而非仅差异基因),这里选择筛选前的差异基因
  genelist=DEG$log2FoldChange  #注意修改 提取logFC，方便进行排序
  names(genelist) <- DEG$gene
  gene_convert <- bitr(DEG$gene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = OrgDb)#注意修改
  DEG = DEG %>% left_join(gene_convert,by=c("gene"="SYMBOL"))%>%
    arrange(desc(log2FoldChange))#按照avg_log2FC进行降序排序
  genelist <- genelist[names(genelist) %in% gene_convert[,1]]
  names(genelist) <- gene_convert[match(names(genelist),gene_convert[,1]),2]#genelist缺失基因过滤(ID转换中缺失的部分)
  genelist <- sort(genelist, decreasing = T)#将genelist按照log2FC值从高到低进行排序：
  return(genelist)
}


setwd("D:/R语言生信分析/数据处理/3.GEO数据下载与ID转换/输出数据/B_vs_C")
#提取差异基因DEG的ID与logFC
DEG <- B_vs_C
DEG <- data.frame(gene=rownames(DEG),log2FoldChange=DEG$logFC) ##修改这里即可

#输入差异分析结果
genelist <- ready_genelist(DEG = DEG)

GSEA_KEGG <- gseKEGG(
  geneList = genelist,
  organism = "hsa", ##人类：hsa
  pvalueCutoff = 1,
  pAdjustMethod = "BH")
GSEA_KEGG <- setReadable(GSEA_KEGG,
                         OrgDb = "org.Hs.eg.db", 
                         keyType = "ENTREZID")#将ENTREZID重转为symbol
GSEA_KEGG_result <- GSEA_KEGG@result
write.csv(GSEA_KEGG_result, file = c('GSEA_KEGG_result-B_vs_C.csv'))


#1.气泡图绘制：
dotplot(
  GSEA_KEGG,
  x = "GeneRatio",#'GeneRatio' and 'Count'
  color = "p.adjust",#pvalue', 'p.adjust' or 'qvalue'
  showCategory = 20,
  font.size = 12,
  title = "",
  orderBy = "x",
  label_format = 60
)
ggsave('GSEA_KEGG dotplot-1.pdf', width=8, height=8)

#分面点图激活和抑制
dotplot(GSEA_KEGG,label_format = 60,split=".sign")+facet_grid(~.sign)
ggsave('GSEA_KEGG dotplot-2.pdf', width=8, height=8)

# 通过改变scales改变Y轴展示
dotplot(GSEA_KEGG,label_format = 60,split=".sign")+facet_wrap(~.sign,scales = "free")
ggsave('GSEA_KEGG dotplot-3.pdf', width=14, height=14)

#2.山峦图绘制：
ridgeplot(
  GSEA_KEGG,
  showCategory = 20,
  fill = "p.adjust", #可以换为其他的
  core_enrichment = TRUE,
  label_format = 60,
  orderBy = "NES",
  decreasing = FALSE
)
ggsave('GSEA_KEGG ridgeplot-p.adjust.pdf', width=8, height=8)

ridgeplot(
  GSEA_KEGG,
  showCategory = 20,
  fill = "NES", #可以换为其他的
  core_enrichment = TRUE,
  label_format = 60,
  orderBy = "NES",
  decreasing = FALSE
)
ggsave('GSEA_KEGG ridgeplot-NES.pdf', width=8, height=8)


###5.2 GSEA可视化:
##选择前10条
terms = GSEA_KEGG_result$ID[1:10]

#单个通路可视化：选择前面的目的通路：Wnt signaling pathway
###显示通路名称和P值
#aim_pathways = "R-HSA-6809371" #这里输入目的通路的ID
# 循环遍历每个通路 ID
for (aim_pathways in terms) {
  p1 = gseaplot2(GSEA_KEGG,
                 geneSetID = aim_pathways,
                 color = "red",       ##color是enrichment score线的颜色
                 rel_heights = c(1.5, 0.5, 1), #子图高度
                 subplots = 1:3, #显示哪些子图
                 base_size=10,#xy轴标题字的大小
                 pvalue_table = T, #是否显示pvalue表
                 title = GSEA_KEGG_result[aim_pathways,][,2],
                 ES_geom = "line") #"dot"将线转换为点
  p1
  ggsave(p1,filename = paste0('GSEA_KEGG_gseaplot-',aim_pathways, '-1.pdf'), height=8, width=8)
  ###不显示P值
  p2 = gseaplot2(GSEA_KEGG,
                 geneSetID = aim_pathways,
                 color = "red",
                 rel_heights = c(1.5, 0.5, 1), #子图高度
                 subplots = 1:3, #显示哪些子图
                 pvalue_table = F, #是否显示pvalue表
                 title = GSEA_KEGG_result[aim_pathways,][,2],
                 ES_geom = "line") #"dot"将线转换为点
  p2
  ggsave(p2,filename = paste0('GSEA_KEGG_gseaplot-',aim_pathways, '-2.pdf'), height=8, width=8)
  ###去掉图片里面的通路名称
  GSEA_KEGG2 <- GSEA_KEGG #备份
  GSEA_KEGG2@result$Description <- ""
  p3 = gseaplot2(GSEA_KEGG2,
                 geneSetID = aim_pathways,
                 color = "red",
                 rel_heights = c(1.5, 0.5, 1), #子图高度
                 subplots = 1:3, #显示哪些子图
                 pvalue_table = T, #是否显示pvalue表
                 title = GSEA_KEGG_result[aim_pathways,][,2],
                 ES_geom = "line") #"dot"将线转换为点
  p3
  ggsave(p3,filename = paste0('GSEA_KEGG_gseaplot-',aim_pathways, '-3.pdf'), height=8, width=8)
}


#多个通路可视化：
p4 <- gseaplot2(GSEA_KEGG,
                geneSetID = terms, #或直接输入基因集ID向量名，如c("hsa04610","hsa00260")
                pvalue_table = F,
                ES_geom = "line")
p4
ggsave(p4, filename = 'GSEA_KEGG_gseaplot多通路-1.pdf', width=13, height=13)

p5 <- gseaplot2(GSEA_KEGG,
                geneSetID = terms, #或直接输入基因集ID向量名，如c("hsa04610","hsa00260")
                color = pal_npg("nrc")(10),
                pvalue_table = TRUE,
                ES_geom = "line")
p5
ggsave(p5, filename = 'GSEA_KEGG_gseaplot多通路-2.pdf', width=13, height=13)

p6 <- gseaplot2(GSEA_KEGG,
                geneSetID = terms, #或直接输入基因集ID向量名，如c("hsa04610","hsa00260")
                color = pal_npg("nrc")(10),
                pvalue_table = F,
                ES_geom = "line")
p6
ggsave(p6, filename = 'GSEA_KEGG_gseaplot多通路-3.pdf', width=13, height=13)

p7 <- gseaplot2(GSEA_KEGG,
                geneSetID = terms, #或直接输入基因集ID向量名，如c("hsa04610","hsa00260")
                pvalue_table = T,
                ES_geom = "line")
p7
ggsave(p7, filename = 'GSEA_KEGG_gseaplot多通路-4.pdf', width=13, height=13)




