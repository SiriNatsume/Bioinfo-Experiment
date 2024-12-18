<div align="center">
     <img alt="img.jpg" height="200" src="https://github.com/SiriNatsume/Bioinfo-Experiment/blob/main/08/img.jpg" width="200"/>
</div>
<h1 align="center">VolcanoGSEA on Python</h1>

## 前言
请您**详细参阅文档对应部分**，使用**已经打包好的 .exe 文件**。该仓库未经修改时可能仅适用于创建者同期 lab8。

## 功能
本仓库实现了
1. 对任意两组基因的 FPKM 数据进行差异化分析
2. 将分析信息输出为火山图
3. 对分析信息执行 GSEA
4. 依据 GSEA 所得 NES 信息绘制前十/后十富集基因集


## 如何使用
1. 从 [release](https://github.com/SiriNatsume/Bioinfo-Experiment/releases) 处下载可执行文件（VolcanoGSEA.exe），解压到任意目录
2. 数据文件格式如下
    - | Hugo_Symbol | TCGA-55-6980-11A-01R-1949-07 | ... |
      |:------------|:-----------------------------|:----|
       | Mar-01      | 155.5                        | ... |
       | ...         | ...                          |     |
3. 准备实验组数据文件，并置于 VolcanoGSEA.exe 所在目录，该文件需保存为 t.csv
4. 准备对照组数据文件，并置于 VolcanoGSEA.exe 所在目录，该文件需保存为 n.csv
5. 双击运行 VolcanoGSEA.exe，等待结果输出
6. 期望输出
    - 根目录/
      - output/
        - de_results.csv （差异分析保存结果）
        - gene_ranked_list.txt（基因排序结果）
        - NES.png（依据 NES 信息绘制前十/后十富集基因**集**结果）
        - volcano.png（依据差异分析绘制前十/后十富集基因结果）
      - gsea_output/
        - prerank/
          - HALLMARK_ADIPOGENESIS.pdf（GSEA 各基因集输出结果）
          - ...
        - gseapy.log（GSEA 日志信息）

## 注意事项
- **请保证对照组与实验组数据文件格式的 Hugo_Symbol 完全一致，同时删去其他所有标签（比如 Eid）**，这是您唯一需要进行的数据处理。
- 如果您的实验内容与创建者不一致，请在阅读源码后确定是否使用。

## 测试环境
- Windows 11 PC
- Windows 11 laptop

## 免责声明
  - 您知悉且接受，由于种种难以预料的因素，最终输出结果可能与标准答案**不相同**。
  - 您知悉且接受，代码会对原始数据进行清洗，但效果**可能不佳**。
  - 您知悉且接受，该仓库地址将被合并在创建者实验报告中**一并提交**。
  - 本仓库旨在分享生物信息学操作流程，**不提倡**抄袭等行为，**不提供**答案。

## 其他
- VolcanoGSEA.exe 内置 hallmark2024 基因集。
- 差异表达分析在本代码中集成实现，分析方法为独立样本 t 检验 + log2 Fold Change + FDR 校正的标准分析流程。
- 最后一份报告了，学点新东西，没想到打包还挺简单的，但是可视化有点麻烦。
- ~~直接运行 python 文件的文档懒得写了，自己看着用吧。~~

## 致谢
- [JetBrains](https://www.jetbrains.com/zh-cn/)
- [ChatGPT](https://www.chatgpt.com)
- Thank you to the professor Xie for the wonderful lectures and to the teaching assistants for their diligent grading of the reports. 

创建者 [@SiriNatsume](https://github.com/SiriNatsume)
祝你愉快 :)
