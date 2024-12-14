<h1 align="center">GSEA on Python</h1>

## 前言
本文档假设您对 python 有初步了解，因此不会涉及 python 环境配置、包管理等内容。如果您没有任何 python 环境配置经验，建议您使用 GSEA 应用而非参照本文档，这往往更快。~~或者您很闲 :D~~
## 基础知识
- 您需要了解以下基础内容来更好地理解并定制化代码。
- 执行 GSEA 前，您需要对原始数据进行差异表达分析，GSEA 仅对分析后的 Log2FC 值进行分析。
- 差异表达分析在本代码中集成实现，分析方法为独立样本 t 检验 + log2 Fold Change + FDR 校正 的标准分析流程，如由特殊需要，可自行更改分析方法。
- 代码使用 GSEApy 的 prerank 方法实现 GSEA，而后通过 matplotlib 另行绘制富集表格。
- 原始数据中可能有多个表达量统计值，您需要选择 **不同** 样本的 **同一类** 统计值进行分析，这也将是您在使用本代码过程中唯一需要输入的原始数据。
## 如何使用
1. 克隆本仓库至本地，不破坏仓库文件目录结构
2. 准备原始数据文件
    - 您需要将原始数据整理成如下格式，以 .csv 格式存放在项目的 data 目录下。
    - | Gene_name | sample1   | sample2   | sample3   | sample4   |
      |:----------|:----------|:----------|:----------|:----------|
      | GeneA     | raw_data1 | raw_data3 | raw_data5 | raw_data7 |
      | GeneB     | raw_data2 | rwa_data4 | raw_data6 | raw_data8 |
    - 首列为 Gene 标识，您可以填写 Eid 、Gene_symbol 、Gene_name 或者任意可分辨 Gene 的内容，这里不作要求。（建议使用 Gene_symbol 直接匹配基因集，否则可能没有富集结果）
    - 后续列为样本表达数据，每一列为一组样品的对应表达原始数据，您可以选择合适的原始数据进行填充，**请保证** 所有样品使用 **同一类型** 的原始数据。
    - data 目录下已经准备了一份 example.csv 文件供参考。
3. 依据代码 01.py import 部分安装所需软件包
4. 更改代码配置部分
    - 您需要对代码 01.py 中标识的配置部分进行修改，接下来将对配置字段进行解释。
    - task：字符串。标注任务名称，每次运行任务生成的文件将以此名称为前缀保存，**重复的任务名称将导致上一次任务的结果被覆盖。**
    - 分组信息：字典。提供样本分组信息，您 **只能** 将样本分为 A、B 两组，在每组 key - value 中，key 为您在 2 中准备的原始数据文件的样本名称，不可重复；value 为样本分组，仅可为 A 或 B。
    - 使用的基因集路径：字符串。提供基因集路径，本仓库已经内置 hallmark2024 基因集于 data 目录下，若您不改变使用基因集，则不需要更改。
5. 运行 01.py，等待结果输出
6. 输出案例
    - 假设任务名称为 example，sample1 设置为 A，sample2 设置为 B。
    - 仓库根目录生成
      - example_top10_positive_nes_table.png
      - example_bottom10_negative_nes_table.png
      - example_gene_ranked_list.txt
      - example_gsea_output 文件夹
## 输出结果
  - 所有在分组信息中配置为 **A** 组的样本将成为 **对照组**。
  - 以输出案例的假设任务为例
    - example_top10_positive_nes_table.png 将显示 **相对于 sample1 (A)** sample2 (B) 的前十基因富集信息，NES 为正数。
    - example_bottom10_negative_nes_table.png 将显示 **相对于 sample2 (B)** sample1 (A) 的前十基因富集信息，NES 为负数。
    - example_gene_ranked_list.txt 将显示 **相对于 sample1 (A)** sample2 (B) 的差异表达分析 log2FC 结果。
    - example_gsea_output 文件夹将包含其他 **相对于 sample1 (A)** sample2 (B) 的 GSEA 分析结果。
- **无论您试图获取何种信息，请时刻注意对照组与实验组对象。**
## 免责声明
  - 您知悉且接受，由于差异表达分析方法不同或其他因素，本仓库代码运行结果可能与 GSEA 应用输出**不完全一致**。
  - 您知悉且接受，代码对原始数据进行了清洗，但效果**可能不佳**。
  - 您知悉且接受，该仓库地址将被合并在创建者实验报告中**一并提交**。
  - 本仓库旨在分享生物信息学操作流程，**不提倡**抄袭等行为，**不提供**答案。
## 其他
- 如果您的数据有生物重复，把它们都放进 .csv 然后分组分在同一组就行。
- 数据清洗理论上应该是有效的，但是 GSEApy 还是报 warning，原因未知。
## 致谢
- [GSEApy](https://gseapy.readthedocs.io/en/latest/introduction.html)
- [阿里云开发者社区-利用 limma 包进行差异表达分析](https://developer.aliyun.com/article/1316090)
- [ChatGPT](https://www.chatgpt.com)
- 如果您希望了解更多 GSEApy 以及差异表达分析（r 语言）方法，可以访问以上链接。

创建者 [@SiriNatsume](https://github.com/SiriNatsume)
祝你愉快 :)