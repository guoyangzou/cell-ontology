# 6 细胞周期分析

细胞周期（cell cycle）是指细胞从一次分裂完成开始到下一次分裂结束所经历的全过程，主要分G1、S、G2和M期。G1期是DNA合成准备期，S期是DNA合成期，G2期是分裂准备期，M期是分裂期。细胞周期调控是机体维持细胞增殖有序性及基因组DNA稳定性的关键，细胞周期失控导致各种疾病，在肿瘤发病中处于极其重要的中心环节，细胞周期的超常快速进行或在DNA复制不完全或有损伤时继续推进，将导致癌变。利用Seurat的AddModuleScore函数对单细胞转录组数据（单细胞表达矩阵）进行细胞周期分析，基于Seurat内置的周期特征蛋白，计算每个细胞内周期特征蛋白的转录表达水平，对每个细胞可能所处的周期状态进行评分，判断细胞是否处于增殖状态。

细胞周期鉴定结果文件路径：BMK\_5\_Cell\_Cycle

**表20 细胞周期鉴定表**

| Cell                  | sample | S.Score | G2M.Score | Phase |
| --------------------- | ------ | ------- | --------- | ----- |
| AAACCCAAGGTGGCTA-1\_1 | M1     | -0.01   | -0.02     | G1    |
| AAACCCACAACCAACT-1\_1 | M1     | -0.0    | 0.02      | G2M   |
| AAACCCACACGGGCTT-1\_1 | M1     | -0.01   | -0.01     | G1    |
| AAACCCACATAGAATG-1\_1 | M1     | -0.01   | 0.01      | G2M   |
| AAACCCAGTCGCCACA-1\_1 | M1     | -0.01   | -0.01     | G1    |

_注：Cell：barcode标记的细胞名。_\
_sample：样本名。_\
_S.Score：计算得到的S期得分。_\
_G2M.Score：计算得到的G2/M期得分_\
_Phase：细胞所处的周期_

All.cell\_cycle.xls.html

![图46 细胞周期TNSE/UMAP图](<../../../.gitbook/assets/image (6).png>)

_注：图中不同的颜色表示细胞所处的不同细胞周期_

_细胞周期鉴定结果文件下载链接：BMK\_5\_Cell\_Cycle/_

_细胞周期鉴定结果文件下载链接_
