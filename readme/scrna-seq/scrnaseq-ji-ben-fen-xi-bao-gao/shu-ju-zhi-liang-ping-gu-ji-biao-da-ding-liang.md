# 数据质量评估及表达定量

## **测序数据质量评估**

高通量测序下机得到的原始图像文件经 CASAVA 碱基识别转化为测序读段（Sequenced Reads），以 FASTQ 格式存储。FASTQ是一种存储生物序列及相应质量值的常用文本格式，格式如下。10xGenomics测序数据每个样本的数据包含I1，R1，R2。I1存储了index信息；R1即read1，28bp为细胞 barcode 和 UMI信息。R2即read2。使用fastqc软件对每个样本的read2数据做质控分析。

![图3 FASTQ格式文件示意图](<../../../.gitbook/assets/fastq (1).png>)

该项目各样品数据产出统计见下表：

**表2 原始数据产出统计**

| ReadSum     | BaseSum        | GC(%) | N(%) | Q20(%) | Q30(%) | sampleID |
| ----------- | -------------- | ----- | ---- | ------ | ------ | -------- |
| 300,000,000 | 90,600,000,000 | 46.99 | 0.0  | 96.23  | 91.01  | M1       |
| 300,000,000 | 90,600,000,000 | 46.68 | 0.0  | 95.71  | 89.81  | M2       |

_注：SampleID：样本名；_

_ReadSum：Raw Data中pair-end Reads总数；_&#x20;

_BaseSum：Raw Data总碱基数；_&#x20;

_GC(%):Raw Data GC含量，即Raw Data中G和C两种碱基占总碱基的百分比；_

_N(%):Raw Data中N占比； Q20(%):Raw Data质量值大于或等于20的碱基所占的百分比；_&#x20;

_Q30(%):Raw Data质量值大于或等于30的碱基所占的百分比。_

**表3 样本信息对照表**

| #ID | read1\_name                          | read2\_name                          |
| --- | ------------------------------------ | ------------------------------------ |
| M1  | Mouse\_demo-SC-020001\_good\_1.fq.gz | Mouse\_demo-SC-020001\_good\_2.fq.gz |
| M2  | Mouse\_demo-SC-020002\_good\_1.fq.gz | Mouse\_demo-SC-020002\_good\_2.fq.gz |

_注：#ID：样本名称；_\
_read1\_name、read2\_name：样本名称对应的双端fq文件名；_\
_样本数量较多时可能会显示不全，建议直接查看BMK\_1\_rawdata/data\_cfg.xls文件；_

测序数据及其质量评估结果文件路径：BMK\_1\_rawData/

测序数据及其质量评估结果文件下载链接

## **测序数据统计**

利用Cell Ranger\[1]进一步对测序数据进行统计分析，CellRanger分析结果网页版报告链接如下:

CellRanger分析结果网页版报告链接：

1. M1.web\_summary.html
2. M2.web\_summary.html

该项目各样品测序数据产出统计见下表：

表4 CellRanger分析序列统计

| sampleID | Number.of.Reads | Valid.Barcodes | Sequencing.Saturation | Q30.Bases.in.Barcode | Q30.Bases.in.RNA.Read | Q30.Bases.in.UMI |
| -------- | --------------- | -------------- | --------------------- | -------------------- | --------------------- | ---------------- |
| M1       | 300,000,000     | 97.4%          | 85.1%                 | 96.1%                | 91.0%                 | 92.5%            |
| M2       | 300,000,000     | 97.5%          | 83.7%                 | 96.1%                | 89.8%                 | 92.2%            |

_注：sampleID：样本ID；_\
_Number of Reads：reads总数；_\
_Valid Barcodes：有效的10X Barcode比例 ；_\
_Sequencing Saturation：测序饱和度；_\
_Q30 Bases in Barcode：Barcode序列中质量值大于或等于30的碱基所占的百分比；_\
_Q30 Bases in RNA Read：reads中质量值大于或等于30的碱基所占的百分比；_\
_Q30 Bases in UMI：UMI序列中质量值大于或等于30的碱基所占的百分比。_

## **比对分析**

Cell Ranger调用STAR\[2]软件将Read2比对到参考基因组上，基于STAR的比对结果，结合参考数据集（gtf／gff文件）里的信息，统计基因组上各个区域的reads覆盖信息，可以得到比对到外显子、内含子、基因间区的比例信息，作为数据质控的参考指标。将reads既比对到已知转录本的外显子上又在同一条链上的作为比对到转录本上的依据，如果该reads比对到已知的单个基因上，将reads称为唯一比对到转录组上，只有比对到转录本上的reads才能作为UMI计数。STAR是一款RNA-Seq数据分析常用的分段比对工具，可以用来发现外显子的连接以及融合现象，其基本工作原理主要分成两步：种子序列的寻找，以及聚类／连接／打分。下图为其原理示意图：

![图4 STAR 比对原理](<../../../.gitbook/assets/image (5).png>)

_注：最大可比对标签（Maximum Mappable Prefix, MMP）的寻找用于发现（a）外显子连接处，（b）错配，（c）多聚 A 尾， 或者接头，或者低质量尾。_

CellRanger分析比对结果统计如下表：

**表5 CellRanger分析比对结果统计**

| sampleID | Reads Mapped to Genome | Reads Mapped Confidently to Genome | Reads Mapped Confidently to Intergenic Regions | Reads Mapped Confidently to Intronic Regions | Reads Mapped Confidently to Exonic Regions | Reads Mapped Antisense to Gene | Reads Mapped Confidently to Transcriptome | Fraction Reads in Cells |
| -------- | ---------------------- | ---------------------------------- | ---------------------------------------------- | -------------------------------------------- | ------------------------------------------ | ------------------------------ | ----------------------------------------- | ----------------------- |
| M1       | 93.7%                  | 91.6%                              | 3.8%                                           | 40.5%                                        | 47.3%                                      | 2.8%                           | 41.6%                                     | 65.5%                   |
| M2       | 92.7%                  | 90.7%                              | 4.0%                                           | 42.8%                                        | 44.0%                                      | 2.7%                           | 38.4%                                     | 62.8%                   |

_注：sampleID：样本ID；_\
_Reads Mapped to Genomes：比对到参考基因组上的Reads在总Reads中占的百分比；_\
_Reads Mapped Confidently to Genome：比对到参考基因组并得到转录本GTF信息支持的Reads在总Reads中占的百分比；_\
_Reads Mapped Confidently to Intergenic Regions：比对到基因间区域的Reads在总Reads中占的百分比；_\
_Reads Mapped Confidently to Intronic Regions：比对到内含子区域的Reads在总Reads中占的百分比；_\
_Reads Mapped Confidently to Exonic Regions：比对到外显子区域的Reads在总Reads中占的百分比；_\
_Reads Mapped Antisense to Gene：比对到基因反义链的Reads在总Reads中占的百分比；_\
_Reads Mapped Confidently to Transcriptome：比对到已知参考转录本的Reads在总Reads中占的百分比；_\
_Fraction Reads in Cells：比对到参考基因且来源于高质量细胞的Reads在总Reads中占的百分比。_

## **细胞鉴定及基因表达定量**

Cell Ranger对每个Barcode的每个基因会去除重复的UMI，统计unique UMI数目作为该细胞内该基因的表达量，通过UMI可以区分一条read是否属于生物学重复还是技术重复，能够有效地去除PCR效应。CellRanger分析细胞统计如下表：

**表6 CellRanger分析细胞信息统计**

| sampleID | Estimated Number of Cells | Mean Reads per Cell | Median UMI Counts per Cell | Median Genes per Cell | Total Genes Detected |
| -------- | ------------------------- | ------------------- | -------------------------- | --------------------- | -------------------- |
| M1       | 6,882                     | 43,592              | 1,292                      | 792                   | 25,556               |
| M2       | 6,567                     | 45,683              | 1,294                      | 740                   | 26,118               |

_注：sampleID：样本ID；_\
_Estimated Number of Cells：检测到的细胞数目；_\
_Mean Reads per Cell：每个细胞平均Reads数目 ；_\
_Median UMI Counts per Cell：每个细胞的UMI中位数；_\
_Median Genes per Cell：每个细胞中基因数目的中位数；_\
_Total Genes Detected：所有细胞的基因总数。_

细胞基因表达量统计结果文件路径：BMK\_2\_cellranger\_analysis/BMK\_1\_summary/\*.raw\_expMat.xls

数据统计及细胞定量的结果文件下载链接

Cell Ranger其它分析结果文件路径：BMK\_2\_cellranger\_analysis
