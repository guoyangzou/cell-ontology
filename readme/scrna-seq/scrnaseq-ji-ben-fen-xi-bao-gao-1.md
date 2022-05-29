# scRNA-seq基本分析报告

## 摘要 <a href="#a9" id="a9"></a>

完成2个样本的单细胞测序，共计得到13449个细胞,平均每个样本原始数据量为90.60 G。

| 项目信息        | 详情              |
| ----------- | --------------- |
| project\_id | scRNA-seq基本分析报告 |
| samples     | M1,M2           |
| groups      | M1\_vs\_M2      |

## 1 背景介绍 <a href="#a10" id="a10"></a>

单细胞转录组测序就是利用高通量测序技术对单个细胞进行RNA-seq，以此了解单细胞水平的基因表达情况。某个细胞在某一生理功能状态下所有转录的mRNA产物的集合，是基因组遗传信息传递和表达的重要步骤和过程。单细胞转录组最大的优势就在于，它可以明确到细胞间的差异性，可以将某个组织中不同的单细胞分别进行研究，获取更多的样本信息，而不同于对整个组织（bulk）进行探究的普通转录组，对于理解单个细胞的多样性有着明显的局限性。每种组织中有多种类型的细胞，每一种细胞类型有着独特的起源和功能，而每一个细胞的谱系和发展的状态又决定了每个细胞如何与周围的细胞及环境互作，那么把高通量测序应用到单个细胞层面，对于我们理解细胞的起源、功能、变异等就有着至关重要的作用。

利用单细胞转录组测序，可以对细胞类型、细胞亚型分类、标记基因进行鉴定，从而实现对细胞群体的划分与细胞群体间基因表达差异的检测，此外该技术还可以预测细胞分化与发育轨迹，在单细胞水平研究细胞相互联系和分子机制，已广泛应用肿瘤异质性、疾病机制、免疫反应、耐药评估以及组织器官生长发育等研究领域，极大拓展了生命科学研究的广度和深度。

## 2 项目流程介绍 <a href="#a13" id="a13"></a>

### **2.1 实验流程**

![图1 单细胞转录组实验流程图](<../../.gitbook/assets/image (1) (1) (1) (1).png>)

#### 1）单细胞悬液的制备和质检

根据组织或细胞样品的来源和特性，选取适合的方法将其在较短时间内制备成单细胞悬液并进行细胞质控，严格把控单细胞悬液的浓度和质量，使用合格的样品进行10X单细胞转录组测序。细胞悬液质量要求：细胞总量＞50万个、细胞活性>85%、细胞浓度700-1200cell/μl、细胞体积>100μl、结团率<15%、细胞直径<40μm，少碎片，不存在逆转录抑制剂和非细胞的核酸分子。

#### 2）单细胞文库构建

在微流控芯片中，将单个细胞、反应所需试剂，与带有细胞标签序列（cellBarcode）的凝胶珠（Gel bead）一起包裹在油滴中，生成 GEMs；在 GEM油滴中，细胞裂解释放出 RNA，在合适的条件下，RNA 与带有 cell Barcode，UMI 的Poly( dT)引物结合，进行互补链延伸，并在延伸链末端加上 3 个 C碱基；CCC 与 TSO的 rGrGrG 互补配对后，以 TSO为模板延伸，完成逆转录反应；随后打破GEMs，回收并通过PCR扩增富集cDNA，进行 cDNA的文库构建。

#### 3）cDNA 回收、扩增和质检

随后打破 GEMs，利用磁珠回收 cDNA；通过 PCR 扩增全长cDNA，并进行纯化。利用qubit 4.0 对 cDNA 浓度进行质控，利用 Aglient 2100 对cDNA 的完整性进行质控。

#### 4）上机测序

选择适量的 cDNA 进行文库构建，文库构建基本包括片段化、末端修复和加 A尾，接头连接和 Index 扩增等基本步骤，利用 qubit 4.0对文库浓度进行质控，Qseq400对文库的片段进行质控；库检合格后，使用二代测序平台进行测序。

### **2.2 分析流程**

![图2 单细胞转录组生信分析流程图](../../.gitbook/assets/企业微信截图\_16536460416984.png)

## 3 数据质量评估及表达定量 <a href="#a26" id="a26"></a>

### **3.1 测序数据质量评估**

高通量测序下机得到的原始图像文件经 CASAVA 碱基识别转化为测序读段（Sequenced Reads），以 FASTQ 格式存储。FASTQ是一种存储生物序列及相应质量值的常用文本格式，格式如下。10xGenomics测序数据每个样本的数据包含I1，R1，R2。I1存储了index信息；R1即read1，28bp为细胞 barcode 和 UMI信息。R2即read2。使用fastqc软件对每个样本的read2数据做质控分析。

![图3 FASTQ格式文件示意图](<../../.gitbook/assets/fastq (1).png>)

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

### **3.2 测序数据统计**

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

### **3.3 比对分析**

Cell Ranger调用STAR\[2]软件将Read2比对到参考基因组上，基于STAR的比对结果，结合参考数据集（gtf／gff文件）里的信息，统计基因组上各个区域的reads覆盖信息，可以得到比对到外显子、内含子、基因间区的比例信息，作为数据质控的参考指标。将reads既比对到已知转录本的外显子上又在同一条链上的作为比对到转录本上的依据，如果该reads比对到已知的单个基因上，将reads称为唯一比对到转录组上，只有比对到转录本上的reads才能作为UMI计数。STAR是一款RNA-Seq数据分析常用的分段比对工具，可以用来发现外显子的连接以及融合现象，其基本工作原理主要分成两步：种子序列的寻找，以及聚类／连接／打分。下图为其原理示意图：

![图4 STAR 比对原理](<../../.gitbook/assets/image (5).png>)

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

### **3.4 细胞鉴定及基因表达定量**

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

## 4 细胞亚群分析 <a href="#a52" id="a52"></a>

利用Cell Ranger对测序数据进行比对及定量，后续利用Seurat\[3]软件做进一步分析，主要包括细胞过滤、数据标准化处理、批次效应校正\[4]、PCA（Principal Component Analysis）分析、t-SNE（t-Distributed Stochastic Neighbor Embedding）分析\[5]、UMAP（Uniform Manifold Approximation and Projection）分析\[6]、细胞聚类、marker gene分析、细胞注释等。

### **4.1 细胞过滤**

理想情况下中的细胞只有1个，但实验过程中也存在空细胞，或者存在2个甚至多个细胞的情况；并且当细胞发生死亡或者裂解时，细胞中将含有大量的线粒体基因，因此需要根据细胞的UMI总数、线粒体基因比例及单个细胞鉴定到的基因数目3个指标进行高质量细胞筛选，本项目过滤参数设置为：

(1) 线粒体基因比例(%)：>= 20，一般认为细胞发生凋亡、破裂会导致线粒体基因偏高，因此需要过滤掉线粒体基因比例过高的细胞；

(2) 单个细胞鉴定到的UMI数目：≥ 100，对单个细胞检出的转录本数目进行过滤；

(3) 单个细胞鉴定到的基因数目：500 \~ 7000，正常情况下一个细胞表达的基因数目处于一定范围内，因此利用基因数目可以判断GEMs油滴中是否包含1个以上细胞、或者细胞是否发生破裂。

同时也会对检出不可信的基因进行过滤，某个基因表达至少在10个细胞中检测到才会被认为是可信检出。

统计过滤前后各个样本细胞的nGene（number of Gene），nUMI（number of UMI）和percent.mt（线粒体基因）占比，结果如下图所示：

![图5 过滤前各样本细胞的小提琴图](<../../.gitbook/assets/image (3) (1) (1).png>)

_注：每个点表示一个细胞。nGene图：过滤前各个样本单个细胞中检测到的基因数量分布情况；nUMI图：过滤前各个样本单个细胞中检测到的nUMI分布情况；percent.mt图：过滤前各个样本单个细胞中线粒体基因表达量的百分比分布情况。_

![图6 过滤后各样本细胞的小提琴图](<../../.gitbook/assets/image (8) (1) (1).png>)

_注：每个点表示一个细胞。nGene图：过滤后各个样本单个细胞中检测到的基因数量分布情况；nUMI图：过滤后各个样本单个细胞中检测到的nUMI分布情况；percent.mt图：过滤后各个样本单个细胞中线粒体基因表达量的百分比分布情况。_

细胞中基因表达数目与UMI数目呈现正相关，细胞中检测到的UMI数目越多，表达的基因数目越多。线粒体基因与UMI数目没有明显的相关性。过滤前后相关性分析散点图如下所示：

![图7 过滤前各样本细胞的散点图](<../../.gitbook/assets/image (2).png>)

_注：每个点表示一个细胞，图上方为相关性系数。左图表示细胞中检测到的UMI数目和线粒体基因含量的关系；右图表示细胞中检测到的UMI数目和表达基因数目的关系。_

![图8 过滤后各样本细胞的散点图](<../../.gitbook/assets/image (14) (1).png>)

_注：每个点表示一个细胞，图上方为相关性系数。左图表示细胞中检测到的UMI数目和线粒体基因含量的关系；右图表示细胞中检测到的UMI数目和表达基因数目的关系。_

**表8 细胞过滤统计表**

| Sample | All\_cells\_number | filter by Mitochondria | filter by UMI | filter by minimal gene expression | filter by maximum gene expression | Final\_cells\_number | percent(%) |
| ------ | ------------------ | ---------------------- | ------------- | --------------------------------- | --------------------------------- | -------------------- | ---------- |
| M1     | 6,882              | 1                      | 0             | 1,190                             | 0                                 | 5,691                | 82.69      |
| M2     | 6,567              | 0                      | 0             | 1,137                             | 2                                 | 5,428                | 82.66      |

_注：Sample:样本名；_\
_All\_cells\_number：过滤前的所有细胞数目；_\
_filter by Mitochondria：根据细胞中线粒体最大比例过滤掉的细胞数目；_\
_filter by UMI：根据细胞中最少的UMI数过滤掉的细胞数目；_\
_filter by minimal gene expression：根据细胞中最少的Gene数过滤掉的细胞数；_\
_filter by maximum gene expression：根据细胞中最多的Gene数过滤掉的细胞数；_\
_Final\_cells\_number：过滤后最终剩余的细胞数；_\
_percent(%)：剩余细胞的百分比(100%)。_

细胞过滤结果文件路径：BMK\_3\_seurat\_analysis/BMK\_1\_CellsFilter/

细胞过滤结果文件下载链接

### **4.2 细胞聚类分析**

多个样本需要对数据合并后再进行后续分析，可以选择是否去除批次效应，来处理样品在不同批次中产生的与生物复杂性无关的差异。利用seurat软件对所有样本数据进行批次效应校正，校正原理是将CCA（anonical Correspondence Analysis）与MNN（Mutual Nearest Neighbors）算法结合起来，认为相同类型和状态的细胞它们之间的基因表达差异是技术偏倚引起的，通过计算并校正技术偏倚引起的基因表达差值，从而实现不同单细胞数据集的整合。

单细胞转录组定量表达数据是一个M\*N的矩阵（行是基因，列是细胞），对这样的矩阵聚类计算量极大。因此在对细胞进行聚类之前，Seurat使用PCA方法对数据进行降维，PCA降维是一种线性降维方法，运用方差分解，将高维的数据映射到低维的空间中；然后基于SNN聚类算法对细胞进行聚类和分群，构建细胞间的聚类关系；最后将降维后的数据传递到t-SNE与UMAP进行可视化展示，细胞之间的基因表达模式越相似，在t-SNE/UMAP图中的距离也越接近。

单样本细胞聚类及后续分析结果文件路径：BMK\_3\_seurat\_analysis/BMK\_2\_SingleAnalysis

单样本细胞聚类及后续分析结果文件下载链接

#### **4.2.1 多样本的Cluster分布统计**

基于多样本的聚类结果,统计不同cluster中各样本的占比情况。如下图所示：

![图9 多样本的Cluster分布柱状图](<../../.gitbook/assets/image (3) (1).png>)

_注：左图横坐标表示相同cluster里面不同样本的细胞数百分比，不同颜色表示不同的样本；右图横坐标表示所有样本在每个cluster里面的细胞总个数分布情况。纵坐标表示不同的cluster。_

**表9 多样本的Cluster统计表**

| cluster   | M1    | M2    | total  |
| --------- | ----- | ----- | ------ |
| cluster0  | 1,449 | 1,476 | 2,925  |
| cluster1  | 1,686 | 1,148 | 2,834  |
| cluster2  | 538   | 609   | 1,147  |
| cluster3  | 198   | 805   | 1,003  |
| cluster4  | 645   | 169   | 814    |
| cluster5  | 500   | 312   | 812    |
| cluster6  | 183   | 239   | 422    |
| cluster7  | 43    | 355   | 398    |
| cluster8  | 248   | 57    | 305    |
| cluster9  | 127   | 61    | 188    |
| cluster10 | 53    | 131   | 184    |
| cluster11 | 21    | 66    | 87     |
| total     | 5,691 | 5,428 | 11,119 |

#### **4.2.2 聚类结果可视化**

将聚类得到的结果使用t-SNE和UMAP进行可视化展示。

![图10 所有样本细胞聚类的t-SNE图](<../../.gitbook/assets/image (17) (1).png>)

_注：tsne\_all.png：所有cluster的t-SNE图，不同颜色表示不同的cluster； tsne\_sample.png：所有样本的t-SNE图，不同颜色表示不同的样本； tsne\_splt.png：不同的cluster在不同样本中的情况。_

多样本的Cluster分布统计结果文件路径：BMK\_3\_seurat\_analysis/BMK\_3\_Integrated/BMK\_3\_Cluster

多样本的Cluster分布统计结果文件路径

### **4.3 细胞亚群注释**

细胞聚类分析是根据细胞之间的相似性，将相似度最高的一群细胞识别为一个亚群，但是得到的细胞亚群并没有生物学意义。细胞注释主要有两种方法，1)基于参考数据库利用软件工具进行自动化；2)基于参考数据库、文献资料确定已知细胞类型的特征表达基因(marker标记)列表和相关通路进行人工注释。SingleR\[7]是基于斯皮尔曼相关性，对scRNA-seq数据实现自动化注释的一个软件，百迈客使用SingleR对细胞类群进行自动注释（只针对人、小鼠）。

注：细胞类型鉴定的过程即复杂又繁琐，自动注释只能提供一个参考，细胞注释要同时结合自动化注释和人工注释，这两种方法并不是独立割裂开来的，需要搭配使用，软件自动化细胞注释方法方便且系统化，完全依赖于参考数据集，注释的结果有时也并不是高置信度注释。在百迈客提供的数据库自动注释的基础上，后续可根据已发表的文章（物种、疾病、组织检索）或者cell marker数据库进行进一步的人工注释校正，可大大提高注释的准确度。

细胞类型注释统计表文件路径：BMK\_3\_seurat\_analysis/BMK\_3\_Integrated/BMK\_6\_cellAnnotation/All.cell\_annotation\_stat.xls

**表10 样本细胞类型注释统计表**

| Cell Type         | M1            | M2            | Total         |
| ----------------- | ------------- | ------------- | ------------- |
| Astrocytes        | 500(8.79%)    | 312(5.75%)    | 812(7.3%)     |
| Endothelial cells | 127(2.23%)    | 61(1.12%)     | 188(1.69%)    |
| Epithelial cells  | 53(0.93%)     | 131(2.41%)    | 184(1.65%)    |
| Fibroblasts       | 893(15.69%)   | 226(4.16%)    | 1,119(10.06%) |
| Neurons           | 1,987(34.91%) | 2,085(38.41%) | 4,072(36.62%) |
| Oligodendrocytes  | 2,131(37.45%) | 2,613(48.14%) | 4,744(42.67%) |
| Total             | 5,691         | 5,428         | 11,119        |

_注：Cell Type：注释得到的细胞类型；_\
_第二列至倒数第二列：该样本中属该细胞类型的细胞数目，括号内为细胞类型占比；_\
_Total：总的细胞数，括号内为细胞类型占比。_

All.cell\_annotation\_stat.xls.html

![图12 细胞类型分布统计图](<../../.gitbook/assets/image (7) (1) (1).png>)

_注：横坐标代表不同样本，纵坐标表示对应的细胞数目百分比，不同颜色表示不同的细胞类型_

**表11 样本细胞类型Cluster统计表**

| Cell Type         | Cluster                          |
| ----------------- | -------------------------------- |
| Astrocytes        | Cluster5(100.0%)                 |
| Endothelial cells | Cluster9(100.0%)                 |
| Epithelial cells  | Cluster10(100.0%)                |
| Fibroblasts       | Cluster4(72.7%); Cluster8(27.3%) |
| Neurons           | Cluster0(71.8%); Cluster2(28.2%) |

_注：Cell Type：注释得到的细胞类型；_\
_Cluster：Cluster细胞数占该细胞类型细胞数的比例。_

All.cluster\_annotation\_result.xls.html

![图13 细胞类型鉴定TNSE/UMAP图](<../../.gitbook/assets/image (20) (1).png>)

多样本整合后的细胞亚群的类型鉴定结果文件路径：BMK\_3\_seurat\_analysis/BMK\_3\_Integrated/BMK\_\*\_cellAnnotation/

多样本整合后的细胞亚群的类型鉴定结果文件下载链接

多样本整合后的细胞亚群的类型鉴定结果文件路径（其它注释数据库）：BMK\_3\_seurat\_analysis/OtherAnnoDatabase/

多样本整合后的细胞亚群的类型鉴定结果文件下载链接（其它注释数据库）

### **4.4 marker基因分析**

#### **4.4.1 marker基因筛选**

基于降维聚类得到的结果，使用非参数检验方法（wilcox秩和检验），鉴定在每个cluster中特异表达的基因，将Fold Change≥ 2 且Thred < 0.1（parameter is: FDR）作为筛选标准，筛选得到每个cluster最显著的差异基因，即marker基因。

细胞亚群间的差异分析结果文件路径：BMK\_3\_seurat\_analysis/BMK\_3\_Integrated/BMK\_4\_MarkerGene

**表12 差异基因数目统计**

| cluster     | cluster0 | cluster1 | cluster2 | cluster3 | cluster4 | cluster5 | cluster6 | cluster7 | cluster8 | cluster9 | cluster10 | cluster11 |
| ----------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | --------- | --------- |
| deg\_number | 41       | 31       | 21       | 26       | 31       | 33       | 15       | 203      | 28       | 24       | 17        | 35        |

**表13 cluster的Marker基因**

| ID                 | symbol | Pvalue | log2FC | pct.1 | pct.2 | Qvalue | clsuterName | cluster0\_count | cluster1\_count | cluster2\_count | cluster3\_count | cluster4\_count | cluster5\_count | cluster6\_count | cluster7\_count | cluster8\_count | cluster9\_count | cluster10\_count | cluster11\_count |
| ------------------ | ------ | ------ | ------ | ----- | ----- | ------ | ----------- | --------------- | --------------- | --------------- | --------------- | --------------- | --------------- | --------------- | --------------- | --------------- | --------------- | ---------------- | ---------------- |
| ENSMUSG00000021268 | Meg3   | 0      | 2.78   | 0.99  | 0.42  | 0      | cluster0    | 31.17           | 0.64            | 16.27           | 2.89            | 1.89            | 1.92            | 2.36            | 1.21            | 2.42            | 3.88            | 3.93             | 1.47             |
| ENSMUSG00000044349 | Snhg11 | 0      | 2.66   | 0.99  | 0.42  | 0      | cluster0    | 36.26           | 0.74            | 23.19           | 3.57            | 2.35            | 2.1             | 2.68            | 1.29            | 2.82            | 4.64            | 3.4              | 1.32             |
| ENSMUSG00000025326 | Ube3a  | 0      | 1.69   | 0.89  | 0.36  | 0      | cluster0    | 4.45            | 0.32            | 1.98            | 0.73            | 0.55            | 0.52            | 0.62            | 0.37            | 0.54            | 0.73            | 0.65             | 0.47             |
| ENSMUSG00000019986 | Ahi1   | 0      | 1.65   | 0.92  | 0.37  | 0      | cluster0    | 4.91            | 0.33            | 2.99            | 0.67            | 0.52            | 0.63            | 1.07            | 0.53            | 0.58            | 0.77            | 1.24             | 0.78             |
| ENSMUSG00000026576 | Atp1b1 | 0      | 1.57   | 0.76  | 0.23  | 0      | cluster0    | 3.55            | 0.18            | 1.68            | 0.31            | 0.29            | 0.46            | 0.32            | 0.3             | 1.84            | 0.55            | 0.31             | 0.3              |

_注:ID:基因ID；_\
_symbol:基因symbol；_\
_Pvalue:显著性p值；_\
_log2FC:差异倍数的log2值；_\
_pct.1:基因在cluster(i) 中有表达的细胞比例；_\
_pct.2:基因在除了cluster(i)中以外所有的cluster中有表达的细胞比例；_\
_Qvalue:校正后的p值；_\
_clusterName：差异基因的clustercluster；_\
_N\_count：差异基因在clusterN的单个细胞reads数平均值。_

All.cluster0.diff\_featuregene.xls.html

通过火山图（Volcano Plot）可查看基因在细胞亚群间的表达水平的差异，以及差异的统计学显著性。各个cluster差异基因火山图展示如下：

![图14 差异表达基因火山图](<../../.gitbook/assets/image (6) (1) (1).png>)

_注：差异表达火山图中的每一个点表示一个基因，横坐标表示某一个基因在cluster中表达量差异倍数的对数值；纵坐标表示错误发现率的负对数值。横坐标绝对值越大，说明表达量在两样品间的表达量倍数差异越大；图中灰色的点代表无差异表达的基因，红色的点代表上调基因，蓝色的点代表下调的基因。_

每个cluster差异基因（top10）的表达变化热图展示。

![图15 各cluster差异基因（top10）的表达变化热图](<../../.gitbook/assets/image (9) (1) (1).png>)

_注：红色表示高表达，蓝色表示低表达，顶部注释条表示cluster名。由于cluster之间的top10的差异基因间可能存在重复，因此实际基因数目可能偏小。热图右边特定选取每个cluster的top2的基因进行标出。_

每个cluster差异基因（top2）的气泡图展示。

![图16 各cluster差异基因（top2）的气泡图](<../../.gitbook/assets/image (18) (1).png>)

_注：颜色越深表示该基因平均表达值越高，点越大，表达比例越大。由于cluster之间的top2的差异基因间可能存在重复，因此实际基因数目可能偏小。_

从各个cluster的差异基因中挑选出差异倍数排名前10的基因进行展示。

![图17 Top10 Marker基因的小提琴图](<../../.gitbook/assets/image (13) (1) (1).png>)

_注：图中纵坐标代表细胞聚类的 cluster ，横坐标代表基因的表达值。通过该图我们可以很直观的看出 marker 基因在不同的单细胞亚群中的表达量高低分布情况_。

![图18 Top10 Marker基因的t-SNE图](<../../.gitbook/assets/image (11) (1).png>)

_注: top10 marker 基因所有细胞类型中表达值在t-SNE聚类结果中的可视化，每张图中紫色标记的细胞即为特异表达该 marker 基因的细胞类型。_

![图19 Top10 Marker基因的UMAP图](<../../.gitbook/assets/image (2) (1) (1).png>)

_注: top10 marker 基因所有细胞类型中表达值在UMAP聚类结果中的可视化，每张图中紫色标记的细胞即为特异表达该 marker 基因的细胞。_

top marker基因展示结果文件路径：BMK\_3\_seurat\_analysis/BMK\_3\_Integrated/BMK\_4\_MarkerGene/BMK\_2\_top10\_marker

top marker基因展示结果文件下载链接

#### **4.4.2 marker基因注释**

将得到的差异基因与NR\[8]，Swiss-Prot\[9]，GO\[10]，COG\[11]，KOG\[12]，Pfam\[13]，KEGG\[14]，Reactome\[15]数据库进行序列比对，得到基因的注释信息。

**表14 差异表达基因注释文件**

| ID                 | symbol      | Pvalue | log2FC | pct.1 | pct.2 | Qvalue | clsuterName | cluster0\_count | cluster1\_count | cluster2\_count | cluster3\_count | cluster4\_count | cluster5\_count | cluster6\_count | cluster7\_count | cluster8\_count | cluster9\_count | cluster10\_count | cluster11\_count | COG\_class | COG\_class\_annotation                 | GO\_annotation                                                                                                                                                                                                                                                                                                                                       | KEGG\_annotation                                                                                                                                                                                        | KEGG\_pathway\_annotation                                                                                                                                                                                                        | KOG\_class | KOG\_class\_annotation                 | Pfam\_annotation                                                                                                                                                                                                                                                       | Swiss-Prot\_annotation                                                                     | eggNOG\_class | eggNOG\_class\_annotation                                    | NR\_annotation                                                                        |
| ------------------ | ----------- | ------ | ------ | ----- | ----- | ------ | ----------- | --------------- | --------------- | --------------- | --------------- | --------------- | --------------- | --------------- | --------------- | --------------- | --------------- | ---------------- | ---------------- | ---------- | -------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------- | -------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------ | ------------- | ------------------------------------------------------------ | ------------------------------------------------------------------------------------- |
| ENSMUSG00000097545 | Mir124a-1hg | 0      | 1.16   | 0.66  | 0.12  | 0      | cluster0    | 1.77            | 0.04            | 1.14            | 0.16            | 0.11            | 0.12            | 0.27            | 0.03            | 0.18            | 0.21            | 0.13             | 0.07             | --         | --                                     | --                                                                                                                                                                                                                                                                                                                                                   | --                                                                                                                                                                                                      | --                                                                                                                                                                                                                               | --         | --                                     | --                                                                                                                                                                                                                                                                     | --                                                                                         | --            | --                                                           | --                                                                                    |
| ENSMUSG00000033061 | Resp18      | 0      | 1.24   | 0.61  | 0.16  | 0      | cluster0    | 2.17            | 0.08            | 1.52            | 0.21            | 0.17            | 0.22            | 0.23            | 0.14            | 0.25            | 0.46            | 0.16             | 0.14             | --         | --                                     | Biological Process: in utero embryonic development (GO:0001701);; Cellular Component: extracellular region (GO:0005576);; Cellular Component: Golgi apparatus (GO:0005794);; Cellular Component: secretory granule (GO:0030141);; Cellular Component: perikaryon (GO:0043204);; Cellular Component: rough endoplasmic reticulum lumen (GO:0048237);; | --                                                                                                                                                                                                      | --                                                                                                                                                                                                                               | --         | --                                     | RESP18 domain                                                                                                                                                                                                                                                          | Regulated endocrine-specific protein 18 OS=Mus musculus OX=10090 GN=Resp18 PE=2 SV=1       | O             | Posttranslational modification, protein turnover, chaperones | regulated endocrine-specific protein 18 precursor \[Mus musculus]                     |
| ENSMUSG00000097451 | Rian        | 0      | 1.5    | 0.85  | 0.2   | 0      | cluster0    | 3.01            | 0.08            | 1.84            | 0.27            | 0.28            | 0.21            | 0.29            | 0.15            | 0.23            | 0.45            | 0.41             | 0.21             | --         | --                                     | --                                                                                                                                                                                                                                                                                                                                                   | --                                                                                                                                                                                                      | --                                                                                                                                                                                                                               | --         | --                                     | --                                                                                                                                                                                                                                                                     | --                                                                                         | --            | --                                                           | --                                                                                    |
| ENSMUSG00000030302 | Atp2b2      | 0      | 1.07   | 0.62  | 0.13  | 0      | cluster0    | 1.52            | 0.06            | 0.63            | 0.13            | 0.09            | 0.33            | 0.3             | 0.08            | 0.08            | 0.26            | 0.13             | 0.05             | \[P]       | Inorganic ion transport and metabolism | --                                                                                                                                                                                                                                                                                                                                                   | K05850\|0\|mmu:11941\|K05850 Ca2+ transporting ATPase, plasma membrane \[EC:3.6.3.8] \| (RefSeq) Atp2b2, D6Abb2e, Gena300, PMCA2, Tmy, dfw, jog, wms, wri; ATPase, Ca++ transporting, plasma membrane 2 | Calcium signaling pathway (ko04020);; cGMP-PKG signaling pathway (ko04022);; cAMP signaling pathway (ko04024);; Adrenergic signaling in cardiomyocytes (ko04261);; Salivary secretion (ko04970);; Pancreatic secretion (ko04972) | \[P]       | Inorganic ion transport and metabolism | Cation transporting ATPase, C-terminus;; E1-E2 ATPase;; Plasma membrane calcium transporter ATPase C terminal;; Cation transport ATPase (P-type);; haloacid dehalogenase-like hydrolase;; Cation transporter/ATPase, N-terminus;; haloacid dehalogenase-like hydrolase | Plasma membrane calcium-transporting ATPase 2 OS=Mus musculus OX=10090 GN=Atp2b2 PE=1 SV=2 | P             | Inorganic ion transport and metabolism                       | ATPase, Ca++ transporting, plasma membrane 2, isoform CRA\_b, partial \[Mus musculus] |
| ENSMUSG00000115783 | Bc1         | 0      | 1.25   | 0.72  | 0.22  | 0      | cluster0    | 2.32            | 0.16            | 1.3             | 0.28            | 0.24            | 0.32            | 0.27            | 0.56            | 0.25            | 0.37            | 0.23             | 0.21             | --         | --                                     | --                                                                                                                                                                                                                                                                                                                                                   | --                                                                                                                                                                                                      | --                                                                                                                                                                                                                               | --         | --                                     | --                                                                                                                                                                                                                                                                     | --                                                                                         | --            | --                                                           | --                                                                                    |

_注:网页版只展示前6行。_\
_ID:基因ID；_\
_symbol:基因symbol；_\
_Pvalue:显著性p值；_\
_log2FC:差异倍数的log2值；_\
_pct.1:基因在cluster(i) 中有表达的细胞比例；_\
_pct.2:基因在除了cluster(i)中以外所有的cluster中有表达的细胞比例；_\
_Qvalue:校正后的p值；_\
_clusterName：差异基因的cluster；_\
_clusterN\_count：差异基因在clusterN的单个细胞reads数平均值；_\
_其余列:COG，GO，KEGG，KOG，Pfam，Swiss-Prot，eggNOG，NR数据库对应的注释信息。_

All.cluster0.diff\_featuregene.annotation.xls.html

差异表达基因注释结果文件路径：BMK\_3\_seurat\_analysis/BMK\_3\_Integrated/BMK\_4\_MarkerGene/BMK\_3\_Anno\_enrichment/All.cluster\*/BMK\_1\_Anno/

功能富集分析及基因注释的结果文件下载链接(展示一个样本)

#### **4.4.3 GO功能富集分析**

GO数据库是GO组织（Gene Ontology Consortium）于2000年构建的一个结构化的标准生物学注释系统，旨在建立基因及其产物知识的标准词汇体系，适用于各个物种。GO注释系统是一个有向无环图，包含三个主要分支，即：生物学过程（Biological Process），分子功能（Molecular Function）和细胞组分（Cellular Component）。对每个cluster的差异基因集，采用ClusterProfiler对基因分别进行生物学过程，分子功能和细胞组分的富集分析。富集分析采用超几何检验方法来寻找与整个基因组背景相比显著富集的GO条目。对富集结果得到的Term采用绘制柱状图气泡图等进行可视化。

![图20 差异表达基因富集柱状图](<../../.gitbook/assets/image (10) (1).png>)

_注：图中横坐标为对应的GO term,纵坐标为-log10(pvalue)。每个柱子上的数字表示富集到该term的基因数。不同的颜色分别代表GO的三个本体：BP、CC、MF_

![图21 差异表达基因富集气泡图](<../../.gitbook/assets/image (15) (1).png>)

_注：图中每一个圆表示一个term，横坐标表示term名称，纵坐标为富集因子（Enrichment Factor），表示差异基因中注释到某term的基因比例与所有基因中注释到该term的基因比例的比值。富集因子越大，表示差异表达基因在该term中的富集水平越显著。圆圈的颜色代表pvalue，pvalue越小，表示差异表达基因在该term中的富集显著性越可靠；圆圈的大小表示term中富集的基因数目，圆圈越大，表示基因越多。_

![图22 差异表达基因富集网络图](<../../.gitbook/assets/image (1) (1).png>)

_注：差异表达基因与GO term的网络图.边的颜色代表不同的term,基因节点的颜色代表差异倍数,term节点越大说明富集到该term的基因数目越多。_

对每个cluster的差异基因进行富集分析，富集到的Term做topGO有向无环图。topGO有向无环图能直观展示差异表达基因富集的GO节点（Term）及其层级关系，是差异表达基因GO富集分析的结果图形化展示，分支代表包含关系，从上至下所定义>的功能描述范围越来越具体。差异表达基因的topGO有向无环图如下:

![图23 GO富集有向无环图](<../../.gitbook/assets/image (4) (1) (1).png>)

_注：对每个GO节点进行富集，最显著的10个节点在图中用方框表示，图中还包含其各层对应关系。每个方框（或椭圆）内给出了该GO节点的内容描述和富集显著性值。不同颜色代表不同的富集显著性，颜色越深，显著性越高。_

功能富集分析之GO分析结果文件路径：BMK\_3\_seurat\_analysis/BMK\_3\_Integrated/BMK\_4\_MarkerGene/BMK\_3\_Anno\_enrichment/All.cluster\*/BMK\_2\_GO\_enrichment/

#### **4.4.4 KEGG功能注释及富集分析**

在生物体内，不同的基因产物相互协调来行使生物学功能，对差异表达基因的通路（Pathway）注释分析有助于进一步解读基因的功能。KEGG（Kyoto Encyclopedia of Genes and Genomes）是系统分析基因功能、基因组信息数据库，它有助于研究者把基因及表达信息作为一个整体网络进行研究。作为有关Pathway的主要公共数据库(Kanehisa,2008），KEGG提供的整合代谢途径(pathway)查询，包括碳水化合物、核苷、氨基酸等的代谢及有机物的生物降解，不仅提供了所有可能的代谢途径，而且对催化各步反应的酶进行了全面的注解，包含有氨基酸序列、PDB库的链接等等，是进行生物体内代谢分析、代谢网络研究的强有力工具。

对差异表达基因KEGG的注释结果按照KEGG中通路类型进行分类，分类图如下图所示：

![图24 差异表达基因KEGG分类图](<../../.gitbook/assets/image (8) (1).png>)

_注：横坐标为注释到该通路下的基因个数及其个数占被注释上的基因总数的比例，纵坐标为KEGG代谢通路的名称。_

差异表达基因的通路注释结果见下图：

![图25 差异表达基因的KEGG通路注释图](<../../.gitbook/assets/image (6) (1).png>)

_注：红色（绿色）框标记的酶与每个cluster差异高（低）表达基因有关，框内的数字代表酶的编号（EC number），而整个通路由多种酶催化的复杂生化反应构成，此通路图中与差异表达基因相关的酶均用颜色标出，根据研究对象间的差异，重点研究某些代谢通路相关基因的差异表达情况，通过通路解释表型差异的根源。_

分析差异表达基因在某一通路上是否发生显著差异（over-presentation）即为差异表达基因的通路富集分析。Pathway显著性富集分析以KEGG数据库中Pathway为单位，应用超几何检验，找出与整个基因组背景相比，在差异表达基因中显著性富集的Pathway。差异表达基因KEGG通路富集分析结果见下图。

![图26 差异表达基因KEGG富集气泡图](<../../.gitbook/assets/image (9) (1).png>)

_注：图中每一个圆表示一个KEGG通路，横坐标表示通路名称，纵坐标为富集因子（Enrichment Factor），表示差异基因中注释到某通路的基因比例与所有基因中注释到该通路的基因比例的比值。富集因子越大，表示差异表达基因在该通路中的富集水平越显著。圆圈的颜色代表pvalue，pvalue越小，表示差异表达基因在该通路中的富集显著性越可靠；圆圈的大小表示通路中富集的基因数目，圆圈越大，表示基因越多。_

![图27 差异表达基因KEGG富集网络图](<../../.gitbook/assets/image (7) (1).png>)

_注：差异表达基因与KEGG通路的网络图.边的颜色代表不同的通路,基因节点的颜色代表差异倍数,通路节点越大说明富集到该通路的基因数目越多。_

功能富集分析之KEGG分析结果文件路径：BMK\_3\_seurat\_analysis/BMK\_3\_Integrated/BMK\_4\_MarkerGene/BMK\_3\_Anno\_enrichment/All.cluster\*/BMK\_3\_KEGG\_enrichment/

#### **4.4.5 Reactome功能富集分析**

Reactome\[15]是一个免费的、开源的信号和代谢分子的关系数据库。Reactome数据库搜集了人类相关的反应和生物学通路（包含 13,827 个人类反应，分为 2,536 条通路，涉及 11,374 种蛋白质）。典型的生物学通路包括：中间代谢、信号传导、转录调控、细胞凋亡和疾病。富集结果如下表所示。

**表15 Reactome通路富集表**

| ReactomeID    | Description                     | GeneRatio | BgRatio  | pvalue               | p.adjust             | qvalue               | geneID                                                                                                                                                                     | geneSymbol                                                | Count |
| ------------- | ------------------------------- | --------- | -------- | -------------------- | -------------------- | -------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------- | ----- |
| R-MMU-382551  | Transport of small molecules    | 9/17      | 653/8880 | 8.46778888659635e-07 | 8.04439944226653e-05 | 5.61548105111126e-05 | ENSMUSG00000002985/ENSMUSG00000026576/ENSMUSG00000030302/ENSMUSG00000027523/ENSMUSG00000032554/ENSMUSG00000017740/ENSMUSG00000007097/ENSMUSG00000033379/ENSMUSG00000040907 | Apoe/Atp1b1/Atp2b2/Gnas/Trf/Slc12a5/Atp1a2/Atp6v0b/Atp1a3 | 9     |
| R-MMU-5578775 | Ion homeostasis                 | 4/17      | 52/8880  | 2.35230031365502e-06 | 8.04756801828705e-05 | 5.6176929102724e-05  | ENSMUSG00000026576/ENSMUSG00000030302/ENSMUSG00000007097/ENSMUSG00000040907                                                                                                | Atp1b1/Atp2b2/Atp1a2/Atp1a3                               | 4     |
| R-MMU-936837  | Ion transport by P-type ATPases | 4/17      | 53/8880  | 2.54133726893275e-06 | 8.04756801828705e-05 | 5.6176929102724e-05  | ENSMUSG00000026576/ENSMUSG00000030302/ENSMUSG00000007097/ENSMUSG00000040907                                                                                                | Atp1b1/Atp2b2/Atp1a2/Atp1a3                               | 4     |
| R-MMU-983712  | Ion channel transport           | 5/17      | 171/8880 | 1.281224051422e-05   | 0.0                  | 0.0                  | ENSMUSG00000026576/ENSMUSG00000030302/ENSMUSG00000007097/ENSMUSG00000033379/ENSMUSG00000040907                                                                             | Atp1b1/Atp2b2/Atp1a2/Atp6v0b/Atp1a3                       | 5     |
| R-MMU-5576891 | Cardiac conduction              | 4/17      | 114/8880 | 5.39013677517281e-05 | 0.0                  | 0.0                  | ENSMUSG00000026576/ENSMUSG00000030302/ENSMUSG00000007097/ENSMUSG00000040907                                                                                                | Atp1b1/Atp2b2/Atp1a2/Atp1a3                               | 4     |

_注：ReactomeID：Reactome编号；_\
_Description：Reactome编号对应的功能描述；_\
_GeneRatio：注释到Reactome编号上的cluster差异高表达基因数与cluster差异高表达基因总数的比值；_\
_BgRatio：注释到Reactome编号上的背景基因数与背景基因总数的比值；_\
_pvalue：显著性检验p值；_\
_padj：校正后的p值；_\
_geneID：注释到Reactome编号上的基因；_\
_geneSymbol：注释到Reactome通路编号上的差异基因名称；_\
_Count：注释到Reactome编号上的基因数。_

All.cluster0\_reactome\_enrich.list.html

对Reactome富集结果得到的通路进行可视化展示，选取最显著的20个通路（如不足20个，则绘制所有）绘制柱状图和气泡图。

![图28 Reactome富集柱状图](<../../.gitbook/assets/image (11).png>)

![图29 Reactome富集气泡图](<../../.gitbook/assets/image (2) (1).png>)

_注：图中每一个圆表示一个Reactome通路，横坐标表示通路名称，纵坐标为GeneRatio。圆圈的颜色代表pvalue，pvalue越小，表示差异表达基因在该通路中的富集显著性越可靠；圆圈的大小表示通路中富集的基因数目，圆圈越大，表示基因越多。_

Reactome分析结果文件路径：BMK\_3\_seurat\_analysis/BMK\_3\_Integrated/BMK\_4\_MarkerGene/BMK\_3\_Anno\_enrichment/All.cluster\*/BMK\_4\_Reactome\_enrichment

#### **4.4.6 差异蛋白调控网络分析**

该分析的目的是研究蛋白-蛋白互作。STRING是收录多个物种预测的和实验验证蛋白质-蛋白质互作的数据库，包括直接的物理互作和间接的功能相关。结合差异表达分析结果和数据库收录的互作关系对，构建差异表达基因互作网络。对于数据库中包含的物种，可直接从数据库中提取出目标基因集的互作关系对构建互作网络；对于数据库中未收录信息的物种，使用BLAST软件，将目的基因与数据库中的蛋白质进行序列比对，寻找同源蛋白，根据同源蛋白的互作关系对构建互作网络。构建完成的蛋白质互作网络可导入Cytoscape软件进行可视化。Cytoscape可视化的差异表达基因蛋白质互作网络如下图：

![图30 差异表达基因蛋白质互作网络图](<../../.gitbook/assets/image (12) (1).png>)

_注：图中的节点为蛋白质，边为互作关系。互作网络中节点(node)的大小与此节点的度(degree)成正比，即与此节点相连的边越多，它的度越大，节点也就越大。节点的颜色代表对应的基因表达量的上下调的趋势，红色代表上调，蓝色代表下调。边(edge)的宽度表示此边连接的两个节点间的互相作用的关系强弱，互相作用的关系越强，边越宽，边的颜色则代表了不同相互作用的类型。没有的组合代表没有互作关系。_

差异表达基因蛋白互作网络结果文件路径：BMK\_3\_seurat\_analysis/BMK\_3\_Integrated/BMK\_4\_MarkerGene/BMK\_4\_PPI/

差异表达基因蛋白互作网络结果文件下载链接（展示一个样本）

#### **4.4.7 转录因子结合位点分析**

转录因子结合位点（Transcription factor binding site，TFBS）是与转录因子结合的DNA片段，长度通常在 5\~20 bp范围内，一个转录因子往往同时调控若干个基因，而它在不同基因上的结合位点具有一定的保守性，又不完全相同。我们用R包TFBStools对差异基因的启动子区域（定义基因上游１kb 左右为潜在的启动子区）上的 TFBS 进行了预测，参考的转录因子motif数据库是JASPAR数据库（http://jaspar.genereg.net/）。预测结果如下表：

**表16 转录因子结合位点预测结果**

| Model\_id  | seqname            | start | end | score | strand | Attributes                                                               | Pvalue |
| ---------- | ------------------ | ----- | --- | ----- | ------ | ------------------------------------------------------------------------ | ------ |
| MA0006.1.1 | ENSMUSG00000000740 | 746   | 751 | 1     | +      | TF=Ahr::Arnt;class=Basic helix-loop-helix factors (bHLH);sequence=TGCGTG | 0      |
| MA0006.1.2 | ENSMUSG00000000740 | 973   | 978 | 0.99  | +      | TF=Ahr::Arnt;class=Basic helix-loop-helix factors (bHLH);sequence=CGCGTG | 0.0    |
| MA0035.1.1 | ENSMUSG00000000740 | 93    | 98  | 0.97  | +      | TF=Gata1;class=Other C4 zinc finger-type factors;sequence=TGATGC         | 0.0    |
| MA0035.1.2 | ENSMUSG00000000740 | 242   | 247 | 0.93  | +      | TF=Gata1;class=Other C4 zinc finger-type factors;sequence=AGATGA         | 0.01   |
| MA0035.1.3 | ENSMUSG00000000740 | 266   | 271 | 0.91  | +      | TF=Gata1;class=Other C4 zinc finger-type factors;sequence=TGATAT         | 0.01   |

_注：Model\_id：转录因子结合位点模体id；_\
_seqname：基因名称；_\
_start：开始位置；_\
_end：结束位置；_\
_score：评分，评分越高，表示该转录因子与输入序列结合的可能性越大；_\
_strand：链方向, + 表示正链，- 表示负链；_\
_Attributes：注释信息。TF：转录因子id；class：转录因子注释；sequence：TFBS序列；_\
_Pvalue：P值。_

allGenes\_TFBS\_predictRes.xls.html

TFBS序列特征示例图如下：

![图31 TFBS序列特征图](<../../.gitbook/assets/image (13) (1).png>)

注：横坐标为 motif 中碱基的相对位置，纵坐标为该位置碱基的序列保守性，而碱基信号的高度代表该碱基在该位置上的相对频率。

转录因子结合位点预测结果文件路径：BMK\_3\_seurat\_analysis/BMK\_3\_Integrated/BMK\_4\_MarkerGene/BMK\_5\_TF\_analysis/TFBS\_Analysis/

转录因子结合位点预测结果文件下载链接（展示一个样本）

## 5 组间差异表达分析 <a href="#a162" id="a162"></a>

为了比较同一种细胞亚型在不同样品之间的差异情况，对同一cluster中的不同分组进行组间差异分析。差异分析结果如下表所示：

### **5.1 差异表达基因分析**

相同Cluster不同样品间的差异分析结果文件路径：BMK\_3\_seurat\_analysis/BMK\_3\_Integrated/BMK\_5\_Group\_Anlysis

**表17 差异结果统计表**

| Groups | deg\_group | cluster0 | cluster1 | cluster2 | cluster3 | cluster4 | cluster5 | cluster6 | cluster7 | cluster8 | cluster9 | cluster10 | cluster11 |
| ------ | ---------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | --------- | --------- |
| group1 | M1\_vs\_M2 | 1        | 5        | 2        | 12       | 4        | 11       | 2        | 103      | 4        | 2        | 4         | 1         |

利用热图展示组间TOP10的差异基因在分组中的表达分布模式。

![图32 top10的差异基因在各样本的表达热图](<../../.gitbook/assets/image (4) (1).png>)

_注：将数据按照各亚群各个分组分开，图片展示了某差异分组(A\_vs\_B)在所有亚群中的差异情况。行表示基因（每个亚群选取top10），列表示各细胞亚群。_

#### **5.1.1 top差异基因展示**

利用小提琴图、t-SNE图和UMAP图展示组间TOP2的差异基因在比较组中的表达分布情况。

![图33 Top2 Marker基因的小提琴图](<../../.gitbook/assets/image (10).png>)

注：展示cluster的top2的 marker 基因在不同样本间的表达分布。不同样本用不同颜色表示，横轴标出样本名，纵轴表示基因在对应cluster中的表达。































