# 组间差异表达分析

为了比较同一种细胞亚型在不同样品之间的差异情况，对同一cluster中的不同分组进行组间差异分析。差异分析结果如下表所示：

## **差异表达基因分析**

相同Cluster不同样品间的差异分析结果文件路径：BMK\_3\_seurat\_analysis/BMK\_3\_Integrated/BMK\_5\_Group\_Anlysis

**表17 差异结果统计表**

| Groups | deg\_group | cluster0 | cluster1 | cluster2 | cluster3 | cluster4 | cluster5 | cluster6 | cluster7 | cluster8 | cluster9 | cluster10 | cluster11 |
| ------ | ---------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | --------- | --------- |
| group1 | M1\_vs\_M2 | 1        | 5        | 2        | 12       | 4        | 11       | 2        | 103      | 4        | 2        | 4         | 1         |

利用热图展示组间TOP10的差异基因在分组中的表达分布模式。

![图32 top10的差异基因在各样本的表达热图](<../../../.gitbook/assets/image (4).png>)

_注：将数据按照各亚群各个分组分开，图片展示了某差异分组(A\_vs\_B)在所有亚群中的差异情况。行表示基因（每个亚群选取top10），列表示各细胞亚群。_

### **top差异基因展示**

利用小提琴图、t-SNE图和UMAP图展示组间TOP2的差异基因在比较组中的表达分布情况。

![图33 Top2 Marker基因的小提琴图](<../../../.gitbook/assets/image (10).png>)

注：展示cluster的top2的 marker 基因在不同样本间的表达分布。不同样本用不同颜色表示，横轴标出样本名，纵轴表示基因在对应cluster中的表达。

