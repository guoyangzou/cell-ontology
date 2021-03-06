# 2 项目流程介绍

## 2.1 **实验流程** <a href="#2.1-shi-yan-liu-cheng" id="2.1-shi-yan-liu-cheng"></a>

![图1 单细胞转录组实验流程图](https://3672080242-files.gitbook.io/\~/files/v0/b/gitbook-x-prod.appspot.com/o/spaces%2F182kdjhopQDf4Z51a1hd%2Fuploads%2FgNERyekzyw6tfVFnQFsX%2Fimage.png?alt=media\&token=4a02a895-4d1d-4695-8ec3-26b94a503100)

#### 1）单细胞悬液的制备和质检 <a href="#1-dan-xi-bao-xuan-ye-de-zhi-bei-he-zhi-jian" id="1-dan-xi-bao-xuan-ye-de-zhi-bei-he-zhi-jian"></a>

根据组织或细胞样品的来源和特性，选取适合的方法将其在较短时间内制备成单细胞悬液并进行细胞质控，严格把控单细胞悬液的浓度和质量，使用合格的样品进行10X单细胞转录组测序。细胞悬液质量要求：细胞总量＞50万个、细胞活性>85%、细胞浓度700-1200cell/μl、细胞体积>100μl、结团率<15%、细胞直径<40μm，少碎片，不存在逆转录抑制剂和非细胞的核酸分子。

#### 2）单细胞文库构建 <a href="#2-dan-xi-bao-wen-ku-gou-jian" id="2-dan-xi-bao-wen-ku-gou-jian"></a>

在微流控芯片中，将单个细胞、反应所需试剂，与带有细胞标签序列（cellBarcode）的凝胶珠（Gel bead）一起包裹在油滴中，生成 GEMs；在 GEM油滴中，细胞裂解释放出 RNA，在合适的条件下，RNA 与带有 cell Barcode，UMI 的Poly( dT)引物结合，进行互补链延伸，并在延伸链末端加上 3 个 C碱基；CCC 与 TSO的 rGrGrG 互补配对后，以 TSO为模板延伸，完成逆转录反应；随后打破GEMs，回收并通过PCR扩增富集cDNA，进行 cDNA的文库构建。

#### 3）cDNA 回收、扩增和质检 <a href="#3cdna-hui-shou-kuo-zeng-he-zhi-jian" id="3cdna-hui-shou-kuo-zeng-he-zhi-jian"></a>

随后打破 GEMs，利用磁珠回收 cDNA；通过 PCR 扩增全长cDNA，并进行纯化。利用qubit 4.0 对 cDNA 浓度进行质控，利用 Aglient 2100 对cDNA 的完整性进行质控。

#### 4）上机测序 <a href="#4-shang-ji-ce-xu" id="4-shang-ji-ce-xu"></a>

选择适量的 cDNA 进行文库构建，文库构建基本包括片段化、末端修复和加 A尾，接头连接和 Index 扩增等基本步骤，利用 qubit 4.0对文库浓度进行质控，Qseq400对文库的片段进行质控；库检合格后，使用二代测序平台进行测序。

## **2.2 分析流程** <a href="#2.2-fen-xi-liu-cheng" id="2.2-fen-xi-liu-cheng"></a>

![单细胞转录组生信分析流程图](../../../.gitbook/assets/企业微信截图\_16536460416984.png)
