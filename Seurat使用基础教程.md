# Seurat使用基础教程

标签（空格分隔）： 单细胞教程

---

[官方原本教程][1]
# 一、安装软件
通过`Rstudio / r-base`默认的`CRAN`平台安装`Seurat`

```R
# Enter commands in R (or R studio, if installed)
install.packages('BiocManager')
BiocManager::install('multtest')
install.packages('Seurat')
library(Seurat)
``` 

如果出现一下警告信息，请输入`y`:

```report
package which is only available in source form, and may need compilation of C/C++/Fortran: 'Seurat'
Do you want to attempt to install these from sources?
y/n:

```
### **注意事项：**

 - 默认安装的都是`Seurat`最新版本，当前最新版本为`Seurat(Version:3.1.4)`，如果`R`版本过低有可能会安装失败或者安装低版本的`Seurat`。
 - 最常用的`Seurat`版本有两种一个是`Seurat 3.0`以后的版本，另一个是`Seurat 2.3.4`，两个版本之间所使用的算法以及代码参数设置都有极大的不同，使用时请注意。
 - 查看`Seurat`版本代码：`sessionInfo(package = "Seurat")`

# 二、下载训练数据

我们从`10X Genomics`官网上获取免费的单细胞`PBMC`测序数据集，我们将`Illumina NextSeq 500`上对`2700`个`PBMC`单细胞进行了测序。

[点击此处下载训练数据][2]

```bash
wget https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
```


# 三、读数原始数据转换成Seurat对象

```R
# 首先载入相关包
library(dplyr) 
library(Seurat)
library(patchwork)
# 使用Read10X函数读取CellRanger输出产生的matrix，并转换为稀疏矩阵。
pbmc.data <- Read10X(data.dir = "../data/pbmc3k/filtered_gene_bc_matrices/hg19/")
# 使用CreateSeuratObject函数将稀疏矩阵转换为Seurat 对象，方便后续的操作。
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc
```

**`Seurat`的所有操作全部建立在`Seurat`对象的基础上，需要熟悉掌握`Seurat`对象的格式。**

```
## An object of class Seurat 
## 13714 features across 2700 samples within 1 assay 
## Active assay: RNA (13714 features)
```

# 四、检测数据并初步过滤
Seurat首先检测线粒体基因在各个细胞中的占比，然后过滤掉捕获率极低的Gene、低质量的细胞、Doublet细胞以及线粒体异常高的细胞。

```R
# 计算线每个细胞的线粒体基因占比，并且保存在pbmc$percent.mt中
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-") 
```

![image_1e5rr7kp61dhhaph9petd51v0q9.png-20.5kB][3]

```
# 可视化数据三个主要质量指标：每个细胞的基因数、每个细胞的UMI数以及每个细胞的线粒体占比
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# 可视化 每个细胞的UMI数与每个细胞线粒体基因占比的关系
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
# 可视化 每个细胞的基因数与每个细胞UMI数之间的关系
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
```

![image_1e5rr8jv1dtc5p57pnug0t4lm.png-121.8kB][4]

```
# 使用subet函数，来过滤基因数过少/多，且线粒体占比过多的的细胞。
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
```
**在这个样本中，我们根据可视化后的数据质量报告：**
 - 过滤掉基因>2500以及<200的细胞
 - 过滤掉线粒体基因占比 >5% 的细胞

# 五、数据标准化
默认使用全局缩放归一化方法：LogNormalize，再将其乘以比例因子(默认为10,000)，存储在pbmc[['RNA']]@data。
```R
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
```

# 六、筛选高变量基因

默认情况下，我们根据建模，计算单细胞中最突出的生物信号，返回2,000个高变异Gene，用于下游分析，包括PCA等。

```R
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(pbmc), 10) # 筛选前10个高变基因
# 将高变基因HVG进行可视化。
plot1 <- VariableFeaturePlot(pbmc) 
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) # 额外标记前10个HVG。
plot1 + plot2 #图层叠加
```
![image_1e5rspfu42rbau61thna0q16n913.png-68.3kB][5]

# 七、二次标准化

 - 数据标准化，是指不同样本之间的数据标准化。为了能够是不同样本数据之间进行比较。
 - 二次标准化，只是同一样本内的不同细胞进行标准化，从而能够实现PCA的特征选择

```R
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes) # 可以对所有基因进行二次标准化
pbmc <- ScaleData(pbmc) # 默认对HVG进行二次标准化
pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt") # 二次标准化时，可以对线粒体基因等指标进行拟合，从而实现对某些特征值的矫正拟合。
```

 - 拟合可以消除批次效应，但是效果并不好，会出现矫正过度的现象。


# 八、降维

然后我们对二次标准化后的数据进行PCA，默认以HVG作为特征输入。

```R
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc),verbose=F)
```

可以通过以下多种方式查看PCA各个基因的结果。

**常规PCA结果输出：**
```R
# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5) # 打印输出PCA各个基因的结果
```
![image_1e5sd5k6h20i126u1j0913fi1a4b1g.png-12.1kB][6]

**散点图展示不同PC的重要程度：**
```R
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca") # 通过气泡图的方式展示PCA结果
```
![image_1e5sd7qms1i5r1ctv1baaaddpg21t.png-65.7kB][7]
**散点图展示每个细胞PC1和PC2坐标：**
```R
DimPlot(pbmc, reduction = "pca") # 展示每个细胞的PCA降维的PC1和PC2的坐标 
```
**热图展示PC1~15贡献最大的基因，并随机挑选500个细胞：**
```R
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE) # 展示PC1~15差异最大的基因，并选500个细胞，用热图展示。
```
![image_1e5sfp85l1o1418lnu5q1q5ff6m2a.png-63.5kB][8]

**确定合适的维数选择(此步骤可以跳过)：**

```R
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15) # 将重要的PC的P值进行可视化，通常情况下此步骤可以跳跃，选择不同数量的PC结果会有很大的不同，且建议使用更多的PC进行下游分析，如果选择过少，则会丢失较多的特征值，对分析产生不利影响。通常默认选择20个，可选择10~50个，结果通常产生不了太大变化。
ElbowPlot(pbmc)
```

# 九、聚类

**基于KNN算法，构建KNN图，并计算聚类结果：**
```R
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)
```
# 十、可视化
**Seurat提供了两种非线性降维技术：**`tSNE/UMAP`

```R
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
pbmc <- RunTSNE(pbmc, dims = 1:10) # 运行tSNE非线性降维方法执行可视化
pbmc <- RunUMAP(pbmc, dims = 1:10) # 运行UMAP非线性降维方法执行可视化
# 两种非线性降维方法执行一种即可UMAP效果更好，通常选择UMAP的方法可视化
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap", label=T) # 作图，默认使用UMAP的方法进行可视化。
```
# 十一、其他功能

## 1、查找差异表达的特征（集群生物标记）

```
Cluster1.markers <- FindMarkers(pbmc, ident.1=1, min.pct=0.25) # 计算Cluster1的Markers
head(Cluster1.markers, n=5)

```

12、版本说明

[训练数据][9]


  [1]: https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html
  [2]: https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
  [3]: http://static.zybuluo.com/czc/d9bsisme08j1egudoqtwmu6u/image_1e5rr7kp61dhhaph9petd51v0q9.png
  [4]: http://static.zybuluo.com/czc/2t0bhjcuo0tnjd0jetpinf1q/image_1e5rr8jv1dtc5p57pnug0t4lm.png
  [5]: http://static.zybuluo.com/czc/g6hpuv4xfbqnswgt60can7id/image_1e5rspfu42rbau61thna0q16n913.png
  [6]: http://static.zybuluo.com/czc/483d8v1plcflyshme4vpm2uy/image_1e5sd5k6h20i126u1j0913fi1a4b1g.png
  [7]: http://static.zybuluo.com/czc/jrt7h320wye6jclxewvji1v4/image_1e5sd7qms1i5r1ctv1baaaddpg21t.png
  [8]: http://static.zybuluo.com/czc/cnan8jc2y65u6rjvlk0mpda6/image_1e5sfp85l1o1418lnu5q1q5ff6m2a.png
  [9]: http://static.zybuluo.com/czc/d9bsisme08j1egudoqtwmu6u/image_1e5rr7kp61dhhaph9petd51v0q9.png