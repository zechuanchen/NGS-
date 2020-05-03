# DoubletFinder识别单细胞数据中的Doublet细胞

标签（空格分隔）： 单细胞教程

---

# 一、文章以及基本原理

## 1、文献与软件作用

**文献链接：[文献链接][1]**
**大致内容**

 - DoubletFinder主要用来检测单细胞测序中的Doublet，也就是指双联体细胞（两个细胞包裹在同一个液滴中的情况）
 - DoubletFinder可以识别不同细胞造成的双峰（所谓双峰指，不同类型的细胞的特征基因表达造成的双峰）
 - DoubletFinder很难识别同类型细胞的双联体以及具有交叉表达谱的细胞。

## 2、算法原理：

 1. 图A：初次降维聚类，不清楚那些细胞是Doublet，那些细胞是Singlet
 2. 图B：根据初次降维聚类结果，从各个初次聚类得到的簇中提取出一部分细胞人工平均细胞基因表达，形成人工Doublet
 3. 图C：将人工得到的Doublet重新加入样本中，进行重聚类。
 4. 图D：人工Doublet与原本数据中的Doublet因为相似性而被聚集在一起，根据重聚类各个簇的人工Doublet含量进行统计
 5. 图E：根据混合高斯模型分离Singlet与Doublet。
 6. 图F：去除掉Doublet

![image_1e7aua60ok5hh8i1bve1dijn39.png-108.8kB][2]


# 二、使用方法（代码教程）

## 1、软件安装与数据载入

**软件安装：**

```
Dependencies<-c("Seurat","Matrix","fields","KernSmooth","ROCR","parallel") #所有的依赖包
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager") #判断是否已经安装过BiocManager，如果没有安装，自行安装
BiocManager::install(setdiff(Dependencies,rownames(installed.packages()))) #自动检测是否已经安装过相关依赖包，如果都按照会自行跳过，未被安装的进行安装
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder') #安装DoubletFinder
```

**数据载入：**
从[上一篇教程][3]的基础上进行后续分析，数据可以从此处下载
```
work_dir<-"/data/zeruo/czc_R_work/github/"
setwd(work_dir)
load("pbmc.RData")
object_ID<-"pbmc" #定义一下做的样本的类型
```

## 2、数据分析流程

**使用DoubletFinder对数据进行分析时最好了解一下参数：**
其中最主要的参数是pN和pK，pN可以默认为0.25，而pK则需要根据数据情况进行计算

    DoubletFinder takes the following arguments:

    seu ~ This is a fully-processed Seurat object (i.e., after NormalizeData, FindVariableGenes, ScaleData, RunPCA, and RunTSNE have all been run).

    PCs ~ The number of statistically-significant principal components, specified as a range (e.g., PCs = 1:10)

    pN ~ This defines the number of generated artificial doublets, expressed as a proportion of the merged real-artificial data. Default is set to 25%, based on observation that DoubletFinder performance is largely pN-invariant (see McGinnis, Murrow and Gartner 2019, Cell Systems).

    pK ~ This defines the PC neighborhood size used to compute pANN, expressed as a proportion of the merged real-artificial data. No default is set, as pK should be adjusted for each scRNA-seq dataset. Optimal pK values should be estimated using the strategy described below.

    nExp ~ This defines the pANN threshold used to make final doublet/singlet predictions. This value can best be estimated from cell loading densities into the 10X/Drop-Seq device, and adjusted according to the estimated proportion of homotypic doublets.

**pK Identification (no ground-truth)：**
首先将Seurat v3的对象转化为 DoubletFinder格式的对象，并计算合理的pK值。

```
sweep.data <- paramSweep_v3(pbmc,PCs=1:20)
sweep.stats <- summarizeSweep(sweep.data, GT = FALSE)
bcmvn= find.pK(sweep.stats)
```

![image_1e7avf0s3vhd11vf1pb7a0ldh3m.png-25.4kB][4]

**Homotypic Doublet Proportion Estimate：**
然后对结果进行预测，输入上面得到的合理的pK值，pK值取bcmvn为最高时候的值

```
homotypic.prop=modelHomotypic(pbmc@meta.data$RNA_snn_res.0.6)
nExp_poi=round(0.075*length(pbmc$orig.ident))
nExp_poi.adj=round(nExp_poi*(1-homotypic.prop))
pbmc=doubletFinder_v3(pbmc, PCs = 1:20, pN = 0.25, pK = as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)])),
                                              nExp = nExp_poi.adj, reuse.pANN = FALSE)
```

**在降维结果中对Doublet进行标明，然后进行剔除**

```
pbmc@meta.data$DF_hi.lo<- pbmc@meta.data[,8]   #注意，这里是按照之前教程的顺序没有添加其他列，所以Doublet排在第8列，如果在其他列，请更改相关代码
pdf(paste0(object_ID,"_DoubletFinder.pdf"),width = 10,height = 10)
P5<-DimPlot(pbmc, group.by ="DF_hi.lo",cols=c("black","gold","red"),reduction = "umap")
show(P5)
dev.off()
#绘图，结果用金色代表Singlet，黑色代表Doublet。添加红色是为了防止出现为添加相关标签的NA值准备，通常结果不会存在红色的点，如果存在说明结果有误。
Doublet<-table(pbmc@meta.data$DF_hi.lo=="Doublet")
pbmc@meta.data[,c(7,8)]=NULL   #同上，删除无关的列
Idents(pbmc ) <- "DF_hi.lo"    #激活相关ID
pbmc<-subset(x = pbmc, idents="Singlet")   #过滤细胞，只保留经过检测的Singlet

```
![image_1e7d69cnc1k6h1opu6441nps187j9.png-99.8kB][5]

**注意事项：**并非所有的Doublet都能够被识别，并非所有被识别的都是准确的，请了解该软件所适用的范围，切勿过度相信软件。

---------
参考文献：
1、McGinnis C S, Murrow L M, Gartner Z J. DoubletFinder: Doublet detection in single-cell RNA sequencing data using artificial nearest neighbors[J]. Cell systems, 2019, 8(4): 329-337. e4.



  [1]: https://www.sciencedirect.com/science/article/pii/S2405471219300730?via=ihub
  [2]: http://static.zybuluo.com/czc/aqeebgxfhueth30zr1uku4on/image_1e7aua60ok5hh8i1bve1dijn39.png
  [3]: https://github.com/zechuanchen/NGS-/blob/master/Seurat%E4%BD%BF%E7%94%A8%E5%9F%BA%E7%A1%80%E6%95%99%E7%A8%8B.md
  [4]: http://static.zybuluo.com/czc/dtuxyslsi23odfbmi47z3csx/image_1e7avf0s3vhd11vf1pb7a0ldh3m.png
  [5]: http://static.zybuluo.com/czc/pkflati7w5smbfskhfvu3un7/image_1e7d69cnc1k6h1opu6441nps187j9.png