# 批次矫正方法之间的比较 Harmony 与Seurat

标签（空格分隔）： 单细胞教程

---


# 一、Harmony算法

## 1、载入并安装相关的包

```
library(Seurat)
library(dplyr)
library(glue)
library(ggplot2)
library(pheatmap)

library(cowplot)
library(harmony)
```

## 2、载入相关数据

这里总共有两个数据，2个都是PBMC但是数据情况应该差距比较大，一个大约有1w个细胞，另一个是10X support的数据，这里不都用10X提供的数据主要是因为10X support之间的数据批次效应比较小。

```
setwd("/home/czc/R_work/github/2_Harmony/")
load("./object.RData")
object_1<-object

load("./pbmc.RData")
object_2<-pbmc

```

![image_1e8k3quieqqffmri0e1mqe1qe39.png-8.9kB][1]

![image_1e8k3rasr1gng1f0u117i1hha38dm.png-8.7kB][2]

**两个数据的情况分别如下:**

![image_1e8k3t7kc154p8ln7c41k24ou413.png-119.7kB][3]

![image_1e8k3u5cvc2gac9n181moaegk1g.png-104.8kB][4]

## 3、数据整和

```
object<-merge(object_1,object_2,add.cell.ids = c("object_1","object_2"))
object@meta.data<-object@meta.data[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt")]

object$sample<-substr(rownames(object@meta.data), 1,8)

```

## 4、常规分析流程

```

pbmc<-object

pbmc <- subset(pbmc, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 20)
pbmc %>% NormalizeData() %>% FindVariableFeatures() %>% 
  ScaleData(assay="RNA")  %>% 
  RunPCA(npcs = 40,assay="RNA") ->  pbmc
pbmc %>% FindNeighbors(assay="RNA",dims = 1:40)%>%
  FindClusters(assay="RNA",resolution = 0.6)->  pbmc
pbmc <- RunUMAP(pbmc, dims = 1:40)
pbmc_raw<-pbmc
save(pbmc_raw,file="pbmc_raw.RData")
```

## 5、Harmony矫正+矫正后的常规流程

```
pbmc <- pbmc %>% 
  RunHarmony("sample", plot_convergence = TRUE)
pbmc <- pbmc %>% 
  RunUMAP(reduction = "harmony", dims = 1:40) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:40) %>% 
  FindClusters(resolution = 0.6) %>% 
  identity()

```

## 6、结果比较

矫正前
![image_1e8k4a9ir2o5avk17ev8in1dmg1t.png-90.2kB][5]
![image_1e8k4apchfrg1p0igstehl14qv2a.png-76.3kB][6]
矫正后
![image_1e8k4bk8k1824usu1h38f8vee92n.png-93.1kB][7]
![image_1e8k4c00djs3890m0e1cp05j734.png-88.6kB][8]


# 二、Seurat V3 自带的CCA+MNN算法

## 1、偷点懒，直接载入上面的中间数据


```
library(Seurat)
library(dplyr)
library(glue)
library(ggplot2)
library(pheatmap)

load("pbmc_raw.RData")
```

## 2、拆分数据并对拆分的数据进行标准化

```
ifnb.list <- SplitObject(pbmc_raw, split.by = "sample")
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
```

## 3、数据进行CCA + fastMNN矫正

```

immune.anchors <- FindIntegrationAnchors(object.list =ifnb.list, dims = 1:40)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:40)
```

## 4、常规流程，不过使用的是矫正后的Assay

```
DefaultAssay(immune.combined) <- "integrated"
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 40, verbose = FALSE)

immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:40)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:40)
immune.combined <- FindClusters(immune.combined, resolution = 0.6)
```

**展示结果**

```
DimPlot(immune.combined, reduction = "umap", group.by = "sample")
DimPlot(immune.combined, reduction = "umap", label = TRUE)
```

![image_1e8k59ra2cb6122c1b5d19el1bo63h.png-203.8kB][9]

![image_1e8k5afgvdhr1l3m1p0c1hg81j4r3u.png-187.5kB][10]



-------

**比较Harmony与Seurat自带的结果发现一件比较可怕的事情，就是Seurat会过度矫正，例如，object_1中并没有Macrophage细胞，Seurat依旧倾向于将一部分细胞认为类似于Macr，以至于靠的非常近。**

![image_1e8k5hcgg8n91j9a1bd5ee41o3b9.png-62.5kB][11]

![image_1e8k5huqgjlm179u16ib1ubrhdqm.png-67.3kB][12]

![image_1e8k5nssv110uj2kbmvcmijp113.png-113.2kB][13]
结果发现只是由于Atypical memory B 细胞与Macrophage都表达少量的SLC11A1，就被强行矫正在一起。主要是由于CCA的过度矫正问题，依旧无法被解决。

  [1]: http://static.zybuluo.com/czc/tewm905jen5gr9mkyhcd8zcz/image_1e8k3quieqqffmri0e1mqe1qe39.png
  [2]: http://static.zybuluo.com/czc/vehfe3wcbutxiaykoxz7j0bs/image_1e8k3rasr1gng1f0u117i1hha38dm.png
  [3]: http://static.zybuluo.com/czc/oue0srm4335ftml5v46rub2q/image_1e8k3t7kc154p8ln7c41k24ou413.png
  [4]: http://static.zybuluo.com/czc/xh1zb8zlwlic4qvqxid94gz6/image_1e8k3u5cvc2gac9n181moaegk1g.png
  [5]: http://static.zybuluo.com/czc/97ox8zoegye9owmw197ue7pu/image_1e8k4a9ir2o5avk17ev8in1dmg1t.png
  [6]: http://static.zybuluo.com/czc/8hkbwh6hrvfvr9a4zi4ipv6c/image_1e8k4apchfrg1p0igstehl14qv2a.png
  [7]: http://static.zybuluo.com/czc/940gaqnazx1un3g284f5cqrt/image_1e8k4bk8k1824usu1h38f8vee92n.png
  [8]: http://static.zybuluo.com/czc/op8vsiv8kwz6bykhpx7l7eu7/image_1e8k4c00djs3890m0e1cp05j734.png
  [9]: http://static.zybuluo.com/czc/owcyfalp08nbriixd09jfpia/image_1e8k59ra2cb6122c1b5d19el1bo63h.png
  [10]: http://static.zybuluo.com/czc/6ruhno52d7qic5fubtigqydl/image_1e8k5afgvdhr1l3m1p0c1hg81j4r3u.png
  [11]: http://static.zybuluo.com/czc/5e5l2k7rpzqtjf72amo2vxna/image_1e8k5hcgg8n91j9a1bd5ee41o3b9.png
  [12]: http://static.zybuluo.com/czc/rsjdprex296ipv4q8o89i60i/image_1e8k5huqgjlm179u16ib1ubrhdqm.png
  [13]: http://static.zybuluo.com/czc/l94y20xhrxre7gwq0tbyf420/image_1e8k5nssv110uj2kbmvcmijp113.png