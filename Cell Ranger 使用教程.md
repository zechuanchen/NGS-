# Cell Ranger 使用教程

标签（空格分隔）： 单细胞教程

---

# 一、CellRanger软件简介

## 1、Chromium 单细胞分析套件

 **Chromium 单细胞分析套件**：可以用于分析和可视化`10X Chromium`平台产生的`3'RNA-seq`和`Barcode` 数据。该套件包括：`Cell Ranger`和`Loupe Browser`。

 - [x] **Cell Ranger**:分析套件，用于样本数据分离，Barcode序列处理以及3'端单细胞数据的定量。
 - [x] **Loupe Browser**: 通过易于使用的可视化软件，在交互平台中，查找单细胞数据重要的基因、细胞类型以及亚型结构。


## 2、10X Genomic scRNA-seq 建库原理与常用术语
目前市场上`10X Genomic`大多数都是使用`V3`版本的试剂盒来建库
![image_1e5alnj061qfm8bmas01ehl135s9.png-44kB][1]
**相比于传统的RNA-seq建库，`scRNA-seq`为了做了两个改进：**

 1. 通过随机的Barcode序列来标记细胞
 2. 通过随机的UMI序列来标记每一条转录本

    这样我们就可以避免bulk RNA-seq中随机打断后无法区分是那个细胞的转录本的问题，而且10X Genomic单细胞测序时并非是随机打断，而是通过核酸酶，切割下3'端的序列，而且只需要检测ploy(A)尾的片段，因为有UMI系列进行标记，所以可以避免PCR造成的扩增偏差，去除掉重复UMI序列的Read，用nUMI代表表达量，即可以实现对于基因长度的标准化。

测序生成的文件：

 - **I1**：8bp Sample Barcode，用以区分样本来源。
 - **R1**：26bp read1，包括16 bp的Chromium Barcode以及10bp的 UMI序列，用以标记细胞和转录本来源
 - **R2**：98bp read2， 转录本序列，用以比对到基因组上。
 - **Barcode**：这里一般指Chromium Barcode，在Drop-seq的油包水中，一个细胞对应一种随机Barcode
 - **UMI**(Unique Molecular Identifier): 标记未扩增前的转录本序列，从而规避PCR扩增偏好，并消除随机打断造成的基因长度偏差。


**注意事项：**

 - [x] **10X Genimic使用Barcode来标记细胞，UMI来标记转录本。**
 - [x] **Cell Ranger版本只能高于试剂盒版本，否则会出现老版本Cell Ranger无法识别新版本试剂盒序列的异常bug,详细情况见下表**
![image_1e5bpi59b1t1o1bll6b4s893mq20.png-73.5kB][2]

## 3、Cell Ranger工作流程
```flow
st=>start: Sequencing
op=>operation: CellRanger mkfastq
op2=>operation: CellRanger count
cond=>condition: Sample combine or No?
e=>end: Expression matrix
op3=>operation: CellRanger aggr

st->op->op2->cond
cond(no)->e
cond(yes)->op3->e
```
 - Cellranger mkfastq：根据Illumina 测序生成的BCL文件，进行解复用，生成FASTQ文件。
 - Cellranger count：根据FASTQ文件，执行全基因组序列比对，实现定量并获得基因表达矩阵
 - Cellranger aggr：汇总多个样本的基因表达矩阵，从而实现多样本数据的整合

Cell Ranger的具体使用也是看具体情况来进行操作的，通常以下述3种情况比较多：

 - 1、一个样本一个流通池：最简单且常规的测序流程

![image_1e5bmb9sdalehif1mtf167n104h9.png-27.5kB][3]

 - 2、一个样本多个流通池：为了加深测序深度

![image_1e5bmbjqr11jq1c5e1qhgdt8o8bm.png-33.9kB][4]

- 3、多个样本一个流通池：多样本最主要的是可以减小批次效应

![image_1e5bmd8gao9dg5t16fhoua13cj13.png-49.4kB][5]



# 二、CellRanger安装
**CellRanger-3.1.0 (July 24, 2019)**
```bash
wget -O cellranger-3.1.0.tar.gz "http://cf.10xgenomics.com/releases/cell-exp/cellranger-3.1.0.tar.gz?Expires=1586356212&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cDovL2NmLjEweGdlbm9taWNzLmNvbS9yZWxlYXNlcy9jZWxsLWV4cC9jZWxscmFuZ2VyLTMuMS4wLnRhci5neiIsIkNvbmRpdGlvbiI6eyJEYXRlTGVzc1RoYW4iOnsiQVdTOkVwb2NoVGltZSI6MTU4NjM1NjIxMn19fV19&Signature=f-i3sfHa0ZRamefq1Uq~DEPKTz5LU6nIqGcrmMxJADqGo0qdcdrmsLJ4wOpWVEllRdaREMtPYETyKIiPWTMigFDBlLQCuT-60na4N9sPCZWMMo15um0h75Sgo~lq-PLRTOpkEsszly89E0RO7rb7UwZVhS~7QUjtXf1NQULRJ0LAOLlEZYIXyScjZgPeBraJmYuXgd7vKhHaPFGkhvj0aH1KjaLYpRHzHcjvNuATT~mnRQTJR~tXEEJ3rfgBorsqYjU1Zax19urnkZrI3XYWVLh13KwpXtf1yRIc~WYNczol62l~BCX6LjYHyG~4qqE8R8w7PuPZ9J~ApXYzalf7CA__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"
```
**References - 3.0.0 (November 19, 2018)**

 - Human reference (GRCh38) dataset required for Cell Ranger.
 - Download – 11 GB – md5sum: edb1dc39a0e379e0f226ed9ee004be3c
 - [Build steps][6]

```bash
wget http://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-GRCh38-3.0.0.tar.gz
```


 - Mouse reference dataset required for Cell Ranger.
 - Download – 9.6 GB – md5sum: 8ce6bc561e2554701fc43871301042e6
 - [Build steps][7]

```bash
wget http://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-mm10-3.0.0.tar.gz
```


**References - 3.1.0 (July 24, 2019)**

 - Human reference (GRCh38) dataset required for Cell Ranger.
 - Download – 9.8 GB – md5sum: 7a7c815b59d9ed965be58a91d1e36c20
 - [Build steps][8]

```bash
wget http://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-GRCh38-and-mm10-3.1.0.tar.gz
```

 - [x] **References 只需要选择合适的,通常选择第一个** 

# 三、CellRanger实操

1、准备数据

```bash
wget http://s3-us-west-2.amazonaws.com/10x.files/samples/cell-exp/3.0.2/5k_pbmc_v3/5k_pbmc_v3_fastqs.tar
wget http://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-GRCh38-3.0.0.tar.gz
tar -zxvf ./*.tar*

```

2、安装CellRanger
```bash
tar /czc/R_work/software/cellranger-3.1.0.tar.gz
PATH=$PATH:/home/czc/R_work/software/cellranger-3.1.0 >> ~/.bashrc
```

3、CellRanger实操

```bash
#!/bin/bash
cat list | while read id  #list 用于存储多个样本的样本名，这里只存储5k_pbmc_v3_fastqs一个样本名就行
do
cellranger count 
    --id ${id}_result #保存的文件夹
    --transcriptome /home/czc/R_work/software/refdata-cellranger-GRCh38-3.0.0 #参考序列 
    --fastqs ./$id #fastq文件（RI和R2必须,I1可有可无 ）
    --sample $id #存放fastq文件的文件名
    --jobmode=local  
    --localcores 4 #使用的本地核心数
    --nopreflight 
done

#这是本人使用CellRanger-2.1.0所用参数，如若使用最新版CellRanger可以适当修本文未注释参数
```


**CellRanger**是一款高度集成的上游分析套件，我们直接从官网上下载相应基因组索引(无需自己构建，若有额外需求，可以根据官方搭建索引的步骤进行修改)。这里我只演示CellRanger的使用方法，通常来说，CellRanger mkfastq无需自己来进行处理，公司负责提供原始的Fastq文件，除非自己实验室负责从自己建库到测序分析全部步骤。另外CellRanger aggr在整合样本的效果并不好，不能很好的消除样本之间的批次效应，虽然最新版的aggr添加了fastmnn算法来优化这些结果，但是我们依旧可以在得到表达矩阵之后再来处理批次效应造成的影响，更加方便利于数据的处理分析。


---

###**附：最新版CellRanger官方教程链接**
[CellRanger mkfastq官方教程](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/mkfastq)

[CellRanger count官方教程](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count)

[Cellranger aggr官方教程](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/aggregate)

---

### **参考文献：**
【1】：https://support.10xgenomics.com/single-cell-gene-expression/software/overview/welcome


  [1]: http://static.zybuluo.com/czc/pijvef5murcpuvb7zh9ni7eo/image_1e5alnj061qfm8bmas01ehl135s9.png
  [2]: http://static.zybuluo.com/czc/5re2izcchba1nygtzh6yphc6/image_1e5bpi59b1t1o1bll6b4s893mq20.png
  [3]: http://static.zybuluo.com/czc/zkh92fpy53xr91lvwj5l7wbx/image_1e5bmb9sdalehif1mtf167n104h9.png
  [4]: http://static.zybuluo.com/czc/u57du6jslmcgrk8vglzce3m9/image_1e5bmbjqr11jq1c5e1qhgdt8o8bm.png
  [5]: http://static.zybuluo.com/czc/n9xedrtrf7ozaqtmqbmx2nhw/image_1e5bmd8gao9dg5t16fhoua13cj13.png
  [6]: https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build#grch38_#%7Bfiles.refdata_GRCh38.version%7D
  [7]: https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build#mm10_#%7Bfiles.refdata_mm10.version%7D
  [8]: https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build#grch38mm10_#%7Bfiles.refdata_GRCh38_and_mm10.version%7D