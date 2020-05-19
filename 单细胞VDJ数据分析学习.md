# 单细胞VDJ数据分析学习

标签（空格分隔）： 单细胞教程

---

# 一、免疫组库基础知识

## 1、什么事免疫组库

**免疫组库** ：淋巴细胞抗原受体多样性的产生。
**适应性免疫的淋巴细胞与固有免疫的非淋巴细胞之间的区别** ：淋巴细胞具有结构多样的抗原受体储备。

![image_1e80htefg1t5r1e848bq1sair6v10.png-446.7kB][1]

**TCR与BCR：** T细胞表面受体 与 B细胞表面受体。

**功能异同：**虽然均能够检测抗原，但是TCR主要是通过识别MHC上来识别抗原，而BCR可以直接与自由存在的可溶性抗原结合。

![image_1e80i1sde1sg1p7t1mpn1g8eim11d.png-171.6kB][2]

## 2、BCR结构 + TCR结构讲解

BCR=> 2条重链(H) + 2条轻链(L) 。 

重链 => 1个可变区(VH) + 3个恒定区(CH1/CH2/CH3)

轻链 => 1个可变区(LH) + 1个恒定区(CL)

CH => IG的恒定区 => 免疫球蛋白重链 => IgA/G/M/D/E + k/λ轻链

VH + CH => IG的可变区 => 决定它与什么抗原结合

![image_1e80if1b61pdf1cimp13448pfs1q.png-416.3kB][3]

可变区 => 4个骨架区(FR1/2/3/4) + 3个互补决定区(CDR1/2/**3**) => FR1/CDR1/FR2/CDR2/FR3/CDR3/FR4

Ig抗原结合部位 => 互补决定区(CDR) => VL + VH => 各有3个HVR 

![image_1e80j1f3aon71bhlu3jtnom727.png-379.9kB][4]

归根到底，淋巴细胞抗原受体多样性的产生是由于胚系基因的重组和重排。

人体4种受体链 => TCRa链、TCRβ链、BCR轻链、BCR重链

TCRa => 5’ L(前导链) + V(可变) + J(连接) + C(恒定)
TCRβ => 5’ L(前导链) + V(可变) + J(连接) + D(多样) + C(恒定)
BCR L => 5’ L(前导链) + V(可变) + J(连接) + C(恒定)
BCR H => 5’ L(前导链) + V(可变) + J(连接) + D(多样) + C(恒定)

![image_1e80je19bedm1asrt781mn64nt2k.png-951.1kB][5]

基因重排发生在B细胞成熟前
体系高高频突变(SHM)发生在B细胞成熟后

重组激活基因：RAG1(gene) 如果被阻断，使VDJ重组失衡、原发性免疫缺陷病


## 3、单细胞免疫组库的目的


![image_1e80jnur2n0b17v8fo11imp1lo031.png-866kB][6]


**目的：**


 1. 揭示克隆性(特异性)、多样性、抗原特异性和细胞环境
 2. 组装并注释V(D)J基因序列
 3. 从单个T细胞识别α和β链序列
 4. 将来自单个B细胞的重链与轻链免疫球蛋白(Ig)序列以全同型分辨率匹配
 5. 同时测定同一个细胞中TCR、B细胞Ig 、细胞表面蛋白与5'基因表达
 6. 同时配对TCRα和β链与TCR-pMHC特异性序列


**应用：**

![image_1e837sfaqm3s6nrkhi17818id9.png-133.1kB][7]

 1. 捕捉肿瘤发生时免疫微环境的变化、寻找免疫治疗的靶点、从而辅助免疫治疗更好的抗击肿瘤
 2. 器官或者骨髓移植时，经常会诱发宿主排斥反应的发生，从而发生慢性移植抗宿主病
 3. 自身免疫免疫性疾病是由于机体对自身抗原发生免疫反应而导致自身组织损害所引起的疾病
 4. 免疫组库在感染性疾病、抗体开发、用药以及疫苗评估等方面均有应用价值。例如通过免疫组库研究，检测感染类疾病过程中的免疫动态变化。在抗体开发方面，获得特征性的BCR序列，缩短抗体开发流程；针对某种疾病用药后的外周血样本进行评估，确认药物是否激发免疫反应及其功效。


## 4、名词解释

 11. Clonetype:通过精准的核苷酸匹配，收集共享一组CDR3生产序列的细胞
 12. Consensus: 对于一个单一的克隆型和链，在所有的细胞和该克隆型之间建立的共识为该链。这个共识是通过从克隆型细胞中重组相应的contig建立的
 13. Contig：通过组装产生的连续序列
 14. Full-length：如果一个contig与V基因的起始部分匹配，继续存在并最终与J基因的末端部分匹配。则其为全长的。
 15. Gem group:当将不同组的gem库合并到一个分析中时，我们在每个读取的条形码上附加一个小整数in silico，以识别读取的来自哪个库。这可以防止条形码冲突，否则会在虚拟双重态的形式中造成混乱。


----
# 二、数据表格说明

## 1、 BCR组装出来的高质量Contig 

![image_1e80r4j751i6o1h1itpbb3k1sl03r.png-24.4kB][8]

 - rownames: barcode
 - clonotype_id : 
 - barcode : 细胞标记序列
 - is_cell: 是否是细胞
 - contig_id: 组装出来的Contig编号
 - high_confidence：组装出来的这个Contig是否可信
 - length：组装出来的Contig长度
 - chain：组装出来的Contig比对到TCR和BCR中的哪个基因上

![image_1e816uca4s5n1l5n1lgh8i11e6i9.png-26.2kB][9]

 - v_gene：组装出来的Contig 比对上的已知V 基因名称
 - d_gene：组装出来的Contig 比对上的已知D 基因名称
 - j_gene：组装出来的Contig 比对上的已知J 基因名称
 - c_gene：组装出来的Contig 比对上的已知C 基因名称
 - full_length：是否全长，是否跨越了V的5'端与J的3'端
 - productive：全长 + V有一个起始密码子，在V-J区间有一个可翻译阅读框的CDR3序列，但是没有终止密码子；如果Contig不能完全跨越V和J则设置为None
 - cdr3：cdr3氨基酸序列
 - cdr3_nt：cdr3的核苷酸序列
 - reads：组装这个Contig使用的序列数
 - umis：组装这个Contig使用的UMI数
 - raw_clonetype_id:克隆类型编号
 - raw_consensus_id:Consensus序列编号

## 2、Contig---clonetype统计结果

![image_1e8191r83vl91akf1j5ntfs1kkpm.png-40.1kB][10]

 - clonetype_id:
 - frequency: 同类型的clonetype数量
 - proportion:数量百分比
 - cdr3s_aa:对应的cdr3区域的氨基酸含量
 - cdr3s_nt:对应的核苷酸序列

----

# 三、数据分析

## 1、CellRanger vdj 测序分析

```
cat list| while read id
do
/home/xmzhang/czc/software/cellranger-3.1.0/cellranger vdj --id=${id}_result --reference=/home/xmzhang/czc/R_work/Zhang/1_2020-5-2/2_VDJ_learn/0_reference/refdata-cellranger-vdj-GRCh38-alts-ensembl-3.1.0 --fastqs=./${id}  --sample=$id --localcores=8  --localmem=100
done
```

## 2、 使用scRepertoire分析TCR和BCR

### 1、加载所需要的R包

[github链接地址][11]

安装包的过程——略 顺便一提，这包安装蛮麻烦的，我github bug了，安不了，本地安装需要自己解决依赖包问题，依赖包详见上述github主页。

```
suppressMessages(library(scRepertoire))
suppressMessages(library(Seurat))
```

### 2、加载TCR的CellRanger结果

```
setwd("./R_work/2020-5-2/1_result/")
tcr <- read.csv("./0_sample/2_TCR/filtered_contig_annotations.csv",stringsAsFactors = F)
tcr$barcode <- gsub("-1", "", tcr$barcode) # gsub的功能见下图
contig_list<-list()
contig_list[[1]]<-tcr
x<-contig_list[[1]]
```
![image_1e895bekav9t1ips95at7l13u99.png-4.4kB][12]

### 3、对TCR数据结果进行最初的整合+数据本身的分类

在TCR数据中有一个问题，就是一条Barcode对应多个TCR序列，包括TCRa与TCRb两条链的信息，甚至有可能出现第三条链，现在的Barcode信息不适合直接与Seurat整合，如果按照[简书上的教程][13]或者是[BioSTAR上的教程][14]都会出现大量数据丢失的情况，而这里使用combineTCR函数可以轻松解决该问题。

另外，即使是针对单个数据，这里依旧需要Run combineTCR函数实现对数据本身的分类，否则格式会报错。

多样本例子教程请[点击这里][15]

```
# 单样本例子
combined <- combineTCR(contig_list, samples = c("PBMC_1"), ID = "Normal", removeNA = T, removeMulti = T,cells = c("T-AB"))
```

### 4、可视化Clonetype -- 可视化每个样本的Unique CloneType

```
quantContig(combined, cloneCall="gene+nt", scale = T)
```

![image_1e895fpmj11ol1rjcpfv1jpa170rm.png-9.9kB][16]

咳咳，可能需要多几个样本的例子看起来比较舒服，不过多样本的例子都有了，就懒得那啥了

### 5、可视化Clonetype -- 可视化每个样本的CloneType的丰度

```
abundanceContig(combined, cloneCall = "gene", scale = F)
```

![image_1e895ivlkmer1d7v1nlnd6q1lav13.png-7.8kB][17]

### 6、可视化Clonetype -- 可视化每个样本的CloneType的CDR3的长度分布情况

```
# 一起看氨基酸的长度
lengthContig(combined, cloneCall="aa", chains = "combined") 
# 分开TCRa和TCRb的核苷酸长度情况
lengthContig(combined, cloneCall="nt", chains = "single") 
```

![image_1e895or5e1f091em2vom6a3ju11g.png-6.2kB][18]

![image_1e895p0u7opq164l1bu0bjq11661t.png-6.8kB][19]


### 7、可视化Clonetype -- 可视化每个样本的CloneType稀有程度分布差异

```
clonalHomeostasis(combined, cloneCall = "gene")
clonalProportion(combined, cloneCall = "gene") 
```

![image_1e895rc411nh6e7886ev1cqvd2a.png-13.3kB][20]

### 8、 可视化CloneType -- 与Seurat联合分析
这里注意一个问题，虽然都是同一批数据，但是很可能出现数据Barcode对不上的情况，就是BCR、TCR和表达谱数据测到的Barcode会有部分对不上，主要是由于测序深度等因素造成的。

**整合**

```
load("./2_T_cells/T_cells.RData")
# 加载已经运行好的表达谱数据
object@meta.data<-object@meta.data[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","seurat_clusters","Cluster")]
# 我这里把之前加的标注信息都删了，就保留了上面几个

combined[[1]]$barcode<-substr(combined[[1]]$barcode,15,30)
# 然后，由于要根据barcode序列整合，我这里就吧之前加在combined上的"sample"和"ID"都删了

seurat <- combineSeurat(combined, object, cloneCall="gene", groupBy = "sample")

```

**整合后的可视化**

```
colorblind_vector <- colorRampPalette(c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6"))

seurat@meta.data$cloneType <- factor(seurat@meta.data$cloneType, levels = c("Hyperexpanded (100 < X <= 500)", "Large (20 < X <= 100)", "Medium (5 < X <= 20)", "Small (1 < X <= 5)", "Single (0 < X <= 1)", NA))

DimPlot(seurat, group.by = "cloneType") +
  scale_color_manual(values = c(rev(colorblind_vector(5))), na.value="grey")

```

![image_1e896l1cr12hi133q1njg189jsos2n.png-39kB][21]

![image_1e896lcj0dge43qvqbe6214j034.png-9.5kB][22]


**分开查阅各个Cluster的CloneType的差异**

```
meta <- data.frame(seurat@meta.data, seurat@active.ident) 
ggplot(meta, aes(x=seurat.active.ident, y=Frequency)) + 
  geom_boxplot(outlier.alpha = 0, aes(fill=seurat.active.ident)) + 
  guides(fill=F) + 
  theme_classic() + 
  theme(axis.title.x = element_blank())
```

![image_1e896ucqr1qo8101h1etm1lf53o3h.png-9.5kB][23]

**将各个CLuster的CloneType与所有数据的CloneType进行比较**

```
seurat$Patient<-"Normal"
alluvialClonotypes(seurat, cloneCall = "gene", compare = "cluster", facet = "Patient")
```

![image_1e8970gbr1ff8maa1jga18qa17jq3u.png-489.1kB][24]

**查看各个CLuster CloneType的多样性**

```
clonalDiversity(seurat, cloneCall = "nt", colorBy = "cluster")
```

![image_1e8970uglm39f4gvhjf4jmv24b.png-11.4kB][25]

**查看各个CLuster CloneType的重叠情况**
```
clonalOverlap(seurat, cloneCall="aa", method="overlap")
```

![image_1e89714a0p0tjgt16g011vb1kqo4o.png-22.8kB][26]

----
参考文献：

【1】：[简书:单细胞数据Seurat整合分析](https://www.jianshu.com/p/b0673107dd15)
【2】：[Cell:肿瘤免疫微环境的类器官建模](https://www.cell.com/cell/fulltext/S0092-8674(18)31513-7?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867418315137%3Fshowall%3Dtrue)
【3】：[10X support:CellRanger多种库类型分析](https://support.10xgenomics.com/single-cell-vdj/software/analysis-of-multiple-libraries/latest/overview)
【4】：[10X support:单细胞免疫分析数据集](https://support.10xgenomics.com/single-cell-vdj/datasets)
【5】：[Biostar VDJ测序数据与Seurat集成](https://www.biostars.org/p/384640/)
【6】：[Seurat将VDJ数据与scRNA-seq相关联](https://www.biostars.org/p/383217/)
【7】：[Cell:通过单细胞测序解释肝癌中的浸润性T细胞](https://sci-hub.tw/10.1016/j.cell.2017.05.035)
【8】：[Cell:同上](https://www.cell.com/cell/fulltext/S0092-8674(17)30596-2?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867417305962%3Fshowall%3Dtrue)
【9】：[知乎：10X Genomics 单细胞测序在免疫组库研究中的应用](https://zhuanlan.zhihu.com/p/44788585)
【10】：[简书：10X Genomics单细胞免疫组库VDJ分析基础](https://www.jianshu.com/p/db4831091a5c)


  [1]: http://static.zybuluo.com/czc/eh3ljtc4rtlcr1708uu7mtyq/image_1e80htefg1t5r1e848bq1sair6v10.png
  [2]: http://static.zybuluo.com/czc/t2epiroz4uruo1tha0gk2o2y/image_1e80i1sde1sg1p7t1mpn1g8eim11d.png
  [3]: http://static.zybuluo.com/czc/4lexc9wndfr768pv0nn57b4q/image_1e80if1b61pdf1cimp13448pfs1q.png
  [4]: http://static.zybuluo.com/czc/wr9etrqodk3o093edylrcfmg/image_1e80j1f3aon71bhlu3jtnom727.png
  [5]: http://static.zybuluo.com/czc/dr2z9l0ej5zi8551czrqpuws/image_1e80je19bedm1asrt781mn64nt2k.png
  [6]: http://static.zybuluo.com/czc/jic5gqyiw80m3sl9tk0d2juv/image_1e80jnur2n0b17v8fo11imp1lo031.png
  [7]: http://static.zybuluo.com/czc/5a9kdp89p7zpt2e95fuwoxse/image_1e837sfaqm3s6nrkhi17818id9.png
  [8]: http://static.zybuluo.com/czc/nxs7sum0mhigim5hkhrkgwxh/image_1e80r4j751i6o1h1itpbb3k1sl03r.png
  [9]: http://static.zybuluo.com/czc/jt1i79yt35aei78lx4t7xyhl/image_1e816uca4s5n1l5n1lgh8i11e6i9.png
  [10]: http://static.zybuluo.com/czc/0qal0zvesbc5ta1ytkz2281p/image_1e8191r83vl91akf1j5ntfs1kkpm.png
  [11]: https://github.com/ncborcherding/scRepertoirege_1e8191r83vl91akf1j5ntfs1kkpm.png
  [12]: http://static.zybuluo.com/czc/nrw8pifakgyiqxxkdge6m3qz/image_1e895bekav9t1ips95at7l13u99.png
  [13]: https://www.jianshu.com/p/b0673107dd15
  [14]: https://links.jianshu.com/go?to=https://www.biostars.org/p/384640/
  [15]: https://ncborcherding.github.io/vignettes/vignette.html
  [16]: http://static.zybuluo.com/czc/r1w8p0exzhhp2jb9b2yku0s5/image_1e895fpmj11ol1rjcpfv1jpa170rm.png
  [17]: http://static.zybuluo.com/czc/rtjg828d51v8iqjoquf6c691/image_1e895ivlkmer1d7v1nlnd6q1lav13.png
  [18]: http://static.zybuluo.com/czc/lpnljg1el8mdunek3wj9j1s7/image_1e895or5e1f091em2vom6a3ju11g.png
  [19]: http://static.zybuluo.com/czc/99r83131bhcvc8ijuccang0y/image_1e895p0u7opq164l1bu0bjq11661t.png
  [20]: http://static.zybuluo.com/czc/n0bnetbhtt4b2jj5k62l8bvx/image_1e895rc411nh6e7886ev1cqvd2a.png
  [21]: http://static.zybuluo.com/czc/kes6bb03v190x499y91ece4s/image_1e896l1cr12hi133q1njg189jsos2n.png
  [22]: http://static.zybuluo.com/czc/2o2u22lcvgvrgoedy5t5x19o/image_1e896lcj0dge43qvqbe6214j034.png
  [23]: http://static.zybuluo.com/czc/zode5g7ei3cc7bwslqliatte/image_1e896ucqr1qo8101h1etm1lf53o3h.png
  [24]: http://static.zybuluo.com/czc/ja6obmsaqy3jj67v4xfozpjg/image_1e8970gbr1ff8maa1jga18qa17jq3u.png
  [25]: http://static.zybuluo.com/czc/ynrlhmuy0cxrhc1ln9mvk8u1/image_1e8970uglm39f4gvhjf4jmv24b.png
  [26]: http://static.zybuluo.com/czc/jd2yzdvco19jxufa6m55pfvk/image_1e89714a0p0tjgt16g011vb1kqo4o.png