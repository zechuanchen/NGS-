# GEO生存曲线数据分析

标签（空格分隔）： 杨思佳

---

参考：
[https://cloud.tencent.com/developer/article/1481945](https://cloud.tencent.com/developer/article/1481945)
[https://www.jianshu.com/p/1537efae5be9](https://www.jianshu.com/p/1537efae5be9) 徐州更
[https://www.jianshu.com/p/8dd7dc1e1719](https://www.jianshu.com/p/8dd7dc1e1719) 生物技能树


这篇主要讲的是怎么从GSE的**矩阵表达文件**到差异基因分析到生存曲线分析，用的gse号是[GSE21050](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE21050)。
# 数据介绍
GEO网站上的数据形式有：

- `Series Matrix File`**（所有样本的表达矩阵文件，经过归一化处理的），里面包括了 实验信息**, **分组信息** 和 **表达信息。是从袁术数据进行批次效应**
- `SOFT formatted family files`**（记录），跟表达矩阵文件差不多，但是存储结构不一样**
- `platforms：GPL570`**，是芯片平台，每个平台都有探针-基因对应关系**
- `TAR of CEL: `**单个样本的原始数据打包**

**
**理解CEL文件和Series Matrix File**
最开始得到的都是CEL文件，CEL文件需要**一系列的步骤才能转换成表达矩阵**，例如去除**批次效应、质控和过滤**等，得到的表达矩阵在上传时会增加元数据信息（处理方法、分组信息），就成为我们下载的`GSEXXXX_series_matrix.txt.gz`. 通过手工解析加R语言简单操作得到了R语言中的数据框(data.frame)， 而GEOquery能够帮助我们完成下载和解析这两个步骤。
三者的优先级为：GEOquery > 手工下载表达量矩阵文件 > 手工下载原始的CEL文件。
![!\[image.png\](https://cdn.nlark.com/yuque/0/2020/png/1239986/1587000345363-107e174a-d06a-45c8-a6ed-67d86f6bb89d.png#align=left&display=inline&height=258&margin=%5Bobject%20Object%5D&name=image.png&originHeight=516&originWidth=986&size=69936&status=done&style=none&width=493)][1]


这篇先介绍简单暴力的**Series Matrix File（所有样本的表达矩阵文件）**
# 数据下载方法
根据文章提供的GSE号，进入NIH相关页面：
[https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE21050](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE21050)


这里介绍两种方法：
数据获取有两种方式，R包**GEOquery**解析和手动下载。其中前面一种最方便，完成了手动数据下载和Bioconductor常见数据结构`ExpressionSet`的构造，关于这个数据结构的具体介绍看Bioconductor的介绍或者视频，简言之，就是用于存放 **实验信息**, **分组信息** 和 **表达信息**, 方便后续调用。
### GEOquery包下载
```r
library(GEOquery)
gset <- getGEO("GSE21050", GSEMatrix =TRUE, getGPL=F)
show(gset)
```
getGEO函数参数：
`getGEO(GEO = NULL, filename = NULL, destdir = tempdir(),`
` GSElimits = NULL, GSEMatrix = TRUE, AnnotGPL = FALSE, getGPL = TRUE,`
` parseCharacteristics = TRUE)`


- `GEO `= 可以是任何
- `destdir `= "." save to current directory, or any diectory u want
- `GSEMatrix `= TRUE download GSEmatrix
- `AnnotGPL `= FALSE 这个应该是也是GPL平台id转换信息，一般也选F 
- `getGPL=` TRUE  if True, download probe id-symbol relationship directly from GEO website，通常用F，因为可以用bioconductor自带包进行probe id转换，用F可以节省很多下载时间



### 直接下载


有时候网速感人，直接用getGEO下载会失败，下载Series Matrix File到当前目录，然后用getGEO函数直接进行读取成为ExpressionSet对象(他会查看当前目录有没有表达矩阵文件)
```r
gset <- getGEO(filename = 'GSE21050', destdir=".", AnnotGPL = F, getGPL = F)
```

# 数据提取与过滤
主要步骤：

- exprs(), pdata(), 分别提取表达矩阵和分组信息。
- 通过GPL平台，寻找对应bioconductor包对应的探针-基因名信息
- 表达矩阵中有些探针没有对应基因名，过滤
- 有些基因有多个探针（因为一个基因有很多transcript），得到rowMeans() 最大的探针位置作为该基因的探针
- 将探针名转换为基因名，得到最后的探针-基因名一对一的表达矩阵做下游分析
### exprs(), pdata(), 分别提取表达矩阵和分组信息
```r
## 表达矩阵提取
exprSet <- exprs(gset[[1]])
## 分组信息等提取
pData <- pData(gset[[1]])
```
### 通过GPL平台，寻找对应bioconductor包对应的探针-基因名信息
[http://www.bio-info-trainee.com/1399.html](http://www.bio-info-trainee.com/1399.html)用这个寻找对应bioconductor包，如果没有的话只能下载GPL平台的注释信息，然后自己进行parse。
```r
## get microarray platform, "GPL570", find corresponding microarray-SYMBOL transform library
## http://www.bio-info-trainee.com/1399.html, according to the website, 
## GPL570 corresponds to hgu133plus2 
gset[[1]]@annotation
library(hgu133plus2.db)

## this contains all id-transform capability
ls("package:hgu133plus2.db")

## dataframe, which has probeid -SYMBOL relationship
## Note that one gene may have multiple probes
ids=toTable(hgu133plus2SYMBOL)
```
### 表达矩阵中有些探针没有对应基因名，过滤
```r
## look at examples of genes which have multiple probe-id
tail(sort(table(ids$symbol)))

## some of probe in expSet does not have corresponding gene SYMBOL
## 对探针进行过滤，把没有对应基因名的探针过滤掉
table(rownames(exp) %in% ids$probe_id)
exp = exp[rownames(exp) %in% ids$probe_id,]
ids=ids[match(rownames(exp),ids$probe_id),]
```
### 有些基因有多个探针（因为一个基因有很多transcript），得到rowMeans() 最大的探针位置作为该基因的探针
```r
## for genes which have multiple probe_id, get the probe_id which has the maximum rowMeans
tmp = by(exp,
         ids$symbol,
         function(x) rownames(x)[which.max(rowMeans(x))])

probes = as.character(tmp)
dim(exp)
exp = exp[rownames(exp) %in% probes,] # 过滤有多个探针的基因
dim(exp)

## changing rownames(exp) from probe_id to SYMBOL
rownames(exp)=ids[match(rownames(exp),ids$probe_id),2]
```
## 差异基因筛选（高低转移）
虽然本文的目的只是画几个基因的生存曲线，不过顺便加上芯片差异分析的code。
### 查看是否符合limma假设
参考：[https://kasperdanielhansen.github.io/genbioconductor/html/limma.html](https://kasperdanielhansen.github.io/genbioconductor/html/limma.html)
一般来说matrix文件下载下来都是经过标准化的，但是这个数据集没有，所以我会用log2转换。因为limma包默认的logFC计算是默认输入矩阵经过logFC转换，然后直接拿两组的group mean相减就可以得到logFC值
而P值计算是改良版t.test，这里t检验的variance estimation来自ebayes (根据整体基因的variance来看单个基因的variance)，所以这是改良版的t.test
-but one where the variance estimation (the denominator of the t-statistics) is moderated by borrowing strength across genes (this is what eBayes() does); this is called a moderated t-statistic


log2转换的理由可以看这里：[https://blog.csdn.net/tuanzide5233/article/details/88542805](https://blog.csdn.net/tuanzide5233/article/details/88542805)
这个数据集只有高低转移的分组，在下面我会挑前几个，画boxplot，查看分布是否满足limma数据假设。
limma假设：数据**符合正态分布**，两组**equal variance（虽然其实不相等的组间差异对limma差异基因检验的影响很小）**。如果两组的variance差异大，可以用limma voom。
![!\[image.png\](https://cdn.nlark.com/yuque/0/2020/png/1239986/1587114248844-0c595701-29d5-4474-97a3-a7cac57bd2b9.png#align=left&display=inline&height=579&margin=%5Bobject%20Object%5D&name=image.png&originHeight=579&originWidth=695&size=33123&status=done&style=none&width=695)][2]
从boxplot可以看到组间差异不大，然而..并不是正态分布emm（可以根据hist 或者 shapiro test 进行判断）
根据[https://stat.ethz.ch/pipermail/bioconductor/2013-November/056214.html](https://stat.ethz.ch/pipermail/bioconductor/2013-November/056214.html)所说，limma对非正态分布也有robustness
>_The limma code is very robust against non-normality.  All the usual_
>_microarray platforms and standard preprocessing procedures produce data_
>_that is normally distributed to a good enough approximation.  Much effort_
>_has been devoted to developing good preprocessing and normalization_
>_algorithms._
以下为boxplot，shapiro.test，plotMDS的代码
```r
test <- exp.clean[,1:6]
test.meta <- data.frame(type = group_list[1:6], variable  = colnames(test))
test.l <- test%>% melt()
final <- merge(test.l, test.meta, by = "variable")
final %>% ggplot(aes(x=variable, y =log2(value+1), fill = type))+  geom_boxplot()
plotMDS(test, labels = test.meta$type)
shapiro.test(log2(test$GSM525806+1)[1:3000])
```
### limma差异分析code
#### 建一个design矩阵
```r
library(limma)
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(exp.clean)
design
```
#### 建一个contrast矩阵
```r
## create contrast matrix
## 这个矩阵声明，我们要把yes and no进行差异分析比较
contrast.matrix<-makeContrasts(yes-no,levels = design)
contrast.matrix
#contrast.matrix<-makeContrasts(Ten-Control,Eleven-Control,Twelve-Control, levels=design)
```
#### 差异分析
```r
##### 差异分析
##step1 线性模型拟合
fit <- lmFit(exp.clean,design)
##step2 根据对比模型进行差值计算  
fit2 <- contrasts.fit(fit, contrast.matrix)  
##step3 贝叶斯检验
fit2 <- eBayes(fit2)  
##step4 生成所有基因的检验结果报告
# coef =1 first pair, n=inf means all gene,
# topTable(fit2, corf =1, adj = "fdr", number =50, sort.by ="logFC")
tempOutput = topTable(fit2, coef=1, n=Inf)

## 根据筛选阈值筛选DGE
dif <- tempOutput[tempOutput[, "P.Value"]<0.01,]
```
开始偷懒：贴一下topTable的output:（随便拿了一个limma结果）
![!\[image.png\](https://cdn.nlark.com/yuque/0/2020/png/1239986/1587265564041-a1aec20b-f237-41f9-897b-097a4897c6a0.png#align=left&display=inline&height=35&margin=%5Bobject%20Object%5D&name=image.png&originHeight=70&originWidth=698&size=12304&status=done&style=none&width=349)][3]
[https://stat.ethz.ch/pipermail/bioconductor/2004-September/006132.html](https://stat.ethz.ch/pipermail/bioconductor/2004-September/006132.html)统计意义的参考
The output from `topTable()` includes

- `logFC`: the log fold-change between cases and controls.
- `t`: the t-statistic used to assess differential expression.
- `P.Value`: the p-value for differential expression; this value is not adjusted for multiple testing.
- `adj.P.Val`: the p-value adjusted for multiple testing. Different adjustment methods are available, the default is Benjamini-Horchberg.
- `B value`: The B-statistic (lods or B) is the log-odds that that gene is 

differentially expressed. Suppose for example that B=1.5. The odds of 
differential expression is exp(1.5)=4.48, i.e, about four and a half to 
one. The probability that the gene is differentially expressed is 
4.48/(1+4.48)=0.82, i.e., the probability is about 82% that this gene is 
differentially expressed. A B-statistic of zero corresponds to a 50-50 
chance that the gene is differentially expressed. The B-statistic is 
automatically adjusted for multiple testing by assuming that 1% of the 
genes, or some other percentage specified by the user, are expected to be 
differentially expressed. If there are no missing values in your data, then 
the moderated t  and B statistics will rank the genes in exactly the same 
order. Even you do have spot weights or missing data, the p-values and 
B-statistics will usually provide a very similar ranking of the genes.
#### 差异分析做图：
火山图，heatmap, MA plot各种都可以画，不赘述。
## 生存（事件）曲线做图
做生存曲线最常用的包是**survival**和**survminer**。
```r
library(survival)
library(survminer)
```
回顾一下，之前从pData <- pData(gset[[1]])中提取了metainfo，这里我们的metainfo包括MFS（metastatus free survival）事件发生与否，以及MFS维持的时间。然后我们将我们感兴趣的基因表达经过二分分为**高和低，然后进行事件曲线做图。**类似于生存曲线做图，我们只需要三个数据就可以做最简单的图：

- 基因表达二分类，或者是别的性别，种族等
- 事件发生与否
- 事件发生时间

绘图矩阵长这样：

![image_1e6bd4vkb1tisn9a3i01vsfcub2a.png-94.7kB][4]

### 创建surv object
第一个为事件发生事件，第二个为兴趣事件（这里我们的兴趣事件是转移）
```r
##with(aml, Surv(time, status))
surv_obj<- Surv(as.numeric(newdata$MFS), newdata$metastasis=="yes")
```
### 决定基因表达量高低cut off阈值
这里有两种常用的办法，

- 一种是人为定义：比如用mean或者median
- 一种是surv包：帮我们找寻p值最小的cut_off点：**用maxstat方法，surv_cutpoint函数**（但是这样出来的图虽然p小，生物学意义可能不大）

根据cut_off就能将我们的基因表达量分为高表达组和地表达组了。


`surv_cutpoint`(data, time = "time", event = "event", variables,
minprop = 0.1, progressbar = TRUE)
`data:` 数据
`time:` time of event
`event:` 需要变成0，1， event是1，因此需要修改一下input
`variable:` 需要函数计算的cutoff point的连续变量
```r
## decide cut-off value as high/low expression for PDE4D 
## 人工选择，选择为mean
cut_off <- mean(newdata$PDE4D)

## 函数选择p值最小的cut_off点，用surv_cutpoint函数
# a <- newdata %>% dplyr::mutate(metastasis.status = ifelse(metastasis =="no", 0, 1  ))
# cut_off_info <-surv_cutpoint(
#   a,
#   time = "MFS",
#   event = "metastasis.status",
#   variables = c("PDE4D")
# )
# 
# cut_off <- cut_off_info$cutpoint$cutpoint

## 根据cut_off点分组
newdata <-newdata %>% dplyr::mutate(PDE4D.status = ifelse(PDE4D <cut_off, "low","high"))

```
用survdiff检验某个变量对事件以及事件发生事件的影响
`survdiff`(formula, data, subset, na.action, rho=0, timefix=TRUE)
`formula:` surv(time, event) ~ variable e.g.  surv(time, event) ~ sex
`na.action:` na data filtering funciton
`rho:` a scalar parameter that controls the type of test。取值有0,1 是不同的统计检验方法
`rho = 0` this is the log-rank or Mantel-Haenszel test, and with `rho = 1` it is equivalent to the Peto & Peto modification of the Gehan-Wilcoxon test.
```r
# survdiff(surv_obj~variable1 + variable2 , data=newdata)
survdiff(surv_obj~PDE4D.status, data=newdata)
```
### 用ggsurvplot画图
检验默认是log rank test
```r
sfit <- survfit(surv_obj~PDE4D.status,data=newdata)
ggsurvplot(sfit, conf.int=F, pval=TRUE, pval.method =T) + ylab("MFS probability")
```

![image_1e6bd5llc17k1vibe2318c91rin2n.png-79.2kB][5]

### 更好看的图
```r
p1<-ggsurvplot(sfit,
              pval = TRUE, conf.int = TRUE,
              risk.table = TRUE, # Add risk table
              risk.table.col = "strata", # Change risk table color by groups
              linetype = "strata", # Change line type by groups
              surv.median.line = "hv", # Specify median survival
              ggtheme = theme_bw(), # Change ggplot2 theme
              palette = c("#E7B800", "#2E9FDF"), title = "Metastasis-Free Survival Curve for PDE4D (GSE21050 n=309)" ,
              ylab ="metastasis-free survival")
p1
```
![image_1e6bd6av81q9u1e29pc9ichnht3h.png-121.1kB][6]
**以后再填芯片原始数据处理以及生存曲线统计介绍的坑。咕咕咕。**


  [1]: http://static.zybuluo.com/czc/lhrwmc36ij7fbp0jkixxast7/image_1e6bd274ie0d1o9s11d01bp411599.png
  [2]: http://static.zybuluo.com/czc/6tk1y6glqg3qadqacefofil6/image_1e6bd32vea411sej1r1rb2f1p7fm.png
  [3]: http://static.zybuluo.com/czc/blwx2c435a5verbqxymhniti/image_1e6bd3laj1indo121kjf9ki1e5113.png
  [4]: http://static.zybuluo.com/czc/szuopg0va28bw51d3j3wuoxv/image_1e6bd4vkb1tisn9a3i01vsfcub2a.png
  [5]: http://static.zybuluo.com/czc/72ytvc5k4dw35wh3f1nkd6c8/image_1e6bd5llc17k1vibe2318c91rin2n.png
  [6]: http://static.zybuluo.com/czc/68ojkr10ndlxy80dmh9sig08/image_1e6bd6av81q9u1e29pc9ichnht3h.png