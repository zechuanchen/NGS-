# TCGA数据分析（上）

标签（空格分隔）： 杨思佳

---

这篇主要是讲怎么从TCGA下载数据（包括mrna fpkm文件，以及对应的clinical data），
下篇讲怎么整合数据，进行生存分析。例子是TCGA-SARC project。
**参考：**

1. [https://mp.weixin.qq.com/s?__biz=MzIyNzk1NjUxOA==&mid=2247483756&idx=1&sn=544cb6cc3a7a7f37a2b4c32731e58d7f&chksm=e8580325df2f8a334167b1dff7f87f5a5735a478471f79506c9ae1a169757f9570f0fa95ce75&scene=21#wechat_redirect](https://mp.weixin.qq.com/s?__biz=MzIyNzk1NjUxOA==&mid=2247483756&idx=1&sn=544cb6cc3a7a7f37a2b4c32731e58d7f&chksm=e8580325df2f8a334167b1dff7f87f5a5735a478471f79506c9ae1a169757f9570f0fa95ce75&scene=21#wechat_redirect)

2. [https://mp.weixin.qq.com/s?__biz=MzIyNzk1NjUxOA==&mid=2247483735&idx=1&sn=04096935a31bba2f7fe80d68e4603f2e&chksm=e858031edf2f8a08fb816df49f40e024ec4de14c7c904f4152374c8fe9bcc7f9f45e1064ed9d&scene=21#wechat_redirect](https://mp.weixin.qq.com/s?__biz=MzIyNzk1NjUxOA==&mid=2247483735&idx=1&sn=04096935a31bba2f7fe80d68e4603f2e&chksm=e858031edf2f8a08fb816df49f40e024ec4de14c7c904f4152374c8fe9bcc7f9f45e1064ed9d&scene=21#wechat_redirect)

# TCGA数据库介绍

基本复制粘贴
TCGA（The Cancer Genome Atlas, 癌症基因组图谱）项目最早始于2005年，由美国政府出资，美国国家癌症研究所(National Cancer Institute)和美国人类基因组研究所(National Human Genome Research Institute)共同监督，旨在应用高通量的基因组分析技术，以帮助人们对癌症有个更好的认知，从而提高对于癌症的预防、诊断和治疗能力。作为目前最大的癌症基因信息数据库，TCGA的全面不仅仅体现在众多癌型上，还体现在多组学数据，包括**基因表达数据(mrna-seq)、miRNA表达数据(mrna-seq)、拷贝数变异(copy number variation)、DNA甲基化(methylation chip)、somatic mutation**。

![!\[image.png\](https://cdn.nlark.com/yuque/0/2020/png/1239986/1588320424984-b42c553c-92f4-42f1-9d7d-ae62e24d8ee6.png#align=left&display=inline&height=265&margin=%5Bobject%20Object%5D&name=image.png&originHeight=529&originWidth=996&size=149572&status=done&style=none&width=498)][1]
TCGA现在的数据均收录在GDC中，在GDC中可以通过**GDC Data Portal **和 **GDC Legacy Archive **这两种方式获得TCGA数据，GDC Data Portal 中的数据是最新经过统一标准整理的，而 GDC Legacy Archive 中的数据是所有未经处理的数据，更全面。**

# 数据下载

数据下载的方法有很多，这里先介绍怎么用**gdc client**进行下载。gdc client下载需要提供**manifest**文件，manifest文件可以在GDC官网先获得，然后可以通过gdc client进行下载。

## 进入官网：

[https://portal.gdc.cancer.gov/](https://portal.gdc.cancer.gov/)

## 点击Respository：

当然，你也可以通过Projects等panel来explore数据。

![!\[image.png\](https://cdn.nlark.com/yuque/0/2020/png/1239986/1588125619269-472dfa3b-76a4-4d3a-9f48-27f04d58ac8a.png#align=left&display=inline&height=209&margin=%5Bobject%20Object%5D&name=image.png&originHeight=417&originWidth=1151&size=208229&status=done&style=none&width=575.5)][2]

## 进入Respository之后，选择cases进行筛选：

可以先从cases进行筛选，选出你感兴趣的数据集。图例选择了TCGA数据库中的TCGA-BRCA。除了TCGA外，GDC还收录了其他的数据库数据。

![!\[image.png\](https://cdn.nlark.com/yuque/0/2020/png/1239986/1588320859765-56146d7a-3201-4a7f-80dc-df51c1f09875.png#align=left&display=inline&height=212&margin=%5Bobject%20Object%5D&name=image.png&originHeight=423&originWidth=310&size=12213&status=done&style=none&width=155)][3]

![!\[image.png\](https://cdn.nlark.com/yuque/0/2020/png/1239986/1588320805371-1824b1e9-277f-4637-aae5-e079736a8772.png#align=left&display=inline&height=385&margin=%5Bobject%20Object%5D&name=image.png&originHeight=769&originWidth=1767&size=177882&status=done&style=none&width=883.5)][4]

这里点击cases，选择TCGA数据库，选择TCGA-SARC

## 可以进一步筛选肿瘤亚型，感兴趣的race，年龄等等




## 下载clinical数据（可以跳过）

#### 点击files，选择clinical数据，下载clinical manifest以及json file

在files里选择clinical, 这里选择bcr xml（其他选项应该是病人资料混合的xml），然后点击右边的manifest，进行clinical manifest的下载，可以改名为clinical.txt

![!\[image.png\](https://cdn.nlark.com/yuque/0/2020/png/1239986/1588323419485-cfec3d39-5024-4556-a73f-6508de003115.png#align=left&display=inline&height=354&margin=%5Bobject%20Object%5D&name=image.png&originHeight=708&originWidth=1864&size=147170&status=done&style=none&width=932)][5]

打开manifest的txt文件是这样的，他提供了md5码，文件大小，文件名称等等。之后我们会用gdc-client 和这个manifest来批量下载xml clinical文件

![!\[image.png\](https://cdn.nlark.com/yuque/0/2020/png/1239986/1588321283548-76b0e9da-1d2b-4885-8e2c-ae2659cd0eb4.png#align=left&display=inline&height=109&margin=%5Bobject%20Object%5D&name=image.png&originHeight=219&originWidth=1248&size=53398&status=done&style=none&width=624)][6]


然后点击右边的json文件，这个文件包含了xml文件（并且这个xml有病人的TCGA.barcode信息） 与病人 uuid的一一对应关系。

![image_1e7ftg0dt1r45uk1nh91soo11ur3h.png-84.4kB][7]

**当然还有另外一种下载临床数据的方法（更方便），后面会讲到。**

## 下载感兴趣的数据（mRNA，methylation等）的manifest文件

例如在file里选择rna-seq，然后选择htseq-count(如果想用Deseq2或者edgeR进行差异分析)，或者HTseq-FPKM（标准化之后的数据，可以直接拿来做生存分析等），HTseq-FPKM-UQ（这个是一个类似于FPKM的标准化方法，但是分母除的是75%分位的表达量，相比于FPKM是除以总表达量）。

![!\[image.png\](https://cdn.nlark.com/yuque/0/2020/png/1239986/1588322072504-48ea875f-c7da-4d92-b86e-00322b18d4fe.png#align=left&display=inline&height=366&margin=%5Bobject%20Object%5D&name=image.png&originHeight=732&originWidth=1204&size=110251&status=done&style=none&width=602)][8]

## 下载clinical文件以及metadata

将上一步的265个mRNA-seq的数据放到cart中，然后进入cart界面，点击下载clinical以及metadata文件（json格式）

![!\[image.png\](https://cdn.nlark.com/yuque/0/2020/png/1239986/1588515450658-cc63bb80-4d56-4462-b80f-e9d4cc623248.png#align=left&display=inline&height=204&margin=%5Bobject%20Object%5D&name=image.png&originHeight=408&originWidth=1562&size=49533&status=done&style=none&width=781)][9]

![!\[image.png\](https://cdn.nlark.com/yuque/0/2020/png/1239986/1588515476265-75300014-90ed-4eba-a102-a3e26626764c.png#align=left&display=inline&height=247&margin=%5Bobject%20Object%5D&name=image.png&originHeight=494&originWidth=1863&size=86128&status=done&style=none&width=931.5)][10]

clinical中有病人所有的临床信息，包括submitter id（病人id）以及case_id（patient uuid）

![!\[image.png\](https://cdn.nlark.com/yuque/0/2020/png/1239986/1588515674098-1b638ae9-6f07-44f5-83a6-18e01d6beae3.png#align=left&display=inline&height=262&margin=%5Bobject%20Object%5D&name=image.png&originHeight=524&originWidth=1034&size=88854&status=done&style=none&width=517)][11]

![!\[image.png\](https://cdn.nlark.com/yuque/0/2020/png/1239986/1588515676321-02b9a638-bd45-4382-a219-0f37ce14ce31.png#align=left&display=inline&height=156&margin=%5Bobject%20Object%5D&name=image.png&originHeight=311&originWidth=901&size=45887&status=done&style=none&width=450.5)][12]

metadata中包括了文件名（FPKM.txt.gz与submitter id 与 case id）

![!\[image.png\](https://cdn.nlark.com/yuque/0/2020/png/1239986/1588515798551-ed4c398b-5391-4d93-97ec-37c034c90de3.png#align=left&display=inline&height=208&margin=%5Bobject%20Object%5D&name=image.png&originHeight=415&originWidth=1120&size=90420&status=done&style=none&width=560)][13]


## 安装GDC下载专用软件：

### 官方介绍：

[https://gdc.cancer.gov/access-data/gdc-data-transfer-tool](https://gdc.cancer.gov/access-data/gdc-data-transfer-tool)

TCGA中存储的测序数据文件，例如FASTQ和BAM文件，特别是全基因组的BAM文件有时可以达到200-300GB，所以需要一个稳定高效的数据下载和上传工具来处理数据库与用户之间的交互。

### 1. 安装gdc-client 软件

该软件支持在Windows, Linux, Mac OS 等不同操作系统上运行，下载下来，安装既可以使用，非常的方便，下载链接

[https://gdc.cancer.gov/access-data/gdc-data-transfer-tool](https://gdc.cancer.gov/access-data/gdc-data-transfer-tool)

根据自己的OS以及python版本下载，下载后解压后是一个名为gdc-client的可执行文件

![!\[image.png\](https://cdn.nlark.com/yuque/0/2020/png/1239986/1588322218552-8f574289-9722-4c46-8446-66c0fcc3adda.png#align=left&display=inline&height=182&margin=%5Bobject%20Object%5D&name=image.png&originHeight=363&originWidth=594&size=29852&status=done&style=none&width=297)][14]


```bash
cd ~/Biosofts 
wget https://gdc.cancer.gov/system/files/authenticated%20user/0/gdc-client_v1.5.0_Ubuntu_x64.zip

# 解压
unzip gdc-client_v1.5.0_Ubuntu_x64.zip  
# 查看是否安装
./gdc-client -h 

```

### 2. GDC网站筛选下载的数据

需要的clinical manifest以及mRNA的manifest文件，之前我们已经下好了。这里也可以直接用patient uuid进行下载。

```bash
## gdc-client 下载的指令很简单：
## gdc-client download -m manifestfile -d destdir

./gdc-client download -m mRNA_gdc_manifest.txt -d mRNA 
## ./gdc-client download -m clinical_gdc_manifest.txt -d clinical

## download by patient uuid
##./gdc-client download a7613c60-93a1-4a4c-9b12-dfd2103b0168
```

`-m: manifest file path`
`-d: destdir`

# TCGA数据下载陷阱：

[https://mp.weixin.qq.com/s?__biz=MzIyNzk1NjUxOA==&mid=2247483841&idx=1&sn=babc0f85ffd85227ef29af919774e89c&chksm=e8580388df2f8a9e2fa99ea7e11bb7c117ceeb9522c21f3b7778239905b8313a7888022db7cb&scene=21#wechat_redirect](https://mp.weixin.qq.com/s?__biz=MzIyNzk1NjUxOA==&mid=2247483841&idx=1&sn=babc0f85ffd85227ef29af919774e89c&chksm=e8580388df2f8a9e2fa99ea7e11bb7c117ceeb9522c21f3b7778239905b8313a7888022db7cb&scene=21#wechat_redirect)

**这篇推文说了关于TCGA数据下载会遇到的陷阱：**

大意如下，GDC Data Portal只提供4个平台的数据，并没有细分。

![!\[image.png\](https://cdn.nlark.com/yuque/0/2020/png/1239986/1588325431471-3262abf4-37f9-4dda-a876-c1393e5c9f65.png#align=left&display=inline&height=75&margin=%5Bobject%20Object%5D&name=image.png&originHeight=150&originWidth=384&size=8925&status=done&style=none&width=192)][15]


当作者选择了miRNA之后，其实miRNA的平台在illumina下面有Hiseq 和 GA 两种平台，而两种平台的数据在作者做了PCA之后，实际上有明显的批次效应。因此在下载数据的时候，要注意有没有混合平台的数据一起下载了。

![!\[image.png\](https://cdn.nlark.com/yuque/0/2020/png/1239986/1588325502337-dacf5424-74c4-42c4-b6c7-0998bb767c36.png#align=left&display=inline&height=350&margin=%5Bobject%20Object%5D&name=image.png&originHeight=700&originWidth=895&size=235589&status=done&style=none&width=447.5)][16]

这时候就可以去GDC Legacy Archive

![!\[image.png\](https://cdn.nlark.com/yuque/0/2020/png/1239986/1588325567574-8e56c802-3cd6-4f59-8c70-4a0fe9f71a9e.png#align=left&display=inline&height=352&margin=%5Bobject%20Object%5D&name=image.png&originHeight=703&originWidth=1885&size=459656&status=done&style=none&width=942.5)][17]

发现legacy平台提供更多的平台选择：

![!\[image.png\](https://cdn.nlark.com/yuque/0/2020/png/1239986/1588325731442-e9684891-0731-48ed-a081-44381b74ddca.png#align=left&display=inline&height=407&margin=%5Bobject%20Object%5D&name=image.png&originHeight=813&originWidth=318&size=37649&status=done&style=none&width=159)][18]

这边我们检查一下sarc的rna-seq是不是多个平台（发现并不是）：

![!\[image.png\](https://cdn.nlark.com/yuque/0/2020/png/1239986/1588325663279-932e0501-15ad-4d86-953c-6ec6d64224e8.png#align=left&display=inline&height=284&margin=%5Bobject%20Object%5D&name=image.png&originHeight=568&originWidth=316&size=22078&status=done&style=none&width=158)][19]

## TCGA.barcode介绍（submitter_id）

介绍一下TCGA.barcode

[https://zhuanlan.zhihu.com/p/48326648](https://zhuanlan.zhihu.com/p/48326648)
![!\[\](https://cdn.nlark.com/yuque/0/2020/jpeg/1239986/1588326139423-260303e2-f94a-4eeb-8e2b-6dbe362fe156.jpeg#align=left&display=inline&height=181&margin=%5Bobject%20Object%5D&originHeight=181&originWidth=454&size=0&status=done&style=none&width=454)][20]

我们将此示图以"-"分割，具体拆开解读一下：

- **TCGA：Project**, 所有TCGA样本名均以这个开头
- **02：Tissue source site**，组织来源编码，

[https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tissue-source-site-codes](https://link.zhihu.com/?target=https%3A//gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tissue-source-site-codes)

- **00001：Participant**, 参与者编号
- **01：Sample**, 这两个数字可以说是最关键、最被大家注意的，其中编号**01~09表示肿瘤**，10~19表示正常对照，所以在TCGA样本名中，这个位置最常见的就是01和11
- **C：Vial**, 在一系列患者组织中的顺序，绝大多数样本该位置编码都是A; 很少数的是B
- **01：Portion**, 同属于一个患者组织的不同部分的顺序编号，同一组织会分割为100-120mg的部分，分别使用
- **D：Analyte**, 分析的分子类型，对应关系如下所示：



![!\[\](https://cdn.nlark.com/yuque/0/2020/jpeg/1239986/1588326272370-6a68b912-7e4c-497d-9376-97cb0def0033.jpeg#align=left&display=inline&height=409&margin=%5Bobject%20Object%5D&originHeight=409&originWidth=627&size=0&status=done&style=none&width=627)][21]


- **0812：Plate**, 在一系列96孔板中的顺序，值大表示制板越晚
- **01：Center**, 测序或鉴定中心编码



# 总结：

## 最后看看我们手上拥有的是什么数据：

本篇主要介绍了通过gdc-client，以及gdc data portal下载：

- 感兴趣的癌症mRNA-seq（FPKM数据）， 以及文件对应的json文件（metadata）。
- 对应的FPKM数据的clinical data



## 数据整合的方向：

mRNA-seq数据包括了每个样本的基因的FPKM信息，需要进行整合，变成一个完整的data table。
metadata有mRNA-seq文件名和submitter_id（其中包含数据处理经过，以及是癌旁/癌组织样本），以及patient uuid
clinical data有submitter_id，以及patient uuid。

因此我们需要通过patient uuid 以及 submitter_id对mRNA-seq数据进行整合，最后将mRNA表达值与病人临床指标联系起来（比如Overall Survival等等）

----



  [1]: http://static.zybuluo.com/czc/wvex4r0vmbqtgelo2c4bcnsk/image_1e7ftai57qd21cuu6jj76jgr9.png
  [2]: http://static.zybuluo.com/czc/r50grgiu8x444dbfl7qi3s72/image_1e7ftb5e2odqkhv1avk1tv51611m.png
  [3]: http://static.zybuluo.com/czc/bukprvqxjthxjqofrtiggo66/image_1e7ftbie5mj0obm10vv1332gjb13.png
  [4]: http://static.zybuluo.com/czc/rw1u5jbnope4xc421go87zeo/image_1e7ftehu51n3e1j0cu1q1f8ets52a.png
  [5]: http://static.zybuluo.com/czc/wsrr0y99a2gdch3xisukcqxb/image_1e7ftf9en1joo3sh13q7h9ndbr2n.png
  [6]: http://static.zybuluo.com/czc/e6rfdrsuqu4h29vwnzgwf69n/image_1e7ftfkk6uj01vvdbfe5ld67p34.png
  [7]: http://static.zybuluo.com/czc/l0krh5dv42uljpb103t8zrga/image_1e7ftg0dt1r45uk1nh91soo11ur3h.png
  [8]: http://static.zybuluo.com/czc/tiliglgrbcpf74m9p87mln03/image_1e7ftgred1a841j2h7qd1eacrq33u.png
  [9]: http://static.zybuluo.com/czc/jgikaxytrsfp25ok2312aubk/image_1e7fthrld3fs7671iqqio412t94b.png
  [10]: http://static.zybuluo.com/czc/mc34rgxfiqkl3rntkk1seiti/image_1e7fti98tjqp1q5m7jm16lpo2d4o.png
  [11]: http://static.zybuluo.com/czc/ls8zhl6codwq10pymuvldnbu/image_1e7ftimgj1g3u10k3oddqht1dc255.png
  [12]: http://static.zybuluo.com/czc/1fav29czt027ovnb6uxk6q98/image_1e7ftj4ga1ju0htp1856voe1nl35i.png
  [13]: http://static.zybuluo.com/czc/hgu9bcxz63c0i6eqfpoizn6i/image_1e7ftjfqs4ic1pn11h41ens12365v.png
  [14]: http://static.zybuluo.com/czc/4m4b88ljkaw3waeadfzq0zjo/image_1e7ftk2u71k5f1f8m1gqa181uimc6c.png
  [15]: http://static.zybuluo.com/czc/hykkbyfk6nqmlfpeh93wpxxl/image_1e7ftl6961a2tgrbq2b1apj1p886p.png
  [16]: http://static.zybuluo.com/czc/mokxd4wmj2rk2r2att4su6pp/image_1e7ftlh3vaddon2l4j192tqt776.png
  [17]: http://static.zybuluo.com/czc/0wxmpzxl3spa3wu4lv9k85kf/image_1e7ftlsmh1cc01evg14fk1q7p17si7j.png
  [18]: http://static.zybuluo.com/czc/d2xx0ut7x6z01yumoaq6f9nn/image_1e7ftm8fk13171uaqbd71km91mqq80.png
  [19]: http://static.zybuluo.com/czc/3fn980h03uwj02bsejw1pp0i/image_1e7ftmnie1jnshge9nm57qqh48d.png
  [20]: http://static.zybuluo.com/czc/tc2233v6vc4w3lc6nz1rr1cq/image_1e7ftnb63b2ai7j1lko7ufvvr8q.png
  [21]: http://static.zybuluo.com/czc/rpqlnk5b8e1ypzfeas50q7b7/image_1e7ftnsojjre10ajmc41cibva797.png