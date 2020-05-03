# ggplot2 绘制散点图

标签（空格分隔）： 绘图

---


大多数人对ggplot2真的是又爱又恨(我大概就是`掘墓人`吧)：爱是因为ggplot2没有啥不能绘制的，恨是因为参数真特么的多，多到无法计数，而且还需要考虑图层之间的顺序以及关系。我这里做一个不断更新的散点图教程，将不断的收集美图进行完善。

# 一、基本绘图

1、获取数据

```
require(ggplot2)   # 
data(diamonds)
set.seed(42)
small <- diamonds[sample(nrow(diamonds), 1000), ]
head(small)
summary(small)
```

![image_1e72j5f8e1j6e19mc1sjdr961krj9.png-20.3kB][1]


## 2、查看数据统计信息

```
summary(small)
```

![image_1e72j67liv9m1ih0l062s1820m.png-37.5kB][2]


## 3、基础绘图
使用克拉作为横坐标；价格作为纵坐标，以散点图的形式展示数据结果
``` 
p1.3 <- ggplot(data = small, mapping = aes(x = carat, y = price)) #输入数据
p1.3 =p1.3 + geom_point() #以散点图的形式展示
p1.3
```
![image_1e72j6ua91sgfp7610bcg1acd713.png-31.2kB][3]


# 二、添加分组
 
## 1、用颜色进行区分分组
根据钻石的颜色来对数据进行区分

``` 
p2.1 <- ggplot(data = small)
p2.1<- p2.1 + geom_point(mapping = aes(x = carat, y = price,color=color)) 
p2.1
#这里是按照cut进行分组，颜色是默认颜色顺序
```

![image_1e72j7kud1tnr8b9l8hf3v1jvp1g.png-44kB][4]


#### 1.1、自定义颜色分组
紧接着p1,添加自己的颜色代码，然后使用颜色代码进行自定义颜色
``` 
mycol <- c("#ca1414","#0baa4f","#c3d316","#FFD121","#11AA4D","#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767")
p2.1.1 <- p2.1 + scale_color_manual(values = mycol[1:length(unique(small$color))])
p2.1.1

```

![image_1e72j865antqgkjbfe1okk9l81t.png-44.1kB][5]


### 2、颜色+形状实现多组区分
使用颜色与形状实现分组展示多组信息
``` 
mycol <- c("#ca1414","#0baa4f","#c3d316","#FFD121","#11AA4D","#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767")
p2.2 <- ggplot(data = small)
p2.2<- p2.2 + geom_point(mapping = aes(x = carat, y = price,color=color, shape=cut)) 
p2.2 <- p2.2 + scale_color_manual(values = mycol[1:length(unique(small$color))])
p2.2
```


![image_1e72j8oi11k5l199qmm430r1bkj2a.png-47.8kB][6]

#### 2.1、自定义散点形状
紧接着p2，添加自己的散点形状代码，从  `0 ~ 25`
详细散点[对照表](http://www.sthda.com/sthda/RDoc/images/points-symbols.png)
``` 
p2.2.1 <- p2.2 + scale_shape_manual(values = c(1:5))  #符号形式是根据数字来定义选择的
p2.2.1
```


![image_1e72j9639j3a7l15vdps71lhr2n.png-55.8kB][7]



# 三、调整散点的大小与透明度

## 1、调整散点的大小

``` 
p3.1 <- ggplot(data = small)
p3.1 <- p3.1 + geom_point(mapping = aes(x = carat, y = price,color=color,size=5)) 
p3.1
```

![image_1e72j9l01nqptq31ao41unkii34.png-50.2kB][8]



## 1.1、连续变量下的散点大小修改

``` 
p3.1.1 <- ggplot(data = small)
p3.1.1 <- p3.1.1 + geom_point(mapping = aes(x = depth, y = price,color=color,size=price)) 
p3.1.1
```

![image_1e72ja32v1p5515ps1q8sgggibf3h.png-64.7kB][9]


## 1.2、分组情况下的自定义散点大小修改

``` 
p3.1.2 <- ggplot(data = small)
p3.1.2 <- p3.1.2 + geom_point(mapping = aes(x = depth, y = price,color=color,size=color))
p3.1.2 <- p3.1.2 + scale_size_manual(values = c(7:1))
p3.1.2
```

![image_1e72jai7n7amoqurjm1lga1kle3u.png-60.3kB][10]


## 2、调整散点的透明度

``` 
p3.2 <- ggplot(data = small)
p3.2 <- p3.2 + geom_point(mapping = aes(x = carat, y = price),alpha=0.1) 
p3.2
#同理，通过scale_alpha_manual可以自定义透明度，但是一般没有实际意义，透明度统一调整才能通过颜色深浅来判断数据的密集程度。
```

![image_1e72jb4bm150d10m31s7i7j3omk4b.png-33.4kB][11]


  [1]: http://static.zybuluo.com/czc/e2rapq4iarggtac06ozwem4s/image_1e72j5f8e1j6e19mc1sjdr961krj9.png
  [2]: http://static.zybuluo.com/czc/rhx1whzvidm1wn64x0uvxpsd/image_1e72j67liv9m1ih0l062s1820m.png
  [3]: http://static.zybuluo.com/czc/ltgw8rbpyyig5q64wxv6bcn8/image_1e72j6ua91sgfp7610bcg1acd713.png
  [4]: http://static.zybuluo.com/czc/6q5wpzrhntm3n1u93p6g17z7/image_1e72j7kud1tnr8b9l8hf3v1jvp1g.png
  [5]: http://static.zybuluo.com/czc/p02ygvc54tmdnaa0pv86da5m/image_1e72j865antqgkjbfe1okk9l81t.png
  [6]: http://static.zybuluo.com/czc/7uatyt3xizeyar6rpl34t9ud/image_1e72j8oi11k5l199qmm430r1bkj2a.png
  [7]: http://static.zybuluo.com/czc/lzfoz37g3d6vo48ysnb0dn66/image_1e72j9639j3a7l15vdps71lhr2n.png
  [8]: http://static.zybuluo.com/czc/bo1ynqi2av6fau6e4rshadlp/image_1e72j9l01nqptq31ao41unkii34.png
  [9]: http://static.zybuluo.com/czc/wym9tfuqj7y84dktucg1wult/image_1e72ja32v1p5515ps1q8sgggibf3h.png
  [10]: http://static.zybuluo.com/czc/94qjas7wfyzag8cbxspgxkgd/image_1e72jai7n7amoqurjm1lga1kle3u.png
  [11]: http://static.zybuluo.com/czc/z8ysocpenxupvta36p7h0izd/image_1e72jb4bm150d10m31s7i7j3omk4b.png