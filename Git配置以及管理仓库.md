# Git配置以及管理仓库

标签（空格分隔）： Github使用教程

---

# 一、注册github账号：

https://github.com/

# 二、Git安装

[下载git Windows版](https://github.com/git-for-windows/git/releases/download/v2.26.0.windows.1/Git-2.26.0-64-bit.exe)

# 三、配置Git

## 1、本地创建ssh key

```bash
ssh-keygen -t rsa -C "user.email@email.com"
```

## 2、生成结果：可以选择默认，在~/生成.ssh文件夹

```bash
Generating public/private rsa key pair.
Enter file in which to save the key (/c/Users/hp/.ssh/id_rsa):
/c/Users/hp/.ssh/id_rsa already exists.
Enter passphrase (empty for no passphrase): #输入github账户的密码
Enter same passphrase again:  #同上
Your identification has been saved in /c/Users/hp/.ssh/id_rsa.
Your public key has been saved in /c/Users/hp/.ssh/id_rsa.pub.  #public key 保存在/c/Users/hp/.ssh/id_rsa.pub

The key fingerprint is:
SHA256:tm+UfmlP5GQ2Jln2qmXM6DGq5WtrdGJZPaoSD0NPq8s user.email@email.com
The key's randomart image is:
```
## 3、然后，可以查看本地生成的公钥 

```bash
cat /c/Users/hp/.ssh/id_rsa.pub

ssh-rsa AAAAB3NzaC1yc2EAAAADAQABAAABAQC93oNWMKRtfEmuCE6KTf97JEDUYzqy3vuPOyuE7NGuImBuu797PdrU7rz7jGEVpjhtmuXN1woiKj7GKWNsVU0QRldMD5D5WkhunmFBajZnoYEK97Fkyt/ZHvpSlBJxzcqW1ToifzQU80+wOlATDTNG7/kE/8EnKMGhz7tIVI895r4/U7UuDnYz0EGOYFOPV0kFipUgqHs/U5LpmN/CLVbFjZccGy0CAyEEF534xKJl7aXEoxTTAqdTjwnjBFpUzPWZrUw6DpxRxIRD4Oy48ln44EJAGupliFP6tINdqhQUkwLhY95c22Y6x+BYfdBvgd5/sM7yiG5JVZ3vsHXr+NMV user.email@email.com
```
## 4、然后回到github上，进入Account Setting。

![image_1e59ogmjd11k4p2s1cnn1khb12141t.png-11.3kB][1]

## 5、配置github ssh公钥：title随便写，本地生成的公钥粘贴到空白区域。

![image_1e59ohlk9173prq51m5i1pon5vd2a.png-31.3kB][2]

## 6、验证是否成功：

```
ssh -T user.name@github.com
```

## 7、设置用户名和邮箱信息以方便上传时自动应用相应设置

```
git config --global user.name "user.name"
 git config --global user.email "user.email@email.com"
```

# 四、本地仓库-远程仓库

## 1、拉取远程仓库到本地
```git-bash
cd d:/Desktop/work
mkdir github
mkdir github/1.Linux
mkdir github/2.NGS
mkdir github/'3.Machine learning'
cd github/1.Linux
git clone https://github.com/user.name/Linux-study.git
cd Linux-study/
```
另外本地使用git来管理文件上传与存储交互
Git的本地仓库管理基本可以分为3个分区之间的传递

**本地文件**<--->**暂存区**<--->**本地仓库**

> * **本地文件**：自己刚刚修改，刚刚删除，刚刚创建的文件被称为本地文件。
> * **暂存区**：临时保存，用于上传或者等待其他文件进行统一上传的。
> * **本地仓库**：本地的github仓库文件，可以用于远程同步到github.com云端的内容

## 2、创建本地文件/修改本地文件/删除本地文件。
```git-bash
touch test.txt #创建本地文件
vim test.txt #使用Vim编辑器打开并进行修改
git rm test.txt #删除已经存在的文件

git status #可以看见Untracked files是红色的
```
## 3、将发生变化的文件添加到暂存区。

```git-bash
git add test.txt #添加发生变化的文件名，可以使用正则表达式、通配符等。
git status #可以看见添加后的文件变绿
```

## 4、将暂存区的文件保存到本地仓库。

```git-bash
git commit -m "change test,txt" #将暂存区的文件放入本地仓库，并标注改动信息“change test.txt”。
git status #暂存区的文件清空了。
```

## 5、将本地的文件同步到云端Github仓库

首先修改本地的仓库配置，使其能够使用自行记录登录账号和密码
```git-bash
vim d:/Desktop/work/github/1.Linux/Linux-study/.git/config   #修改配置文件，在默认链接方式中添加自己github的账号和密码，主要是为了方便上传同步。
```

#### **修改方式如下：**

1、**原始信息：**

![image_1e59niuu21iq1nq2ci67lnijn13.png-8.1kB][3]

2、**修改后的信息**

![image_1e59nkeb9rum1e1o1na1k8a1ipb1g.png-9.8kB][4]

3、**同步本地仓库到远程仓库**
```git-bash
git push #完活
```

----

参考文献：
【1】：https://www.bilibili.com/video/BV1Xx411m7kn?p=1


  [1]: http://static.zybuluo.com/czc/8o2s19hy8dry28d9hrw0pvhz/image_1e59ogmjd11k4p2s1cnn1khb12141t.png
  [2]: http://static.zybuluo.com/czc/uigvlnswte0ypizt8x2pd0hx/image_1e59ohlk9173prq51m5i1pon5vd2a.png
  [3]: http://static.zybuluo.com/czc/c063cmc8se5ydl31ejy89i7f/image_1e59niuu21iq1nq2ci67lnijn13.png
  [4]: http://static.zybuluo.com/czc/iipk9e9lwv0dbxy9qnpmv177/image_1e59nkeb9rum1e1o1na1k8a1ipb1g.png