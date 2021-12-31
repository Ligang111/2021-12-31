
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("corrplot")
#install.packages("circlize")

setwd("C://Users//Desktop//TGIF1")        #设置工作目录
inputFile="normalize.txt"                                  #输入文件
gene="TGIF1"                                           #基因名字

#读取输入文件，对数据进行处理
library(limma)
rt=read.table(inputFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
rt=data[rowMeans(data)>0.5,]

#提取相关基因，计算基因间相关系数
gene=read.table("gene.txt",header=F)
data=t(rt[as.vector(gene[,]),])
cor1=cor(data)

#引用圈图可视化R包
library(circlize)
library(corrplot)
options(stringsAsFactors=F)

#设置图形颜色
pdf("circos.pdf")
col = c(rgb(1,0,0,seq(1,0,length=32)),rgb(0,1,0,seq(0,1,length=32)))
cor1[cor1==1]=0
par(mar=c(2,2,2,4))
c1 = ifelse(c(cor1)>=0,rgb(1,0,0,abs(cor1)),rgb(0,1,0,abs(cor1)))
col1 = matrix(c1,nc=ncol(data))

#绘制圈图
circos.par(gap.degree =c(3,rep(2, nrow(cor1)-1)),start.degree = 180)
chordDiagram(cor1,grid.col=rainbow(ncol(data)),transparency = 0.5,col=col1,symmetric = T)
par(xpd=T)
colorlegend(col, vertical = T,labels=c(1,0,-1),xlim=c(1.1,1.3),ylim=c(-0.4,0.4))       #绘制图例
dev.off()
circos.clear()
