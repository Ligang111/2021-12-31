
#install.packages("beeswarm")

setwd("C://Users//Administrator//Desktop//家林 处理又出的 6月29号//TGIF1//4.GEO数据库验证//4.diff//GSE50161 好")         #工作目录
inputFile="singleGene.txt"                                      #输入文件
yMin=4                      #纵轴最小值
yMax=14                     #纵轴最大值
ySeg=yMax*0.94

library(limma)
library(beeswarm)

rt=read.table(inputFile,sep="\t",header=T,row.names=1,check.names=F)
geneName=colnames(rt)[1]
labels=c("Normal","Tumor")
colnames(rt)=c("expression","Type")
wilcoxTest<-wilcox.test(expression ~ Type, data=rt)
wilcoxP=wilcoxTest$p.value
pvalue=signif(wilcoxP,4)
pval=0
if(pvalue<0.001){
     pval=signif(pvalue,4)
     pval=format(pval, scientific = TRUE)
}else{
     pval=round(pvalue,3)
}

outFile=paste(geneName,".pdf",sep="")
pdf(file=outFile,width=7,height=5)
par(mar = c(4,7,3,3))
boxplot(expression ~ Type, data = rt,names=labels,
     ylab = paste(geneName," expression",sep=""),
     cex.main=1.5, cex.lab=1.3, cex.axis=1.2,ylim=c(yMin,yMax),outline = FALSE)
beeswarm(expression ~ Type, data = rt, col = c("blue","red"),lwd=0.1,
     pch = 16, add = TRUE, corral="wrap")
segments(1,ySeg, 2,ySeg);segments(1,ySeg, 1,ySeg*0.96);segments(2,ySeg, 2,ySeg*0.96)
text(1.5,ySeg*1.05,labels=paste("p=",pval,sep=""),cex=1.2)
dev.off()

