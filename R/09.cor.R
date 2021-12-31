

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

setwd("C://Users//Desktop//TGIF1/")        #设置工作目录
inputFile="normalize.txt"                               #输入文件
gene="TGIF1"                                        #基因或lncRNA名字
corFilter=0.4                                           #相关系数过滤值
pFilter=0.01                                           #统计学p值过滤值

library(limma)
rt=read.table(inputFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
rt=data[rowMeans(data)>0.5,]

x=as.numeric(rt[gene,])
gene1=unlist(strsplit(gene,"\\|",))[1]
outputFile=paste(gene1,".cor.xls",sep="")
outTab=data.frame()

for(j in rownames(rt)){
        y=as.numeric(rt[j,])
		    gene2=unlist(strsplit(j,"\\|",))[1]
		    corT=cor.test(x,y)

				z=lm(y~x)
				cor=corT$estimate
				cor=round(cor,3)
				pvalue=corT$p.value
				if(pvalue<0.001){
				  pval=signif(pvalue,4)
				  pval=format(pval, scientific = TRUE)
				}else{
				  pval=round(pvalue,3)}

        #输出相关性图片
				if((abs(cor)>corFilter) & (pvalue<pFilter)){
				  pdfFile=paste(gene1,"_",gene2,".cor.pdf",sep="")
					pdf(file=pdfFile,width =6,height = 6)
					plot(x,y, type="p",pch=16,col="blue",main=paste("Cor=",cor," (p-value=",pval,")",sep=""),
					    cex=1, cex.lab=1, cex.main=1,cex.axis=1,
					    xlab=paste(gene1,"expression"),
					    ylab=paste(gene2,"expression") )
					lines(x,fitted(z),col=2)
					dev.off()
					outTab=rbind(outTab,cbind(gene1,gene2,cor,pvalue))
				}
}
write.table(file=outputFile,outTab,sep="\t",quote=F,row.names=F)
