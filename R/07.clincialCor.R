

#install.packages("beeswarm")

library(beeswarm)
setwd("C://Users//Desktop//TGIF1/")                    #修改工作目录
file="singleGeneCliData.txt"                                                       #输入文件
rt=read.table(file,sep="\t",header=T,check.names=F,row.names=1)                    #读取表达数据文件
gene="TGIF1"                                                                   #基因名字

#临床相关性分析，输出图形结果
for(clinical in colnames(rt[,4:ncol(rt)])){
      #定义颜色
	    xlabel=vector()
			tab1=table(rt[,clinical])
			labelNum=length(tab1)
			dotCol=c(2,3)
			if(labelNum==3){
				dotCol=c(2,3,4)
			}
			if(labelNum==4){
				dotCol=c(2,3,4,5)
			}
			if(labelNum>4){
				dotCol=rainbow(labelNum)
			}
			for(i in 1:labelNum){
			  xlabel=c(xlabel,names(tab1[i]) )
			}
	    
	    #相关性检验
	    i=gene
		  rt1=rbind(expression=rt[,i],clinical=rt[,clinical])
		  rt1=as.matrix(t(rt1))
		  if(labelNum==2){
		    rtTest<-wilcox.test(expression ~ clinical, data=rt1)
		  }else{
		    rtTest<-kruskal.test(expression ~ clinical, data = rt1)}
		  pValue=rtTest$p.value
		  pval=0
		  if(pValue<0.001){
			  pval="<0.001"
			}else{
			   pval=paste0("=",sprintf("%.03f",pValue))
		  }

      #可视化
      if(pValue<0.05){
					b = boxplot(expression ~ clinical, data = rt1,outline = FALSE, plot=F)
					yMin=min(b$stats)
					yMax = max(b$stats/5+b$stats)
					n = ncol(b$stats)
					outPdf=paste0(i,".",clinical,".pdf")
					width=ifelse(clinical=="Histology",14,7)
					pdf(file=outPdf,width = width,height = 5)
					par(mar = c(4.5,6,3,3))
					boxplot(expression ~ clinical, data = rt1,names=xlabel,
						     ylab = "Gene expression",main=paste0(i," (p",pval,")"),xlab=clinical,
						     cex.main=1.4, cex.lab=1.4, cex.axis=1.3,ylim=c(yMin,yMax),outline = FALSE)
				  beeswarm(expression ~ clinical, data = rt1, col =dotCol, lwd=0.1,
				         pch = 16, add = TRUE, corral="wrap")
				  dev.off()
		  }
}
