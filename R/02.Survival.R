

library(survival)
library("survminer")

setwd("C://Users/Desktop//TGIF1")                #工作目录（需修改）
rt=read.table("singleGeneSurData.txt",header=T,sep="\t",check.names=F,row.names=1)     #读取输入文件
gene="TGIF1"                                                                   #基因名字

a=ifelse(rt[,gene]<=median(rt[,gene]),"low","high")
diff=survdiff(Surv(futime, fustat) ~a,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
fit=survfit(Surv(futime, fustat) ~ a, data = rt)
if(pValue<0.001){
		      pValue="<0.001"
		  }else{
		      pValue=paste0("=",round(pValue,3))
		  }
surPlot=ggsurvplot(fit, 
			           data=rt,
			           conf.int=TRUE,
			           pval=paste0("p",pValue),
			           pval.size=5,
			           risk.table=T,
			           legend.labs=c("high","low"),
			           legend.title=paste0(gene," level"),
			           xlab="Time(years)",
			           break.time.by = 1,
			           risk.table.title="",
			           palette=c("red", "blue"),
			           risk.table.height=.25)          
pdf(file=paste(gene,".survival.pdf",sep=""), width = 6, height = 5.5,onefile = FALSE)
print(surPlot)
dev.off()

summary(fit)                 #查看五年生存率