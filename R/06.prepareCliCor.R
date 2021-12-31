
setwd("C://Users//Desktop//TGIF1")                            #修改工作目录
expFile="rocSigExp.txt"                                                            #表达数据文件
clinicalFile="clinical.txt"                                                       #临床数据文件
gene="TGIF1"

exp=read.table(expFile,sep="\t",header=T,check.names=F,row.names=1)                #读取表达数据文件
cli=read.table(clinicalFile,sep="\t",header=T,check.names=F,row.names=1)           #读取临床数据文件
samSample=intersect(row.names(exp),row.names(cli))
exp=exp[samSample,]
cli=cli[samSample,]
selectCol=c("futime","fustat",gene)
outTab=cbind(exp[,selectCol],cli)
outTab=cbind(id=row.names(outTab),outTab)
write.table(outTab,file="singleGeneCliData.txt",sep="\t",row.names=F,quote=F)
