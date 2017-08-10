setwd("D:\\data")
library(rjson)
fromJSON(json_str,file="clinical.project-TCGA-LUAD.2016-07-04T05-49-08.389307.json",method = "c",unexpected.escape = "error")
jsondata<-fromJSON(file="clinical.project-TCGA-LUAD.2016-07-04T05-49-08.389307.json",method = "C",unexpected.escape = "error")
write.csv(jsondata,file="cptcgaluad.csv",row.names=FALSE)
write.table(jsondata,file="cptcgaluad.txt",append=FALSE,quote=F,sep="",eol="\n",na="NA",dec=".",row.names=F,col.names=F,qmethod=c("escape","double"),fileEncoding = "" )
library(xlsx)
install.packages("mygene")
install.packages("")
library(rJava)
source("https://bioconductor.org/biocLite.R")
biocLite("mygene")
biocLite("GenomicFeatures")
library(ggplot2)
install.packages("mygene")
library(mygene)
install.packages("GenomicFeatures")
biocLite("mygene")
library(mygene)


gene<-getGene("1017",fields="all")
length(gene)
gene$name
gene$pathway
gene$taxid
gene$uniprot
gene$refseq
humangene<-getGenes(geneids="all",fields=c("entrezgene","name"),species="human")
human-gene<-getGenes(geneids="all",fields=c("symbol","name","entrezgene"),return.as = c("DataFrame"),species="human")

gene$go
GO_MF<-gene$go$MF
GO_MF[1:10]
write.csv(GO_MF)
write(GO_MF,file="D:\\data\\gomf",ncolumns=if(is.character(GO_MF)) 1 else 5,append=FALSE,sep=" ")
write(append=FALSE,sep=" ")

setwd("D:\\data")
genenames<-read.table("D:\\data\\Annotation_AllGenes_Clean.txt",colClasses=c("Character","NULL","NULL","NULL","NULL","NULL",header=TRUE))   #读入文件第一列数据
rm(list)
rm( list = ls ( all = TRUE)) 
gc()

for (i in 1:10) {
  gene[i]<-getGene("i",fields="all")
  write.table(gene[i],file="test1.txt")
}

#2016-7-12
#读取第一列数据并写出
gtmp<- read.table("Annotation_AllGenes_Clean.txt",header=TRUE)   #读入文件第一列数据
genenames<-gtmp[,1]
write.table(genenames,file="test1.txt")
colnames(genenames)<-c("Genes")

#getGene
testgene1<-getGene("1017",fields="all")
testgene1$name
testgene1$symbol
testgene1$entrezgene
GO_MF<-testgene1$go$GO_MF
GO_CC<-testgene1$go$GO_CC
GO_BP<-testgene1$go$GO_BP
KEGG_PATHWAY<-testgene1$kegg
REACTOME_PATHWAY<-testgene1$reactome

write.table(c("testgene1$name","GO_MF","GO_CC","GO_BP","KEGG_PATHWAY","REACTOME_PATHWAY"),file="test2.txt",append=FALSE,quote=TRUE,sep=" ",eol="\n",na="Unknow",dec=".",row.names=FALSE,col.names=FALSE,qmethod=c("escape","double"),fileEncoding="")

  
genesname<-read.table("hgenes.txt",header=F)
genesnameid<-getGenes("genesname",fields=c("symbol","name"))

genes<-getGenes()  

#2016-7-13
#########################
#创建进度条
pb<-txtProgressBar(min=0,max=59956,style = 1)
#设置进度条
Sys.sleep(0.0001)
setTxtProgressBar(pb,i)


genesname<-read.table("hgenes.txt",header=F)
i=1
for (i in 1:100){
 ghnames<- getGenes("1017",fields=c("name","entrezgene"))
 
}

write.table(genesname,file="ghnames.txt")
write.table(ghnames,file="ghnames.txt",append=F,sep=" ")

rm(list=ls())
i=1
for (i in 1:100){
ghnamestest1<- getGene("i",fields=c("name","entrezgene"))
}
#2016-7-14
############################
#biostar上给出的一个例子
library(mygene)
xli<-c('DDX26B','CCDC83',  'MAST3', 'RPL11', 'ZDHHC20',  'LUC7L3',  'SNORD49A',  'CTSH', 'ACOT8')
queryMany(xli, scopes="symbol", fields=c("uniprot", "ensembl.gene", "reporter"), species="human")


###########重新确定试##########
setwd("D:\\data")
gtmp<- read.table("Annotation_AllGenes_Clean.txt",header=TRUE)   #读入文件第一列数据
genenames<-gtmp[,1]
xli<-c(genenames[,1])
mmm<-queryMany(genenames,scopes="symbol",fields=c("ensembl.gene","entrezgene"),species="human")

write.table(mmm,file="tsts.txt")  #没有name信息。下一步将增加name信息在fields里。
#第四列是Gene_id 
di4lie<-mmm[,4]
write.table(di4lie,file="di4lie.txt")
llobs<-getGenes(di4lie,fields = c("name","GO","pathway"))  #写出来有pathway，没有GO的东西。
write.table(llobs,file="llobs.txt")
write.csv(llobs,file="llob.csv")# 写出来格式窜了，不能用这种格式
####查找GO格式
gene1017<-getGene("1017",fields="all")
gene1017$GO      #NULL
write.table(gene1017,file="gene1017.txt")
gene1017$go   
#这样子行，但会提示些提醒
#Warning messages:
#1: In .HTMLsearch(query) : Unrecognized search field: title
#2: In .HTMLsearch(query) : Unrecognized search field: keyword
#3: In .HTMLsearch(query) : Unrecognized search field: alias
###############new 提取#######
genetes1<-getGene("1017",fields=c("name","BP","CC","MF","kegg","reactome"))  #不行，后面的都不认识
genetes2<-getGene("1017",fields=c("name","go$BP","go$CC","go$MF","pathway$kegg","pathway$reactome")) #不行，后面的都不认识
genetes3<-getGene("1017",fields=c("name","go","pathway"))#这样的就有了，不过是list的。

akb1<-getGenes(di4lie,fields=c("name","go","pathway"))
write.table(akb1,file="akb1.txt")
genehgname<-c(akb1[,1-4,6,8-9])#出错 data[c(1,3)],所以如果改成genehgname2<-akb1[,c(1,2,3,4,6,8,9)]应该没问题
write.table(genehgname,file="genehnametest1.txt")#出错
genehgname<-akb1[,1:5]#不知道怎么调顺序。。。只能1到5
genehgname<-cbind(akb1[6],akb1[4],akb1[1],akb1[2],akb1[3],akb1[8],akb1[9])#可以，并且调整编号顺序，可以提前或置后相应的列
write.table(genehgname,file="genehgname2.txt",quote=F,row.names=F,sep="\t")###去除引号，行名，以制表符做分隔
genehgname2<-akb1[,c(6,4,1,2,3,8,9)]  #果然也可以，调编号顺序就能改列顺序






#############################7-16，将全部pathway信息都要######
genehgname3<-akb1[,c(6,4,1,3,2,8,9,7,10,11,12,13,14,15)] 
write.table(genehgname3,file="genehgname3.txt",quote=F,row.names=F,sep="\t")

#######练习######
x<-rnorm(20,mean=5,sd=2)
y<-rnorm(20,mean=10,sd=1)
plot(x,y)
barplot(x)
pie(y)
boxplot(x,y)
mx<-matrix(x,nrow=4,ncol=5)
stars(mx)
