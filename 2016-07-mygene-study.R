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
genenames<-read.table("D:\\data\\Annotation_AllGenes_Clean.txt",colClasses=c("Character","NULL","NULL","NULL","NULL","NULL",header=TRUE))   #�����ļ���һ������
rm(list)
rm( list = ls ( all = TRUE)) 
gc()

for (i in 1:10) {
  gene[i]<-getGene("i",fields="all")
  write.table(gene[i],file="test1.txt")
}

#2016-7-12
#��ȡ��һ�����ݲ�д��
gtmp<- read.table("Annotation_AllGenes_Clean.txt",header=TRUE)   #�����ļ���һ������
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
#����������
pb<-txtProgressBar(min=0,max=59956,style = 1)
#���ý�����
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
#biostar�ϸ�����һ������
library(mygene)
xli<-c('DDX26B','CCDC83',  'MAST3', 'RPL11', 'ZDHHC20',  'LUC7L3',  'SNORD49A',  'CTSH', 'ACOT8')
queryMany(xli, scopes="symbol", fields=c("uniprot", "ensembl.gene", "reporter"), species="human")


###########����ȷ����##########
setwd("D:\\data")
gtmp<- read.table("Annotation_AllGenes_Clean.txt",header=TRUE)   #�����ļ���һ������
genenames<-gtmp[,1]
xli<-c(genenames[,1])
mmm<-queryMany(genenames,scopes="symbol",fields=c("ensembl.gene","entrezgene"),species="human")

write.table(mmm,file="tsts.txt")  #û��name��Ϣ����һ��������name��Ϣ��fields�
#��������Gene_id 
di4lie<-mmm[,4]
write.table(di4lie,file="di4lie.txt")
llobs<-getGenes(di4lie,fields = c("name","GO","pathway"))  #д������pathway��û��GO�Ķ�����
write.table(llobs,file="llobs.txt")
write.csv(llobs,file="llob.csv")# д������ʽ���ˣ����������ָ�ʽ
####����GO��ʽ
gene1017<-getGene("1017",fields="all")
gene1017$GO      #NULL
write.table(gene1017,file="gene1017.txt")
gene1017$go   
#�������У�������ʾЩ����
#Warning messages:
#1: In .HTMLsearch(query) : Unrecognized search field: title
#2: In .HTMLsearch(query) : Unrecognized search field: keyword
#3: In .HTMLsearch(query) : Unrecognized search field: alias
###############new ��ȡ#######
genetes1<-getGene("1017",fields=c("name","BP","CC","MF","kegg","reactome"))  #���У�����Ķ�����ʶ
genetes2<-getGene("1017",fields=c("name","go$BP","go$CC","go$MF","pathway$kegg","pathway$reactome")) #���У�����Ķ�����ʶ
genetes3<-getGene("1017",fields=c("name","go","pathway"))#�����ľ����ˣ�������list�ġ�

akb1<-getGenes(di4lie,fields=c("name","go","pathway"))
write.table(akb1,file="akb1.txt")
genehgname<-c(akb1[,1-4,6,8-9])#���� data[c(1,3)],��������ĳ�genehgname2<-akb1[,c(1,2,3,4,6,8,9)]Ӧ��û����
write.table(genehgname,file="genehnametest1.txt")#����
genehgname<-akb1[,1:5]#��֪����ô��˳�򡣡���ֻ��1��5
genehgname<-cbind(akb1[6],akb1[4],akb1[1],akb1[2],akb1[3],akb1[8],akb1[9])#���ԣ����ҵ������˳�򣬿�����ǰ���ú���Ӧ����
write.table(genehgname,file="genehgname2.txt",quote=F,row.names=F,sep="\t")###ȥ�����ţ����������Ʊ������ָ�
genehgname2<-akb1[,c(6,4,1,2,3,8,9)]  #��ȻҲ���ԣ������˳����ܸ���˳��






#############################7-16����ȫ��pathway��Ϣ��Ҫ######
genehgname3<-akb1[,c(6,4,1,3,2,8,9,7,10,11,12,13,14,15)] 
write.table(genehgname3,file="genehgname3.txt",quote=F,row.names=F,sep="\t")

#######��ϰ######
x<-rnorm(20,mean=5,sd=2)
y<-rnorm(20,mean=10,sd=1)
plot(x,y)
barplot(x)
pie(y)
boxplot(x,y)
mx<-matrix(x,nrow=4,ncol=5)
stars(mx)