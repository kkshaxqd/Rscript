#FFPE 谢主任12对肿瘤样本分析  R语言处理

#1  KEGG   首先是ID转换 

if("org.Hs.eg.db" %in% rownames(installed.packages()) == FALSE) {source("http://bioconductor.org/biocLite.R");biocLite("org.Hs.eg.db")}
suppressMessages(library(org.Hs.eg.db))  #我比较喜欢这样加载包
install.packages("annotate") #NO use


source("https://bioconductor.org/biocLite.R")
biocLite("clusterProfiler")

source("https://bioconductor.org/biocLite.R")
biocLite("annotate")    
library(annotate) #一般都是这样加载包   #OK

temp<-gene_data2
write.table(temp,file="D:/data/17-08/xiezhurenSNPINDEL/18-geneshow-great6.txt",quote=F,sep="\t",eol="\n",na="NA",row.names = F,col.names = T)
tempG<-read.table("D:/data/17-08/xiezhurenSNPINDEL/18-geneshow-great6.txt",header=F)



#fangfa1
library(org.Hs.eg.db)
symbols <- org.Hs.egSYMBOL[as.character(temp)]  #many cannot find

#way2
library(annotate)
geneID<-getSYMBOL(as.character(temp), data='org.Hs.eg')

#way3
source("https://bioconductor.org/biocLite.R")
biocLite("mygene")
library(mygene)
?mygene
a<-queryMany(tempG, scopes=c("symbol", "reporter", "accession"), fields="all", species='human') #get all info of these genes
b<-a$entrezgene
write.table(b,file="D:/data/17-08/xiezhurenSNPINDEL/18-geneshow-id.txt",quote=F,sep="\t",eol="\n",na="NA",row.names = F,col.names = T)

a<-queryMany(temp, scopes=c("symbol", "reporter", "accession"), fields=c("go","pathway"), species='human') 

write.table(a,file="D:/data/17-08/xiezhurenSNPINDEL/18-geneshow-all-info.txt",quote=F,sep="\t",eol="\n",na="NA",row.names = F,col.names = T)


#GO PATHWAY  CLUSTERPROFILER
require(DOSE)
install.packages("DOSE") #NO USE
source("https://bioconductor.org/biocLite.R")
biocLite("DOSE")
library(DOSE)
library(clusterProfiler)
install.packages("ggplot2")
remove.packages(c("ggplot2", "data.table"))
remove.packages("S4Vectors")
install.packages('ggplot2', dep = TRUE)
install.packages('data.table', dep = TRUE)
install.packages('S4Vectors', dep = TRUE)
gene=as.character(temp)

remove.packages(c("ggplot2", "data.table","DOSE","clusterProfiler"))

install.packages("ggplot2")  #only 2.1 install

install.packages("devtools")
devtools::install_github("hadley/ggplot2")  #install 2.2

source("https://bioconductor.org/biocLite.R")
biocLite("clusterProfiler")
biocLite()  #升级

remove.packages("scales")
install.packages("scales",dep=T)
biocLite("org.Hs.eg.db") 

library(DOSE)
library(clusterProfiler)
library(org.Hs.eg.db)

ego <- enrichGO(gene=temp,'org.Hs.eg.db',ont="MF",pvalueCutoff=0.01,readable=TRUE)

ekk <- enrichKEGG(gene=temp,pvalueCutoff=0.01)

write.csv(summary(ekk),"KEGG-enrich.csv",row.names =F)

write.csv(summary(ego),"GO-enrich.csv",row.names =F)

remove.packages("BiocInstaller")
source("https://bioconductor.org/biocLite.R")
?biocLite

data(geneList)
de <- names(geneList)[1:100]

temp<-read.table("D:/data/17-08/xiezhurenSNPINDEL/18-geneshow-id.txt",header=F)

enrichDAVID(temp, idType = "ENTREZ_GENE_ID", listType = "Gene",
            minGSSize = 10, maxGSSize = 500, annotation = "GOTERM_BP_FAT",
            pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2,
            species = NA, david.user)
biocLite("RDAVIDWebService") 
library(rJava)
install.packages("rJava")

out18path<-a$pathway.kegg


#109gene
tempG109<-read.table("D:/data/17-08/xiezhurenSNPINDEL/109gene.txt",header=F)

a<-queryMany(tempG109, scopes=c("symbol", "reporter", "accession"),  species='human') #get all info of these genes
name<-a$`_id`
write.table(name,file="D:/data/17-08/xiezhurenSNPINDEL/109gene-id.txt",quote=F,sep="\t",eol="\n",na="NA",row.names = F,col.names = T)

#ego <- enrichGO(gene=tempG109,'org.Hs.eg.db',ont="CC",pvalueCutoff=0.01,readable=TRUE)
ego_MF <- enrichGO(gene =tempG109,
                   OrgDb=org.Hs.eg.db,
                   ont = "BP",
                   pAdjustMethod = "MF",
                   minGSSize = 1,
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.01,
                   readable = TRUE)

kk <- enrichKEGG(gene = tempG109, organism ="human",
                 pvalueCutoff = 0.01,
                 qvalueCutoff = 0.01,
                 minGSSize = 1,
                 #readable = TRUE ,
                 use_internal_data =FALSE)

write.csv(summary(ekk),"KEGG-enrich.csv",row.names =F)

write.csv(summary(ego),"GO-enrich.csv",row.names =F)



#差异分析
m<-read.table("clipboard", header = T, sep = '\t') 

chayi.pop<-matrix(data=c(149,115,1139,1009),nrow=2)
fisher.test(chayi.pop,alternative="greater")


##wego注释

#构建wego识别文件
GENEdata<-read.table("D:/data/17-08/xiezhurenSNPINDEL/zuihougene447.txt",header=T,sep='\t')
ycmgene<-read.table("clipboard",header=F,sep='\t')
ycmdfgene<-ycmgene
wxhdfgene<-data.frame(GENEdata[GENEdata$WXH.DF==1,]$GENE)
dljdfgene<-data.frame(GENEdata[GENEdata$DLJ.DF ==1,]$GENE)
ycjdfgene<-data.frame(GENEdata[GENEdata$YCJ.DF ==1,]$GENE)
lqqdfgene<-data.frame(GENEdata[GENEdata$LQQ.DF==1,]$GENE)
fjydfgene<-data.frame(GENEdata[GENEdata$FJY.DF==1,]$GENE)
ctxdfgene<-data.frame(GENEdata[GENEdata$CTX.DF==1,]$GENE)
cxyffgene<-data.frame(GENEdata[GENEdata$CXY.FF==1,]$GENE)
nlffgene<-data.frame(GENEdata[GENEdata$NL.FF==1,]$GENE)
zllffgene<-data.frame(GENEdata[GENEdata$ZLL.FF==1,]$GENE)
zcyffgene<-data.frame(GENEdata[GENEdata$ZCY.FF==1,]$GENE)
cjlffgene<-data.frame(GENEdata[GENEdata$CJL.FF==1,]$GENE)



