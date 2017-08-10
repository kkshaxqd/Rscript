#learn bioconductor
setwd("D:\\data\\17-06\\tex\\QC")

rm(list=ls())
source("http://www.bioconductor.org/biocLite.R")
biocLite("limma")
vignettes()
x <- 1:10
names(x) <- letters[1:10]
class(x)
x[1:3]
x
x[c("a", "b")]
.Machine$integer.max
x <- list(1:3, letters[1:3], is.numeric)
x
names(x) <- c("numbers", "letters", "function")
x[1:2]


###R markdown##########
install.packages("rmarkdown")
qdata <- read.table("result.txt")
install.packages("knitr")
install.packages("rmarkdown", repos="https://cloud.r-project.org")
install.packages("ggplot2")
install.packages('DT')
install.packages("gridExtra")
library(DT)
datatable(head(iris))

#test1-test
library("ggplot2")
qdata <- read.table("N170057.GeneRegionMutationCount.txt")
pdata<-data.frame(qdata$V2,qdata$V1)
colnames(pdata) <- c("GeneFuction","Number")
pdata<-pdata[order(pdata$Number,decreasing = TRUE, na.last = NA),]
pdata$GeneFuction <-factor(pdata$GeneFuction,levels=pdata$GeneFuction  )
ggplot(pdata)+geom_bar(aes(x=GeneFuction,y=Number,fill=factor(GeneFuction)),stat= 'identity', position = 'stack')+theme_bw()+theme(panel.grid=element_blank())+labs(title="GeneFuction")+theme(plot.title = element_text(hjust = 0.5))+theme(axis.text.x = element_text(size = 10,color = "black", face = "bold", vjust = 1, hjust = 1,angle=90),plot.margin=unit(c(0.2,0.2,0.2,0.2),"inches"))+geom_text(aes(x=GeneFuction,y=Number,label=Number),position=position_dodge(1.2),vjust = 0)

pdata$GeneFuction <-factor(pdata$GeneFuction,levels=pdata$GeneFuction  )   #  ggplot2输出图默认是按因子排序的，这个步骤使输出图按因子来输出
myLabel = as.vector(pdata$GeneFuction)   ## 转成向量，否则图例的标签可能与实际顺序不一致
myLabel = paste(myLabel, "(", round(pdata$Number / sum(pdata$Number) * 100, 2), "%)", sep = "")   
ggplot(pdata)+geom_bar(aes(x=GeneFuction,y=Number,fill=factor(GeneFuction)),stat= 'identity', position = 'stack')+theme_bw()+
  theme(panel.grid=element_blank())+labs(x="",title="Gene Fuction")+theme(plot.title = element_text(hjust = 0.5),legend.position="top",legend.title = element_blank())+
  scale_fill_discrete(breaks = pdata$GeneFuction, labels = myLabel)+
  theme(axis.text.x = element_text(size = 10,color = "black", face = "bold",angle=90))+
  ylim(0,max(pdata$Number)+2000)+
  geom_text(aes(x=GeneFuction,y=Number,label=Number),position=position_dodge(0.8),vjust = 0)




#dt$obj = factor(dt$obj, levels=c('D','B','C','A','E')) ## 设置柱条的顺序
Pdata <- read.table("N170044.SNVtypeCount.txt")

ldata<-read.table("N170057.average_Coverage.txt",sep="\t")
colnames(ldata) <- c("DifferentTargetCoverage","Ratio")
rname<-c(">=1X",">=5X",">=10X",">=20X",">=50X",">=100X",">=200X")
ldata[,1]<-rname
ldata<-ldata[order(ldata$Ratio, decreasing = TRUE, na.last = NA),]
#ldata$DifferentTargetCoverage<-factor(ldata$DifferentTargetCoverage,levels=rname)
myLabel = as.vector(ldata$DifferentTargetCoverage)   ## 转成向量，否则图例的标签可能与实际顺序不一致
myLabel = paste(myLabel, "(", round((ldata$Ratio) * 100, 2), "%) ", sep = "")   ## 用 round() 对结果保留两位小数

ggplot(ldata)+geom_bar(aes(x="",y=Ratio,fill=factor(DifferentTargetCoverage)),stat= 'identity', position = 'stack',width = 0.2)+coord_polar(theta = "y") +
  theme_bw()+theme(panel.grid=element_blank())+theme(panel.border=element_blank()) +
  labs(title="Fraction of target covered")+theme(plot.title = element_text(hjust = 0.5),legend.title = element_blank(),legend.position="top")+
  scale_fill_discrete(breaks = ldata$DifferentTargetCoverage, labels = myLabel)   ## 将原来的图例标签换成现在的myLabel
####################晕####也许。。。这个做成柱形图更好看一些。。
colnames(ldata) <- c("DifferentTargetCoverage","Ratio")
rname<-c(">=1X",">=5X",">=10X",">=20X",">=50X",">=100X",">=200X")
ldata[,1]<-rname
#ldata$DifferentTargetCoverage= as.vector(ldata$DifferentTargetCoverage) 
ldata$Ratio=paste(round((ldata$Ratio) * 100, 2),"%",sep="")
ldata<-ldata[order(ldata$Ratio,decreasing = TRUE, na.last = NA),]
ldata$DifferentTargetCoverage<-factor(ldata$DifferentTargetCoverage,levels=ldata$DifferentTargetCoverage )

ggplot(ldata)+geom_bar(aes(x=DifferentTargetCoverage,y=Ratio,fill=factor(DifferentTargetCoverage)),stat= 'identity', position = 'stack')+
  theme_bw()+theme(panel.grid=element_blank())+labs(x="",title="Fraction of target covered")+theme(plot.title = element_text(hjust = 0.5),legend.title = element_blank(),legend.position="top")+
  theme(plot.margin=unit(c(0.2,0.2,0.2,0.2),"inches"))+
  geom_text(aes(x=DifferentTargetCoverage,y=Ratio,label=Ratio),position=position_dodge(1.2),vjust = 0)


#######################################
mdata <- read.table("N170057.SNVtypeCount.txt",sep="\t")
ndata<-data.frame(mdata$V2,mdata$V1)
colnames(ndata) <- c("ExonicFunction","Number")
ndata<-ndata[order(ndata$Number, decreasing=TRUE, na.last = NA),]  ## 用 order() 让数据框的数据按 Number 列数据从大到小排序,order返回的是数据大小的位置

#ndata<-sort(ndata)
knitr::kable(
  ndata, format="html",table.attr = "class=\"table table-bordered\"", caption = 'SNV突变类型质控结果表'
)
ndata$ExonicFunction<-factor(ndata$ExonicFunction,levels=ndata$ExonicFunction )
myLabel = as.vector(ndata$ExonicFunction)   ## 转成向量，否则图例的标签可能与实际顺序不一致
myLabel = paste(myLabel, "(", round(ndata$Number / sum(ndata$Number) * 100, 2), "%)", sep = "")   
##环图
ggplot(ndata)+geom_bar(aes(x="",y=Number,fill=factor(ExonicFunction)),stat= 'identity', position = 'stack',width = 0.2)+
  coord_polar(theta = "y")+
  theme_bw()+theme(panel.grid=element_blank())+labs(x="",y="",title="Exonic Function")+
  theme(plot.title = element_text(hjust = 0.5),legend.position="top",legend.title = element_blank())+
  scale_fill_discrete(breaks = ndata$ExonicFunction, labels = myLabel)+
  theme(axis.text.x = element_blank())+
  geom_text(aes(y = (ndata$Number/sum(ndata$Number)) , x = ndata$Number/sum(ndata$Number), label = myLabel), size = 5)   ## 在图中加上百分比：x 调节标签到圆心的距离, y 调节标签的左右位置

##柱状图
ggplot(ndata)+geom_bar(aes(x=ExonicFunction,y=Number,fill=factor(ExonicFunction)),stat= 'identity', position = 'stack')+
  theme_bw()+theme(panel.grid=element_blank())+labs(x="",title="ExonicFunction")+
  theme(plot.title = element_text(hjust = 0.5),legend.position="top",legend.title = element_blank())+
  scale_fill_discrete(breaks = ndata$ExonicFunction, labels = myLabel)+
  geom_text(aes(x=ExonicFunction,y=Number,label=Number),position=position_dodge(1.2),vjust = 0)+
  theme(axis.text.x = element_text(size = 10,color = "black", face="bold",angle=90))


#############表1测试##
sampleid="N170057"

hdata <- read.table("N170057.reads_bases_GC.txt",sep=":")
hdata[1,2]
colnames(hdata) <- c("Sample",sampleid)
rname=c("Total reads","Total bases","GC content(%)")
hdata$Sample<-rname
hdata[1,2]<-round(hdata[1,2],0)   #ceiling(x)
hdata[2,2]<-round(hdata[2,2],0)
hdata[3,2]<-round(hdata[3,2],2)  #signif(x, digits = 6)
#读入Q20 Q30
hsdata<-read.table("result_Q20Q30.txt",sep = "\t",fill = TRUE)
q20<-as.numeric(as.character(hsdata[3,6])) 
q30<-as.numeric(as.character(hsdata[3,7]))
newhang1<-c("Q20",q20)  #原来这样就行
newhang2<-c("Q30",q30)
hdata<-rbind(hdata,newhang1)
hdata<-rbind(hdata,newhang2)
# 读入捕获特异性
hhdata<-read.table("Capture_specificity.txt")
ebot<-round(hhdata[1,1],0)                         #Effective bases on target
teb<-round(hhdata[2,1],0)              #Total effective bases
cs<-round((hhdata[3,1])*100,2)                #Capture specificity (%)
newhang3<-c("Effective bases on target",ebot)
newhang4<-c("Total effective bases",teb)
newhang5<-c("Capture specificity (%)",cs)
hdata<-rbind(hdata,newhang3)
hdata<-rbind(hdata,newhang4)
hdata<-rbind(hdata,newhang5)
#读入平均深度
hjdata<-read.table("N170057.average_Depth.txt",sep="\t")
asdot<-round(hjdata$V2,2)                        #Average sequencing depth on target
newhang6<-c("Average sequencing depth on target",asdot)
hdata<-rbind(hdata,newhang6)
#读入覆盖度情况
hldata<-read.table("N170057.average_Coverage.txt",sep="\t")
hldata$V1<-paste(hldata$V1,"(%)",sep="" )
hldata$V2<-round((hldata$V2)*100,2)
colnames(hldata) <- c("Sample",sampleid)
hdata<-rbind(hdata,hldata)


library(ggplot2)
dt = data.frame(A = c(2, 7, 4, 10, 1), B = c('B','A','C','D','E'))

myLabel = as.vector(dt$B)   ## 转成向量，否则图例的标签可能与实际顺序不一致
myLabel = paste(myLabel, "(", round(dt$A / sum(dt$A) * 100, 2), "%)        ", sep = "")   ## 用 round() 对结果保留两位小数

p = ggplot(dt, aes(x = "", y = A, fill = B)) + 
  geom_bar(stat = "identity", width = 1) + 
  geom_rect(colour = "grey30", show.legend = FALSE) +
  coord_polar(theta = "y") + 
  labs(x = "", y = "", title = "") + 
  theme(axis.ticks = element_blank()) + 
  theme(legend.title = element_blank(), legend.position = "top") + 
  scale_fill_discrete(breaks = dt$B, labels = myLabel)   ## 将原来的图例标签换成现在的myLabel
p
####################

##```{r eval=TRUE}

#qdata <- read.table("result.txt")

#library(DT)

#DT::datatable(qdata, class='cell-border stripe')
###```

mydata <- data.frame(a=1:50, b=rnorm(50))

mytable <- cbind(sites=c("site 1","site 2","site 3","site 4"),mydata[10:13,])

rownames(mytable) <- c("a","b","c","d")  重命名

#```{r }
qdata <- read.table("result.txt")
library(DT)
library(ggplot2)
####表格
datatable(qdata)
###3绘图
ggplot(qdata)+geom_bar(aes(x=V2,y=V1,fill=factor(V2)),stat= 'identity', position = 'stack')+theme_bw()+labs(title="Quaity result")+theme(axis.text.x = element_text(size = 8,color = "black", face = "bold", vjust = 1, hjust = 1, angle = 90))+geom_text(aes(x =V2,y=V1,label=V1),position=position_dodge(.9))

#```
#############学习CNVPanelizer
source("https://bioconductor.org/biocLite.R")
biocLite("CNVPanelizer")
biocLite("GenomicRanges")
library(devtools)
install_github("biostuff/CNVPanelizer")
biocLite("RCurl")
biocLite("Rsamtools")
biocLite("Biobase")

library(CNVPanelizer)
?CNVPanelizer
#example
data(sampleReadCounts)
data(referenceReadCounts)
## Gene names should be same size as row columns
geneNames <- row.names(referenceReadCounts)
ampliconNames <- NULL
normalizedReadCounts <- CombinedNormalizedCounts(sampleReadCounts,
                                                 referenceReadCounts,
                                                 ampliconNames = ampliconNames)
# After normalization data sets need to be splitted again to perform bootstrap
samplesNormalizedReadCounts = normalizedReadCounts["samples"][[1]]
referenceNormalizedReadCounts = normalizedReadCounts["reference"][[1]]
#Values above 10000 should be used
replicates <- 10
# Perform the bootstrap based analysis
bootList <- BootList(geneNames,
                     samplesNormalizedReadCounts,
                     referenceNormalizedReadCounts,
                     replicates = replicates)
background <- Background(geneNames,
                         samplesNormalizedReadCounts,
                         referenceNormalizedReadCounts,
                         bootList,
                         replicates = replicates,
                         significanceLevel = 0.1)

#BedToGenomicRanges BedToGenomicRanges
bedFilepath <- file.path("D:\\myscript\\mydatabase\\03_2742.bed")
ampliconColumn <- 4
genomicRangesFromBed <- BedToGenomicRanges(bedFilepath, ampliconColumn)

genomicRangesFromBed[1]


backgroundNoise <- Background(geneNames,
                              samplesNormalizedReadCounts,
                              referenceNormalizedReadCounts,
                              bootList,
                              replicates = replicates)


reportTables <- ReportTables(geneNames,
                             samplesNormalizedReadCounts,
                             referenceNormalizedReadCounts,
                             bootList,
                             backgroundNoise)
PlotBootstrapDistributions(bootList, reportTables, save = FALSE)



#####shiny学习33
install.packages("shiny")

#####shiny支持中文办法##
install.packages("devtools")
devtools::install_github('yihui/shiny@bugfix/native-encoding')  #用以下语句对shiny包打上其他语种的编码补丁错误，目录被换了
devtools::install_github('rstudio/shiny')
devtools::install_github('rstudio/htmltools') #想解决DT包中文乱码问题，未安装成功。可能要R.3.3.3 #莫名其妙安装成功了
#########rstudio生成pdf中中文显示问题

#@1.需要安装Rstudio的最新版本

#2.安装rticles包

#library(devtools)

#devtools::install_github("rstudio/rticles")

#3.新建rmarkdown文件-From Template--Ctex Document

#这里的ctex模版是文章模版，如果需要presentation格式，需要自己修改下latex模版。

但报错了

install.packages("tufte", type = "source") #试试这个模板

#####
library(shiny)
runExample("01_hello")
runExample("02_text")


