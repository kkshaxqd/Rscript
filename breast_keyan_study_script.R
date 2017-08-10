#Rscriptforbreastcancerstudy
library(ggplot2)
library(gridExtra)  #用来把多图放到一起
rm(p2)
setwd("D:\\data\\17-05")   

peiduidata<-read.table(file="peiduidata-chuli.txt",header=T,sep="\t")

#genecoord<-paste(peiduidata$Chr,peiduidata[,2],peiduidata[,3],peiduidata[,4],peiduidata[,5],sep=":", collapse = NULL)

#newdata<-cbind(genecoord,peiduidata[,6],peiduidata[,7],peiduidata[,9],peiduidata[,10])

tiff("testbreastcancer3peiduiananly.png",width = 2160, height = 1080)

P1<-ggplot(peiduidata)+geom_bar(aes(x =SAMPLE,y=Chr.Start.End.Ref.Alt,fill=factor(ExonicFunc.refGene)),stat="identity")+theme_bw() +theme_classic()+ labs(title="SOMATIC SNV RESILT")+theme(axis.text.y = element_text(size = 8,color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 90))+coord_flip() 
P2<-ggplot(peiduidata)+geom_bar(aes(x =SAMPLE,y=Gene.refGene,fill=factor(ExonicFunc.refGene)),stat="identity")+theme_bw() +theme_classic()+ labs(title="SOMATIC SNV RESILT")+theme(axis.text.y = element_text(size = 8,color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 90))+coord_flip() 
grid.arrange(P1, P2, ncol=2)
dev.off()



breast12somticmutationdata<-read.table("clipboard", header = T, sep = '\t')
bk<-summary(breast12somticmutationdata)
write.table(bk,file="summaryresults.txt",quote=F,sep="\t",eol="\n",na="NA",row.names = F,col.names = T)
#####给12样本的体细胞突变基因绘制图
data3<-read.table("clipboard", header = T, sep = '\t')
data3<-na.omit(data3) ##删除缺失行
tiff("12samplesomaticmutationgenes2.png",width = 2160, height = 1080)
ggplot(data3)+geom_bar(aes(x=GENE,y=COUT,fill=factor(SAMPLE)),stat= 'identity', position = 'stack')+theme_bw() +labs(title="12 SAMPLE SOMATIC MUTATION GENES RESILT")+theme(axis.text.x = element_text(size = 8,color = "black", face = "bold", vjust = 1, hjust = 1, angle = 90))
dev.off()



library(gridExtra)  #用来把多图放到一起

data31<-data3[1:600,]
data32<-data3[601:1200,]
data33<-data3[1201:1800,]
data34<-data3[1801:2400,]
tiff("12samplesomaticmutationgenes3.png",width = 2160, height = 2160)
p1<-ggplot(data31)+geom_bar(aes(x=GENE,y=COUT,fill=factor(SAMPLE)),stat= 'identity', position = 'stack')+theme_bw() +labs(title="12 SAMPLE SOMATIC MUTATION GENES RESILT1")+theme(axis.text.x = element_text(size = 8,color = "black", face = "bold", vjust = 0.5, hjust = 1, angle = 90))+ theme(panel.grid =element_blank())+ ylim(limits=c(0,12))
p2<-ggplot(data32)+geom_bar(aes(x=GENE,y=COUT,fill=factor(SAMPLE)),stat= 'identity', position = 'stack')+theme_bw() +labs(title="12 SAMPLE SOMATIC MUTATION GENES RESILT2")+theme(axis.text.x = element_text(size = 8,color = "black", face = "bold", vjust = 0.5, hjust = 1, angle = 90))+ theme(panel.grid =element_blank())+ ylim(limits=c(0,12))
p3<-ggplot(data33)+geom_bar(aes(x=GENE,y=COUT,fill=factor(SAMPLE)),stat= 'identity', position = 'stack')+theme_bw() +labs(title="12 SAMPLE SOMATIC MUTATION GENES RESILT3")+theme(axis.text.x = element_text(size = 8,color = "black", face = "bold", vjust = 0.5, hjust = 1, angle = 90))+ theme(panel.grid =element_blank())+ ylim(limits=c(0,12))
p4<-ggplot(data34)+geom_bar(aes(x=GENE,y=COUT,fill=factor(SAMPLE)),stat= 'identity', position = 'stack')+theme_bw() +labs(title="12 SAMPLE SOMATIC MUTATION GENES RESILT4")+theme(axis.text.x = element_text(size = 8,color = "black", face = "bold", vjust = 0.5, hjust = 1, angle = 90))+ theme(panel.grid =element_blank())+ ylim(limits=c(0,12))

grid.arrange(p1,p2,p3,p4,nrow=2,ncol=2)
dev.off()


#数据统计
DATASTEMP<-read.table("clipboard", header = T,sep='\t')

data.frame(summary(DATASTEMP$Gene.refGene))

nedata<-data.frame(paste(DATASTEMP$SAMPLE,DATASTEMP$Gene.refGene,sep="\t"))
tiff("datazhanshi-12sample2-2.png",width = 2160, height = 1080)
p1<-ggplot(DATASTEMP)+geom_bar(aes(x=Alignment.location.reduced.software,y=DATASTEMP[,-1],fill=factor(SAMPLE)),stat= 'identity', position = 'stack')+theme_bw()+theme_classic()+labs(title="12 SAMPLE SOMATIC MUTATION GENES RESILT")+theme(axis.text.x = element_text(size =12,color = "black", face = "bold", vjust = 0.5, hjust = 1, angle = 90))
p1+ theme(axis.text.y = element_blank())  #去除刻度标签
p1+ theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) #去除刻度线

dev.off()
install.packages("vcd")
library(vcd)
attach(Arthritis)
counts <- table(data$gene, data[,2:13])
spine(data2, main="Spinogram Example")
detach(Arthritis)


#数据分析17-08-02
rm(list=ls())
DATASTEMP<-read.table("clipboard", header = T,sep='\t')










