library(ggplot2)
install.packages("gridExtra")
library(gridExtra)  #用来把多图放到一起
cnv88sampletongjidelamp<-read.table("clipboard", header = T, sep = '\t')
colnames(cnv88sampletongjidelamp)[2]="CNVnumber" #重命名列名
tiff("cnv88sampletongji2.png",width = 1080, height = 720)
ggplot(cnv88sampletongjidelamp)+geom_bar(aes(x=Sample,y=CNVnumber,fill=Type), stat="identity")+theme_bw() +theme_classic()+labs(title="XHMM 88 Samples CNV RESILT")+theme(axis.text.x = element_text(size = 12,color = "black", face = "bold", vjust = 1, hjust = 1, angle = 60))+scale_y_continuous(breaks=seq(-30, 30, 10))
dev.off()
#5kbd
cnv88sampletongjidelamp5kb<-read.table("clipboard", header = T, sep = '\t')
colnames(cnv88sampletongjidelamp5kb)[2]="CNVnumber" #重命名列名
tiff("cnv88sampletongji3.png",width = 1080, height = 720)
ggplot(cnv88sampletongjidelamp5kb)+geom_bar(aes(x=Sample,y=CNVnumber,fill=Type), stat="identity")+theme_bw() +theme_classic()+labs(title="5KB XHMM 88 Samples CNV RESILT")+theme(axis.text.x = element_text(size = 12,color = "black", face = "bold", vjust = 1, hjust = 1, angle = 60))+scale_y_continuous(breaks=seq(-30, 30, 10))
dev.off()
#20kb
cnv88sampletongjidelamp20kb<-read.table("clipboard", header = T, sep = '\t')
colnames(cnv88sampletongjidelamp20kb)[2]="CNVnumber" #重命名列名
tiff("cnv88sampletongji4.png",width = 1080, height = 720)
ggplot(cnv88sampletongjidelamp20kb)+geom_bar(aes(x=Sample,y=CNVnumber,fill=Type), stat="identity")+theme_bw() +theme_classic()+labs(title="20KB XHMM 88 Samples CNV RESILT")+theme(axis.text.x = element_text(size = 12,color = "black", face = "bold", vjust = 1, hjust = 1, angle = 60))+scale_y_continuous(breaks=seq(-30, 30, 10))
dev.off()
#50kb
cnv88sampletongjidelamp50kb<-read.table("clipboard", header = T, sep = '\t')
colnames(cnv88sampletongjidelamp50kb)[2]="CNVnumber" #重命名列名
tiff("cnv88sampletongji5.png",width = 1080, height = 720)
ggplot(cnv88sampletongjidelamp50kb)+geom_bar(aes(x=Sample,y=CNVnumber,fill=Type), stat="identity")+theme_bw() +theme_classic()+labs(title="50KB XHMM 88 Samples CNV RESILT")+theme(axis.text.x = element_text(size = 12,color = "black", face = "bold", vjust = 1, hjust = 1, angle = 60))+scale_y_continuous(breaks=seq(-10, 10, 5))
dev.off()
#100kb
cnv88sampletongjidelamp100kb<-read.table("clipboard", header = T, sep = '\t')
colnames(cnv88sampletongjidelamp100kb)[2]="CNVnumber" #重命名列名
tiff("cnv88sampletongji6.png",width = 1080, height = 720)
ggplot(cnv88sampletongjidelamp100kb)+geom_bar(aes(x=Sample,y=CNVnumber,fill=Type), stat="identity")+theme_bw() +theme_classic()+labs(title="100KB XHMM 88 Samples CNV RESILT")+theme(axis.text.x = element_text(size = 12,color = "black", face = "bold", vjust = 1, hjust = 1, angle = 60))+scale_y_continuous(breaks=seq(-5, 5, 5))
dev.off()


#cnv
cnvtamp<-read.table("clipboard", header = T, sep = '\t')
tiff("cnvtamp.png",width = 1080, height = 720)
ggplot(cnvtamp)+geom_bar(aes(x=GENE,y=CNV,fill=SAMPLE), stat="identity")+theme_bw() +theme_classic()+labs(title="CNV RESILT")+theme(axis.text.x = element_text(size = 8,color = "black", face = "bold", vjust = 1, hjust = 1, angle = 90))+scale_y_continuous(breaks=seq(-3, 10, 1))
dev.off()

summary(cnvtamp$GENE)
