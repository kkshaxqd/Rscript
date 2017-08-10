install.packages("xlsx")
library(xlsx)
workbook<-"D:\\data\\16-10\\测序深度和覆盖度计算过程\\5830dui.depth.final.xlsx"
stuscore<-read.xlsx(workbook,1)#1表示sheet1
zhuanzhistuscore<-t(stuscore)
aaaa<-zhuanzhistuscore[2:89,1:80]
bbb<-zhuanzhistuscore[1,1:80]
boxplot(zhuanzhistuscore[2:89,1:80])
sample<-zhuanzhistuscore[2:89,0]
gene<-zhuanzhistuscore[1,1:80]
boxplot(sample~gene,data=zhuanzhistuscore,main="58+30对样本深度")

par(mfrow=c(1,2))
neirong1<-c(5,95)
neirong2<-c(5,10)
lbls<-c("")
rm(list=ls())
#######R可视化学习###### 
####edgelist####
install.packages("igraph")
install.packages("network")
install.packages("sna")
install.packages("ndtv")
setwd("D:\\data\\16-12\\12月学习\\Polnet2015\\Data")
nodes <- read.csv("Dataset1-Media-Example-NODES.csv", header=T, as.is=T)
links <- read.csv("Dataset1-Media-Example-EDGES.csv", header=T, as.is=T)
head(nodes)
head(links)
nrow(nodes);length(unique(nodes$id))
nrow(links);nrow(unique(links[,c("from","to")]))
links<-aggregate(links[,3],links[,-3],sum)
links<-links[order(links$from,links$to),]
colnames(links)[4]<-"weight"
rownames(links)<-NULL

####matrix#####
nodes2 <- read.csv("Dataset2-Media-User-Example-NODES.csv", header=T, as.is=T)
links2 <- read.csv("Dataset2-Media-User-Example-EDGES.csv", header=T, row.names=1)
head(nodes2)
head(links2)
links2 <- as.matrix(links2)
dim(links2)
dim(nodes2)

####network visualization:first steps with igraph
library(igraph)
net <- graph.data.frame(links, nodes, directed=T)
net
E(net) # The edges of the "net" object
V(net) # The vertices of the "net" object
E(net)$type # Edge attribute "type"
V(net)$media # Vertex attribute "media"
# You can also manipulate the network matrix directly:
net[1,]
net[5,7]
plot(net) # not a pretty picture!
net <- simplify(net, remove.multiple = F, remove.loops = T)
plot(net, edge.arrow.size=.4,vertex.label=NA)
help(igraph)
plot(x=1:10, y=rep(5,10), pch=19, cex=3, col="dark red")
points(x=1:10, y=rep(6, 10), pch=19, cex=3, col="557799")
points(x=1:10, y=rep(4, 10), pch=19, cex=3, col=rgb(.25, .5, .3))
plot(x=1:5, y=rep(5,5), pch=19, cex=12, col=rgb(.25, .5, .3, alpha=.5), xlim=c(0,6))
par(bg="white")
col.tr <- grDevices::adjustcolor("557799", alpha=0.7)
plot(x=1:5, y=rep(5,5), pch=19, cex=12, col=col.tr, xlim=c(0,6))

####R中颜色研究
colors() # List all named colors
grep("blue", colors(), value=T) # Colors that have "blue" in the name
pal1 <- heat.colors(5, alpha=1) # 5 colors from the heat palette, opaque
pal2 <- rainbow(5, alpha=.5) # 5 colors from the heat palette, transparent
plot(x=1:10, y=1:10, pch=19, cex=5, col=pal1)
plot(x=1:10, y=1:10, pch=19, cex=5, col=pal2)
palf <- colorRampPalette(c("gray80", "skyblue"))
plot(x=10:1, y=1:10, pch=19, cex=5, col=palf(10))
# If you don't have R ColorBrewer already, you will need to install it:
install.packages("RColorBrewer")
library(RColorBrewer)
display.brewer.all()
display.brewer.pal(8, "Set3")
display.brewer.pal(8, "Spectral")
display.brewer.pal(8, "Blues")
pal3 <- brewer.pal(10, "Set3")
plot(x=10:1, y=10:1, pch=19, cex=4, col=pal3)

#A brief detour II: Fonts in R plots R中字体研究
install.packages("extrafont")
library(extrafont)
# Import system fonts - may take a while, so DO NOT run this during the workshop.
font_import()
fonts() # See what font families are available to you now.
loadfonts(device = "win") # use device = "pdf" for pdf plot output.
plot(net, vertex.size=60)
plot(net, vertex.size=20,edge.arrow.size=.4,vertex.color="skyblue",vertex.label.family="Arial Black" )

# First you may have to let R know where to find ghostscript on your machine:
Sys.setenv(R_GSCMD = "C:/Program Files/gs/gs9.10/bin/gswin64c.exe")
## pdf() will send all the plots we output before dev.off() to a pdf file:
pdf(file="ArialBlack.pdf")
plot(net, vertex.size=30, vertex.color,vertex.label.family="Shonar Bangla" )
dev.off()
embed_fonts("ArialBlack.pdf", outfile="ArialBlack_embed.pdf")

#plotting networks
# Plot with curved edges (edge.curved=.1) and reduce arrow size:
plot(net, edge.arrow.size=.4, edge.curved=.1)

####2017-03-29
require(ggplot2)
datause<-read.table(file = "clipboard", header = TRUE) 
ggplot(datause)+geom_bar(aes(x=Inputsize))
ggplot(datause)+geom_bar(aes(x=pipline,y=Time,fill=factor(MACHINE)),stat="identity")
ggplot(datause)+geom_bar(aes(x=pipline, fill=Time))+windrose()
ggplot(datause)+geom_bar(aes(x=factor(1), fill=Time))#+coord_polar(theta="y")

ggplot(datause)+geom_bar(aes(x =pipline,y=Inputsize,fill=factor(MACHINE)),stat="identity")
ggplot(datause)+geom_bar(aes(x =pipline,y=Time,fill=factor(MACHINE)),stat="identity")+coord_flip()
ggplot(datause)+geom_bar(aes(x =pipline,y=COST,fill=factor(MACHINE)),stat="identity")
par(mfrow=c(1,3))                               
library(ggplot2)
library(gridExtra)  #用来把多图放到一起

opts(panel.background=theme_blank(),axis.title.x=theme_blank(), axis.title.y=theme_blank(),axis.text.x=theme_blank(),axis.text.y=theme_blank(),axis.ticks=theme_blank())+coord_flip()

setwd("D:\\data\\17-04")                               
png("test1.png")

p1<-ggplot(datause)+geom_bar(aes(x =pipline,y=Inputsize,fill=factor(MACHINE)),stat="identity")+opts(panel.background=theme_blank(),axis.title.x=theme_blank(), axis.title.y=theme_blank(),axis.text.x=theme_blank(),axis.text.y=theme_blank(),axis.ticks=theme_blank())
p2<-ggplot(datause)+geom_bar(aes(x =pipline,y=Time,fill=factor(MACHINE)),stat="identity")+coord_flip()+opts(panel.background=theme_blank(),axis.title.x=theme_blank(), axis.title.y=theme_blank(),axis.text.x=theme_blank(),axis.text.y=theme_blank(),axis.ticks=theme_blank())
p3<-ggplot(datause)+geom_bar(aes(x =pipline,y=COST,fill=factor(MACHINE)),stat="identity")+opts(panel.background=theme_blank(),axis.title.x=theme_blank(), axis.title.y=theme_blank(),axis.text.x=theme_blank(),axis.text.y=theme_blank(),axis.ticks=theme_blank())


p1<-ggplot(datause)+geom_bar(aes(x =pipline,fill=factor(MACHINE)),stat="identity")+facet_grid(Diet ~ Inputsize)
p2<-ggplot(datause)+geom_bar(aes(x =pipline,y=Time,fill=factor(MACHINE)),stat="identity")
p3<-ggplot(datause)+geom_bar(aes(x =pipline,y=COST,fill=factor(MACHINE)),stat="identity")

grid.arrange(p1, p2, p3, ncol=1,nrow=3,widths=3)
dev.off()

###so最后画图命令为###
library(ggplot2)
library(gridExtra)  #用来把多图放到一起
datause<-read.table(file = "clipboard", header = TRUE) 
setwd("D:\\data\\17-04")                               
tiff("test3.png")
p1<-ggplot(datause)+geom_bar(aes(x =pipline,y=Inputsize,fill=factor(MACHINE)),stat="identity")+geom_text(aes(x =pipline,y=Inputsize,label=Inputsize),position=position_dodge(.9))
p2<-ggplot(datause)+geom_bar(aes(x =pipline,y=Time,fill=factor(MACHINE)),stat="identity") +geom_text(aes(x =pipline,y=Time,label=Time),position=position_dodge(.9))+coord_flip()
p3<-ggplot(datause)+geom_bar(aes(x =pipline,y=COST,fill=factor(MACHINE)),stat="identity") +geom_text(aes(x =pipline,y=COST,label=COST),position=position_dodge(.9))
grid.arrange(p1, p2, p3, ncol=1,nrow=3,widths=3)
dev.off()                             



######ggplot2教学分享#####
library(ggplot2)
dt
lmx <- lm(dt$price ~ dt$carat)
gs1 <- geom_line(aes(y = lmx$fitted.values), size = 3, color = "red")
gs2 <- geom_smooth(aes(group = 1), method = "lm", se = FALSE, size = 1.5)
p + gs1 + gs2 
p + gs1 + gs2 + scale_x_log10()


library(ggplot2)
attach(iris)
p <- ggplot(data=iris,aes(x = Sepal.Length,y = Sepal.Width))
p + geom_point(aes(colour = Species)) + stat_smooth() + 
  labs(title = "Iris of Sepal.length \n According to the Sepal.Width") +
  theme_classic() + theme_bw() +annotate("text",x=7,y=4,parse = T,label = "x[1]==x[2]",size=6, family="serif",fontface="italic", colour="darkred") +facet_wrap()
Species
example(facet_grid)



setwd("D:\\data\\17-04")                               
tiff("test211234.png")
ggplot(data = iris, aes(x = Sepal.Length, y = Sepal.Width, color=Species)) + geom_point(aes(shape=Species), size=3)+ stat_smooth() 
dev.off() 


set.seed(1234)
x <- sample(c(1,2,4,6,7), size = 1000, replace = TRUE,prob = c(0.1,0.2,0.2,0.3,0.2))
ggplot(data = data.frame(x = x), mapping = aes(x = factor(x), y = ..count..))+ geom_bar(stat = 'count', fill = 'steelblue', colour = 'darkred')


x <- rep(1:5, each = 3)
y <- rep(c('A','B','C'),times = 5)
set.seed(1234)
z <- round(runif(min = 10, max = 20, n = 15)) 
df <- data.frame(x= x, y = y, z = z)
ggplot(data = df, mapping = aes(x = factor(x), y = z,fill = y)) + geom_bar(stat = 'identity', position = 'dodge')

with(mpg,table(class,year))
p <- ggplot(data=mpg,aes(x=class,fill=factor(year)))
p + geom_bar(position='dodge')
p + geom_bar(position='stack')
p + geom_bar(position='fill')
p + geom_bar(position='identity')

y=c(1.1,1.8,2.5,3.6,3.1,2.7,1.9,-0.1,-3.5,3.0)
x=2001:2010
data=data.frame(x,y)
p=ggplot(data,aes(x,y,fill=y))
p+geom_bar(stat="identity")+ 
  geom_abline(intercept = 0, slope = 0,size=1,colour='gray')+
  geom_text(aes(label=y),hjust=0.5, vjust=-0.5 )+
  scale_y_continuous(limits=c(-3.8,4.2))+
  labs(x='年份', y='GDP增长率%')+
  ggtitle("美国GDP增长率")
example(opt)

a <- ggplot(mpg, aes(x=hwy))
a + stat_bin(aes(fill=..count.., color=-1*..ndensity..), binwidth = 1)

x <- c('A','B','C','D','E','F','G')
y <-c('xx','yy','yy','xx','xx','xx','yy')
z <- c(10,33,12,9,16,23,11) 
df<- data.frame(x = x, y = y, z = z)
ggplot(data = df, mapping = aes(x= x, y = z, fill = y)) + geom_bar(stat = 'identity')

set.seed(12)
x <- 1980 + 1:35
y <- round(100*rnorm(35))
df <- data.frame(x = x,y = y)
# 判断y是否为正值
df <- transform(df,judge = ifelse(y>0,"YES","NO"))
# 去除图例用theme()主题函数
ggplot(df,aes(x = x,y = y,fill = judge))+
  geom_bar(stat = "identity")+
  theme(legend.position= "")+
  xlab("Year")+
  scale_fill_manual(values = c("darkred","blue"))

ggplot(data = df, mapping = aes(x = x, y = y, fill = judge))+ 
  geom_bar(stat = 'identity', position = 'identity')+ 
  scale_fill_manual(values = c('blue','red'), guide = FALSE)+ 
  xlab('Year')

example(guides)

source("https://bioconductor.org/biocLite.R")
biocLite("ChIPseeker")
require(ChIPseeker)
files <- getSampleFiles()
peak <- readPeakFile(files[[4]])
peak[1:2,]



####batterman
f1 <- function(x) {
  y1 <- 3*sqrt(1-(x/7)^2)
  y2 <- -3*sqrt(1-(x/7)^2)
  y <- c(y1,y2)
  d <- data.frame(x=x,y=y)
  d <- d[d$y > -3*sqrt(33)/7,]
  return(d)
}
x1 <- c(seq(3, 7, 0.001), seq(-7, -3, 0.001))
d1 <- f1(x1)
p1 <- ggplot(d1,aes(x,y)) + geom_point(color="red") +xlab("") + ylab("") + theme_bw()
p1
x2 <- seq(-4,4, 0.001)
y2 <- abs(x2/2)-(3*sqrt(33)-7)*x2^2/112-3 + sqrt(1-(abs(abs(x2)-2)-1)^2)
d2 <- data.frame(x2=x2, y2=y2)
p2 <- p1 + geom_point(data=d2, aes(x=x2,y=y2), color="yellow")
p2
x3 <- c(seq(0.75,1,0.001), seq(-1,-0.75,0.001))
y3 <- 9-8*abs(x3)
d3 <- data.frame(x3=x3, y3=y3)
p3 <- p2+geom_point(data=d3, aes(x=x3,y=y3), color="green")
p3
x4 <- c(seq(0.5,0.75,0.001), seq(-0.75,-0.5,0.001))
y4 <- 3*abs(x4)+0.75
d4 <- data.frame(x4=x4,y4=y4)
p4 <- p3+geom_point(data=d4, aes(x=x4,y=y4), color="steelblue")
p4
x5 <- seq(-0.5,0.5,0.001)
y5 <- rep(2.25,length(x5))
d5 <- data.frame(x5=x5,y5=y5)
p5 <- p4+geom_point(data=d5, aes(x=x5,y=y5))
p5
x6 <- c(seq(-3,-1,0.001), seq(1,3,0.001))
y6 <- 6 * sqrt(10)/7 +
  (1.5 - 0.5 * abs(x6)) * sqrt(abs(abs(x6)-1)/(abs(x6)-1)) -
  6 * sqrt(10) * sqrt(4-(abs(x6)-1)^2)/14
d6 <- data.frame(x6=x6,y6=y6)
d6 <- d6[-which(is.na(d6[,2])), ]
p6 <- p5+geom_point(data=d6,aes(x=x6,y=y6), colour="blue")
p6


ggsave(filename="test7.png",plot=p6,width=8,height=8,units="cm",dpi=300)

######################


get_gene_name_vect<- function(gtf_data){
  gene_info_vect <- gtf_data[,9]
  gene_info_vect[gtf_data[,3]!='gene'] = ''
  gene_name_vect <- sapply(gene_info_vect,function(x){
    strsplit(strsplit(x,';')[[1]][3], '\"')[[1]][2]
  })
  names(gene_name_vect) = NULL
  return(gene_name_vect)
}

find_target_data <- function(gene, all_data, gene_name_vect){
  target_row <- which(gene_name_vect == toupper(gene))
  i <- 1
  while( is.na(gene_name_vect[target_row +i])){
    i <- i+ 1
  }
  target_data<-all_data[target_row:(target_row+i-1),]
  colnames(target_data) =c('chr','db','record','start','end','tmp1','strand','tmp3','tmp4')
  return(target_data)
}

plot_gene_structure <- function (target_data, gene, rect_width=0.8, intron_width=0.2) {
  gene_start <- target_data$start[1]
  gene_end <- target_data$end[1]
  transcript_num <- length(which(target_data$record == 'transcript'))
  tmp_colors <- c('green', 'red', 'black', 'blue', 'yellow', 'blue','grey','grey')
  names(tmp_colors) <- c('gene', 'CDS', 'transcript', 'exon','start_codon','stop_codon','three_prime_utr','five_prime_utr')
  rect_sep <- c(intron_width,rep(rect_width,7))
  names(rect_sep) <-c('transcript','gene', 'exon', 'CDS','start_codon','stop_codon','three_prime_utr','five_prime_utr')
  par(mar=c(6,6,2,6), bty='n', new=F)
  plot(gene_start:gene_start, 1, type = 'n', xlab='', ylab ='',ylim = c(0,transcript_num*2+1), xlim = c(gene_start,gene_end), yaxt='n',xaxt='n')
  title(main = gene,sub = paste("chr",target_data$chr,": ",gene_start,"-",gene_end,sep=""))
  if (target_data$strand[1]=='-'){
    strand <- 1
  }else { strand <-2 }
  j <- 0
  y_lab <- c(j)
  rect_ybottom <- j-rect_sep[1]
  rect_ytop <- j+rect_sep[1]
  #rect(gene_start, rect_ybottom, gene_end, rect_ytop, col = tmp_colors['gene'], border = F)
  arrows(gene_start, j+rect_sep[1], gene_end, j+rect_sep[1], code= strand, col = tmp_colors['gene'], lwd = 5, angle=30)
  for (rows in (2:dim(target_data)[1])){
    if (target_data$record[rows] == 'transcript') {
      j <- j + rect_sep[2] * 2.5
      y_lab <- c(y_lab,j)
    }
    rect_ybottom <- j-rect_sep[target_data$record[rows]]
    rect_ytop <- j+rect_sep[target_data$record[rows]]
    #print(tmp_colors[target_data$record[rows]])
    rect(target_data$start[rows], rect_ybottom, target_data$end[rows], rect_ytop, col=tmp_colors[target_data$record[rows]], border = NA)
  }
  
  axis(1,c(gene_start,gene_start+(gene_end-gene_start)/5,gene_start+(gene_end-gene_start)*2/5,gene_start+(gene_end-gene_start)*3/5,gene_start+(gene_end-gene_start)*4/5,gene_end))
  
  axis(2,y_lab, labels=c('gene',paste('transcipt', seq(1:transcript_num), sep=' ')),las=1)
  legend(x=gene_end, transcript_num*(rect_sep[2] * 2.5)+1,horiz=F,
         box.lty=0,xpd=T,
         c('gene','non_CDS_exon','CDS_exon','UTR_exon','intron'),
         fill=c('green','blue','red','grey','black'),
         border = F,bty='n',
         cex=0.8)
}


setwd('D:\\data\\17-04')
file <- 'Homo_sapiens.GRCh38.87.chr.gtf'

library(data.table)
#比较耗时
all_data<-as.data.frame(fread(file, header=F, sep='\t', stringsAsFactors=F, skip =5))
#load('all.Rdata')

gene='tp53'
#gene='ANxa1'

#比较耗时
gene_name_vect <- get_gene_name_vect(all_data)

target_data <- find_target_data (gene, all_data, gene_name_vect)

png(paste(gene,'.png'),width = 1200, height = 680)
plot_gene_structure(target_data, gene)
dev.off()

rm(list=ls())



library(ggplot2)

data <- read.table(file = "transcript.txt", sep = "\t", stringsAsFactors = F)

xlab_name <- paste(data[1,1], paste(data[1,3],data[1,4],sep = "-"), sep = ":")
title <- "ANXA1"
data <- data[-1,]
data <- data[data[,2]!="transcript",]

num <- length(unique(as.factor(data[,1])))
color <- rainbow(num)

data_min <- min(c(data[,3], data[,4]))
data$new_start <- data[,3] - data_min
data$new_end <- data[,4] - data_min
data$new_color <- color[as.numeric(data[,1])]

p <- ggplot()+geom_line()+
  annotate("segment", x = data$new_start, xend = data$new_end, y = data[,1], yend = data[,1], color = data$new_color,size = 2)+
  labs(x = xlab_name, y = "number", title = title)
head(all_data)

library(ggplot2)

data <- read.table(file = "transcript.txt", sep = "\t", stringsAsFactors = F)

xlab_name <- paste(data[1,1], paste(data[1,3],data[1,4],sep = "-"), sep = ":")
title <- "ANXA1"
data <- data[-1,]
data <- data[data[,2]!="transcript",]

num <- length(unique(as.factor(data[,1])))
color <- rainbow(num)

data_min <- min(c(data[,3], data[,4]))
data$new_start <- data[,3] - data_min
data$new_end <- data[,4] - data_min
data$new_color <- color[as.numeric(data[,1])]

png("test12.png")
p <- ggplot()+geom_line(x = data$new_start, xend = data$new_end, y = data[,1], yend = data[,1])+
  annotate("segment", x = data$new_start, xend = data$new_end, y = data[,1], yend = data[,1], color = data$new_color,size = 2)+
  labs(x = xlab_name, y = "number", title = title)+theme_bw()+theme(panel.grid=element_blank())
dev.off()

ggsave(filename="test13.png",plot=p,width=18,height=12,units="cm",dpi=700)



#####CNV可视化绘制###
#gene	start	stop	cnv
#ABR	11628	12352	-1
#AGPS	8952	9216	3
datause<-read.table(file = "clipboard", header = TRUE) 
gene<-datause[,1]
gene_start <- datause[,2]
gene_end <- datause[,3]
cnv_stat<-datause[,4]
genesite<-c(gene_start,gene_end)
xlim = c(gene_start,gene_end)
xlab_name <- paste(gene, paste(gene_start,gene_end,sep = "-"), sep = ":")

  labs(x = xlab_name, y = "number", title = title)+theme_bw()+theme(panel.grid=element_blank())
p <- ggplot(datause,aes(x=genesite,cnv_stat))+geom_point()+geom_line(linetype = 2)+labs(x = xlab_name, y = "CNVstat", title = title)+theme_bw()+theme(panel.grid=element_blank())
p

