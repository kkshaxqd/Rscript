rm(list=ls())
library(mygene)
gtmp<- read.table("humangene_list.txt") 
mmm<-queryMany(gtmp,scopes="symbol",fields=c("symbol","go","pathway"),species="human")
akb35<-mmm[,c(1,2,4,3,12,10,14)]
write.table(akb35,file="humangeneanontion.txt",sep="\t",quote=F,na="NA",row.names=F)
###########新思路，既然整体放一起很难用perl处理，那么我就分开一个个处理，最后合到一起。###
########首先处理go.bp的#######
gobp<-mmm[,c(1,2)]
write.table(gobp,file="humangenegobp.txt",sep="\t",quote=F,na="NA",row.names=F)


###循环######
pb<-txtProgressBar(min=0,max=74843,style=3)
test<-akb35$go.BP[1,]
test<-unlist(test)
for (i in 1:10) {
  test<-unlist(akb35$go.BP[i])
  bql<-paste(test,sep="\t")
}
#####将test内容输出到bbs.txt里#####
sink("bbs.txt")# 将输出重定向到bbs.txt
test
sink()# 结束重定向
dev.off()
sink()
writeLines(test,"bbs")    ##一个GO号会成为一行
writeLines()

ske1<-sapply(akb35$go.BP,unlist)
tnf<-paste(ske1[[1]],sep="\t")
rm(i,test,ske1)

for (i in 1:10) {
  test<-sapply(akb35$go.BP[10],unlist)
  dfn<-sapply(test,cbind,sep="\t")
}
test10<-unlist(akb35$go.BP[10])


#####v0###########


pb<-txtProgressBar(min=0,max=10,style=3)    #看完成进度
for (i in 1:10) {
  yname[i]<-akb35[i,1]
  for(j in 2:7){
    xtemp[i]<-unlist(akb35[i,j], recursive = TRUE) 
    k<-length(xtemp$id)
    n=1
    if(n<=k){
      nty[n]<- paste(xtemp$id[n],xtemp$term[n],xtemp$name[n],xtemp$pubmed[n],xtemp$evidence[n],sep="|")
      n<-n+1
    }
    else nty[n]<-NULL
  }
  lbg[i]<-paste(nty[n],sep="\t")
  gname[i]<-paste(yname[i],lbg[i],sep="\t")
  bgname<-paste(gname[i], sep="\n")
  Sys.sleep(0.00001)
  setTxtProgressBar(pb, i)
}
write.table(bgname,file="humangeneanonotion160725.txt",sep="\t",quote=F,na="NA",row.names=F)

###########v1######################

pb<-txtProgressBar(min=0,max=10,style=3)    #看完成进度,以10个为例
for (i in 1:10) {
  yname<-akb35[i,1]                  #得到基因名
  for(j in 2:7){
    xtemp<-unlist(akb35[i,j], recursive = TRUE)  #将列表数据展开，如果是确定的数字，比如xtemp<-unlist(akb35[1,2])就行
    k<-length(xtemp$id)  #得到id向量的长度
    n=1
    nty<-c()
    lbg<-c()
    gname<-c()
    if(n<=k){
      nty[n]<- paste(xtemp$id[n],xtemp$term[n],xtemp$name[n],xtemp$pubmed[n],xtemp$evidence[n],sep="|") #核心，转化成想要的格式
      n<-n+1
    }
    else nty[n]<-NULL
  }
  lbg[i]<-paste(nty[n],sep="\t")             #将全部放到一行里
  gname[i]<-paste(yname,lbg[i],sep="\t")      #将与基因名放到一行里
  bgname<-paste(gname[i], sep="\n")             #将不同的行放到bgname里
  Sys.sleep(0.00001)
  setTxtProgressBar(pb, i)
}
write.table(bgname,file="humangeneanonotion160725.txt",sep="\t",quote=F,na="NA",row.names=F)










