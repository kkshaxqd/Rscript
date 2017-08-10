rm(list=ls())
library(mygene)
gtmp<- read.table("humangene_list.txt") 
mmm<-queryMany(gtmp,scopes="symbol",fields=c("symbol","go","pathway"),species="human")
akb35<-mmm[,c(1,2,4,3,12,10,14)]
write.table(akb35,file="humangeneanontion.txt",sep="\t",quote=F,na="NA",row.names=F)
###########��˼·����Ȼ�����һ�������perl��������ô�Ҿͷֿ�һ�������������ϵ�һ��###
########���ȴ���go.bp��#######
gobp<-mmm[,c(1,2)]
write.table(gobp,file="humangenegobp.txt",sep="\t",quote=F,na="NA",row.names=F)


###ѭ��######
pb<-txtProgressBar(min=0,max=74843,style=3)
test<-akb35$go.BP[1,]
test<-unlist(test)
for (i in 1:10) {
  test<-unlist(akb35$go.BP[i])
  bql<-paste(test,sep="\t")
}
#####��test���������bbs.txt��#####
sink("bbs.txt")# ������ض���bbs.txt
test
sink()# �����ض���
dev.off()
sink()
writeLines(test,"bbs")    ##һ��GO�Ż��Ϊһ��
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


pb<-txtProgressBar(min=0,max=10,style=3)    #����ɽ���
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

pb<-txtProgressBar(min=0,max=10,style=3)    #����ɽ���,��10��Ϊ��
for (i in 1:10) {
  yname<-akb35[i,1]                  #�õ�������
  for(j in 2:7){
    xtemp<-unlist(akb35[i,j], recursive = TRUE)  #���б�����չ���������ȷ�������֣�����xtemp<-unlist(akb35[1,2])����
    k<-length(xtemp$id)  #�õ�id�����ĳ���
    n=1
    nty<-c()
    lbg<-c()
    gname<-c()
    if(n<=k){
      nty[n]<- paste(xtemp$id[n],xtemp$term[n],xtemp$name[n],xtemp$pubmed[n],xtemp$evidence[n],sep="|") #���ģ�ת������Ҫ�ĸ�ʽ
      n<-n+1
    }
    else nty[n]<-NULL
  }
  lbg[i]<-paste(nty[n],sep="\t")             #��ȫ���ŵ�һ����
  gname[i]<-paste(yname,lbg[i],sep="\t")      #����������ŵ�һ����
  bgname<-paste(gname[i], sep="\n")             #����ͬ���зŵ�bgname��
  Sys.sleep(0.00001)
  setTxtProgressBar(pb, i)
}
write.table(bgname,file="humangeneanonotion160725.txt",sep="\t",quote=F,na="NA",row.names=F)









