rm(list=ls())
bbq<-read.table("genehgnameX.txt")  #读不成功
gtmp<- read.table("Annotation_AllGenes_Clean.txt",header=TRUE)  
genenames<-gtmp[,1]
library(mygene)
mmm<-queryMany(genenames,scopes="symbol",fields=c("ensembl.gene","entrezgene"),species="human")
di4lie<-mmm[,4]
llobs<-getGenes(di4lie,fields = c("symbol","go","pathway")) 
akb21<-llobs[,c(6,4,5,3,8,7,11)]
akb23<-cbind(mmm[,1],akb21)  #不行。。
akb25<-cbind(di5file,akb21)  #显示第一行就。。
akb27<-cbind

rm(akb25)


akb29<-cbind(llobs[,c(6,4,5,3,8,7,11)],mmm[,1])

##########new#########
di5file<-mmm[,1]
llobs21<-getGenes(di5file,fields = c("symbol","go","pathway")) #not found
######################
di6lie<-mmm[,4]
llobs6<-getGenes(di6lie,fields = c("symbol","go","pathway"))

##########new############
gtmp<- read.table("Annotation_AllGenes_Clean.txt",header=TRUE)  
genenames<-gtmp[,1]
library(mygene)
mmm<-queryMany(genenames,scopes="symbol",fields=c("ensembl.gene","entrezgene"),species="human",returnall=TRUE)
####直接就有了。。####
akb12<-mmm[,c(2,5,7,6,9,10,12)]
write.table(akb12,file="genehgname20160720.txt",quote=F,row.names=F,sep="\t")
####现在开始研究list格式数据转换吧##
unlist(akb12[1,2])  ##可以提取
listtest<-unlist(akb12[1,3])   #listtest是数据框.然后怎么把GO号与term对应
listtest$id[1]
listtest$term[1]
unlist(akb12[5,5])



unlist(akb12[1,2])

##################################编程处理提取list内容。然而并没成功，发现R对于变量与循环的处理能力好差。###
gtmp<- read.table("Annotation_AllGenes_Clean.txt",header=TRUE) 
genenames<-gtmp[,1]
mmm<-queryMany(genenames,scopes="symbol",fields=c("symbol","go","pathway"),species="human")
write.table(mmm,file="genetemp.txt",sep="\t")
akb32<-mmm[,c(2,6,5,3,4,7,10,11)]
#####################################
relation_matrix<-matrix(0,74451,8)
relation_matrix<-data.matrix(akb32)

pb<-txtProgressBar(min=0,max=100,style=3)
for(i in 1:100)
{
  xtemp<-unlist(akb32[i,3],recursive = TRUE)
  k<-length(xtemp$id) 
  while(k>0){
    goid<-xtemp$id[k]
    goneirong<-xtemp$term[k]
    gopubmed<-xtemp$pubmed[k]
    goevdience<-xtemp$evidence[k]
    k<-k-1
    rgobp<-paste(goid,goneirong,gopubmed,goevdience,sep="|")
  }
  gobp<-paste(rgobp,sep="\n")
  Sys.sleep(0.00001)
  setTxtProgressBar(pb, i)
} 



xtemp<-unlist(akb32[1,3])
yname<-akb32[1,2]
k<-length(xtemp$id)
while(k>0)
{
  nty[[k]]<-paste(xtemp$id[[k]],xtemp$term[[k]],xtemp$name[[k]],xtemp$pubmed[[k]],xtemp$evidence[[k]],sep="|")
  lbg<-cbind(lpg,nty[[k]]) 
  k<-k-1
}
gname<-paste(yname,lbg,sep="    ")





pathwayname<-xtemp$name[l]


pb<-txtProgressBar(min=0,max=10,style=3)  
for (i in 1:10) {
  for(j in 2:7){
    xtemp<-unlist(akb35[i,j], recursive = TRUE)         
    yname<-akb35[i,1]
    k<-length(xtemp$id)
    while(n<k){
      nty<- paste(xtemp$id[n],xtemp$term[n],xtemp$name[n],xtemp$pubmed[n],xtemp$evidence[n],sep="|")
      lbg<-paste("nty",sep="\t")
      n<-n+1
    }
    gname<-paste(yname,lbg,sep="\t")
    bgname<-paste("gname", sep="\n")
  }
  Sys.sleep(0.00001)
  setTxtProgressBar(pb, i)
}

rm(i,j,k,pb,xtemp,yname)


vector<-c(1:10)
varj<-c(3:8)
for (i in vector) {
  for(j in varj){
    xtemp<-apply(akb32[i,j],1,unlist)
  }
  yname<-akb32[i,2]
  
}







for (i in 1:10) {
  for(j in 2:8){
    xtemp<-unlist(akb32[i,j],recursive = T)      
    yname<-akb32[i,2]
    k<-length(xtemp[[id]])
    while(n<k){
      nty<- paste(xtemp$id[[n]],xtemp$term[[n]],xtemp$name[[n]],xtemp$pubmed[[n]],xtemp$evidence[[n]],sep="|")
      lbg<-paste("nty",sep="\t")
      n<-n+1
    }
    gname<-paste(yname,lbg,sep="\t")
    bgname<-paste("gname", sep="\n")
  }
}

write.table(bgname,file="newgeneannotionX.txt")   ################另外，R中for语句执行效率太低。也更糟糕。。
#####################不行，R中竟然不识别这种xtemp$id[n]这种变量，只识别确定的数xtemp$id[1]，这点太糟糕了。。###
rm(i,k,pb,xtemp,goevdience,gobp,goid,goneirong,gopubmed,rgobp)

rm(akp,gname,j,lpg,lbg,nty,yname,akh,bgname)

rm(see,bee)
#############
b=1
m=3
unlist(akb12[b,m])   
#######这样可以###
############
x<-c(1:10)
y<-c(2:7)
i=1
j=2
for (i in 1:10) {
  for (j in 2:7) {
    zzz<-unlist(akb12[i])  
  }
}
m=c(2:7)

i<-i+1
i=1
while(i<10){
  for (j in 2:7) {
    zzz<-unlist(akb12[]) ;
  }
  i<-i+1}
#########################################
install.packages("rlist")
library(rlist)

assssssssssssss<-list.all(akb12,query,na.rm=F)



rm(i)
akh<-list() 
for (i in 3:8) {
  akh<-unlist(akb32[1,8])
  akp<-apply(akh, 1, paste)
}

see<-getGenes(c(1:90000))
write.table(see,file="allhumangenes.txt") 
see<-getGenes(c(120000:200000))
write.table(see,file="allhumangenes.txt") 
bee<-see[see$notfound!=TRUE]

  