# G4A1<-read.table("clipboard",head=T)
# E4B1<-read.table("clipboard",head=T)
# H4C1<-read.table("clipboard",head=T)
# A4D1<-read.table("clipboard",head=T)
# F4E1<-read.table("clipboard",head=T)
# D4F1<-read.table("clipboard",head=T)
# E3G1<-read.table("clipboard",head=T)
# G3H1<-read.table("clipboard",head=T)
# C4A2<-read.table("clipboard",head=T)
# C3B2<-read.table("clipboard",head=T)
# F3C2<-read.table("clipboard",head=T)
# D3D2<-read.table("clipboard",head=T)
# H3E2<-read.table("clipboard",head=T)
# B3F2<-read.table("clipboard",head=T)
# B4G2<-read.table("clipboard",head=T)
# A3H2<-read.table("clipboard",head=T)
# gene<-read.table("clipboard",head=T)  #读取剪切板上的数据。
#GENEDATA<-rbind(gene,G4A1,E4B1,H4C1,A4D1,F4E1,D4F1,E3G1,G3H1,C4A2,C3B2,F3C2,D3D2,H3E2,B3F2,B4G2,A3H2)
rm(list=ls())
gene<-read.table("clipboard",head=T)

library(mygene)
entreID<-queryMany(gene,scopes="reporter", species="human")


help("htmlParse")
#####爬取网页 COSMIC INSERT DATA xml的网页
library(XML)
library(RCurl)

##### <dt class=\"l30\">Nucleotides Inserted:</dt> <dd class=\"l50\">gccgtgtgtggctgcagcccaggccac </dd>
###http://juntai.me/2016/03/30/R中的正则表达式  方法来自于此###
u<-"https://grch37-cancer.sanger.ac.uk/cosmic/mutation/overview?id=1155808"
res<-getURL(u)
nucleotides_inserted<-regexpr(pattern=">[atcg]+" , perl=TRUE,res)##pipei
xulie<-regmatches(x=res, nucleotides_inserted)  ##tiqu

link_region<-regexpr(pattern="link-region\\\">\\S+:\\d+..\\d+<" , perl=TRUE,res)
weidian<-regmatches(x=res, link_region)

mutation_id<-regexpr(pattern = ">COSM\\d+<",perl=TRUE,res)
bianhao<-regmatches(x=res, mutation_id)

cancer_zhonglei<-regexpr(pattern = "/browse/tissue#sn=\\D+&in=t",perl=T,res)
laiyuanzuzhi<-regmatches(x=res,cancer_zhonglei)

cds_mutation<-regexpr(pattern = "c.\\S+\\s\\(Insertion\\)",perl=T,res)
cds_xinxi<-regmatches(x=res,cds_mutation)

aaa<-paste(mutation_id,weidian,cds_xinxi,nucleotides_inserted,laiyuanzuzhi,sep="\t")
bbb<-paste(aaa,sep="\n")
###(new_apol1,file="1000G-APOL1.txt",sep="\t",quote=F,na="NA",row.names=F,col.names =T )
write.table(bbb,file="cosmic-insert-data.txt",sep="\t",quote=F,na="NA",row.names=F,col.names =T)  

testshujuji<-c("11134","12324","1231414")



i<-c("COSM20959")

xunhuan<-sub(pattern = "COSM",replacement = "id=",x=i)

u<-paste('https://grch37-cancer.sanger.ac.uk/cosmic/mutation/overview?',xunhuan,sep="")

res<-getURL(u)
######the way read data by lines cankao###
con <- file("e:/data.txt", "r")
line=readLines(con,n=1)
while( length(line) != 0 ) {
  print(line)
  line=readLines(con,n=1)
}
close(con)
#如果需要将一行的文字劈成多段，再进行处理，可以用strsplit函数，除此之外，还有一些常用的字符串处理函数，记录如下：
#substr(),nchar(), grep(), regexpr(), sub(), gsub()
################zhengshi-run#########
setwd("D:\\data\\17-01")
con<-file("Ins-cosmid-3.txt","r")
cosmicid<-readLines(con,n=1)
bbb=c()

i=1
pb<-txtProgressBar(min=0,max=1010,style=3)
while(length(cosmicid)!=0)
{
  cosmicid=readLines(con,n=1)
  print(i)
  print(cosmicid)
  aaa=c()
  nucleotides_inserted=c()
  link_region=c()
  mutation_id=c()
  cancer_zhonglei=c()
  cds_mutation=c()
  
  xunhuan<-sub(pattern = "COSM",replacement = "id=",x=cosmicid)
  u<-paste('https://grch37-cancer.sanger.ac.uk/cosmic/mutation/overview?',xunhuan,sep="")
  res<-getURL(u)
  
  nucleotides_inserted<-regexpr(pattern="class=\"l50\">\\w+\\s" , perl=TRUE,res)##pipei
  xulie<-regmatches(x=res, nucleotides_inserted)  ##tiqu
  
  link_region<-regexpr(pattern="link-region\\\">\\S+:\\d+..\\d+<" , perl=TRUE,res)
  weidian<-regmatches(x=res, link_region)
  
  mutation_id<-regexpr(pattern = ">COSM\\d+<",perl=TRUE,res)
  bianhao<-regmatches(x=res, mutation_id)
  
  cancer_zhonglei<-regexpr(pattern = "/browse/tissue#sn=\\D+&in=t",perl=T,res)
  laiyuanzuzhi<-regmatches(x=res,cancer_zhonglei)
  
  cds_mutation<-regexpr(pattern = "c.\\S+\\s\\(Insertion\\)",perl=T,res)
  cds_xinxi<-regmatches(x=res,cds_mutation)

  aaa<-paste(bianhao,weidian,cds_xinxi,xulie,laiyuanzuzhi,sep="\t")
  bbb<-paste(aaa,sep="\n")
  write.table(bbb,file="cosmic-insert-data-1.txt",append=T,sep="\t",quote=F,na="NA",row.names=F,col.names =F) ##use append can direct add the data to file
  i=i+1
  Sys.sleep(20)
  setTxtProgressBar(pb, i)
}

###################################第一天爬取后调整######
con<-file("Ins-cosmid-3.txt","r")
cosmicid<-readLines(con,n=1)
bbb=c()

i=1
pb<-txtProgressBar(min=0,max=1010,style=3)
while(length(cosmicid)!=0)
{
  cosmicid=readLines(con,n=1)
  print(i)
  print(cosmicid)
  aaa=c()
  nucleotides_inserted=c()
  link_region=c()
  mutation_id=c()
  cancer_zhonglei=c()
  cds_mutation=c()
  u=c()
  res=c()
  
  xunhuan<-sub(pattern = "COSM",replacement = "id=",x=cosmicid)
  u<-paste('https://grch37-cancer.sanger.ac.uk/cosmic/mutation/overview?',xunhuan,sep="")
  res<-getURL(u)
  
  nucleotides_inserted<-regexpr(pattern="class=\"l50\">\\w+\\s" , perl=TRUE,res)##pipei
  xulie<-regmatches(x=res, nucleotides_inserted)  ##tiqu
  
  link_region<-regexpr(pattern="link-region\\\">\\S+:\\d+..\\d+<" , perl=TRUE,res)
  weidian<-regmatches(x=res, link_region)
  
  mutation_id<-regexpr(pattern = ">COSM\\d+<",perl=TRUE,res)
  bianhao<-regmatches(x=res, mutation_id)
  
  cancer_zhonglei<-regexpr(pattern = "/browse/tissue#sn=\\D+&in=t",perl=T,res)
  laiyuanzuzhi<-regmatches(x=res,cancer_zhonglei)
  
  cds_mutation<-regexpr(pattern = "c.\\S+\\s\\(Insertion\\)",perl=T,res)
  cds_xinxi<-regmatches(x=res,cds_mutation)
  
  aaa<-paste(bianhao,weidian,cds_xinxi,xulie,laiyuanzuzhi,sep="\t")
  bbb<-paste(aaa,sep="\n")
  write.table(bbb,file="cosmic-insert-data-1-1.txt",append=T,sep="\t",quote=F,na="NA",row.names=F,col.names =F) ##use append can direct add the data to file
  i=i+1
  Sys.sleep(0.0010)
  setTxtProgressBar(pb, i)
}

##TissueHistogram( {"link":1,"legend":["Gene name( Frequency )","Samples with mutation","All samples","url"],"max":47,"data":[["Lung",47,"https://grch37-cancer.sanger.ac.uk/cosmic/browse/tissue#sn=lung&in=t"],["Ovary",6,"https://grch37-cancer.sanger.ac.uk/cosmic/browse/tissue#sn=ovary&in=t"],["Breast",1,"https://grch37-cancer.sanger.ac.uk/cosmic/browse/tissue#sn=breast&in=t"]],"x":"Tissue","x_legend":0,"total":3,"y":"Total number of samples","title":"Tissue Distribution"} ) } );</script> 
COSM20959
i<-c("COSM20959")
xunhuan<-sub(pattern = "COSM",replacement = "id=",x=i)
u<-paste('https://grch37-cancer.sanger.ac.uk/cosmic/mutation/overview?',xunhuan,sep="")
res<-getURL(u)
cancer_zhonglei<-regexpr(pattern = "Tissue[\"a-zA-Z0-9\\t,\\(\\)\\[\\]:\\-\\.\\/\\#=&]</script>",perl=T,res)
laiyuanzuzhi<-regmatches(x=res,cancer_zhonglei) ###这样怎么也匹配不了全部的任意字符，所以把res的数据写入文本，后面用perl来处理比较好。
