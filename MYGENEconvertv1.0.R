################edit by 张桥石############
##########20160728###########
############将mygene包获得的人类全基因及注释转化为方便观看的形式############
rm(list=ls())
library(mygene)
gtmp<- read.table("humangene_list.txt") 
mmm<-queryMany(gtmp,scopes="symbol",fields=c("symbol","go","pathway"),species="human")  #每次运行这个，出来的bp,cc,mf,pathway顺序会有些不同，根据这个相应调整就行。
akb37<-mmm[,c(1,14,13,15,8,4,6)]
write.table(akb37,file="humangeneanontion.txt",sep="\t",quote=F,na="NA",row.names=F) 
####################最新版命令##########
gobptemp<-data.frame()
gocctemp<-data.frame()
gomftemp<-data.frame()
keggtemp<-data.frame()
reactometemp<-data.frame()
syname<-c()
# nbgname<-c("GENE","GO.BP","GO.CC","GO.MF","PATHWAY.KEGG","PATHWAY.REACTOME") 这样子反而不行呢。。会重复几遍。。输出。
nbgname<-c()
pb<-txtProgressBar(min=0,max=74844,style=3)    #看完成进度,以10个为例
for (i in 1:74844) {
  syname[i]<-akb37[i,1]   #得到基因名
  sgobp<-c("GO.BP")
  sgocc<-c("GO.CC")
  sgomf<-c("GO.MF")
  skegg<-c("P.KEGG")
  sreactome<-c("P.REACTOME")
  
                    ###################paste函数是一直往里面加进去，不删除原来的，只是不断的加加加????
                    ###################所以每次循环完都重新归原位
  sgkl<-c()  
  if(is.null(gobptemp<-unlist(akb37[i,2]))){  ###gobp
    sgobp<- paste(sgobp,"unknow",sep="|")
   }
    else{ 
      k<-length(gobptemp$id) 
        n=1
        while(n<=k){
        sgobpnty<-paste(gobptemp$id[n],gobptemp$term[n],sep=",") 
        sgobp<-paste(sgobp,sgobpnty,sep="|")
        n=n+1
      }
    }
  
  if(is.null(gocctemp<-unlist(akb37[i,3]))){  ####gocc
    sgocc<- paste(sgocc,"unknow",sep="|")
  }
  else{ 
    k<-length(gocctemp$id) 
    n=1
    while(n<=k){
      sgoccnty<-paste(gocctemp$id[n],gocctemp$term[n],sep=",")
      sgocc<-paste(sgocc,sgoccnty,sep="|")
      n=n+1
    }
  }
  if(is.null(gomftemp<-unlist(akb37[i,4]))){  ####gomf
    sgomf<- paste(sgomf,"unknow",sep="|")
  }
  else{ 
    k<-length(gomftemp$id) 
    n=1
    while(n<=k){
      sgomfnty<-paste(gomftemp$id[n],gomftemp$term[n],sep=",")
      sgomf<-paste(sgomf,sgomfnty,sep="|") 
      n=n+1
    }
  }
  if(is.null(keggtemp<-unlist(akb37[i,5]))){  ###kegg
    skegg<- paste(skegg,"unknow",sep="|")
  }
  else{ 
    k<-length(keggtemp$id) 
    n=1
    while(n<=k){
      skeggnty<-paste(keggtemp$id[n],keggtemp$name[n],sep=",")
      skegg<-paste(skegg,skeggnty,sep="|")
      n=n+1
    }
  }
  if(is.null(reactometemp<-unlist(akb37[i,6]))){     ###reactome
    sreactome<-paste(sreactome,"unknow",sep="|")
  }
  else{ 
    k<-length(reactometemp$id) 
    n=1
    while(n<=k){
      sreactomenty<-paste(reactometemp$id[n],reactometemp$name[n],sep=",")
      sreactome<-paste(sreactome,sreactomenty,sep="|")
      n=n+1
    }
  }
  sgkl<-paste(syname[i],sgobp,sgocc,sgomf,skegg,sreactome,sep="\t")
  nbgname<-paste(nbgname,sgkl,sep="\n")   #将不同的行放到nbgname里
  Sys.sleep(0.00001)
  setTxtProgressBar(pb, i)
}

write.table(nbgname,file="human_gene_annotion_clean2016072817.txt",sep="\t",quote=F,na="NA",row.names=F,col.names=T)