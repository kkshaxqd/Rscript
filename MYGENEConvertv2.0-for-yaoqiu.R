################edit by ����ʯ############
##########20160729###########
############��mygene����õ�����ȫ����ע��ת��Ϊ����ۿ�����ʽ############
rm(list=ls())
library(mygene)
gtmp<- read.table("humangene_list.txt") 
mmm<-queryMany(gtmp,scopes="symbol",fields=c("symbol","go","pathway"),species="human")  #ÿ�����������������bp,cc,mf,pathway˳�����Щ��ͬ�����������Ӧ�������С�
akb37<-mmm[,c(1,14,13,15,8,4,6)]
write.table(akb37,file="humangeneanontion.txt",sep="\t",quote=F,na="NA",row.names=F) 
####################���°�����##########
gobptemp<-data.frame()
gocctemp<-data.frame()
gomftemp<-data.frame()
keggtemp<-data.frame()
reactometemp<-data.frame()
syname<-c()
# nbgname<-c("GENE","GO.BP","GO.CC","GO.MF","PATHWAY.KEGG","PATHWAY.REACTOME") �����ӷ��������ء������ظ����顣�������
hbgname<-c()
pb<-txtProgressBar(min=0,max=74484,style=3)    #����ɽ���,��10��Ϊ��
for (i in 1:74484) {
  syname[i]<-akb37[i,1]   #�õ�������
  sgobp<-c()
  sgocc<-c()
  sgomf<-c()
  skegg<-c()
  sreactome<-c()
  a<-c()
  b<-c()
  ###################paste������һֱ������ӽ�ȥ����ɾ��ԭ���ģ�ֻ�ǲ��ϵļӼӼ�????
  ###################����ÿ��ѭ���궼���¹�ԭλ
  sgkl<-c()  
  if(is.null(gobptemp<-unlist(akb37[i,2]))){  ###gobp
    sgobp<- paste(sgobp,"unknow",sep=" ")
  }
  else{ 
    k<-length(gobptemp$id) 
    n=1
    while(n<=k){
      a<-gobptemp$id[n]
      b<-gobptemp$id[n+1]
      sgobpnty<-paste(gobptemp$id[n],gobptemp$term[n],sep=",") 
      if (is.na(b)|a==b){} else{sgobp<-paste(sgobp,sgobpnty,sep="|")}
      n=n+1
    }
  }
  
  if(is.null(gocctemp<-unlist(akb37[i,3]))){  ####gocc
    sgocc<- paste(sgocc,"unknow",sep=" ")
  }
  else{ 
    k<-length(gocctemp$id) 
    n=1
    while(n<=k){
      a<-gocctemp$id[n]
      b<-gocctemp$id[n+1]
      sgoccnty<-paste(gocctemp$id[n],gocctemp$term[n],sep=",")
     if (is.na(b)|a==b){} else{sgocc<-paste(sgocc,sgoccnty,sep="|")}
      n=n+1
    }
  }
  if(is.null(gomftemp<-unlist(akb37[i,4]))){  ####gomf
    sgomf<- paste(sgomf,"unknow",sep=" ")
  }
  else{ 
    k<-length(gomftemp$id) 
    n=1
    while(n<=k){
      a<-gomftemp$id[n]
      b<-gomftemp$id[n+1]
      sgomfnty<-paste(gomftemp$id[n],gomftemp$term[n],sep=",")
      if (is.na(b)|a==b){} else{sgomf<-paste(sgomf,sgomfnty,sep="|") }
      n=n+1
    }
  }
  if(is.null(keggtemp<-unlist(akb37[i,5]))){  ###kegg
    skegg<- paste(skegg,"unknow",sep=" ")
  }
  else{ 
    k<-length(keggtemp$id) 
    n=1
    while(n<=k){
      a<-keggtemp$id[n]
      b<-keggtemp$id[n+1]
      skeggnty<-paste(keggtemp$id[n],keggtemp$name[n],sep=",")
      if (is.na(b)|a==b){} else{skegg<-paste(skegg,skeggnty,sep="|")}
      n=n+1
    }
  }
  if(is.null(reactometemp<-unlist(akb37[i,6]))){     ###reactome
    sreactome<-paste(sreactome,"unknow",sep=" ")
  }
  else{ 
    k<-length(reactometemp$id) 
    n=1
    while(n<=k){
      sreactomenty<-paste(reactometemp$id[n],reactometemp$name[n],sep=",")
      a<-reactometemp$id[n]
      b<-reactometemp$id[n+1]
      if (is.na(b)|a==b){} else{sreactome<-paste(sreactome,sreactomenty,sep="|")}
      n=n+1
    }
  }
  sgkl<-paste(syname[i],sgobp,sgocc,sgomf,skegg,sreactome,sep="\t")
  hbgname<-paste(hbgname,sgkl,sep="\n")   #����ͬ���зŵ�nbgname��
  Sys.sleep(0.00001)
  setTxtProgressBar(pb, i)
}

write.table(hbgname,file="human_gene_annotion_clean20160730new.txt",sep="\t",quote=F,na="NA",row.names=F,col.names=T)