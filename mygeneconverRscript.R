rm(list=ls())
library(mygene)
gtmp<- read.table("humangene_list.txt") 
mmm<-queryMany(gtmp,scopes="symbol",fields=c("symbol","go","pathway"),species="human")  #ÿ�����������������bp,cc,mf,pathway˳�����Щ��ͬ�����������Ӧ�������С�
akb37<-mmm[,c(1,14,13,15,8,4,6)]
write.table(akb37,file="humangeneanontion.txt",sep="\t",quote=F,na="NA",row.names=F) 
gobp<-akb37[,c(1,2)]
write.table(gobp,file="humangenegobp.txt",sep="\t",quote=F,na="NA",row.names=F)
#####�޸ĸ�ʽ����######
yname<-c()
xtemp<-array()
pb<-txtProgressBar(min=0,max=10,style=3)    #����ɽ���,��10��Ϊ��
for (i in 1:10) {
  yname[i]<-akb37[i,1]#�õ�������
  for(j in 2:7){
    xtemp[j]<-unlist(akb37[i,j], recursive = TRUE)  #���б�����չ���������ȷ�������֣�����xtemp<-unlist(akb35[1,2])����
    k<-length(xtemp[j]$id)  #�õ�id�����ĳ���
    n=1
    nty<-c()
    lbg<-c()
    gname<-c()
    while (n<=k){
      nty[n]<- paste(xtemp[j]$id[n],xtemp[j]$term[n],xtemp[j]$name[n],xtemp[j]$pubmed[n],xtemp[j]$evidence[n],sep="|")#���ģ�ת������Ҫ�ĸ�ʽ
      n<-n+1
    }
    lbg[j]<-cbind(nty[1]:nty[n]) #��ȫ���ŵ�һ����
    gname[i]<-cbind(yname[i],lbg[1]:lbg[j])  #����������ŵ�һ����
  }           
  bgname<-rbind(gname[1]:gname[i])             #����ͬ���зŵ�bgname��
  Sys.sleep(0.00001)
  setTxtProgressBar(pb, i)
}
write.table(bgname,file="humangeneanonotion160725.txt",sep="\t",quote=F,na="NA",row.names=F)


rm(i,j,lbg,n,nty,pb,xtemp,yname,k,gname)
rm(bgname)
use warning()
 

xtemp<-list()
i=2
j=3
n=4
xtemp<-unlist(akb35[i,j]) 



##############����ܵõ���Ҫ����ʽ���������ǲ�̫�ÿ�###############
yname<-c()
xtemp<-data.frame()
gname<-c()
lbg<-0
gkl<-0
pb<-txtProgressBar(min=0,max=10,style=3)    #����ɽ���,��10��Ϊ��
for (i in 1:10) {
  yname[i]<-akb35[i,1]#�õ�������
  for(j in 2:7){
    xtemp<-unlist(akb35[i,j], recursive = TRUE)  #���б�����չ���������ȷ�������֣�����xtemp<-unlist(akb35[1,2])����
    k<-length(xtemp$id)  #�õ�id�����ĳ���
    n=1
    while (n<=k){
      nty<- paste(xtemp$id[n],xtemp$term[n],xtemp$name[n],xtemp$pubmed[n],xtemp$evidence[n],sep="|")#���ģ�ת������Ҫ�ĸ�ʽ
      lbg<-paste(lbg,nty,sep="//") #��һ��ע����Ŀ����ȫ���ŵ�һ����
      n<-n+1
    }
    gkl<-paste(gkl,lbg,sep="\t")  #��ȫ��ע�����ݷŵ�һ����
  }  
  gname[i]<-paste(yname[i],gkl,sep=" ")
  bgname<-rbind(gname[1:i])             #����ͬ���зŵ�bgname��
  Sys.sleep(0.00001)
  setTxtProgressBar(pb, i)
}
write.table(bgname,file="humangeneanonotion160725.txt",sep="\t",quote=F,na="NA",row.names=F)
######################################################

















####################�����Ľ�##########################
xtemp<-data.frame()
gname<-c()
yname<-c()


pb<-txtProgressBar(min=0,max=74844,style=3)    #����ɽ���,��10��Ϊ��
for (i in 1:74844) {
  yname[i]<-akb37[i,1]#�õ�������
  gkl<-c()                                  #�����ÿ��ѭ�����ܵõ�һ��ȫ��ע����Ϣ�������ǰ���������ظ�д  
  for(j in 2:7){
    xtemp<-unlist(akb37[i,j], recursive = TRUE)  #���б�����չ���������ȷ�������֣�����xtemp<-unlist(akb35[1,2])����
    k<-length(xtemp$id)  #�õ�id�����ĳ���
    n=1
    lbg<-c()
    while (n<=k){
      nty<- paste(xtemp$id[n],xtemp$term[n],xtemp$name[n],xtemp$pubmed[n],xtemp$evidence[n],sep="|")#���ģ�ת������Ҫ�ĸ�ʽ
      lbg<-paste(lbg,nty,sep="/") #��һ��ע����Ŀ����ȫ���ŵ�һ����
      n<-n+1
    }
    gkl<-paste(gkl,lbg,sep="\t")  #��ȫ��ע�����ݷŵ�һ����
  }  
  gname[i]<-paste(yname[i],gkl,sep=" ")
  bgname<-paste(gname[1:i],sep="\n")             #����ͬ���зŵ�bgname��
  Sys.sleep(0.00001)
  setTxtProgressBar(pb, i)
}
write.table(bgname,file="humangeneanonotion2016072716.txt",sep="\t",quote=F,na="NA",row.names=F)
rm(pb)




rm(i,j,k,n,pb,lbg,nty,gkl)

###########################�򻯰汾################
sxtemp<-data.frame()
sgname<-c()
syname<-c()


pb<-txtProgressBar(min=0,max=74844,style=3)    #����ɽ���,��10��Ϊ��
for (i in 1:74844) {
  syname[i]<-akb37[i,1]#�õ�������
  sgkl<-c()                                  #�����ÿ��ѭ�����ܵõ�һ��ȫ��ע����Ϣ�������ǰ���������ظ�д  
  for(j in 2:7){
    sxtemp<-unlist(akb37[i,j], recursive = TRUE)  #���б�����չ���������ȷ�������֣�����xtemp<-unlist(akb35[1,2])����
    k<-length(sxtemp$id)  #�õ�id�����ĳ���
    n=1
    slbg<-c()
    while (n<=k){
      snty<- paste(sxtemp$id[n],sxtemp$term[n],sxtemp$name[n],sep=",")#���ģ�ת������Ҫ�ĸ�ʽ
      slbg<-paste(slbg,snty,sep="|") #��һ��ע����Ŀ����ȫ���ŵ�һ����
      n<-n+1
    }
    sgkl<-paste(sgkl,slbg,sep="\t")  #��ȫ��ע�����ݷŵ�һ����
  }  
  sgname[i]<-paste(syname[i],sgkl,sep="\t")
  sbgname<-paste(sgname[1:i],sep="\n")             #����ͬ���зŵ�bgname��
  Sys.sleep(0.00001)
  setTxtProgressBar(pb, i)
}
write.table(sbgname,file="humangeneanonotion20160728.txt",sep="\t",quote=F,na="NA",row.names=F,col.names=T)

##############################����ĵĸ����һЩ
rm(sgobp,sgocc,sgomf,skegg,sreactome,gobptemp,gocctemp,gomftemp,keggtemp,reactometemp,sgkl,nbgname,syname)
gobptemp<-data.frame()
gocctemp<-data.frame()
gomftemp<-data.frame()
keggtemp<-data.frame()
reactometemp<-data.frame()
syname<-c()
# nbgname<-c("GENE","GO.BP","GO.CC","GO.MF","PATHWAY.KEGG","PATHWAY.REACTOME") �����ӷ��������ء������ظ����顣�������
nbgname<-c()
pb<-txtProgressBar(min=0,max=74844,style=3)    #����ɽ���,��10��Ϊ��
for (i in 1:74844) {
  syname[i]<-akb37[i,1]   #�õ�������
  sgobp<-c()
  sgocc<-c()
  sgomf<-c()
  skegg<-c()
  sreactome<-c()
  
                    ###################paste������һֱ������ӽ�ȥ����ɾ��ԭ���ģ�ֻ�ǲ��ϵļӼӼ�????
                    ###################����ÿ��ѭ���궼���¹�ԭλ
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
  nbgname<-paste(nbgname,sgkl,sep="\n")   #����ͬ���зŵ�nbgname��
  Sys.sleep(0.00001)
  setTxtProgressBar(pb, i)
}

write.table(nbgname,file="human_gene_annotion_clean20160729.txt",sep="\t",quote=F,na="NA",row.names=F,col.names=T)
rm(sbgname)
rm(sgobpnty,sgoccnty,sgomfnty,skeggnty,sreactomenty)


i=19
asaas<-unlist(akb37[i,2])


##################20160729#############################
rm(sgobp,sgocc,sgomf,skegg,sreactome,gobptemp,gocctemp,gomftemp,keggtemp,reactometemp,sgkl,syname,a,b,hbgname)
gobptemp<-data.frame()
gocctemp<-data.frame()
gomftemp<-data.frame()
keggtemp<-data.frame()
reactometemp<-data.frame()
syname<-c()
# nbgname<-c("GENE","GO.BP","GO.CC","GO.MF","PATHWAY.KEGG","PATHWAY.REACTOME") �����ӷ��������ء������ظ����顣�������
ahbgname<-c()
pb<-txtProgressBar(min=0,max=74844,style=3)    #����ɽ���,��10��Ϊ��
for (i in 0:74844) {
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
  ahbgname<-paste(ahbgname,sgkl,sep="\n")   #����ͬ���зŵ�nbgname��
  Sys.sleep(0.00001)
  setTxtProgressBar(pb, i)
}

write.table(ahbgname,file="human_gene_annotion_clean2016073008.txt",sep="\t",quote=F,na="NA",row.names=F,col.names=T)
rm(sbgname,i,k,m,n,pb,asaas,yname,hbgname,nbgname)
rm(sgobpnty,sgoccnty,sgomfnty,skeggnty,sreactomenty)

