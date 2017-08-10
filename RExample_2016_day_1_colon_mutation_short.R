
rm(list=ls())

setwd("D:\\LiYingWork\\Course_2016\\2016-Dragon-Xian\\Rcode_2016\\dragon_star_data\\TCGA_colon_cancer_data")

setwd("E:\\2016-Dragon-Xian\\Rcode_2016\\dragon_star_data\\TCGA_colon_cancer_data")



getwd()

rm(list = ls())

list.files()

dir(pattern="*clinic")

load("COAD_clinical_data.RData")
ls()

mode(COAD_clinical_data)

names(COAD_clinical_data)

data_clinical <- COAD_clinical_data
head(data_clinical)

mode(data_clinical)

names(data_clinical)

data_clinical$pathologic_stage <- as.character(data_clinical$pathologic_stage)

all_stages <- unique(sort(data_clinical$pathologic_stage))
all_stages 



all_I_II_stages <- all_stages[c(2,3,4,5,6,7)]


data_stage_I_II <- list()

for (i in 1:length(all_I_II_stages)) {
  id <- which(data_clinical$pathologic_stage == all_I_II_stages[i])
  data_stage_I_II <- rbind(data_stage_I_II,data_clinical[id,])
}

data_stage_I_II
head(data_stage_I_II)

write.csv(data_stage_I_II,'data_stage_I_II.csv')

list.files()




##### »­Í¼
stage_num = table(data_clinical$pathologic_stage)
stage_num

stage_num = stage_num[-1]

plot(c(1:length(stage_num)),stage_num,xaxt = 'n',main = "stage_num",col = 'red',xlab ='stage',ylab = 'number of pateint')
mtext(names(stage_num),at = c(1:length(stage_num)),side = 1,las=2)

lines(
  c(1:length(stage_num)),stage_num,main = "stage_num",col = 'red',xlab = 'stage',ylab =
    'number of pateint'
)
#mtext(names(stage_num),at = c(1:length(stage_num)),side = 1,las=2)

barplot(
  stage_num, xaxt = 'n', main = 'stage_num', xlab =
    'stage',ylab = 'number',
  col = topo.colors(13)
);
mtext(names(stage_num),at = c(1:length(stage_num)),side = 1,las=2)


barplot(
  stage_num, legend = rownames(stage_num),xaxt = 'n', main = 'stage_num', xlab =
    'stage',ylab = 'number',
  col = topo.colors(13), beside = TRUE
);

ss = prop.table(stage_num)

pie(ss,main = 'colon cancer stage distribution')


############################


load("COAD_clinical_data.RData")
ls()

Type_Stages=list()

Type_Stages[["Stage I"]]=c("Stage I","Stage IA","Stage IB","Stage IC")
Type_Stages[["Stage II"]]=c("Stage II","Stage IIA","Stage IIB","Stage IIC")
Type_Stages[["Stage III"]]=c("Stage III","Stage IIIA","Stage IIIB","Stage IIIC")
Type_Stages[["Stage IV"]]=c("Stage IV","Stage IVA","Stage IVB","Stage IVC")





data_clinical <- COAD_clinical_data
names(data_clinical)

data_clinical$pathologic_stage <- as.character(data_clinical$pathologic_stage)

data_clinical$pathologic_stage

all_stages <- unique(sort(data_clinical$pathologic_stage))

Samples_stage=list()



for (i in names(Type_Stages)) {
  aa=intersect(all_stages,Type_Stages[[i]])
  NN=length(aa) 
  if (NN>0)
  {
   for (j in 1:NN) {
   id <- which(data_clinical$pathologic_stage == aa[j])
   Samples_stage[[i]] <- rbind( Samples_stage[[i]],data_clinical[id,])
   }
    }
}



Ages_stage=list()


for (i in names(Samples_stage)){
  
 a=as.numeric(Samples_stage[[i]][["age_at_initial_pathologic_diagnosis"]])

names(a)=Samples_stage[[i]][["X"]]

Ages_stage[[i]]=a
}
  


boxplot(Ages_stage, col =c("red","blue","green","pink"),ylab='age',main="COAD_Stage_Age")  

pdf("COAD_Stage_Age.pdf")
boxplot(Ages_stage, col =c("red","blue","green","pink"),ylab='age',main="COAD_Stage_Age")  
dev.off()

#######################  
##  age and Satge correlation

N=length(Ages_stage)

ss_stage=list()

nn_stage=unlist(lapply(Ages_stage,length))

for (i in 1: length(nn_stage))
{
  ss_stage[[i]]=rep(i,nn_stage[i])
}
  


########################

hist(unlist(Ages_stage))

shapiro.test(unlist(Ages_stage)) 

ks.test(unlist(Ages_stage),y="pnorm") 


range(unlist(Ages_stage))

which.max(unlist(Ages_stage))

which.min(unlist(Ages_stage))

