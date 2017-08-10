
setwd("D:\\LiYingWork\\Course_2016\\2016-Dragon-Xian\\Rcode_2016\\dragon_star_data\\Basis")

setwd("E:\\2016-Dragon-Xian\\Rcode_2016\\dragon_star_data\\Basis")  


### 安装package

#install.packages("RODBC")

install.packages("XLConnect")

library("RODBC")
detach(package:RODBC)


library()

list.files()

list.dirs()

dir(pattern="*.Rdata")

apropos("test")

x=scan()


set.seed(0)
sample(1:300,10)  

sample(1:300,10)  

sample(1:300,10,replace=T)







path=getwd()


## vector
Dragon_name <- c("li","shi","wu","guan","zhang","feng")
Dragon_name

print(Dragon_name)

Dragon_id = c(1,2,3,4,5,6)
Dragon_id

Dragon_cloth = c(FALSE,TRUE,FALSE,TRUE,FALSE,FALSE)

Dragon_cloth


Dragon_cloth[1]


Dragon_id[2:5]


Dragon_name[-2]

Dragon_name[-2:-3]

Dragon_name[c(FALSE,TRUE,FALSE,TRUE,FALSE,TRUE)]


###factor
Dragon_gender <-
  factor(c("female","male","male","male","male","female"))

Dragon_size <-
  factor(c("L","XL","XL","XXL","XXXL","L"),levels = c("S","M","L","XL","XXL","XXXL"))

Dragon_ages <-c(30,18,20,19,39,30)
Gragon_degree<-c('Doctor','Doctor','Doctor','Doctor','Master','Master')


###list
Dragon = list(
  name = Dragon_name,id = Dragon_id,cloth = Dragon_cloth,gender = Dragon_gender,size =
    Dragon_size,age=Dragon_ages,degree=Gragon_degree
)



Dragon

names(Dragon)


length(Dragon)

tt_Dragon=table(Dragon)
ftable(tt_Dragon,row.vars='gender',col.vars="degree")


Dragon = list(Dragon_name,Dragon_id,Dragon_cloth,Dragon_gender,Dragon_size)


print(Dragon)

names(Dragon) = c("name","id","cloth","gender","size")


Dragon = list()
Dragon[[1]] = Dragon_name
Dragon[[2]] = Dragon_id
Dragon[[3]] = Dragon_cloth
Dragon[[4]] = Dragon_gender
Dragon[[5]] = Dragon_size
names(Dragon) = c("name","id","cloth","gender","size")

Dragon

Dragon[1]
Dragon[[1]]
Dragon[[1]][4]



Dragon$size

Dragon['id']
Dragon["id"]
Dragon[1:2]

Dragon[c("name","id")]

Dragon$note = c(9,9)
Dragon

### data.frame

FraDragon = data.frame(
  name = Dragon_name,id = Dragon_id,cloth = Dragon_cloth,gender = Dragon_gender,size =
    Dragon_size
)

names(FraDragon)

rownames(FraDragon) = FraDragon$name

colnames(FraDragon)

FraDragon = FraDragon[,-1]
FraDragon

FraDragon["li",]


FraDragon[,"size"]


##
FraDragon = as.data.frame(Dragon[1:5])

write.csv(Dragon[1:5],file = 'Dragon.csv')

##### 读取文件
#install.packages("RODBC")
library(RODBC)
a = odbcConnectExcel("Dragon_1.xls")

Dragon = sqlFetch(a,"Dragon")
close(a)


############
if (!require(XLConnect)) {
  install.packages('XLConnect')
  library(XLConnect)
}



connect<-loadWorkbook('Dragon_1.xls')

a=readWorksheet(connect,'Dragon')

b=readWorksheet(connect,'Dragon_c')




#############

mode(Dragon)
Dragon

####
Dragon = read.csv("Dragon.csv")

###

sink("a.txt")
Dragon
sink()

###
Dragon = read.table("Dragon.csv")

Dragon
mode(Dragon)

Dragon = read.table("a.txt")
Dragon

###################NA

x = c(1,2,4,NA)
is.na(x)
x[!is.na(x)]

sum(is.na(x))

#####  重复
Dragon[7:19,] = Dragon[2,]
Dragon

duplicated(Dragon)

sum(duplicated(Dragon))

Dragon[!duplicated(Dragon),]

unique(Dragon)
Dragon[]

save(Dragon,file = "Dragon.Rdata")
rm(list = ls())
load("Dragon.Rdata")
ls()

x = matrix(c(1:9),ncol = 3,nrow = 3)
x
rownames(x) = c('A','B','C')
colnames(x) = c('C','D','E')
x
x['A','C']
ncol(x)
nrow(x)
dim(x)


sink("b.txt")
for (i in 1:5) {
  print("Hello world!")
  print(i)
}
sink()
list.files()
bb = read.table("b.txt")
bb

a = c(2:10);
for (i in c(3,8,9)) {
  print(a[i])
}



a = c(2:10);
class(a)
mode(a)

class(Dragon)
mode(Dragon)

score <- 60
if (score < 60)
{
  print("考试不及格")
}else
{
  print("考试及格")
}

score <- 95
if (score < 60)
{
  print("考试不及格")
}else if (score < 80)
{
  print("考试中等")
}else if (score < 90)
{
  print("考试良好")
}else
  print("考试优秀")

a = 0
while (a <= 5) {
  print("Enjoy!")
  a = a + 1
}

#a = 0
#while (a <= 5) {
#  print("Enjoy!")

#}

a = 5;
class(a);
a = c(1:10);
class(a);
a = matrix(c(1:9), 3, 3)
class(a)
a = data.frame();
class(a)
a = list();
class(a)
a = 'ggg'
class(a)

a = 5;
mode(a);
a = c(1:10);
mode(a);
a = matrix(c(1:9), 3, 3)
mode(a)
a = data.frame();
mode(a)
a = list();
mode(a)
a = 'ggg'
mode(a)



fun_static_data = function(x) {
  L = data.frame(
    mean = mean(x),median = median(x),sd = sd(x),min = min(x),max = max(x),
    Q1 = quantile(x,probs = 0.25),Q3 = quantile(x,probs = 0.75),IQR =
      IQR(x)
  );
  return(L);
}

x = runif(20)
static_x = fun_static_data(x)
static_x
summary(x)

x = rnorm(20)
static_x = fun_static_data(x)
static_x
summary(x)

## standardization

y=scale(x)
mean(y)
sd(y)


########################
y=runif(1000)

x = rnorm(5000)

par(mfrow=c(1,2))
qqnorm(x)
qqnorm(y)




x = rnorm(5000)
y=runif(1000)
z = rnorm(5000,mean=2,sd=10)

par(mfrow=c(1,3))
qqplot(x,y)
qqplot(x,z)
qqplot(y,z)



shapiro.test(x)
shapiro.test(y)
shapiro.test(z)

ks.test(x)
ks.test(y)
ks.test(z)


ks.test(rnorm(5000),runif(1000))
ks.test(rnorm(5000),rnorm(5000,mean=2,sd=10))
ks.test(runif(1000),rnorm(5000,mean=2,sd=10))
ks.test(rnorm(200),rnorm(5000))

ks.test(rnorm(5000),"punif")
ks.test(rnorm(5000),'pnorm')
ks.test(runif(100),"pnorm")
ks.test(runif(100),"punif")

################## 列联表+test

tt_Dragon=table(Dragon)
a= ftable(tt_Dragon,row.vars='gender',col.vars="size")

fisher.test(a)

#################


x=rnorm(2000)
par(mfrow=c(1,1))
boxplot(x,outline=TRUE,col="red")

X=list()
X[["China"]]=x
X[["USA"]]=x-0.5
boxplot(X,outline=TRUE,col=c("red","blue"))


