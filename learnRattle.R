#学习R数据挖掘工具rattle
setwd("D:\\myscript\\")
#need R >= 3.4

rm(list=ls())

tempdata<-read.table("outcnvfreq.txt",header=T)

m<-read.table("clipboard", header = T, sep = '\t') 

install.packages("RGtk2")
install.packages("rattle")

remove.packages("RGtk2")

library("rattle")
rattle( )  #not work
rattle( useGtkBuilder = TRUE)

install.packages("D:/myscript/R脚本库/rattle_5.0.10.zip", repos = NULL, type = "win.binary")

install.packages("D:\\myscript\\R脚本库\\rattle_5.0.10\\rattle\\")   #error
install.packages("rattle",contriburl = "http://rattle.togaware.com/bin/windows/contrib/3.3/rattle_5.0.10.zip",dependencies = TRUE) #error



library("rattle")
rattle( ) 
rattle( useGtkBuilder = TRUE)  #but it still doesn't work

###++++++++++++++##正式可以使用的反而是下面这样====================
install.packages("rattle", repos="http://rattle.togaware.com")  # install this ,the it can wor看！！！ 
library("rattle")
rattle( ) 
