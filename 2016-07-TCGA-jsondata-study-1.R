setwd("D:\\data")
library(rjson)
fromJSON(json_str,file="clinical.project-TCGA-LUAD.2016-07-04T05-49-08.389307.json",method = "c",unexpected.escape = "error")
jsondata<-fromJSON(file="clinical.project-TCGA-LUAD.2016-07-04T05-49-08.389307.json",method = "C",unexpected.escape = "error")
