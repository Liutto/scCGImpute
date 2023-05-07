rm(list = ls())
gc()
library(Rtsne)
library(ggplot2)
library(parallel)
library(ClusterR)
library(tidyr)
library(dplyr)
library(gridExtra)
library(grid)
library(kernlab)
library(RColorBrewer)
library(ggeasy)   #easy_add_legend_title
library(patchwork)  #p1+p2
library(ggsignif)  #箱线图p值
library(ggpubr)#箱线图p值
source("D:/software/RStudio/single_cell/scCGImpute/code/CalculateARI.R")
source("D:/software/RStudio/single_cell/scCGImpute/code/MatrixFusion.R")
source("D:/software/RStudio/single_cell/scCGImpute/code/CalculateNMI.R")
source("D:/software/RStudio/single_cell/scCGImpute/code/CalculatePurity.R")

########################################################
#箱线图
###Ziegenhain
corercc <-function(ercc,methodss){
    corpeason<-c()
    for (i in 1:dim(methodss)[2]) {
      corpeason<-c(corpeason,cor(ercc,methodss[,i]))
    }
    return(corpeason)
  }

library("xlsx")
ercc<-read.xlsx("D:/software/RStudio/single_cell/datasets/Ziegenhain/ERCC.xlsx",1)
ercc = as.matrix(ercc[,c(2,4)])
erccloc = order(ercc[,1])
ercc = ercc[erccloc,]
ercc = as.numeric(ercc[,2])
dataercclong<-c()
datasets<-c("SmartSeq2","CELseq2","MARSseq","SCRBseq","SmartSeq")
for (i in 1:length(datasets)) {
  
   eval(parse(text = paste0(datasets[i],"s_scCGImpute<-read.csv('D:/software/RStudio/single_cell/datasets/Ziegenhain/",datasets[i],"/s_scCGImpute/scCGImpute_count.csv',header = T)")))
   eval(parse(text = paste0(datasets[i],"rawdata<-read.csv('D:/software/RStudio/single_cell/datasets/Ziegenhain/",datasets[i],".csv',header = T)")))
   eval(parse(text = paste0("erccrawloc = grep('ERCC',",datasets[i],"rawdata[,1])")))
   eval(parse(text = paste0("erccraw = ",datasets[i],"rawdata[erccrawloc,]")))
   erccrawloc = order(erccraw[,1])
   erccraw = erccraw[erccrawloc,]
   
   eval(parse(text = paste0("erccscCGImpute = ",datasets[i],"s_scCGImpute[grep('ERCC',",datasets[i],"s_scCGImpute[,1]),]")))
   erccscCGImpute = erccscCGImpute[erccrawloc,]
   erccscCGImpute  = apply(erccscCGImpute[,-1],2,as.numeric)
   
   erccrawcor = corercc(log(ercc),log(erccraw[,-1]+1))
   erccscCGImpute = corercc(log(ercc),log(erccscCGImpute+1))
   
   dataercc=cbind(erccrawcor,erccscCGImpute)
   dataercc=data.frame(dataercc)
   colnames(dataercc)<-c("raw","scCGImpute")
   dataercc <- dataercc %>% 
     gather(key = 'Ziegenhain',value = 'correlation') #命名两个列名
   eval(parse(text = paste0("methods<-rep('",datasets[i],"',length(erccrawcor))")))
   dataercctmp = cbind(methods,dataercc)
   dataercclong = rbind(dataercclong,dataercctmp)
   #内容没有 格式转换
   # InsectSprays %>% {
   eval(parse(text = paste0("p",datasets[i],"<- ggplot(dataercc,mapping = aes(x = '',y = correlation , fill = factor(Ziegenhain,levels = c('raw','scCGImpute')) ))+
     geom_boxplot(outlier.shape = 21,outlier.colour = 'red',outlier.fill ='blue',
                                      col = 'black',width=1.2)+#,fill = brewer.pal(7,'Paired')
     theme_bw()+
     # theme(axis.text.x = element_text(angle = 45,size = 10, color = 'black',hjust = 1,vjust = 1))+##设置x轴字体大小
     labs(title = '",datasets[i],"',x = ' ')+
     easy_add_legend_title('Methods')+
      stat_compare_means(method = 't.test')")))
   
}

pSmartSeq2+pCELseq2+pSmartSeq+pSCRBseq+pMARSseq+guide_area()+ plot_layout(guides = "collect")
library(ggpubr)
#分组绘图
p<-ggboxplot(dataercclong,x="methods",y="correlation",fill ="Ziegenhain",palette = 'mycol')
ercclong<-p+stat_compare_means(aes(group = Ziegenhain),method = 't.test')+
  easy_add_legend_title('Methods')
write.csv(dataercclong,"D:/software/RStudio/single_cell/datasets/Ziegenhain/result/dataercclong.csv")
ggsave(ercclong,filename = "D:/software/RStudio/single_cell/datasets/Ziegenhain/result/ercclonggroup.pdf",width = 12,height = 9)



########################################################
###聚类结果评价
#Blakeley_gold
x <- read.table("D:/software/RStudio/single_cell/datasets/Blakeley_gold/GSE66507_human_blastocyst_rnaseq_counts.txt",
                header = TRUE, check.names = FALSE, row.names = 1)
x <- as.matrix(x)
dim(x)
x1<-x[rowSums(x)>0.001,]
dim(x1)
x2<-x1[rowSums(x1!=0)>=3,]
dim(x2)

write.csv(x2,"D:/software/RStudio/single_cell/datasets/Blakeley_gold/Blakeley_rawdata.csv")

rawdata<-read.csv("D:/software/RStudio/single_cell/datasets/Blakeley_gold/Blakeley_rawdata.csv")
scImpute<-read.csv("D:/software/RStudio/single_cell/datasets/Blakeley_gold/scImpute/scimpute_count.csv")
saver<-read.csv("D:/software/RStudio/single_cell/datasets/Blakeley_gold/Blakeley_saver.csv")
magic<-read.csv("D:/software/RStudio/single_cell/datasets/Blakeley_gold/Blakeley_magic.csv")
deepimpute<-read.csv("D:/software/RStudio/single_cell/datasets/Blakeley_gold/deepimpute.csv")#445.37313175201416 seconds
VIPER<-read.csv("D:/software/RStudio/single_cell/datasets/Blakeley_gold/VIPER.csv")  #53.90018 secs
scRMD<-read.csv("D:/software/RStudio/single_cell/datasets/Blakeley_gold/scRMD.csv")  #2.131514 secs
scRecover<-read.csv("D:/software/RStudio/single_cell/datasets/Blakeley_gold/outDir_scRecover/scRecover+scImpute.csv") #8.084504 mins
scCGImpute<-read.csv("D:/software/RStudio/single_cell/datasets/Blakeley_gold/s_scCGImpute/scCGImpute_count.csv")

TEloc<-grep("TE",colnames(rawdata[,-1]))
PEloc<-grep("PE",colnames(rawdata[,-1]))
EPIloc<-grep("EPI",colnames(rawdata[,-1]))
label<-rep(0,30)
label[TEloc]=1
label[PEloc]=2
label[EPIloc]=3

methodss<-c("rawdata","saver","scImpute","deepimpute","VIPER","scRMD","scRecover","scCGImpute")

dataend<-c(rep(0,length(methodss)*3))
for (i in 1:10) {
  print(i)
  for (methodssi in 1:length(methodss)) {
    print(methodssi)
        eval(parse(text = paste0("totalCounts_by_cell = colSums(",methodss[methodssi],"[,-1])")))
      eval(parse(text = paste0("data = sweep(",methodss[methodssi],"[,-1], MARGIN = 2, 10^6/totalCounts_by_cell, FUN = '*')")))
      count_lnorm = log10(data + 1.01)
      pcax=prcomp(t(count_lnorm))
      specc_res = specc(pcax$x[,1:2], centers = 3, kernel = 'rbfdot')
      eval(parse(text = paste0(methodss[methodssi],"ARI<-CalculateARI(label,specc_res@.Data)")))
      eval(parse(text = paste0(methodss[methodssi],"NMI<-CalculateNMI(label,specc_res@.Data)")))
      eval(parse(text = paste0(methodss[methodssi],"Purity<-CalculatePurity(label,specc_res@.Data)")))
      
  }
  assessmethod<-c("ARI","NMI","Purity")
  name<-methodss
  for (j in 1:length(assessmethod)) {
    eval(parse(text = paste0(assessmethod[j],"<-c()")))
    for (i in 1:length(methodss)) {
      eval(parse(text = paste0(assessmethod[j],"<-cbind(",assessmethod[j],",",methodss[i],assessmethod[j],")")))
    }
    eval(parse(text = paste0("filll<-rep('",assessmethod[j],"',length(methodss))")))
    eval(parse(text = paste0(assessmethod[j],"bar<-rbind(t(name),t(filll),as.numeric(",assessmethod[j],"))")))
  } 
  databar<-cbind(ARIbar,NMIbar,Puritybar)
  databar<-t(databar)
  
  dataend<-dataend+as.numeric(databar[,3])
}
dataend = dataend/10
databar[,3]=dataend
databar<-data.frame(databar)
databar[,3] = as.numeric(databar[,3])
# write.csv(databar,"E:/Blakeley_gold_cluster.csv")
colnames(databar)<-c("variable","methods","value")
databar<-data.frame(databar)
# write.csv(databar,"E:/Blakeley_gold.csv")
#######绘图#########
# RColorBrewer::display.brewer.all()
Blakeley_gold1<-ggplot(databar, aes(x=methods, y=value, fill=variable)) +
  geom_bar(stat="identity",position=position_dodge(),
           color="black", width=.8) +
  theme_bw()+
  labs(title = 'Blakeley  30cells × 22251genes',x = "")+
  scale_fill_brewer(palette = "Paired")+
  theme(axis.text.x = element_text(size = 14, color = "black"))+##设置x轴字体大小
  theme(axis.text.y = element_text(size = 14, color = "black"))+##设置y轴字体大小
  theme(title=element_text(size=13))#设置标题字体大小
Blakeley_gold1
Blakeley_goldd2<-ggplot(databar, aes(x=reorder(variable,-value), y=value, fill=methods)) +
  geom_bar(stat="identity",position = "stack",
           color="black", width=.8) +
  theme_bw()+
  labs(title = 'Blakeley  30cells × 22251genes',x = "")+
  scale_fill_brewer(palette = "Paired")+
  theme(axis.text.x = element_text(size = 10, color = "black",angle = 45, hjust = 1,
                                   vjust = 1))+    ##设置x轴字体大小
  theme(axis.text.y = element_text(size = 14, color = "black"))+   ##设置y轴字体大小
  theme(title=element_text(size=13))#设置标题字体大小
Blakeley_goldd2
write.csv(databar,"D:/software/RStudio/single_cell/datasets_unify/Blakeley3_nor+pca+specc_databar.csv")
ggsave(Blakeley_gold1,filename = "D:/software/RStudio/single_cell/datasets_unify/Blakeley3_nor+pca+specc_cluster1.pdf",width = 12,height = 9)
ggsave(Blakeley_goldd2,filename = "D:/software/RStudio/single_cell/datasets_unify/Blakeley3_nor+pca+specc_cluster2.pdf",width = 12,height = 9)
