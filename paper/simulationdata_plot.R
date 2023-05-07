############
#Splat simulation parameters
library("splatter")
library("scater")
library("ggplot2")
library("pheatmap")
##########
#dropout.mid -> 0-4.5
#生成数据
dropoutrate<-c()
for (value in seq(from=0,by=0.5,length=10)){
  data_dir="D:/software/RStudio/single_cell/datasets/simulate_splat_dropouttest/"
  params.groups <- newSplatParams(batchCells = 300, nGenes = 1000)
  eval(parse(text = paste0("sim1 <- splatSimulateGroups(params.groups, group.prob = c(0.3, 0.3 , 0.4), de.prob =  c(0.1, 0.1, 0.1),dropout.type = 'batch',dropout.mid = ",value,",verbose = FALSE)")))
  
  dropoutdata1<-sim1@assays@data@listData[["TrueCounts"]]
  refdata<-sim1@assays@data@listData[["TrueCounts"]]
  dropoutdata1[which(sim1@assays@data@listData[["Dropout"]]==TRUE)]=0
  simulate_splat_labels = sim1@colData@listData[["Group"]]
  simulate_splat_labels<-as.character(simulate_splat_labels)  #将factor转换成字符串
  write.csv(simulate_splat_labels,"D:/software/RStudio/single_cell/datasets/simulate_splat_dropouttest/simulate_splat_labels.csv")
  eval(parse(text = paste0("write.csv(dropoutdata1,paste0(data_dir,'simulate_splat_",value,".csv'))")))
  eval(parse(text = paste0("write.csv(refdata,paste0(data_dir,'simulate_splat_nodropout_",value,".csv'))")))
  dropoutratet=length(which(sim1@assays@data@listData[["Dropout"]]==TRUE))/300000
  dropoutrate<-c(dropoutrate,dropoutratet)
}

#聚类评价 折线图
for (value in seq(from=0,by=0.5,length=10)) {
  print(value)
  eval(parse(text = paste0("data_dir='D:/software/RStudio/single_cell/datasets/simulate_splat_dropouttest/simulate_splat_",value,"/'")))
  method<-c("saver","deepimpute","VIPER","scRMD","ALRA")#
  eval(parse(text = paste0("rawdata<-read.csv(paste0(data_dir,'simulate_splat_",value,".csv'))")))
  eval(parse(text = paste0("refdata<-read.csv(paste0(data_dir,'simulate_splat_nodropout_",value,".csv'))")))
  scImpute<-read.csv(paste0(data_dir,"simulate_splat_scImpute/scimpute_count.csv"))
  # scCGImpute<-read.csv(paste0(data_dir,"simulate_splat_scCGImpute/scCGImpute_count.csv"))
  scCGImpute<-read.csv(paste0(data_dir,"s_scCGImpute/scCGImpute_count.csv"))
  scRecover<-read.csv(paste0(data_dir,"outDir_scRecover/scRecover+scImpute.csv"))
  for (i in 1:length(method)) {
    eval(parse(text = paste0(method[i],"<-read.csv(paste0(data_dir,'simulate_splat_",method[i],".csv'))")))
  }
  simulate_splat_labels<-read.csv("D:/software/RStudio/single_cell/datasets/simulate_splat_dropouttest/simulate_splat_labels.csv")
  simulate_splat_labels<-simulate_splat_labels[,-1]
  simulate_splat_labelschar<-simulate_splat_labels
  for (i in 1:3) {
    eval(parse(text = paste0("simulate_splat_labels[which(simulate_splat_labels=='Group",i,"')]=i")))
  }
  simulate_splat_labels = as.numeric(unlist(simulate_splat_labels))   #将列表转换为向量。它通过保留所有组件简化了生成向量的过程
  source("D:/software/RStudio/single_cell/scCGImpute/code/CalculateARI.R")
  source("D:/software/RStudio/single_cell/scCGImpute/code/MatrixFusion.R")
  source("D:/software/RStudio/single_cell/scCGImpute/code/CalculateNMI.R")
  source("D:/software/RStudio/single_cell/scCGImpute/code/CalculatePurity.R")
  
  # methodss<-c("rawdata","saver","scImpute","scCGImpute","deepimpute","VIPER","scRMD","refdata","ALRA","scRecover","s_scCGImpute","s_nolimit_scCGImpute","s_nolimitnoiter_scCGImpute","s_nolimit10_scCGImpute","s10_scCGImpute")
  methodss<-c("rawdata","saver","scImpute","deepimpute","VIPER","scRMD","refdata","scRecover","scCGImpute","ALRA")
  
  dataend<-c(rep(0,length(methodss)*3))
  for (methodssi in 1:length(methodss)) {
    eval(parse(text = paste0(methodss[methodssi],value,"ARI<-0")))
    eval(parse(text = paste0(methodss[methodssi],value,"NMI<-0")))
    eval(parse(text = paste0(methodss[methodssi],value,"Purity<-0")))
  }
  for (i in 1:10) {
    print(i)
    gc()
    for (methodssi in 1:length(methodss)) {
      print(methodssi)
      gc()
      eval(parse(text = paste0("spec_res",methodss[methodssi],"=specc(t(",methodss[methodssi],"[,-1]), centers = 3, kernel = 'rbfdot')")))
      eval(parse(text = paste0(methodss[methodssi],"ARI<-CalculateARI(simulate_splat_labels,spec_res",methodss[methodssi],"@.Data)")))
      eval(parse(text = paste0(methodss[methodssi],"NMI<-CalculateNMI(simulate_splat_labels,spec_res",methodss[methodssi],"@.Data)")))
      eval(parse(text = paste0(methodss[methodssi],"Purity<-CalculatePurity(simulate_splat_labels,spec_res",methodss[methodssi],"@.Data)")))
      
      eval(parse(text = paste0(methodss[methodssi],value,"ARI<-",methodss[methodssi],value,"ARI+",methodss[methodssi],"ARI")))
      eval(parse(text = paste0(methodss[methodssi],value,"NMI<-",methodss[methodssi],value,"NMI+",methodss[methodssi],"NMI")))
      eval(parse(text = paste0(methodss[methodssi],value,"Purity<-",methodss[methodssi],value,"Purity+",methodss[methodssi],"Purity")))
      
    }   
    
  }
}
assessmethod<-c("ARI","NMI","Purity")
for (j in 1:length(assessmethod)) {
  for (i in 1:length(methodss)) {
    eval(parse(text = paste0(methodss[i],"<-c()")))
    for (value in seq(from=0,by=0.5,length=10)) {
      eval(parse(text = paste0(methodss[i],"<-cbind(",methodss[i],",",methodss[i],value,assessmethod[j],")")))
    }
    eval(parse(text = paste0("name<-rep('",methodss[i],"',length(seq(from=0,by=0.5,length=10)))")))
    filll=seq(from=0,by=0.5,length=10)
    for (q in 1:length(seq(from=0,by=0.5,length=10))) {
      filll[q]=paste0(filll[q],"(",round(dropoutrate[q],3)*100,"%)")
    }
    eval(parse(text = paste0(methodss[i],"bar<-rbind(t(name),t(filll),as.numeric(",methodss[i],"))")))
  }
  
  eval(parse(text = paste0(assessmethod[j],"bar<-",methodss[1],"bar")))
  for (i in 2:length(methodss)) {
    eval(parse(text = paste0(assessmethod[j],"bar<-cbind(",assessmethod[j],"bar,",methodss[i],"bar)")))
  }
  eval(parse(text = paste0(assessmethod[j],"bar<-t(",assessmethod[j],"bar)")))
  eval(parse(text = paste0(assessmethod[j],"bar<-data.frame(",assessmethod[j],"bar)")))
  eval(parse(text = paste0(assessmethod[j],"bar[,3] = as.numeric(",assessmethod[j],"bar[,3])/10")))
  
  eval(parse(text = paste0("colnames(",assessmethod[j],"bar)<-c('methods','dropout.rate','value')")))
  
}
##折线图
p1<-ggplot(ARIbar,aes(x=dropout.rate,y=value,group = factor(methods,levels = c("scCGImpute","rawdata","refdata","saver","scImpute","deepimpute","VIPER","scRMD","scRecover","ALRA")),color=factor(methods,levels = c("scCGImpute","rawdata","refdata","saver","scImpute","deepimpute","VIPER","scRMD","scRecover","ALRA")),shape=factor(methods,levels = c("scCGImpute","rawdata","refdata","saver","scImpute","deepimpute","VIPER","scRMD","scRecover","ALRA"))))+
  geom_point()+
  geom_line()+
  xlab("dropout rate")+#横坐标名称
  ylab("ARI value")+#纵坐标名称
  easy_add_legend_title("Methods")+
  theme_bw() #去掉背景灰色
p2<-ggplot(NMIbar,aes(x=dropout.rate,y=value,group = factor(methods,levels = c("scCGImpute","rawdata","refdata","saver","scImpute","deepimpute","VIPER","scRMD","scRecover","ALRA")),color=factor(methods,levels = c("scCGImpute","rawdata","refdata","saver","scImpute","deepimpute","VIPER","scRMD","scRecover","ALRA")),shape=factor(methods,levels = c("scCGImpute","rawdata","refdata","saver","scImpute","deepimpute","VIPER","scRMD","scRecover","ALRA"))))+
  geom_point()+
  geom_line()+
  xlab("dropout rate")+#横坐标名称
  ylab("NMI value")+#纵坐标名称
  easy_add_legend_title("Methods")+
  theme_bw() #去掉背景灰色
p3<-ggplot(Puritybar,aes(x=dropout.rate,y=value,group = factor(methods,levels = c("scCGImpute","rawdata","refdata","saver","scImpute","deepimpute","VIPER","scRMD","scRecover","ALRA")),color=factor(methods,levels = c("scCGImpute","rawdata","refdata","saver","scImpute","deepimpute","VIPER","scRMD","scRecover","ALRA")),shape=factor(methods,levels = c("scCGImpute","rawdata","refdata","saver","scImpute","deepimpute","VIPER","scRMD","scRecover","ALRA"))))+
  geom_point()+
  geom_line()+
  xlab("dropout rate")+#横坐标名称
  ylab("Purity value")+#纵坐标名称
  easy_add_legend_title("Methods")+
  theme_bw() #去掉背景灰色
p1/p2/p3+ plot_layout(guides = "collect")
ggsave(p1,filename = "D:/software/RStudio/single_cell/datasets/simulate_splat_dropouttest/ARIresult+alras_scCGImpute.pdf",width = 12,height = 9)
ggsave(p2,filename = "D:/software/RStudio/single_cell/datasets/simulate_splat_dropouttest/NMIresult+alras_scCGImpute.pdf",width = 12,height = 9)
ggsave(p3,filename = "D:/software/RStudio/single_cell/datasets/simulate_splat_dropouttest/Purityresult+alras_scCGImpute.pdf",width = 12,height = 9)
data_dir = "D:/software/RStudio/single_cell/datasets/simulate_splat_dropouttest/"
write.csv(ARIbar,paste0(data_dir,"ARIresult+alra_s_scCGImpute.csv"))
write.csv(NMIbar,paste0(data_dir,"NMIresults+alra_scCGImpute.csv"))
write.csv(Puritybar,paste0(data_dir,"Purityresults+alra_scCGImpute.csv"))

########
#相关性
corp<-function(ref,methodss){
  corpeason<-c()
  for (i in 2:dim(ref)[2]) {
    corpeason<-c(corpeason,cor(ref[,i],methodss[,i]))
  }
  return(corpeason)
}
corpgene<-function(ref,methodss){
  corpeason<-c()
  for (i in 1:dim(ref)[1]) {
    corpeason<-c(corpeason,cor(c(t(ref[i,2:dim(ref)[2]])),c(t(methodss[i,2:dim(ref)[2]]))))
  }
  return(corpeason)
}
methodss<-c("refdata","scCGImpute")
datalong<-c()
datagenelong<-c()
for (value in seq(from=0,by=0.5,length=10)) {
  print(value)
  eval(parse(text = paste0("data_dir='D:/software/RStudio/single_cell/datasets/simulate_splat_dropouttest/simulate_splat_",value,"/'")))
  eval(parse(text = paste0("rawdata<-read.csv(paste0(data_dir,'simulate_splat_",value,".csv'))")))
  eval(parse(text = paste0("refdata<-read.csv(paste0(data_dir,'simulate_splat_nodropout_",value,".csv'))")))
  scCGImpute<-read.csv(paste0(data_dir,"scCGImpute/scCGImpute_count.csv"))

  for (i in 1:length(method)) {
    eval(parse(text = paste0(method[i],"<-read.csv(paste0(data_dir,'simulate_splat_",method[i],".csv'))")))
  }
  
  for (methodssi in 1:length(methodss)) {
    print(methodss[methodssi])
    eval(parse(text = paste0("cor",methodss[methodssi],"=corp(rawdata,",methodss[methodssi],")")))
    eval(parse(text = paste0("cor",methodss[methodssi],"gene=corpgene(rawdata,",methodss[methodssi],")")))
    
  }
  data=cbind(correfdata,corscCGImpute)
  data=data.frame(data)
  colnames(data)<-c("refdata","scCGImpute")
  data <- data %>% 
  gather(key = 'methods',value = 'correlation') #命名两个列名
  datagene=cbind(correfdatagene,corscCGImputegene)
  datagene=data.frame(datagene)
  colnames(datagene)<-c("refdata","scCGImpute")
  datagene <- datagene %>% 
    gather(key = 'methods',value = 'correlation') #命名两个列名
 
  eval(parse(text = paste0("values<-rep(paste0('",value,"(",round(dropoutrate[value/0.5+1],3)*100,"%)'),length(correfdata)*2)")))
  datatmp = cbind(values,data)
  datalong = rbind(datalong,datatmp)
  eval(parse(text = paste0("values<-rep('",value,"(",round(dropoutrate[value/0.5+1],3)*100,"%)',length(correfdatagene)*2)")))
  datatmp = cbind(values,datagene)
  datagenelong = rbind(datagenelong,datatmp)
  
  p <- ggplot(data,mapping = aes(x=' ',correlation,fill = factor(methods,levels = c("refdata","scCGImpute"))))
  pgene <- ggplot(datagene,mapping = aes(x=' ',correlation,fill = factor(methods,levels = c("refdata","scCGImpute"))))
  
  eval(parse(text = paste0("cor1 <- p + geom_boxplot(outlier.shape = 21,width=1.2,outlier.colour = 'red',outlier.fill = 'blue',
                           col = 'black')+
    labs(title = 'Cell drop.mid=",value,"',x = '')+
    theme(plot.title = element_text(size = 20))+
    theme_bw()+
    easy_add_legend_title('Methods')")))
  eval(parse(text = paste0("cor1",value,"<-cor1")))
  
  eval(parse(text = paste0("cor2 <- pgene + geom_boxplot(outlier.shape = 21,width=1.2,outlier.colour = 'red',outlier.fill = 'blue',
                           col = 'black')+
    labs(title = 'Gene drop.mid=",value,"',x='')+
    theme(plot.title = element_text(size = 20))+
    theme_bw()+
    easy_add_legend_title('Methods')")))
  eval(parse(text = paste0("cor2",value,"<-cor2")))
  cor1+cor2+ plot_layout(guides = "collect")
}
pd<-ggboxplot(datalong,x="values",y="correlation",fill ="methods",palette = 'mycol')
dataercclong<-pd+stat_compare_means(aes(group = methods),method = 't.test',label = 'p.signif')+
  easy_add_legend_title('Methods')
dataercclong
pdg<-ggboxplot(datagenelong,x="values",y="correlation",fill ="methods",palette = 'mycol')
datafeneercclong<-pdg+stat_compare_means(aes(group = methods),method = 't.test',label = 'p.signif')+
  easy_add_legend_title('Methods')
datafeneercclong




