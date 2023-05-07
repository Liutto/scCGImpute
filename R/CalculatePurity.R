# ClusterPurity <- function(clusters, classes) {
#   sum(apply(table(classes, clusters), 2, max)) / length(clusters)
# }
CalculatePurity<-function(true_lab,pre_lab){
  N=length(true_lab)
  v=unique(true_lab)
  u=unique(pre_lab)
  sumnum=0
  for (i in 1:length(v)) {
    maxnum=0
    for (j in 1:length(u)) {
      num=length(intersect(which(pre_lab==u[i]),which(true_lab==v[j])))
      if(num>maxnum){
        maxnum=num
      }
    }
    sumnum=sumnum+maxnum
  }
  return(sumnum/N)
}