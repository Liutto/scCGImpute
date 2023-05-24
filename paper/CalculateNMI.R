CalculateNMI<-function(true_lab,pre_lab){
  N=length(true_lab)
  v=unique(true_lab)
  u=sort(unique(pre_lab))
  #联合概率
  jointprobability<-matrix(0,nrow = length(u),ncol = length(v))
  for (i in 1:length(u)) {
    for (j in 1:length(v)) {
      jointprobability[i,j]=length(intersect(which(pre_lab==u[i]),which(true_lab==v[j])))/N
    }
  }
  # print(jointprobability)
  
  eps=1.4e-100
  #
  A=rowSums(jointprobability)
  B=colSums(jointprobability)
  Hx=-sum(mapply('*',A,log2(A+eps)))
  Hy=-sum(mapply('*',B,log2(B+eps)))
  # print(Hx)
  # print(Hy)
  C=-sum(mapply('*',jointprobability,log2(jointprobability+eps)))
  # print(C)
  # print("MI")
  # print(Hx+Hy-C)
  return(2*(Hx+Hy-C)/(Hx+Hy))
}