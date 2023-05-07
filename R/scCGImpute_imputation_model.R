find_hv_genes = function(count, I, J){
  count_nzero = lapply(1:I, function(i) setdiff(count[i, ], log10(1.01)))
  mu = sapply(count_nzero, mean) #默认对列求均值
  mu[is.na(mu)] = 0
  sd = sapply(count_nzero, sd)  #计算标准差
  sd[is.na(sd)] = 0 
  cv = sd/mu    #标准差除以样本平均值,得到的是变异系数 变量值平均水平高，其离散程度的测度值越大，反之越小
  cv[is.na(cv)] = 0
  high_var_genes = which(mu >= 1 & cv >= quantile(cv, 0.25))  #quantile求百分位比
  if(length(high_var_genes) < 500){ 
    high_var_genes = 1:I}
  count_hv = count[high_var_genes, ]
  return(count_hv)
}

find_neighbors = function(count_hv, J, Kcluster = NULL, 
                          ncores, cell_labels = NULL){

    print("dimension reduction ...")
    if(J < 5000){
      var_thre = 0.4
      pca = prcomp(t(count_hv))
    }else{
      var_thre = 0.6
      pca = rpca(t(count_hv), k = 1000, center = TRUE, scale = FALSE) 
    }
  eigs = (pca$sdev)^2
  var_cum = cumsum(eigs)/sum(eigs)
  if(max(var_cum) <= var_thre){
    npc = length(var_cum)
  }else{
    npc = which.max(var_cum > var_thre)
    npc = max(npc, Kcluster)
  }
  
    if (npc < 3){ npc = 3 }
    mat_pcs = t(pca$x[, 1:npc]) # columns are cells
    clnum<-detectCores()  
    cl <- makeCluster(getOption("cl.cores", clnum)) 
    
    print("calculating cell distances ...")
    dist_cells_list = parLapply(cl,1:J, function(id1){
      d = sapply(1:id1, function(id2){
        sse = sum((mat_pcs[, id1] - mat_pcs[, id2])^2)
        sqrt(sse)
      })
      return(c(d, rep(0, J-id1)))
    })
    dist_cells = matrix(0, nrow = J, ncol = J)
    for(cellid in 1:J){dist_cells[cellid, ] = dist_cells_list[[cellid]]}
    dist_cells = dist_cells + t(dist_cells)
    
    stopCluster(cl)
    min_dist = sapply(1:J, function(i){
      min(dist_cells[i, -i])
    })
    iqr = quantile(min_dist, 0.75) - quantile(min_dist, 0.25)
    outliers = which(min_dist > 1.5 * iqr + quantile(min_dist, 0.75))
    
    ## clustering
    non_out = setdiff(1:J, outliers)
    spec_res = specc(t(mat_pcs[, non_out]), centers = Kcluster, kernel = "rbfdot")
    nbs = rep(NA, J)
    nbs[non_out] = spec_res
    
    return(list(dist_cells = dist_cells, clust = nbs))
}

find_va_genes = function(parslist, subcount){   
  point = log10(1.01)
  valid_genes = which( (rowSums(subcount) > point * ncol(subcount)) &
                         complete.cases(parslist) )   #complete.cases返回一个逻辑向量，指示哪些情况是完整的，即没有缺失值
  if(length(valid_genes) == 0) return(valid_genes)
  # find out genes that violate assumption
  mu = parslist[, "mu"]
  sgene1 = which(mu <= log10(1+1.01))
  # sgene2 = which(mu <= log10(10+1.01) & mu - parslist[,5] > log10(1.01))
  
  dcheck1 = dgamma(mu+1, shape = parslist[, "alpha"], rate = parslist[, "beta"])
  dcheck2 = dnorm(mu+1, mean = parslist[, "mu"], sd = parslist[, "sigma"])
  sgene3 = which(dcheck1 >= dcheck2 & mu <= 1)
  sgene = union(sgene1, sgene3)
  valid_genes = setdiff(valid_genes, sgene)#求两个向量中不同的元素（只取x中不同的元素）
  return(valid_genes)
}
calculate_weight =function (x, paramt) {
    pz1 = paramt[1] * dgamma(x, shape = paramt[2], rate = paramt[3])
    pz2 = (1 - paramt[1]) * dnorm(x, mean = paramt[4], sd = paramt[5])
    pz = pz1/(pz1 + pz2)#pz是pz1的条件概率
    pz[pz1 == 0] = 0
    return(cbind(pz, 1 - pz))
}
dmix =function (x, pars) {
    pars[1] * dgamma(x, shape = pars[2], rate = pars[3]) + (1 - 
                                                              pars[1]) * dnorm(x, mean = pars[4], sd = pars[5])
}
### root-finding equation
fn = function(alpha, target){
  log(alpha) - digamma(alpha) - target
}

### update parameters in gamma distribution
update_gmm_pars = function(x, wt){
  tp_s = sum(wt)  #条件概率的和
  tp_t = sum(wt * x) 
  tp_u = sum(wt * log(x))
  tp_v = -tp_u / tp_s - log(tp_s / tp_t)
  if (tp_v <= 0){
    alpha = 20
  }else{
    alpha0 = (3 - tp_v + sqrt((tp_v - 3)^2 + 24 * tp_v)) / 12 / tp_v
    if (alpha0 >= 20){alpha = 20
    }else{
      alpha = uniroot(fn, c(0.9, 1.1) * alpha0, target = tp_v, 
                      extendInt = "yes")$root
    }
  }
  ## need to solve log(x) - digamma(x) = tp_v
  ## We use this approximation to compute the initial value
  beta = tp_s / tp_t * alpha   #gamma分布  E(X)=alpha/beta
  return(c(alpha, beta))
}

### estimate parameters in the mixture distribution
get_mix = function(xdata, point){
  # print("get_mix start...")
  inits = rep(0, 5)
  inits[1] = sum(xdata == point)/length(xdata)
  if (inits[1] == 0) {inits[1] = 0.01}
  inits[2:3] = c(0.5, 1)
  xdata_rm = xdata[xdata > point]
  inits[4:5] = c(mean(xdata_rm), sd(xdata_rm))
  if (is.na(inits[5])) {inits[5] = 0}
  paramt = inits
  eps = 10
  iter = 0
  loglik_old = 0
  
  while(eps > 0.5) {
    wt = calculate_weight(xdata, paramt) #gamma分布和normal分布的条件概率，第一列是gamma分布
    paramt[1] = sum(wt[, 1])/nrow(wt)    #计算权重
    paramt[4] = sum(wt[, 2] * xdata)/sum(wt[, 2])#计算normal分布的均值
    paramt[5] = sqrt(sum(wt[, 2] * (xdata - paramt[4])^2)/sum(wt[, 2]))#计算normal分布的方差
    paramt[2:3] = update_gmm_pars(x=xdata, wt=wt[,1])  #计算gamma分布的参数
    
    loglik = sum(log10(dmix(xdata, paramt)))
    eps = (loglik - loglik_old)^2
    loglik_old = loglik
    iter = iter + 1
    if (iter > 100) 
      break
  }
  return(paramt)
}
get_mix_parameters =function (count, point = log10(1.01), path, ncores = 8) 
  {
    print("get_mix_parameters start!")
    count = as.matrix(count)
    null_genes = which(abs(rowSums(count) - point * ncol(count)) < 1e-10)
    no_cores <- detectCores(no_cores) - 1  # 检测核数，并设置使用的数量
    cl = makeCluster(no_cores)
    registerDoParallel(cl)
    parslist =foreach(ii=1:nrow(count))%do%{
      if (ii %in% null_genes) {
        return(rep(NA, 5))
      }
      xdata = count[ii, ]
      paramt = try(get_mix(xdata, point), silent = TRUE)
      if (class(paramt) == "try-error"){
        paramt = rep(NA, 5)
      }
      paramt
    }
    stopCluster(cl)
    save(parslist, file = path)
    parslist = Reduce(rbind, parslist)#Reduce使用一个二元函数依次组合给定向量的元素和可能给定的初始值
    colnames(parslist) = c("rate", "alpha", "beta", "mu", "sigma")
    saveRDS(parslist, file = path)
    print("get_mix_parameters end!")
    return(0)
}
rmix =function (pars, n) {
  n1 = round(n * pars[1])
  n2 = n - n1
  x1 = rgamma(n1, shape = pars[2], rate = pars[3])
  x2 = rnorm(n2, mean = pars[4], sd = pars[5])
  return(c(x1, x2))
}


imputation_model = function(count, point, drop_thre = 0.5, Kcluster = Kcluster, 
                             out_dir, ncores){
  count = as.matrix(count)
  I = nrow(count)
  J = ncol(count)
  count_imp = count
  tryerror = 0

  #pca  cell
  # find highly variable genes -> 表达量差别最大的基因
  count_hv = find_hv_genes(count, I, J)  #genes×cells
  print("searching candidate neighbors ... ")
  if(Kcluster == 1){
    clust = rep(1, J)
    if(J < 5000){
      var_thre = 0.4
      pca = prcomp(t(count_hv))#主成分分析   t() ->cells×genes
    }else{
      var_thre = 0.6
      pca = rpca(t(count_hv), k = 1000, center = TRUE, scale = FALSE) 
      #rpca 随机主成分分析 利用随机奇异值分解快速计算主成分分析
    }
    eigs = (pca$sdev)^2  #sdev主成分的标准差(即协方差/相关矩阵特征值的平方根)
    var_cum = cumsum(eigs)/sum(eigs)
    if(max(var_cum) <= var_thre){ 
      npc = length(var_cum)
    }else{
      npc = which.max(var_cum > var_thre)
      npc = max(npc, Kcluster)
    }
    if (npc < 3){ npc = 3 } #至少选择3列的主成分
    mat_pcs = t(pca$x[, 1:npc]) # t()之后columns are cells
    
    clnum<-detectCores()  ##detectCores函数可以告诉你你的CPU可使用的核数
    cl <- makeCluster(getOption("cl.cores", clnum)) ##设置参与并行的CPU核数目，这里我们使用了所有的CPU核，也就是我们刚才得到的clnum，具体到这个案例

    dist_cells_list = parLapply(cl,1:J, function(id1){
      d = sapply(1:id1, function(id2){
        sse = sum((mat_pcs[, id1] - mat_pcs[, id2])^2)
        sqrt(sse)
      })
      return(c(d, rep(0, J-id1)))
    })
    dist_cells = matrix(0, nrow = J, ncol = J)
    View(dist_cells_list)
    for(cellid in 1:J){dist_cells[cellid, ] = dist_cells_list[[cellid]]}
    dist_cells = dist_cells + t(dist_cells)

    stopCluster(cl)
  }else{
    print("inferring cell similarities ...")
    set.seed(Kcluster)
    neighbors_res = find_neighbors(count_hv = count_hv, J = J, 
                                   Kcluster = Kcluster, ncores = ncores)
    dist_cells = neighbors_res$dist_cells
    clust = neighbors_res$clust
  }
  saveRDS(clust, file = paste0(out_dir, "clust.rds"))
  # mixture model
  nclust = sum(!is.na(unique(clust)))#unique删除重复元素

  for(cc in 1:nclust){
    print(paste("estimating dropout probability for type", cc, "..."))
    paste0(out_dir, "pars", cc, ".rds")
    get_mix_parameters(count = count[, which(clust == cc), drop = FALSE],
                       point = log10(1.01),
                       path = paste0(out_dir, "pars", cc, ".rds"), ncores = ncores)

    print("get_mix_parameters end!")
    cells = which(clust == cc)
    if(length(cells) <= 1) { next }
    parslist = readRDS(paste0(out_dir, "pars", cc, ".rds"))

    print("searching for valid genes ...")
    valid_genes = find_va_genes(parslist, subcount = count[, cells])#分别对每一类型细胞进行处理
    print("searching for valid genes ...end")
    if(length(valid_genes) <= 10){ next }
    subcount = count[valid_genes, cells, drop = FALSE]
    Ic = length(valid_genes)
    Jc = ncol(subcount)  #该类细胞的数目
    parslist = parslist[valid_genes, , drop = FALSE]

    print("droprate...")

    #缺失率   一个基因一个基因地计算
    droprate = t(sapply(1:Ic, function(i) {
      wt = calculate_weight(subcount[i, ], parslist[i, ])
      return(wt[, 1])
    }))
    #大于均值的为TRUE
    mucheck = sweep(subcount, MARGIN = 1, parslist[, "mu"], FUN = ">")
    droprate[mucheck & droprate > drop_thre] = 0    #将均值较大但dropout值>阈值的dropout置为0

    # dropouts
    setA = lapply(1:Ic, function(geneid){
      which(droprate[geneid,] > drop_thre)
    })
    # non-dropouts
    setB = lapply(1:Ic, function(geneid){
      which(droprate[geneid,] <= drop_thre)
    })

    print(paste("imputing dropout values for type", cc, "..."))
    # imputation
    gc() 
    #以平均值填充
    print("cell imputation")
    for(geneid in 1:Ic){
         tmpeffic =which(droprate[geneid,]<drop_thre)
         tmp = subcount[geneid,tmpeffic]
         tmp[which(tmp<=log10(1.01))]==0
         sumsubcout=sum(tmp);
         average=sumsubcout/(length(tmpeffic));

         tmpdropthre=which(droprate[geneid,]>=drop_thre)  #dropout大于阈值
         tmpaverage=which(subcount[geneid,]<average)    #小于平均值
         tmpintersect=intersect(tmpdropthre,tmpaverage)   #取交集

         subcount[geneid,tmpintersect]=average;
    
    }
    #randomForest regression
    print("randomForest regression start")
    pearsongene<-matrix(0,nrow = dim(subcount)[1],ncol = dim(subcount)[1])
    pearsongene = cor(t(subcount))
    gc()
  
    for(geneid in 1:Ic) {
      geneloc<-setdiff(which(abs(pearsongene[,geneid])>=0.5),geneid)
      cellid_drop = setA[[geneid]]
      cellid_obs = setB[[geneid]]
      if(length(cellid_drop) ==0 ||length(cellid_obs)<=1||length(geneloc) <= 1){
        next}
      if(dim(t(subcount[geneloc,cellid_obs]))[1]<dim(t(subcount[geneloc,cellid_obs]))[2]){
        next}
      rf.train<-randomForest(t(subcount[geneloc,cellid_obs]),subcount[geneid,cellid_obs])
      rf.test <- predict(rf.train,newdata=t(subcount[geneloc,cellid_drop]))

      rf.test[rf.test<0]=0
      subcount[geneid,cellid_drop]=rf.test        
    }
    count_imp[valid_genes, cells] = subcount
     }

  count_imp[count_imp < point] = point
  return(count_imp = count_imp)
}

