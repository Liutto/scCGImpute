#' Title
#'
#' @param count_matrix A matrix specifying the scRNA-seq data
#' @param out_dir A character specifying the path of the output directory
#' @param Kcluster  An integer specifying the number of cell subpopulations
#'
#' @export
#' @import foreach
#' @import parallel
#' @import penalized
#' @import survival
#' @import doParallel
#' @import foreach
#' @import iterators
#' @import kernlab
#' @import randomForest
#' @import rpca
scCGImpute <-
  function (count_matrix, out_dir, Kcluster = NULL)
  {

    print("reading in raw count matrix ...")
    dir.create(out_dir, recursive = TRUE)

    raw_count = as.matrix(count_matrix)
    zeroproportion = length(which(raw_count ==0))/(dim(raw_count)[1]*dim(raw_count)[2])
    print(paste("number of genes in raw count matrix", nrow(raw_count)))
    print(paste("number of cells in raw count matrix", ncol(raw_count)))
    totalCounts_by_cell = colSums(raw_count)
    totalCounts_by_cell[totalCounts_by_cell == 0] = 1
    raw_count = sweep(raw_count, MARGIN = 2, 10^6/totalCounts_by_cell, FUN = "*")
    if (min(raw_count) < 0) {
      stop("smallest read count cannot be negative!")
    }
    count_lnorm = log10(raw_count + 1.01)

    print("reading finished!")
    genenames = rownames(count_lnorm)
    cellnames = colnames(count_lnorm)
    count = as.matrix(count_lnorm)

    print("imputation starts ...")
    count_imp = imputation_model(count,
                         point = log10(1.01), drop_thre = 0.5,
                         Kcluster = Kcluster,
                         out_dir = out_dir, ncores = ncores)

    # outliers = res_imp$outlier
    count_imp = 10^count_imp - 1.01
    rownames(count_imp) = genenames
    colnames(count_imp) = cellnames
    print("writing imputed count matrix ...")
    count_imp = sweep(count_imp, MARGIN = 2, totalCounts_by_cell/10^6,
                      FUN = "*")
    count_imp = round(count_imp, digits = 2)
    imp_zeropro = length(which(count_imp==0))/(dim(count_imp)[1]*dim(count_imp)[2])
    write.csv(count_imp, file = paste0(out_dir, "scCGImpute_count.csv"))

    print(paste0("The original ratio of zero :",zeroproportion))
    print(paste0("The ratio of zero after imputaion:",imp_zeropro))
    print(paste0("reduce zero's percentage of ",round(zeroproportion - imp_zeropro,2)*100,"% "))
  }

