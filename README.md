# scCGImpute: An imputation method for single-cell RNA se-quencing data based on similarities between cells and rela-tionships among genes
Tiantian Liu,Yuanyuan Li

## Introduction
Single-cell RNA sequencing (scRNA-seq) has become a powerful technique to investigate cellular heterogeneity and complexity in various fields by revealing the gene expression status of individual cells. Despite the undeniable benefits of scRNA-seq, it is not immune to its inherent limitations, such as sparsity and noise, which would hinder downstream analysis.We introduce scCGImpute to address the challenges of sparsity in scRNA-seq data through imputation. 

## Quick start
scCGImpute(

  count_matrix = count_matrix, #scRNA-seq data
  
  out_dir = out_dir,           #full path to output directory
  
  Kcluster = 3                 #the number of cell type
  
  )    
