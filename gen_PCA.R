#!/usr/bin/env Rscript
.libPaths(c("/data/wangx53",.libPaths()))
library(flashpcaR)
oncopca=flashpca("/data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/euro_pruned", ndim=10,verbose = T)
icogspca=flashpca("/data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/euro_pruned", ndim=10,verbose = T)
save(oncopca,icogspca,file="../result/gen_PCA.RData")
