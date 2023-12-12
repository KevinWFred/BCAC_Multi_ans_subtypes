#!/usr/bin/env Rscript
.libPaths(c("/data/wangx53",.libPaths()))
library(flashpcaR)
icogspca=flashpca("/data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/merged1_pruned", ndim=10,verbose = T)
save(icogspca,file="../result/icogs_flashpca.RData")
oncopca=flashpca("/data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/merged1_pruned", ndim=10,verbose = T)
save(oncopca,file="../result/onco_flashpca.RData")

oncopca=flashpca("/data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/euro_pruned", ndim=10,verbose = T)
icogspca=flashpca("/data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/euro_pruned", ndim=10,verbose = T)
save(oncopca,icogspca,file="../result/gen_PCA.RData")
tmp=icogspca$projection
colnames(tmp)=paste0("pc",1:10)
tmp1=unlist(strsplit(rownames(tmp),":"))
tmp1=unlist(strsplit(tmp1,"_"))
rownames(tmp)=tmp1[seq(1,length(tmp1),4)]

bcacpca=flashpca("/data/DCEGLeiSongData/Kevin/tmp/BCAC_pruned", ndim=10,verbose = T)
save(oncopca,icogspca,bcacpca,file="../result/gen_PCA.RData")


         