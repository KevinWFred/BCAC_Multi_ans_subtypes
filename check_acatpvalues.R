#!/usr/bin/env Rscript

oldres=read.csv("../result/allnovelsnps.csv")
load("../result/metascoreinfo4_newQC.RData")
metascoreinfo4_new=metascoreinfo4
load("../result/metascoreinfo4.RData")

load("../result/compute_metapvalues_newQC.RData")
intrinsicpvalues_new=intrinsicpvalues
scoretestallpvaues_new=scoretestallpvaues
allpvalues_new=allpvalues
load("../result/compute_metapvalues_new.RData")
idx=match(oldres$ID,names(allpvalues))
oldres$acat=allpvalues[idx]
idx=match(oldres$ID,names(allpvalues_new))
oldres$acat_new=allpvalues_new[idx]
idx=match(oldres$ID,rownames(intrinsicpvalues))
oldres$intr=intrinsicpvalues[idx,1]
idx=match(oldres$ID,rownames(intrinsicpvalues_new))
oldres$intr_new=intrinsicpvalues_new[idx,1]
idx=match(oldres$ID,rownames(scoretestallpvaues))
oldres$score=scoretestallpvaues[idx,1]
idx=match(oldres$ID,rownames(scoretestallpvaues_new))
oldres$score_new=scoretestallpvaues_new[idx,1]
idx=match(oldres$ID,rownames(metascoreinfo4$avail))
oldres=cbind(oldres,metascoreinfo4$avail[idx,])
idx=match(oldres$ID,rownames(metascoreinfo4_new$avail))
oldres=cbind(oldres,metascoreinfo4_new$avail[idx,])
# idx=match(oldres$ID,rownames(metascoreinfo4_new$allfreq))
# oldres=cbind(oldres,metascoreinfo4_new$allfreq[idx,])
euro_onco_info=as.data.frame(fread("/data/BB_Bioinformatics/Kevin/BCAC/data/imp_infoscore/data/euro_onco/myinfoscore.txt"))
euro_icogs_info=as.data.frame(fread("/data/BB_Bioinformatics/Kevin/BCAC/data/imp_infoscore/data/euro_icogs/myinfoscore.txt"))
asian_icogs_info=as.data.frame(fread("/data/BB_Bioinformatics/Kevin/BCAC/data/imp_infoscore/data/asian_icogs/myinfoscore.txt"))
asian_onco_info=as.data.frame(fread("/data/BB_Bioinformatics/Kevin/BCAC/data/imp_infoscore/data/asian_onco/myinfoscore.txt"))
african_onco_info=as.data.frame(fread("/data/BB_Bioinformatics/Kevin/BCAC/data/imp_infoscore/data/african_onco/myinfoscore.txt"))
idx=match(oldres$ID,euro_onco_info$ID)
oldres$euro_onco_info=euro_onco_info$R2[idx]
idx=match(oldres$ID,euro_icogs_info$ID)
oldres$euro_icogs_info=euro_icogs_info$R2[idx]
idx=match(oldres$ID,asian_icogs_info$ID)
oldres$asian_icogs_info=asian_icogs_info$R2[idx]
idx=match(oldres$ID,asian_onco_info$ID)
oldres$asian_onco_info=asian_onco_info$R2[idx]
idx=match(oldres$ID,african_onco_info$ID)
oldres$african_onco_info=african_onco_info$R2[idx]

#check HWE
euro_onco_hwe=as.data.frame(fread("../result/imp_QC/onco/euro_hwe.snplist",header=F))
euro_icogs_hwe=as.data.frame(fread("../result/imp_QC/icogs/euro_hwe.snplist",header=F))
asian_icogs_hwe=as.data.frame(fread("../result/imp_QC/icogs/asian_hwe.snplist",header=F))
asian_onco_hwe=as.data.frame(fread("../result/imp_QC/onco/asian_hwe.snplist",header=F))
african_onco_hwe=as.data.frame(fread("../result/imp_QC/onco/african_hwe.snplist",header=F))
idx=match(oldres$ID,euro_onco_hwe$V1)
oldres$euro_onco_hwe=(!is.na(euro_onco_hwe$V1[idx]))
idx=match(oldres$ID,euro_icogs_hwe$V1)
oldres$euro_icogs_hwe=(!is.na(euro_icogs_hwe$V1[idx]))
idx=match(oldres$ID,asian_icogs_hwe$V1)
oldres$asian_icogs_hwe=(!is.na(asian_icogs_hwe$V1[idx]))
idx=match(oldres$ID,asian_onco_hwe$V1)
oldres$asian_onco_hwe=(!is.na(asian_onco_hwe$V1[idx]))
idx=match(oldres$ID,african_onco_hwe$V1)
oldres$african_onco_hwe=(!is.na(african_onco_hwe$V1[idx]))
colnames(oldres)[14:23]=c("euro_onco", "euro_icogs", "asian_icogs", "asian_onco", "african_onco",
                          "euro_onco_new", "euro_icogs_new", "asian_icogs_new", "asian_onco_new", "african_onco_new")
write.csv(oldres,file="../result/check_ascatpvaluesRES.csv",row.names = F,quote=F)
