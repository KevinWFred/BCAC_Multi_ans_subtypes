#!/usr/bin/env Rscript
#compare p-values in previous results
#new results
load("../result/compute_metapvalues_population.RData")
length(collect_pvalues_euro$acatpvalues) #10142249
eurexistingtable=as.data.frame(fread("/data/BB_Bioinformatics/RQ/Dataset/bc_summary_gwas.txt"))
dim(eurexistingtable)
# [1] 10760767       47
tmp1=paste0(eurexistingtable$chr.iCOGs,":",eurexistingtable$Position.iCOGs)
library(rtracklayer)
chain=import.chain("../../tools/liftover/hg19ToHg38.over.chain")

knownvar$hg38position=tmp$start
compareresults=function(oldres=eurexistingtable,newres=collect_pvalues_euro$acatpvalues)
{
  tmp=oldres[,1]
  tmp=unlist(strsplit(tmp,"_"))
  oldchr=tmp[seq(1,length(tmp),4)]
  oldchr[which(oldchr=="23")]="X"
  oldpos=tmp[seq(2,length(tmp),4)]
  gr_olddat=GRanges(seqnames = paste0("chr",oldchr),ranges=IRanges(start=as.integer(oldpos),width=1))
  tmp=liftOver(gr_olddat,chain)
  tmp=as.data.frame(tmp)
  idx=1:nrow(oldres)
  idx=idx[!idx %in% tmp$group]
  #can't find hg38 position
  print(length(idx))
  oldres=oldres[tmp$group,]
  oldres$hg38pos=tmp$start
  newres=data.frame(ID=names(newres),newP=newres)
  tmp=unlist(strsplit(newres$ID,":"))
  newres$chr=tmp[seq(1,length(tmp),4)]
  newres$chr=gsub("chr","",newres$chr)
  newres$hg38pos=as.integer(tmp[seq(2,length(tmp),4)])
  newres$ref=tmp[seq(3,length(tmp),4)]
  newres$alt=tmp[seq(4,length(tmp),4)]
  tmp1=paste0(oldres$chr.iCOGs,":",oldres$hg38pos,":",oldres$Baseline.Gwas,":",oldres$Effect.Gwas)
  tmp2=paste0(oldres$chr.iCOGs,":",oldres$hg38pos,":",oldres$Effect.Gwas,":",oldres$Baseline.Gwas)
  tmp3=paste0(newres$chr,":",newres$hg38pos,":",newres$ref,":",newres$alt)
  tmp4=intersect(tmp3,tmp1)
  tmp5=intersect(tmp4,tmp1)
  newres$oldP=NA
  idx1=match(tmp4,tmp1)
  idx2=match(tmp4,tmp3)
  newres$oldP[idx2]=oldres$p.meta[idx1]
  idx1=match(tmp5,tmp1)
  idx2=match(tmp5,tmp3)
  newres$oldP[idx2]=oldres$p.meta[idx1]
  cor(-log10(newres$oldP),-log10(newres$newP),use="complete")
  plot(-log10(newres$oldP),-log10(newres$newP))
  
  
}
gr_tmp=GRanges(seqnames = "chr1",ranges=IRanges(start=1581586,width=10))
tmp1=liftOver(gr_tmp,chain)
