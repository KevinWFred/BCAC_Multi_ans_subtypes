#!/usr/bin/env Rscript
.libPaths(c("/data/wangx53",.libPaths()))
library(data.table)
plink="/usr/local/apps/plink/1.9.0-beta4.4/plink"
plink2="/usr/local/apps/plink/2.3-alpha/plink2"

#for EUR
metafile="../result/meta_euro_lr.meta"
metadat=as.data.frame(fread(metafile))

#CT
runCT=function(sumstatfile="../result/six10k_sumstats.txt",prefix_tun="../result/six10k_train",prefix_val="../result/six10k_test",outprefix="six10k")
{
  #CT method
  #parameters for clumping
  pthr=1
  r2thr=0.1
  kbpthr=500
  
  cmd=paste0(plink," --bfile ",prefix_tun," --clump ",sumstatfile," --clump-p1 ",
             pthr," --clump-r2 ",r2thr," --clump-kb ",kbpthr," --clump-snp-field rsid --clump-field p --out ../result/",outprefix)
  system(cmd)
  tmp=read.table(paste0("../result/",outprefix,".clumped"),header=T)
  write.table(tmp$SNP,file=paste0("../result/",outprefix,".clumpedsnp"),row.names=F,col.names = F,quote=F)
  sumstat=as.data.frame(fread(sumstatfile))
  tmp=data.frame(SNP=sumstat$rsid,A1=sumstat$a1,beta=sumstat$beta)
  write.table(tmp,file=paste0("../result/",outprefix,".score"),row.names=F,col.names=T,sep=" ",quote=F)
  tmp=data.frame(SNP=sumstat$rsid,P=sumstat$p)
  write.table(tmp,file=paste0("../result/",outprefix,".pvalue"),row.names=F,col.names=T,sep=" ",quote=F)
  cmd=paste0(plink2," --bfile ",prefix_tun," --score ../result/",outprefix,".score ",
             "--q-score-range /data/BB_Bioinformatics/Kevin/PRS_EASLC/result/range_list ../result/",outprefix,".pvalue --extract ../result/",outprefix,".clumpedsnp ",
             "--out ../result/",outprefix,"_tun")
  system(cmd)
  cmd=paste0(plink2," --bfile ",prefix_val," --score ../result/",outprefix,".score ",
             "--q-score-range /data/BB_Bioinformatics/Kevin/PRS_EASLC/result/range_list ../result/",outprefix,".pvalue --extract ../result/",outprefix,".clumpedsnp ",
             "--out ../result/",outprefix,"_val")
  system(cmd)
  pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-05,5E-03,5E-02,5E-01,1) 
  famtun=read.table(paste0(prefix_tun,".fam"))
  famval=read.table(paste0(prefix_val,".fam"))
  auc_tun=rep(0,length(pthres))
  for (i in 1:length(pthres))
  {
    prs=read.table(paste0("../result/",outprefix,"_tun.p_value_",i,".sscore"))
    if (any(!is.na(prs$V6)))
    {
      #all(prs$V1==famtun$V1)
      pheno.prs=data.frame(y=famtun$V6,prs=prs$V6)
      model1 <- glm(I(y==2)~prs, data=pheno.prs,family = "binomial")
      predicted1 <- predict(model1,pheno.prs, type="response")
      auc_tun[i]=as.numeric(pROC::auc(pheno.prs$y,predicted1,quiet=T))
    }
  }
  
  idx_optimal=which.max(auc_tun)
  prs=read.table(paste0("../result/",outprefix,"_val.p_value_",idx_optimal,".sscore"))
  pheno.prs=data.frame(y=famval$V6,prs=prs$V6)
  model1 <- glm(I(y==2)~prs, data=pheno.prs,family = "binomial")
  predicted1 <- predict(model1,pheno.prs, type="response")
  CTauc_val=as.numeric(pROC::auc(pheno.prs$y,predicted1,quiet=T)) #0.595
  return(CTauc_val)
}
