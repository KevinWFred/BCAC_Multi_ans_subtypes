#!/usr/bin/env Rscript
.libPaths(c("/data/wangx53",.libPaths()))
setwd("/data/BB_Bioinformatics/Kevin/BCAC/code")

library(data.table)
plink="/usr/local/apps/plink/1.9.0-beta4.4/plink"
plink2="/usr/local/apps/plink/2.3-alpha/plink2"

#for EUR
metafile="../result/PRS1/euro_training_sumstats.txt"
metadat=as.data.frame(fread(metafile))

#CT
runCT=function(sumstatfile="../result/PRS1/euro_training_sumstats.txt",prefix_tun="../result/PRS1/euro_onco_tuning",outprefix="../result/PRS1/euro/euro")
{
  #CT method
  #parameters for clumping
  pthr=1
  r2thr=0.1
  kbpthr=500
  
  cmd=paste0(plink," --bfile ",prefix_tun," --clump ",sumstatfile," --clump-p1 ",
             pthr," --clump-r2 ",r2thr," --clump-kb ",kbpthr," --clump-snp-field rsid --clump-field p --out ",outprefix," --memory 128000")
  system(cmd)
  clumpsnp=read.table(paste0(outprefix,".clumped"),header=T)
  write.table(clumpsnp$SNP,file=paste0(outprefix,".clumpedsnp"),row.names=F,col.names = F,quote=F)
  sumstat=as.data.frame(fread(sumstatfile))
  tmp=data.frame(SNP=sumstat$rsid,A1=sumstat$a1,beta=sumstat$beta)
  write.table(tmp,file=paste0(outprefix,".score"),row.names=F,col.names=T,sep=" ",quote=F)
  tmp=data.frame(SNP=sumstat$rsid,P=sumstat$p)
  write.table(tmp,file=paste0(outprefix,".pvalue"),row.names=F,col.names=T,sep=" ",quote=F)
  cmd=paste0(plink2," --bfile ",prefix_tun," --score ",outprefix,".score ",
             "--q-score-range /data/BB_Bioinformatics/Kevin/PRS_EASLC/result/range_list ",outprefix,".pvalue --extract ",outprefix,".clumpedsnp ",
             "--out ",outprefix,"_tun"," --memory 128000")
  system(cmd)
  prefix_val1="../result/PRS1/african_onco_validation"
  score_val1=paste0(outprefix,"_validation_onco_african")
  cmd=paste0(plink2," --bfile ",prefix_val1," --score ",outprefix,".score ",
             "--q-score-range /data/BB_Bioinformatics/Kevin/PRS_EASLC/result/range_list ",outprefix,".pvalue --extract ",outprefix,".clumpedsnp ",
             "--out ",score_val1," --memory 128000")
  system(cmd)
  prefix_val2="../result/PRS1/african_icogs_validation"
  score_val2=paste0(outprefix,"_validation_icogs_african")
  cmd=paste0(plink2," --bfile ",prefix_val2," --score ",outprefix,".score ",
             "--q-score-range /data/BB_Bioinformatics/Kevin/PRS_EASLC/result/range_list ",outprefix,".pvalue --extract ",outprefix,".clumpedsnp ",
             "--out ",score_val2," --memory 128000")
  system(cmd)
  prefix_val3="../result/PRS1/asian_onco_validation"
  score_val3=paste0(outprefix,"_validation_onco_asian")
  cmd=paste0(plink2," --bfile ",prefix_val3," --score ",outprefix,".score ",
             "--q-score-range /data/BB_Bioinformatics/Kevin/PRS_EASLC/result/range_list ",outprefix,".pvalue --extract ",outprefix,".clumpedsnp ",
             "--out ",score_val3," --memory 128000")
  system(cmd)
  prefix_val4="../result/PRS1/euro_onco_validation"
  score_val4=paste0(outprefix,"_validation_onco_euro")
  cmd=paste0(plink2," --bfile ",prefix_val4," --score ",outprefix,".score ",
             "--q-score-range /data/BB_Bioinformatics/Kevin/PRS_EASLC/result/range_list ",outprefix,".pvalue --extract ",outprefix,".clumpedsnp ",
             "--out ",score_val4," --memory 128000")
  system(cmd)
  prefix_val5="../result/PRS1/hispanic_onco_validation"
  score_val5=paste0(outprefix,"_validation_onco_hispanic")
  cmd=paste0(plink2," --bfile ",prefix_val5," --score ",outprefix,".score ",
             "--q-score-range /data/BB_Bioinformatics/Kevin/PRS_EASLC/result/range_list ",outprefix,".pvalue --extract ",outprefix,".clumpedsnp ",
             "--out ",score_val5," --memory 128000")
  system(cmd)
  
  pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-05,5E-03,5E-02,5E-01,1) 
  famtun=read.table(paste0(prefix_tun,".fam"))
  auc_tun=rep(0,length(pthres))
  for (i in 1:length(pthres))
  {
    prs=read.table(paste0(outprefix,"_tun.p_value_",i,".sscore"))
    if (any(!is.na(prs$V6)))
    {
      #all(prs$V1==famtun$V1)
      pheno.prs=data.frame(y=famtun$V6,prs=prs$V6)
      model1 <- glm(I(y==2)~prs, data=pheno.prs,family = "binomial")
      predicted1 <- predict(model1,pheno.prs, type="response")
      auc_tun[i]=as.numeric(pROC::auc(pheno.prs$y,predicted1,quiet=T))
    }
  }
  
  idx_optimal=which.max(auc_tun) #euro:3,asian 5
  pvalue_opitmal=pthres[idx_optimal]
  print(paste0("pvalue_optimal: ",pvalue_opitmal))
  print(paste0("number of snps: ",sum(clumpsnp$P<pvalue_opitmal)))
  
  #add tuning PRS on the other tuning data(for example tuning eur PRS on tuning asian data),used for weighted PRS
  other_prefix_tun="../result/PRS1/euro_onco_tuning"
  other_score_tun="../result/PRS1/asian/asian_tun_euro_tun" #asian PRS on euro onco
  if (prefix_tun=="../result/PRS1/euro_onco_tuning")
  {
    other_prefix_tun="../result/PRS1/asian_onco_tuning"
    other_score_tun="../result/PRS1/euro/euro_tun_asian_tun"
  }
  cmd=paste0(plink2," --bfile ",other_prefix_tun," --score ",outprefix,".score ",
             "--q-score-range /data/BB_Bioinformatics/Kevin/PRS_EASLC/result/range_list ",outprefix,".pvalue --extract ",outprefix,".clumpedsnp ",
             "--out ",other_score_tun," --memory 128000")
  system(cmd)
  validationprefix=c(prefix_val1,prefix_val2,prefix_val3,prefix_val4,prefix_val5)
  validationscores=c(score_val1,score_val2,score_val3,score_val4,score_val5)
  CTauc_val=data.frame(matrix(NA,nrow=4,ncol=6))
  colnames(CTauc_val)=c("onco_african","icogs_african","onco_asian","onco_euro","onco_hispanic","african")
  rownames(CTauc_val)=c("AUC","totaln","casen","controln")
  afrprs=afrfam=NULL
  for (i in 1:5)
  {
    prs=read.table(paste0(validationscores[i],".p_value_",idx_optimal,".sscore"))
    famval=read.table(paste0(validationprefix[i],".fam"))
    if (i %in% c(1,2))
    {
      afrprs=rbind(afrprs,prs)
      afrfam=rbind(afrfam,famval)
    }
    all(prs$V1==famval$V1)
    pheno.prs=data.frame(y=famval$V6,prs=prs$V6)
    model1 <- glm(I(y==2)~prs, data=pheno.prs,family = "binomial")
    predicted1 <- predict(model1,pheno.prs, type="response")
    CTauc_val[1,i]=as.numeric(pROC::auc(pheno.prs$y,predicted1,quiet=T))
    CTauc_val[2,i]=nrow(famval)
    CTauc_val[3,i]=sum(famval$V6==2)
    CTauc_val[4,i]=sum(famval$V6==1)
  }
  all(afrprs$V1==afrfam$V1)
  pheno.prs=data.frame(y=afrfam$V6,prs=afrprs$V6)
  model1 <- glm(I(y==2)~prs, data=pheno.prs,family = "binomial")
  predicted1 <- predict(model1,pheno.prs, type="response")
  CTauc_val[1,6]=as.numeric(pROC::auc(pheno.prs$y,predicted1,quiet=T))
  return(CTauc_val)
}

#euro
#5e-6, 387 SNPs
#          onco_african icogs_african   onco_asian    onco_euro onco_hispanic african
# AUC         0.5597353     0.5522557    0.5973346 6.274929e-01     0.5676048 0.5574493
# totaln   5569.0000000  1758.0000000 5406.0000000 2.771000e+04  2413.0000000
# casen    3481.0000000   941.0000000 2663.0000000 1.480900e+04  1195.0000000
# controln 2088.0000000   817.0000000 2743.0000000 1.290100e+04  1218.0000000
#asian
#5e-5, 127 SNPs
#          onco_african icogs_african   onco_asian    onco_euro onco_hispanic afircan
# AUC         0.5233717     0.5236688    0.5600331 5.339368e-01     0.5219497 0.524
# totaln   5569.0000000  1758.0000000 5406.0000000 2.771000e+04  2413.0000000
# casen    3481.0000000   941.0000000 2663.0000000 1.480900e+04  1195.0000000
# controln 2088.0000000   817.0000000 2743.0000000 1.290100e+04  1218.0000000

Weighted_CTprs=function()
{
  famtun1=as.data.frame(fread("../result/PRS1/euro_onco_tuning.fam"))
  famtun2=as.data.frame(fread("../result/PRS1/asian_onco_tuning.fam"))
  famtun=rbind(famtun1,famtun2)
  idx_optimal1=3
  idx_optimal2=5
  prs1_tun1=read.table(paste0("../result/PRS1/euro/euro_tun.p_value_",idx_optimal1,".sscore"))
  prs1_tun2=read.table(paste0("../result/PRS1/euro/euro_tun_asian_tun.p_value_",idx_optimal1,".sscore"))
  prs1_tun=rbind(prs1_tun1,prs1_tun2) #euro first,asian second
  prs2_tun1=read.table(paste0("../result/PRS1/asian/asian_tun.p_value_",idx_optimal2,".sscore"))
  prs2_tun2=read.table(paste0("../result/PRS1/asian/asian_tun_euro_tun.p_value_",idx_optimal2,".sscore"))
  prs2_tun=rbind(prs2_tun2,prs2_tun1) #euro first,asian second
  all(prs1_tun$V2==famtun$V2)
  all(prs2_tun$V2==famtun$V2)
  dat_tun=cbind.data.frame(y=famtun$V6-1,prs1=prs1_tun$V6,prs2=prs2_tun$V6)
  model1 <- glm(y~prs1+prs2, data=dat_tun,family = "binomial")
  
  #print(summary(model1))
  w1=summary(model1)$coefficients[2,1]
  w2=summary(model1)$coefficients[3,1]
  
  #apply on validations
  outprefix1="../result/PRS1/euro/euro"
  score_val1_1=paste0(outprefix1,"_validation_onco_african")
  score_val2_1=paste0(outprefix1,"_validation_icogs_african")
  score_val3_1=paste0(outprefix1,"_validation_onco_asian")
  score_val4_1=paste0(outprefix1,"_validation_onco_euro")
  score_val5_1=paste0(outprefix1,"_validation_onco_hispanic")
  outprefix2="../result/PRS1/asian/asian"
  score_val1_2=paste0(outprefix2,"_validation_onco_african")
  score_val2_2=paste0(outprefix2,"_validation_icogs_african")
  score_val3_2=paste0(outprefix2,"_validation_onco_asian")
  score_val4_2=paste0(outprefix2,"_validation_onco_euro")
  score_val5_2=paste0(outprefix2,"_validation_onco_hispanic")
  
  prefix_val1="../result/PRS1/african_onco_validation"
  prefix_val2="../result/PRS1/african_icogs_validation"
  prefix_val3="../result/PRS1/asian_onco_validation"
  prefix_val4="../result/PRS1/euro_onco_validation"
  prefix_val5="../result/PRS1/hispanic_onco_validation"
  
  validationprefix=c(prefix_val1,prefix_val2,prefix_val3,prefix_val4,prefix_val5)
  validationscores1=c(score_val1_1,score_val2_1,score_val3_1,score_val4_1,score_val5_1) #euro prs
  validationscores2=c(score_val1_2,score_val2_2,score_val3_2,score_val4_2,score_val5_2) #asian prs
  Weighted_CTauc_val=data.frame(matrix(NA,nrow=4,ncol=6))
  colnames(Weighted_CTauc_val)=c("onco_african","icogs_african","onco_asian","onco_euro","onco_hispanic","african")
  rownames(Weighted_CTauc_val)=c("AUC","totaln","casen","controln")
  afrprs=afrfam=NULL
  for (i in 1:5)
  {
    prs1=read.table(paste0(validationscores1[i],".p_value_",idx_optimal1,".sscore"))
    prs2=read.table(paste0(validationscores2[i],".p_value_",idx_optimal2,".sscore"))
    prs=prs1$V6*w1+prs2$V6*w2
    famval=read.table(paste0(validationprefix[i],".fam"))
    all(prs1$V1==famval$V1)
    all(prs2$V1==famval$V1)
    if (i %in% c(1,2))
    {
      afrprs=c(afrprs,prs)
      afrfam=rbind(afrfam,famval)
    }
    
    pheno.prs=data.frame(y=famval$V6,prs=prs)
    model1 <- glm(I(y==2)~prs, data=pheno.prs,family = "binomial")
    predicted1 <- predict(model1,pheno.prs, type="response")
    Weighted_CTauc_val[1,i]=as.numeric(pROC::auc(pheno.prs$y,predicted1,quiet=T))
    Weighted_CTauc_val[2,i]=nrow(famval)
    Weighted_CTauc_val[3,i]=sum(famval$V6==2)
    Weighted_CTauc_val[4,i]=sum(famval$V6==1)
  }

  pheno.prs=data.frame(y=afrfam$V6,prs=afrprs)
  model1 <- glm(I(y==2)~prs, data=pheno.prs,family = "binomial")
  predicted1 <- predict(model1,pheno.prs, type="response")
  Weighted_CTauc_val[1,6]=as.numeric(pROC::auc(pheno.prs$y,predicted1,quiet=T))
  
  return(Weighted_CTauc_val)
}
#          onco_african icogs_african   onco_asian    onco_euro onco_hispanic   african
# AUC         0.5590266     0.5518596    0.5967572 6.269937e-01     0.5671613 0.5568423


#LDpred2
library(bigsnpr)
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)

# 1. Read in the phenotype and covariate files
# read in phenotype and covariates
library(data.table)
library(magrittr)
run_ldpred2=function(tmpdir0="euro",sumstatfile="../result/PRS1/euro_training_sumstats.txt",prefix_tun="../result/PRS1/euro_onco_tuning",outprefix="../result/PRS1/euro/euro",opfiltertuning=T)
{
  print(sumstatfile)
  print(prefix_tun)
  #print(prefix_val)
  #2. Obtain HapMap3 SNPs
  #LDpred2 authors recommend restricting the analysis to only the HapMap3 SNPs
  #load HapMap3 SNPs
  # info <- readRDS(runonce::download_file(
  #   "https://ndownloader.figshare.com/files/25503788",
  #   fname = "map_hm3_ldpred2.rds"))
  info =readRDS("/data/BB_Bioinformatics/Kevin/tools/ldpred2/map_hm3_ldpred2.rds")
  #write.table(info$rsid,file="../result/PRS1/ldpred2_hw3snp.txt",row.names = F,col.names = F,quote=F)
  #3. Load and transform the summary statistic file
  #Load summary statistic file
  # Read in the summary statistic file
  sumdat=as.data.frame(fread(sumstatfile)) #8242741,6753434
  print("sumdat dim:")
  print(dim(sumdat))
  # LDpred 2 require the header to follow the exact naming. a1 :effect allelle
  sumdat1=data.frame(chr=sumdat$chr,pos=sumdat$pos,rsid=sumdat$rsid,a0=sumdat$a0,a1=sumdat$a1,n_eff=sumdat$n_eff,beta_se=sumdat$beta_se,p=sumdat$p,beta=sumdat$beta)
  sumdat1=sumdat1[sumdat1$beta!=0,]
  fwrite(sumdat1,file=paste0(sumstatfile,".ldpred"),row.names=F,sep="\t")
  sumstats <- bigreadr::fread2(paste0(sumstatfile,".ldpred")) 
  
  # Filter out hapmap SNPs
  sumstats <- sumstats[sumstats$rsid%in% info$rsid,] #1012780,811530
  write.table(sumstats$rsid,file=paste0(outprefix,"_sumstats_hw3snp.txt"),row.names = F,col.names = F,quote=F)
  print("sumstat in HM3:")
  print(nrow(sumstats))
  #3. Calculate the LD matrix
  # Get maximum amount of cores
  NCORES <- nb_cores()
  # Open a temporary file
  tmpdir <- tempfile(tmpdir = paste0("tmp_data_",tmpdir0))
  #on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)
 
  #to remove mono snps in tuning data (which result NAs in computing corr)
  if (opfiltertuning)
  {
    cmd=paste0(plink," --bfile ",prefix_tun," --maf 0.01 --extract ../result/PRS1/ldpred2_hw3snp.txt"," --make-bed --out ",prefix_tun,"_ldpred_maf01 --memory 128000 --threads 8")
    system(cmd)
    cmd=paste0(plink," --bfile ",prefix_tun,"_ldpred_maf01 --extract ",outprefix,"_sumstats_hw3snp.txt --make-bed --out ",prefix_tun,"_ldpred_maf01 --memory 128000 --threads 8")
    system(cmd)
  }
  
  # preprocess the bed file (only need to do once for each data set)
  if (opfiltertuning)
  {
    file.remove(paste0(prefix_tun,"_ldpred_maf01.bk"))
    #if (!file.exists(paste0(prefix_tun,"_ldpred_maf01.rds")))
    snp_readBed(paste0(prefix_tun,"_ldpred_maf01.bed"))
    # now attach the genotype object
    obj.bigSNP <- snp_attach(paste0(prefix_tun,"_ldpred_maf01.rds"))
  }else
  {
    #if (!file.exists(paste0(prefix_tun,".rds")))
    snp_readBed(paste0(prefix_tun,".bed"))
    obj.bigSNP <- snp_attach(paste0(prefix_tun,".rds"))
  }
  
  # extract the SNP information from the genotype
  map <- obj.bigSNP$map[-3]
  names(map) <- c("chr", "rsid", "pos", "a1", "a0")
  table(map$rsid %in% sumstats$rsid)
  # perform SNP matching
  info_snp <- snp_match(sumstats, map)
  table(info_snp$rsid==map$rsid)
  write.table(info_snp,file=paste0(outprefix,"_LDpred_","info_snp.txt"),row.names=F,sep="\t",quote=F)
  # Assign the genotype to a variable for easier downstream analysis
  genotype <- obj.bigSNP$genotypes
  genotype1 = snp_fastImputeSimple(genotype)
  # Rename the data structures
  CHR <- map$chr
  POS <- map$pos
  # get the CM information from 1000 Genome
  # will download the 1000G file to the current directory (".")
  POS2 <- snp_asGeneticPos(CHR, POS, dir = ".")
  
  # Initialize variables for storing the LD score and LD matrix
  corr <- NULL
  ld <- NULL
  
  # calculate LD
  for (chr in 1:22) {
    print(chr)
    # Extract SNPs that are included in the chromosome
    ind.chr <- which(info_snp$chr == chr)
    ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
    # Calculate the LD
    corr0 <- snp_cor(
      genotype,
      ind.col = ind.chr2,
      ncores = NCORES,
      infos.pos = POS2[ind.chr2],
      size = 3 / 1000
    )
    if (chr == 1) {
      ld <- Matrix::colSums(corr0^2)
      corr <- as_SFBM(corr0, tmpdir)
    } else {
      ld0=Matrix::colSums(corr0^2)
      ld <- c(ld, Matrix::colSums(corr0^2))
      corr$add_columns(corr0, nrow(corr))
    }
    if (sum(is.na(ld))>0) stop(chr)
  }
  save(ld,corr,file=paste0(outprefix,"_LDpred_","ld.RData"))
  #4. Perform LD score regression
  df_beta <- info_snp[,c("beta", "beta_se", "n_eff", "_NUM_ID_")]
  rownames(df_beta)=info_snp$rsid.ss
  ldsc <- snp_ldsc(   ld,
                      length(ld),
                      chi2 = (df_beta$beta / df_beta$beta_se)^2,
                      sample_size = df_beta$n_eff,
                      blocks = NULL)
  # ldsc <- snp_ldsc2(corr,df_beta)
  h2_est <- ldsc[["h2"]] #asn:0.058
  print(paste0("h2_est:",h2_est))
  if (ldsc[['h2']] < 0) print('h2 negative')
  
  
  #6 grid model
  # Prepare data for grid model
  p_seq <- signif(seq_log(1e-5, 1, length.out = 21), 2)
  h2_seq <- round(h2_est * c(0.3,0.7, 1, 1.4), 4)
  grid.param <-
    expand.grid(p = p_seq,
                h2 = h2_seq,
                sparse = c(FALSE, TRUE))
  # Get adjusted beta from grid model
  set.seed(1000) # to get the same result every time
  beta_grid <-
    snp_ldpred2_grid(corr, df_beta, grid.param, ncores = 1)
  
  
  famtun=read.table(paste0(prefix_tun,".fam"))
  pred_grid <- big_prodMat( genotype1, 
                            beta_grid, 
                            ind.col = info_snp$`_NUM_ID_`)
  rownames(pred_grid)=famtun$V2
  
  #to get AUC on tunning
  auc_tun=rep(0,ncol(pred_grid))
  
  for (i in 1:ncol(pred_grid))
  {
    if (any(!is.na(pred_grid[,i])))
    {
      pheno.prs=cbind.data.frame(y=famtun$V6,prs=pred_grid[,i])
      
      model1 <- glm(I(y==2)~prs, data=pheno.prs,family = "binomial")
      predicted1 <- predict(model1,pheno.prs, type="response")
      auc_tun[i]= as.numeric(pROC::auc(pheno.prs$y,predicted1,quiet=T))
    }
  }
  idx_optimal=which.max(auc_tun)
  param_optimal=grid.param[idx_optimal,]
  ldpred_beta=data.frame(snp=info_snp$rsid.ss,effect=info_snp$a1,beta=beta_grid[,idx_optimal])
  betafile=paste0(outprefix,"_LDpred_weights.txt")
  write.table(ldpred_beta,file=betafile,row.names = F,sep="\t",quote=F)
  
  #save PRS on tuning data
  #euro:../result/PRS1/euro/euro_tuning_ldpredprs.txt
  #asian:../result/PRS1/asian/asian_tuning_ldpredprs.txt
  tun.prs=cbind.data.frame(famtun,prs=pred_grid[,idx_optimal])
  write.table(tun.prs,file=paste0(outprefix,"_tuning_ldpredprs.txt"),row.names = F,sep="\t",quote=F)
  save(ld,corr,beta_grid,pred_grid,ldsc,df_beta,grid.param,auc_tun,idx_optimal,param_optimal,file=paste0(outprefix,"LDpred_pred.RData"))
  #add tuning PRS on the other tuning data(for example tuning eur PRS on tuning asian data),used for weighted PRS
  other_prefix_tun="../result/PRS1/euro_onco_tuning" #asian PRS on euro onco
  if (prefix_tun=="../result/PRS1/euro_onco_tuning")
  {
    other_prefix_tun="../result/PRS1/asian_onco_tuning"
  }
  #euro prs on onco asian tuning: ../result/PRS1/euro/euro_asian_onco_tuning_ldpred.score
  #asian prs on onco euro tuning: ../result/PRS1/asian/asian_euro_onco_tuning_ldpred.score
  other_prs_tun=get_PRS_genotype(betafile=betafile,betaprefix=paste0(outprefix,"_ldpredweight"),prefix=other_prefix_tun,outprefix=outprefix)
  
  #validation
  prefix_val1="../result/PRS1/african_onco_validation"
  prefix_val2="../result/PRS1/african_icogs_validation"
  prefix_val3="../result/PRS1/asian_onco_validation"
  prefix_val4="../result/PRS1/euro_onco_validation"
  prefix_val5="../result/PRS1/hispanic_onco_validation"
 
  validationprefix=c(prefix_val1,prefix_val2,prefix_val3,prefix_val4,prefix_val5)
  
  LDpredauc_val=data.frame(matrix(NA,nrow=4,ncol=6))
  colnames(LDpredauc_val)=c("onco_african","icogs_african","onco_asian","onco_euro","onco_hispanic","african")
  rownames(LDpredauc_val)=c("AUC","totaln","casen","controln")
  afrprs=afrfam=NULL
  for (i in 1:5)
  {
    prs=get_PRS_genotype(betafile=betafile,betaprefix=paste0(outprefix,"_ldpredweight"),prefix=validationprefix[i],outprefix=outprefix)

    #prs=read.table(paste0(validationscores[i],".p_value_",idx_optimal,".sscore"))
    famval=read.table(paste0(validationprefix[i],".fam"))
    if (i %in% c(1,2))
    {
      afrprs=rbind(afrprs,prs)
      afrfam=rbind(afrfam,famval)
    }
    all(prs$IID==famval$V1)
    pheno.prs=data.frame(y=famval$V6,prs=prs$SCORESUM)
    model1 <- glm(I(y==2)~prs, data=pheno.prs,family = "binomial")
    predicted1 <- predict(model1,pheno.prs, type="response")
    LDpredauc_val[1,i]=as.numeric(pROC::auc(pheno.prs$y,predicted1,quiet=T))
    LDpredauc_val[2,i]=nrow(famval)
    LDpredauc_val[3,i]=sum(famval$V6==2)
    LDpredauc_val[4,i]=sum(famval$V6==1)
  }
  all(afrprs$V1==afrfam$V1)
  pheno.prs=data.frame(y=afrfam$V6,prs=afrprs$SCORESUM)
  model1 <- glm(I(y==2)~prs, data=pheno.prs,family = "binomial")
  predicted1 <- predict(model1,pheno.prs, type="response")
  LDpredauc_val[1,6]=as.numeric(pROC::auc(pheno.prs$y,predicted1,quiet=T))
  write.table(LDpredauc_val,file=paste0(outprefix,"_ldpred_valauc.txt"),row.names=F,sep="\t",quote=F)
  # 
  
  return(LDpredauc_val)
}

#euro:
#         onco_african icogs_african   onco_asian    onco_euro onco_hispanic   african
# AUC         0.5680855     0.5579379    0.6117766 6.582231e-01     0.6067224 0.5628277
# totaln   5569.0000000  1758.0000000 5406.0000000 2.771000e+04  2413.0000000        NA
# casen    3481.0000000   941.0000000 2663.0000000 1.480900e+04  1195.0000000        NA
# controln 2088.0000000   817.0000000 2743.0000000 1.290100e+04  1218.0000000        NA

#asian
#          onco_african icogs_african   onco_asian   onco_euro onco_hispanic   african
# AUC         0.5505593     0.5346977    0.5874594 5.75309e-01     0.5638608 0.5518476
# totaln   5569.0000000  1758.0000000 5406.0000000 2.77100e+04  2413.0000000        NA
# casen    3481.0000000   941.0000000 2663.0000000 1.48090e+04  1195.0000000        NA
# controln 2088.0000000   817.0000000 2743.0000000 1.29010e+04  1218.0000000        NA
get_PRS_genotype=function(betafile="../result/PRS1/euro/euro_LDpred_weights.txt",betaprefix="../result/PRS1/euro_ldpredweight",prefix="../result/PRS1/euro_onco_tuning",outprefix="../result/PRS1/euro/euro")
{
  beta=as.data.frame(fread(betafile))
  bim=as.data.frame(fread(paste0(prefix,".bim")))
  comsnps=intersect(bim$V2,beta$snp)
  idx1=match(comsnps,bim$V2)
  idx2=match(comsnps,beta$snp)
  bim=bim[idx1,]
  beta=beta[idx2,]
  idx1=which(beta$effect==bim$V5)
  idx2=which(beta$effect==bim$V6)
  a2=bim$V6
  a2[idx2]=bim$V5[idx2]
  tmp=data.frame(snp=beta$snp,a2=a2)
  a2file=paste0(betaprefix,".a2")
  write.table(tmp,file=a2file,col.names = F,row.names = F,sep="\t",quote=F)
  snpfile=paste0(betaprefix,".snp")
  write.table(tmp$snp,file=snpfile,col.names = F,row.names = F,sep="\t",quote=F)
  ldpredprefix=paste0(outprefix,"_",basename(prefix),"_ldpred")
  cmd=paste0(plink2," --bfile ",prefix, " --ref-allele ",a2file," --extract ",snpfile," --make-bed --out ",ldpredprefix)
  system(cmd)
  cmd=paste0(plink," --bfile ",ldpredprefix," --allow-no-sex --score ",betafile," header sum --out ",ldpredprefix," --memory 64000 --threads 8")
  system(cmd)
  PRS=read.table(paste0(ldpredprefix,".profile"),header=T)
  system(paste0("rm ",ldpredprefix,".bed"))
  system(paste0("rm ",ldpredprefix,".bim"))
  system(paste0("rm ",ldpredprefix,".fam"))
  scorefile=paste0(ldpredprefix,".score")
  write.table(PRS,file=scorefile,row.names = F,sep="\t",quote=F)
  return(PRS)
}

Weighted_ldpredprs=function()
{
  famtun1=as.data.frame(fread("../result/PRS1/euro_onco_tuning.fam"))
  famtun2=as.data.frame(fread("../result/PRS1/asian_onco_tuning.fam"))
  famtun=rbind(famtun1,famtun2)
  prs1_tun1=read.table("../result/PRS1/euro/euro_tuning_ldpredprs.txt",header=T)
  prs1_tun2=read.table("../result/PRS1/euro/euro_asian_onco_tuning_ldpred.score",header = T)
  prs1_tun=rbind(data.frame(IID=prs1_tun1$V2,prs=prs1_tun1$prs),data.frame(IID=prs1_tun2$IID,prs=prs1_tun2$SCORESUM)) #euro first,asian second
  prs2_tun1=read.table("../result/PRS1/asian/asian_tuning_ldpredprs.txt",header=T)
  prs2_tun2=read.table("../result/PRS1/asian/asian_euro_onco_tuning_ldpred.score",header = T)
  prs2_tun=rbind(data.frame(IID=prs2_tun2$IID,prs=prs2_tun2$SCORESUM),data.frame(IID=prs2_tun1$V2,prs=prs2_tun1$prs)) #euro first,asian second
  all(prs1_tun$IID==famtun$V2)
  all(prs2_tun$IID==famtun$V2)
  dat_tun=cbind.data.frame(y=famtun$V6-1,prs1=prs1_tun$prs,prs2=prs2_tun$prs)
  model1 <- glm(y~prs1+prs2, data=dat_tun,family = "binomial")
  
  #print(summary(model1))
  w1=summary(model1)$coefficients[2,1]
  w2=summary(model1)$coefficients[3,1]
  
  #apply on validations
  prefix_val1="../result/PRS1/african_onco_validation"
  prefix_val2="../result/PRS1/african_icogs_validation"
  prefix_val3="../result/PRS1/asian_onco_validation"
  prefix_val4="../result/PRS1/euro_onco_validation"
  prefix_val5="../result/PRS1/hispanic_onco_validation"
  
  validationprefix=c(prefix_val1,prefix_val2,prefix_val3,prefix_val4,prefix_val5)
  #euro prs
  validationprefix1=paste0("../result/PRS1/euro/euro_",basename(validationprefix),"_ldpred")
  #asian prs
  validationprefix2=paste0("../result/PRS1/asian/asian_",basename(validationprefix),"_ldpred")
  
  Weighted_LDpredauc_val=data.frame(matrix(NA,nrow=4,ncol=6))
  colnames(Weighted_LDpredauc_val)=c("onco_african","icogs_african","onco_asian","onco_euro","onco_hispanic","african")
  rownames(Weighted_LDpredauc_val)=c("AUC","totaln","casen","controln")
  afrprs=afrfam=NULL
  for (i in 1:5)
  {
    prs1=read.table(paste0(validationprefix1[i],".score"),header=T)
    prs2=read.table(paste0(validationprefix2[i],".score"),header=T)
    prs=prs1$SCORESUM*w1+prs2$SCORESUM*w2
    famval=read.table(paste0(validationprefix[i],".fam"))
    all(prs1$V1==famval$V1)
    all(prs2$V1==famval$V1)
    if (i %in% c(1,2))
    {
      afrprs=c(afrprs,prs)
      afrfam=rbind(afrfam,famval)
    }
    
    pheno.prs=data.frame(y=famval$V6,prs=prs)
    model1 <- glm(I(y==2)~prs, data=pheno.prs,family = "binomial")
    predicted1 <- predict(model1,pheno.prs, type="response")
    Weighted_LDpredauc_val[1,i]=as.numeric(pROC::auc(pheno.prs$y,predicted1,quiet=T))
    Weighted_LDpredauc_val[2,i]=nrow(famval)
    Weighted_LDpredauc_val[3,i]=sum(famval$V6==2)
    Weighted_LDpredauc_val[4,i]=sum(famval$V6==1)
  }
  
  pheno.prs=data.frame(y=afrfam$V6,prs=afrprs)
  model1 <- glm(I(y==2)~prs, data=pheno.prs,family = "binomial")
  predicted1 <- predict(model1,pheno.prs, type="response")
  Weighted_LDpredauc_val[1,6]=as.numeric(pROC::auc(pheno.prs$y,predicted1,quiet=T))
  
  return(Weighted_LDpredauc_val)
}
#          onco_african icogs_african   onco_asian    onco_euro onco_hispanic   african
# AUC         0.5755086     0.5609777    0.6187005 6.530507e-01     0.6060151 0.5717066
