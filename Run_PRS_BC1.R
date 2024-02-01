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
runCT=function(sumstatfile="../result/PRS1/euro_training_sumstats.txt",prefix_tun="../result/PRS1/euro_onco_tuning",outprefix="PRS1/euro/euro")
{
  #CT method
  #parameters for clumping
  pthr=1
  r2thr=0.1
  kbpthr=500
  
  cmd=paste0(plink," --bfile ",prefix_tun," --clump ",sumstatfile," --clump-p1 ",
             pthr," --clump-r2 ",r2thr," --clump-kb ",kbpthr," --clump-snp-field rsid --clump-field p --out ../result/",outprefix," --memory 128000")
  system(cmd)
  clumpsnp=read.table(paste0("../result/",outprefix,".clumped"),header=T)
  write.table(clumpsnp$SNP,file=paste0("../result/",outprefix,".clumpedsnp"),row.names=F,col.names = F,quote=F)
  sumstat=as.data.frame(fread(sumstatfile))
  tmp=data.frame(SNP=sumstat$rsid,A1=sumstat$a1,beta=sumstat$beta)
  write.table(tmp,file=paste0("../result/",outprefix,".score"),row.names=F,col.names=T,sep=" ",quote=F)
  tmp=data.frame(SNP=sumstat$rsid,P=sumstat$p)
  write.table(tmp,file=paste0("../result/",outprefix,".pvalue"),row.names=F,col.names=T,sep=" ",quote=F)
  cmd=paste0(plink2," --bfile ",prefix_tun," --score ../result/",outprefix,".score ",
             "--q-score-range /data/BB_Bioinformatics/Kevin/PRS_EASLC/result/range_list ../result/",outprefix,".pvalue --extract ../result/",outprefix,".clumpedsnp ",
             "--out ../result/",outprefix,"_tun"," --memory 128000")
  system(cmd)
  prefix_val1="../result/PRS1/african_onco_validation"
  score_val1=paste0("../result/",outprefix,"_validation_onco_african")
  cmd=paste0(plink2," --bfile ",prefix_val1," --score ../result/",outprefix,".score ",
             "--q-score-range /data/BB_Bioinformatics/Kevin/PRS_EASLC/result/range_list ../result/",outprefix,".pvalue --extract ../result/",outprefix,".clumpedsnp ",
             "--out ",score_val1," --memory 128000")
  system(cmd)
  prefix_val2="../result/PRS1/african_icogs_validation"
  score_val2=paste0("../result/",outprefix,"_validation_icogs_african")
  cmd=paste0(plink2," --bfile ",prefix_val2," --score ../result/",outprefix,".score ",
             "--q-score-range /data/BB_Bioinformatics/Kevin/PRS_EASLC/result/range_list ../result/",outprefix,".pvalue --extract ../result/",outprefix,".clumpedsnp ",
             "--out ",score_val2," --memory 128000")
  system(cmd)
  prefix_val3="../result/PRS1/asian_onco_validation"
  score_val3=paste0("../result/",outprefix,"_validation_onco_asian")
  cmd=paste0(plink2," --bfile ",prefix_val3," --score ../result/",outprefix,".score ",
             "--q-score-range /data/BB_Bioinformatics/Kevin/PRS_EASLC/result/range_list ../result/",outprefix,".pvalue --extract ../result/",outprefix,".clumpedsnp ",
             "--out ",score_val3," --memory 128000")
  system(cmd)
  prefix_val4="../result/PRS1/euro_onco_validation"
  score_val4=paste0("../result/",outprefix,"_validation_onco_euro")
  cmd=paste0(plink2," --bfile ",prefix_val4," --score ../result/",outprefix,".score ",
             "--q-score-range /data/BB_Bioinformatics/Kevin/PRS_EASLC/result/range_list ../result/",outprefix,".pvalue --extract ../result/",outprefix,".clumpedsnp ",
             "--out ",score_val4," --memory 128000")
  system(cmd)
  prefix_val5="../result/PRS1/hispanic_onco_validation"
  score_val5=paste0("../result/",outprefix,"_validation_onco_hispanic")
  cmd=paste0(plink2," --bfile ",prefix_val5," --score ../result/",outprefix,".score ",
             "--q-score-range /data/BB_Bioinformatics/Kevin/PRS_EASLC/result/range_list ../result/",outprefix,".pvalue --extract ../result/",outprefix,".clumpedsnp ",
             "--out ",score_val5," --memory 128000")
  system(cmd)
  
  pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-05,5E-03,5E-02,5E-01,1) 
  famtun=read.table(paste0(prefix_tun,".fam"))
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
  pvalue_opitmal=pthres[idx_optimal]
  print(paste0("pvalue_optimal: ",pvalue_opitmal))
  print(paste0("number of snps: ",sum(clumpsnp$P<pvalue_opitmal)))
  validationprefix=c(prefix_val1,prefix_val2,prefix_val3,prefix_val4,prefix_val5)
  validationscores=c(score_val1,score_val2,score_val3,score_val4,score_val5)
  CTauc_val=data.frame(matrix(NA,nrow=4,ncol=5))
  colnames(CTauc_val)=c("onco_african","icogs_african","onco_asian","onco_euro","onco_hispanic")
  rownames(CTauc_val)=c("AUC","totaln","casen","controln")
  for (i in 1:5)
  {
    prs=read.table(paste0(validationscores[i],".p_value_",idx_optimal,".sscore"))
    famval=read.table(paste0(validationprefix[i],".fam"))
    all(prs$V1==famval$V1)
    pheno.prs=data.frame(y=famval$V6,prs=prs$V6)
    model1 <- glm(I(y==2)~prs, data=pheno.prs,family = "binomial")
    predicted1 <- predict(model1,pheno.prs, type="response")
    CTauc_val[1,i]=as.numeric(pROC::auc(pheno.prs$y,predicted1,quiet=T))
    CTauc_val[2,i]=nrow(famval)
    CTauc_val[3,i]=sum(famval$V6==2)
    CTauc_val[4,i]=sum(famval$V6==1)
  }
  
  return(CTauc_val)
}

#euro
#5e-6, 387 SNPs
#          onco_african icogs_african   onco_asian    onco_euro onco_hispanic
# AUC         0.5597353     0.5522557    0.5973346 6.274929e-01     0.5676048
# totaln   5569.0000000  1758.0000000 5406.0000000 2.771000e+04  2413.0000000
# casen    3481.0000000   941.0000000 2663.0000000 1.480900e+04  1195.0000000
# controln 2088.0000000   817.0000000 2743.0000000 1.290100e+04  1218.0000000
#asian
#5e-5, 127 SNPs
#          onco_african icogs_african   onco_asian    onco_euro onco_hispanic
# AUC         0.5233717     0.5236688    0.5600331 5.339368e-01     0.5219497
# totaln   5569.0000000  1758.0000000 5406.0000000 2.771000e+04  2413.0000000
# casen    3481.0000000   941.0000000 2663.0000000 1.480900e+04  1195.0000000
# controln 2088.0000000   817.0000000 2743.0000000 1.290100e+04  1218.0000000

#LDpred2
library(bigsnpr)
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)

# 1. Read in the phenotype and covariate files
# read in phenotype and covariates
library(data.table)
library(magrittr)
run_ldpred2=function(tmpdir=3,sumstatfile="../result/PRS1/euro_training_sumstats.txt",prefix_tun="../result/PRS1/euro_onco_tuning",prefix_val="../result/six10k_test",outprefix="PRS1/euro/euro",opfiltertuning=F)
{
  print(sumstatfile)
  print(prefix_tun)
  print(prefix_val)
  #2. Obtain HapMap3 SNPs
  #LDpred2 authors recommend restricting the analysis to only the HapMap3 SNPs
  #load HapMap3 SNPs
  # info <- readRDS(runonce::download_file(
  #   "https://ndownloader.figshare.com/files/25503788",
  #   fname = "map_hm3_ldpred2.rds"))
  info =readRDS("/data/BB_Bioinformatics/Kevin/tools/ldpred2/map_hm3_ldpred2.rds")
  #3. Load and transform the summary statistic file
  #Load summary statistic file
  # Read in the summary statistic file
  sumdat=as.data.frame(fread(sumstatfile)) #5678649
  print("sumdat dim:")
  print(dim(sumdat))
  # LDpred 2 require the header to follow the exact naming. a1 :effect allelle
  sumdat1=data.frame(chr=sumdat$chr,pos=sumdat$pos,rsid=sumdat$rsid,a0=sumdat$a0,a1=sumdat$a1,n_eff=sumdat$n_eff,beta_se=sumdat$beta_se,p=sumdat$p,beta=sumdat$beta)
  fwrite(sumdat1,file=paste0(sumstatfile,".ldpred"),row.names=F,sep="\t")
  sumstats <- bigreadr::fread2(paste0(sumstatfile,".ldpred")) 
  
  # Filter out hapmap SNPs
  sumstats <- sumstats[sumstats$rsid%in% info$rsid,] #689779
  print("sumstat in HM3:")
  print(nrow(sumstats))
  #3. Calculate the LD matrix
  # Get maximum amount of cores
  NCORES <- nb_cores()
  # Open a temporary file
  tmp <- tempfile(tmpdir = paste0("tmp-data",tmpdir))
  #on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)
  # Initialize variables for storing the LD score and LD matrix
  corr <- NULL
  ld <- NULL
  
  #to remove mono snps in tuning data (which result NAs in computing corr)
  if (opfiltertuning)
  {
    cmd=paste0(plink," --bfile ",prefix_tun," --maf 0.01 --make-bed --out ",prefix_tun,"_maf05 --memory 128000 --threads 8")
    system(cmd)
  }
  
  # preprocess the bed file (only need to do once for each data set)
  if (opfiltertuning)
  {
    if (!file.exists(paste0(prefix_tun,"_maf05.rds")))
      snp_readBed(paste0(prefix_tun,"_maf05.bed"))
  }else
  {
    if (!file.exists(paste0(prefix_tun,".rds")))
      snp_readBed(paste0(prefix_tun,".bed"))
  }
  
  # now attach the genotype object
  
  obj.bigSNP <- snp_attach(paste0(prefix_tun,".rds"))
  # extract the SNP information from the genotype
  map <- obj.bigSNP$map[-3]
  names(map) <- c("chr", "rsid", "pos", "a1", "a0")
  # perform SNP matching
  info_snp <- snp_match(sumstats, map)
  write.table(info_snp,file=paste0("../result/",outprefix,"_LDpred_","info_snp.txt"),row.names=F,sep="\t",quote=F)
  # Assign the genotype to a variable for easier downstream analysis
  genotype <- obj.bigSNP$genotypes
  genotype1 = snp_fastImputeSimple(genotype)
  # Rename the data structures
  CHR <- map$chr
  POS <- map$pos
  # get the CM information from 1000 Genome
  # will download the 1000G file to the current directory (".")
  POS2 <- snp_asGeneticPos(CHR, POS, dir = ".")
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
      corr <- as_SFBM(corr0, tmp)
    } else {
      ld0=Matrix::colSums(corr0^2)
      ld <- c(ld, Matrix::colSums(corr0^2))
      corr$add_columns(corr0, nrow(corr))
    }
    if (sum(is.na(ld))>0) stop(chr)
  }
  save(ld,corr,file=paste0("../result/",outprefix,"_LDpred_","ld.RData"))
  #4. Perform LD score regression
  df_beta <- info_snp[,c("beta", "beta_se", "n_eff", "_NUM_ID_")]
  rownames(df_beta)=info_snp$rsid.ss
  ldsc <- snp_ldsc(   ld,
                      length(ld),
                      chi2 = (df_beta$beta / df_beta$beta_se)^2,
                      sample_size = df_beta$n_eff,
                      blocks = NULL)
  # ldsc <- snp_ldsc2(corr,df_beta)
  h2_est <- ldsc[["h2"]] #0.1
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
  
  #validation
  if(! file.exists(paste0(prefix_val,".rds")))
    snp_readBed(paste0(prefix_val,".bed"))
  val.obj.bigSNP <- snp_attach(paste0(prefix_val,".rds"))
  genotype_val <- val.obj.bigSNP$genotypes
  genotype1_val = snp_fastImputeSimple(genotype_val)
  famval=read.table(paste0(prefix_val,".fam"))
  pred_grid_val <- big_prodMat( genotype1_val, 
                                beta_grid, 
                                ind.col = info_snp$`_NUM_ID_`)
  rownames(pred_grid_val)=famval$V2
  
  save(ld,corr,beta_grid,pred_grid_val,beta_grid,pred_grid,ldsc,df_beta,grid.param,file=paste0("../result/LDpred_",outprefix,"_pred.RData"))
  
  #to get AUC
  famtun=read.table(paste0(prefix_tun,".fam"))
  famval=read.table(paste0(prefix_val,".fam"))
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
  # 
  #auc on validation set
  
  pheno.prs=cbind.data.frame(y=famval$V6,prs=pred_grid_val[,idx_optimal])
  model1 <- glm(I(y==2)~prs, data=pheno.prs,family = "binomial")
  predicted1 <- predict(model1,pheno.prs, type="response")
  LDpredauc_val=as.numeric(pROC::auc(pheno.prs$y,predicted1,quiet=T))#0.593
  
  return(LDpredauc_val)
}

