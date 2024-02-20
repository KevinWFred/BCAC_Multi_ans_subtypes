#!/usr/bin/env Rscript
#step4, find the best super learner based on tuning data and compute final PRS on validation data
.libPaths(c("/data/BB_Bioinformatics/Kevin/tools/Rpackages",.libPaths()))

args = commandArgs(trailingOnly=TRUE) #config input file

library(data.table)
library(caret)
library(CTSLEB)
library(dplyr)
library(SuperLearner)
library(ranger)
library(glmnet)

setwd("/data/BB_Bioinformatics/Kevin/BCAC/code")
set.seed(1000)

helper_PreparePlinkFile <- function(plink_list,
                                    results_dir)
  
{
  
  temp.dir <- paste0(results_dir,"temp/")
  scores <- plink_list[[1]]
  score_file <- paste0(temp.dir,"score_file")
  # write.table(scores,
  #             file = score_file ,
  #             row.names = F,
  #             col.names = F,
  #             quote=F)
  
  p_values <- plink_list[[2]]
  p_value_file <- paste0(temp.dir,"p_value_file")
  
  unique_infor <- plink_list[[3]]
  
  q_range <- plink_list[[4]]
  q_range_file <- paste0(temp.dir,"q_range_file")
  
  # write.table(q_range,
  #             file = q_range_file,
  #             row.names = F,
  #             col.names = F,
  #             quote=F)
  assign("q_range", q_range, envir = .GlobalEnv)
  assign("scores", scores, envir = .GlobalEnv)
  assign("p_values", p_values, envir = .GlobalEnv)
  assign("unique_infor", unique_infor, envir = .GlobalEnv)
  assign("score_file", score_file, envir = .GlobalEnv)
  assign("p_value_file", p_value_file, envir = .GlobalEnv)
  assign("q_range_file", q_range_file, envir = .GlobalEnv)
  
  names <-c("score_file","p_value_file","q_range_file")
  values <- list(score_file, p_value_file, q_range_file)
  list <- setNames(values, names)
  
  return(list)
}
PreparePlinkFileEBayes <- function(
                                   snp_list,
                                   clump_info,
                                   post_clump_info,
                                   post_beta,
                                   results_dir){
  
  print("Executing PreparePlinkFileEBayes() ... ")
  unique_id <- post_clump_info$SNP
  names(unique_id) <- "SNP"
  
  #number of ancestry
  n_ans <- ncol(post_beta)
  
  n_col = length(snp_list)
  n_row <- nrow(post_clump_info)
  post_beta <- as.matrix(post_beta)
  post_beta[is.na(post_beta)] <- 0
  
  beta_mat <- matrix(rep(post_beta,n_col),nrow =n_row,ncol =n_col*n_ans)
  names <- rep("c",n_col*n_ans)
  temp <- 0
  for(ldx in 1:n_col){
    ld <- snp_list[[ldx]]
    names(ld) = "SNP"
    idx <- which(clump_info$SNP%in%ld$SNP==F)
    beta_mat[idx,(1:n_ans)+temp] = 0
    names[(1:n_ans)+temp] <- paste0(names(snp_list[[ldx]]),
                                    "_",
                                    colnames(post_beta))
    temp <- temp + n_ans
  }
  colnames(beta_mat) <- names
  this_scores <- data.frame(SNP = unique_id,
                            A1 = post_clump_info$A1,
                            beta_mat)
  temp_dir <- paste0(results_dir,"temp/")
  score_eb_file <- paste0(temp_dir,"score_eb_file")
  this_p_values <- data.frame(SNP = unique_id,
                              P = post_clump_info$P)
  p_values_eb_file <-paste0(temp_dir,"p_values_eb_file")
  names <- c("scores_eb",
             "p_values_eb",
             "score_eb_file",
             "p_values_eb_file")
  values <- list(this_scores,
                 this_p_values,
                 score_eb_file,
                 p_values_eb_file)
  
  print("PreparePlinkFileEBayes() complete ...")
  result <- setNames(values, names)
  
  return(result)
}


helper_prefix <- function(out_prefix,
                          ebayes=FALSE){
  if (is.null(out_prefix)) {
    prefix <- paste0("")
    eb_prefix <- paste0("eb_")
    
    if (ebayes) {
      return(eb_prefix)
    }else{
      return(prefix)
    }
    
  } else {
    prefix <- paste0(out_prefix, "_")
    eb_prefix <- paste0(out_prefix, "_eb_")
    
    if (ebayes) {
      return(eb_prefix)
    }else{
      return(prefix)
    }
  }
}

#get prs_mat based on results
helper_CombinePRS <- function(scores,
                              pthres,
                              prs_p_other_)
{
  
  print("combining PRS results ...")
  prs_list <- list()
  temp <- 1
  names <- colnames(scores[3:ncol(scores)])
  for(k1 in 1:length(pthres)){
    for(k2 in 1:length(pthres)){
      prs_file <- paste0(prs_p_other_,k1,".p_tar_",k2,".sscore")
      if (file.exists(prs_file))
      {
        #print(paste0("reading :", prs_file))
        prs_temp <- fread(prs_file)
        prs_list[[temp]] <- prs_temp[,6:ncol(prs_temp)]
        colnames(prs_list[[temp]]) <- paste0(names,"_",
                                             "p_other_",
                                             pthres[k1],
                                             "_p_tar_",
                                             pthres[k2])
        temp <- temp + 1
      }
    }
  }
  prs_mat <- as.data.frame(cbind(prs_temp[,1:2],bind_cols(prs_list)))
  print("prs_mat complete... ")
  return(prs_mat)
}

configfile=args[1] #CT_SLEB_never.config
#configfile="CT_SLEB_EAS.config"
config=read.table(configfile)
if (!file.exists(configfile)) stop("no input file.")
#if (nrow(config)!=9 || ncol(config)!=2) stop("input file format is not right")

#Specify the directory for the summary statistics
EUR_sumstats_file <- config[1,2] #"../result/EUnever_CTSLEB.sumdat" #  reference population
Target_sumstats_file <- config[2,2]  #"../result/EASnever_CTSLEB.sumdat" #  target population

#Specify the directory to the reference data for clumping purpose
EUR_ref_plinkfile <- config[3,2] #"../../tools/CTSLEB/EUR/chr_all_FLCCA"
Target_ref_plinkfile <- config[4,2] #"../result/FLCCA_1kg"#"../../tools/CTSLEB/EAS/chr_all"
#Specify the tuning and validation data directory
Target_test_plinkfile <- config[5,2] #"../result/FLCCA_1kg"
#output folder
out_dir = config[6,2] # "../result/CTSLEB_Never_swarm/"
outprefix <- config[7,2] #"EASnever_CTSLEB"
plink19_exec=config[8,2] #"/usr/local/apps/plink/1.9/plink"
plink2_exec=config[9,2] #"/usr/local/apps/plink/2.3-alpha/plink2"
#famtrainfile=config[10,2] # /data/BB_Bioinformatics/Kevin/PRS_EASLC/result/EAS_never_train1.fam
#famtestfile=config[11,2] # /data/BB_Bioinformatics/Kevin/PRS_EASLC/result/EAS_never_test1.fam

#used for split tuning/validation data
famtrain=read.table(paste0(Target_ref_plinkfile,".fam"))
famtest=read.table(paste0(Target_test_plinkfile,".fam"))

sum_EUR <- fread(EUR_sumstats_file,header=T)
sum_Target <- fread(Target_sumstats_file,header=T)
sum_com <- AlignSum(sum_target = sum_Target,
                    sum_ref = sum_EUR)
# #Specify the directory to the reference data for clumping purpose
# EUR_ref_plinkfile <- "../../tools/CTSLEB/EUR/chr_all_FLCCA"
# Target_ref_plinkfile <- "../result/FLCCA_1kg_train"#"../../tools/CTSLEB/EAS/chr_all"
# #Specify the tuning and validation data directory
# Target_test_plinkfile <- "../result/FLCCA_1kg_train"

# PRS_farm <- SetParamsFarm(plink19_exec = plink19_exec,
#                           plink2_exec = plink2_exec,threads = 2,mem=64000,
#                           wc_base_vec = 100, r2_vec = c(0.01,0.1,0.5,0.8))
PRS_farm <- SetParamsFarm(plink19_exec = plink19_exec,
                          plink2_exec = plink2_exec,threads = 6,mem=64000)
out.prefix <- paste0(outprefix, "_")
wc_base_vec = PRS_farm$wc_base_vec
r2_vec = PRS_farm$r2_vec
temp=1
snp_list <-list()
temp.dir <- paste0(out_dir,"temp/")
for(r_ind in 1:length(r2_vec)){
  wc_vec <- wc_base_vec/r2_vec[r_ind]
  for(w_ind in 1:length(wc_vec)){
    pthr <-1
    r2thr <- r2_vec[r_ind]
    kbpthr <- wc_vec[w_ind]
    ref_outfile <- paste0(temp.dir, out.prefix, "ref_CT_rind_",r_ind,"_wcind_",w_ind)
    target_outfile <- paste0(temp.dir, out.prefix,"target_CT_rind_",r_ind,"_wcind_",w_ind)
    LD_ref <-fread(paste0(ref_outfile,
                          ".clumped"))[,3,drop=F]
    LD_target <-fread(paste0(target_outfile,
                             ".clumped"))[,3,drop=F]
    # print(paste0("binding rows for ",
    #              ref_outfile,
    #              ".clumped and ",
    #              target_outfile,
    #              ".clumped"))
    LD  <- rbind(LD_ref,LD_target)
    snp_list[[temp]] <- LD
    # print(paste0("creating snp list for clump_r2_",
    #              r2thr,
    #              "_ws_",
    #              kbpthr))
    names(snp_list[[temp]]) <- paste0("clump_r2_",
                                      r2thr,
                                      "_ws_",
                                      kbpthr)
    temp <- temp + 1
  }
}
plink_list <- PreparePlinkFile(params_farm = PRS_farm,
                               snp_list = snp_list,
                               sum_com = sum_com,
                               results_dir = out_dir)
file_list <- helper_PreparePlinkFile(plink_list = plink_list,
                                     results_dir = out_dir)
plink_list <- c(plink_list,file_list)
assign("plink_list", plink_list, envir = .GlobalEnv)

scores <- plink_list[[1]]
#prs_p_other_="../result/CTSLEB_Never/temp/EASnever_CTSLEB_train_prs_p_other_"
prs_p_other_ <- paste0(temp.dir, out.prefix, "prs_p_other_")
prs_mat <- helper_CombinePRS(scores = scores,
                             pthres = pthres,
                             prs_p_other_ = prs_p_other_)
prs_tune=prs_mat[prs_mat[,1] %in% famtrain[,1],]
#prs_tune <- prs_mat[1:ceiling(nrow(prs_mat)/2),]
#n.total.prs <- length(pthres)^2*length(r2_vec)*length(wc_base_vec)
n.total.prs=ncol(prs_mat)-2
prs_auc_vec_test <- rep(0,n.total.prs)

library(pROC)
# phenotype=read_excel("../data/FLCCA_Phenotype.xlsx")
# phenotype=phenotype[phenotype$HISTOLOGY %in% c("ADENO","ADENO_OTHER","CONTROL"),]
# phenotype$y=0
# phenotype$y[which(phenotype$CASECONTROL=="CASE")]=1
# phenotype1=phenotype[match(prs_tune[,1],phenotype$GWAS_ID),]
phenoicogs=read.table("../data/concept_750_zhang_icogs_pheno_v15_02_age.txt",header=T,sep="\t")
phenoicogs$ID=phenoicogs$SG_ID
idx=which(colnames(phenoicogs) %in% paste0("pc",1:10))
phenoicogs[,idx]=NULL
tmp=read.table("../result/imp_icogs/merged1.eigenvec")
colnames(tmp)=c("ID",paste0("pc",1:20))
tmp1=unlist(strsplit(tmp$ID,"_"))
tmp$ID=tmp1[seq(1,length(tmp1),2)]
idx=match(phenoicogs$ID,tmp$ID)
phenoicogs=cbind(phenoicogs,tmp[idx,2:11])
phenoicogs$y=NA
phenoicogs$y[which(phenoicogs$Behaviour1==1)]=1
phenoicogs$y[which(is.na(phenoicogs$Behaviour1))]=0
phenoonco=read.table("../data/concept_750_zhang_onco_pheno_v15_02_corrected_age.txt",header=T,sep="\t")
phenoonco$ID=phenoonco$Onc_ID
colnames(phenoonco)=gsub("PC_","pc",colnames(phenoonco))
phenoonco$y=NA
phenoonco$y[which(phenoonco$Behaviour1==1)]=1
phenoonco$y[which(is.na(phenoonco$Behaviour1))]=0
idx=match(prs_tune[,1],phenoonco$ID)
phenotype1=phenoonco[idx,]
for(p_ind in 1:n.total.prs){
  pheno.prs=cbind.data.frame(phenotype1,prs=prs_tune[,(2+p_ind)])
  #model1 <- glm(y~prs+AGE+EV1+EV2+EV3+EV4+EV5+EV6+EV7+EV8+EV9+EV10, data=pheno.prs,family = "binomial")
  model1 <- glm(y~prs, data=pheno.prs,family = "binomial")
  predicted1 <- predict(model1,pheno.prs, type="response")
  prs_auc_vec_test[p_ind]= as.numeric(auc(pheno.prs$y,predicted1,quiet=T))
  
}
max_ind <- which.max(prs_auc_vec_test)
print(colnames(prs_tune)[max_ind+2])

best_snps <- colnames(prs_tune)[max_ind+2]

scores <- plink_list[[1]]
clump_info <- plink_list[[3]]
#number of snps
print(paste0("nrow(clump_info):",nrow(clump_info)))
best_snp_set <- GetSNPSet(snp_ind = best_snps,
                          scores = scores,
                          clump_info = clump_info)

clump_info_post <- EBayesPostMean(x = clump_info,
                                  y = best_snp_set)
post_coef_mat <- cbind(clump_info_post$BETA_EB_target,
                       clump_info_post$BETA_EB_ref)
colnames(post_coef_mat) <- c("EB_target","EB_ref")

plinklist_eb <- PreparePlinkFileEBayes(
                                       snp_list = snp_list,
                                       clump_info = clump_info,
                                       post_clump_info = clump_info_post,
                                       post_beta = post_coef_mat,
                                       results_dir = out_dir)


scores_eb <- plinklist_eb[[1]]
colnames(scores_eb)
# [1] "SNP"                              "A1"                              
# [3] "clump_r2_0.01_ws_5000_EB_target"  "clump_r2_0.01_ws_5000_EB_ref"    
# [5] "clump_r2_0.01_ws_10000_EB_target" "clump_r2_0.01_ws_10000_EB_ref"   
# [7] "clump_r2_0.05_ws_1000_EB_target"  "clump_r2_0.05_ws_1000_EB_ref"    
# [9] "clump_r2_0.05_ws_2000_EB_target"  "clump_r2_0.05_ws_2000_EB_ref"    
# [11] "clump_r2_0.1_ws_500_EB_target"    "clump_r2_0.1_ws_500_EB_ref"      
# [13] "clump_r2_0.1_ws_1000_EB_target"   "clump_r2_0.1_ws_1000_EB_ref"     
# [15] "clump_r2_0.2_ws_250_EB_target"    "clump_r2_0.2_ws_250_EB_ref"      
# [17] "clump_r2_0.2_ws_500_EB_target"    "clump_r2_0.2_ws_500_EB_ref"      
# [19] "clump_r2_0.5_ws_100_EB_target"    "clump_r2_0.5_ws_100_EB_ref"      
# [21] "clump_r2_0.5_ws_200_EB_target"    "clump_r2_0.5_ws_200_EB_ref"      
# [23] "clump_r2_0.8_ws_62.5_EB_target"   "clump_r2_0.8_ws_62.5_EB_ref"     
# [25] "clump_r2_0.8_ws_125_EB_target"    "clump_r2_0.8_ws_125_EB_ref" 
outfile_prefix <- helper_prefix(out_prefix = outprefix,
                                ebayes = TRUE)
temp.dir <- paste0(out_dir,"temp/")
this_prs_p_other_ <- paste0(temp.dir, outfile_prefix, "prs_p_other_")
#testbfile including tuning and validation
CalculateEBEffectSize_onall <- function(bfile,
                                         testbfile, #for testing PRS
                                  snp_ind,
                                  plink_list,
                                  memory = 8000,
                                  threads = 2,
                                  out_prefix,
                                  results_dir,
                                  params_farm = as.null()){
  print("Executing CalculateEBEffectSize() ... ")
  scores <- plink_list[[1]]
  clump_info <- plink_list[[3]]
  
  best_snp_set <- GetSNPSet(snp_ind = snp_ind,
                            scores = scores,
                            clump_info = clump_info)
  
  clump_info_post <- EBayesPostMean(x = clump_info,
                                    y = best_snp_set)
  post_coef_mat <- cbind(clump_info_post$BETA_EB_target,
                         clump_info_post$BETA_EB_ref)
  colnames(post_coef_mat) <- c("EB_target","EB_ref")
  
  plinklist_eb <- PreparePlinkFileEBayes(snp_list = snp_list,
                                         clump_info = clump_info,
                                         post_clump_info = clump_info_post,
                                         post_beta = post_coef_mat,
                                         results_dir = results_dir)
  
  assign("best_snps_set", best_snp_set, envir = .GlobalEnv)
  assign("unique_infor_post", clump_info_post, envir = .GlobalEnv)
  assign("plink_list_eb", plinklist_eb, envir = .GlobalEnv)
  
  scores_eb <- plinklist_eb[[1]]
  score_eb_file <- as.character(unlist(plinklist_eb["score_eb_file"]))
  write.table(scores_eb,
              file = score_eb_file,
              row.names = F,
              col.names = F,
              quote=F)
  p_values_eb <- plinklist_eb[[2]]
  
  ebayes_prs <- PRSscoreEBayes(bfile = testbfile, #changed to compute score on testbfile
                               eb_plink_list = plinklist_eb,
                               plink_list = plink_list,
                               results_dir = results_dir,
                               out_prefix = out_prefix,
                               params_farm = params_farm)
  
  assign("scores_eb", scores_eb, envir = .GlobalEnv)
  assign("score_eb_file", score_eb_file, envir = .GlobalEnv)
  assign("p_values_eb", p_values_eb, envir = .GlobalEnv)
  
  return(ebayes_prs)
  
}

this_prs_p_other_ <- helper_ebscore_loop(params_farm = params_farm, 
                                         plink2_exec = plink2_exec, bfile = bfile, eb_plink_list = eb_plink_list, 
                                         plink_list = plink_list, pthres = pthres, threads = threads, 
                                         memory = memory, results_dir = results_dir, out_prefix = out_prefix)
scores_eb <- eb_plink_list[[1]]
this_prs_mat <- helper_CombinePRS(scores = scores_eb, pthres = pthres, 
                                  prs_p_other_ = this_prs_p_other_)
print("prs_mat_eb object created")

prs_mat_eb <- CalculateEBEffectSize_onall(bfile = Target_ref_plinkfile,
                                          testbfile = Target_test_plinkfile,
                                          snp_ind = best_snps,
                                          plink_list = plink_list,
                                          out_prefix = outprefix,
                                          results_dir = out_dir,
                                          params_farm = PRS_farm)



PRSTrainValSplit <- function(x,
                             n = 0.50,famtrain=NULL,famtest=NULL,corcutoff=0.98) {
  mat_eb <- x
  if (is.null(famtrain) || is.null(famtest))
  {
    n.test <- ceiling(dim(mat_eb)[1]*n)
    super_tune <- as.data.frame(mat_eb[1:n.test,-c(1:2),drop=F])
    rownames(super_tune)=mat_eb[1:n.test,1]
    super_validate <- as.data.frame(mat_eb[(n.test+1):nrow(mat_eb),-c(1:2),drop=F])
    rowna.es(super_validate)=mat_et[(n.test+1):nrow(mat_eb),1]
  }else
  {
    idx=match(famtrain$V1,mat_eb[,1])
    super_tune <- as.data.frame(mat_eb[idx,-c(1:2),drop=F])
    rownames(super_tune)=famtrain$V1
    idx=match(famtest$V1,mat_eb[,1])
    super_validate <- as.data.frame(mat_eb[idx,-c(1:2),drop=F])
    rownames(super_validate)=famtest$V1
  }
  
  #select/prune tuning prs based on outcomes
  alltunauc=rep(NA,ncol(super_tune))
  #phenotype1=phenotype[match(rownames(super_tune),phenotype$GWAS_ID),]
  idx=match(rownames(super_tune),phenoonco$ID)
  phenotype1=phenoonco[idx,]
  
  for(p_ind in 1:length(alltunauc)){
    pheno.prs=cbind.data.frame(phenotype1,prs=super_tune[,p_ind])
    #model1 <- glm(y~prs+AGE+EV1+EV2+EV3+EV4+EV5+EV6+EV7+EV8+EV9+EV10, data=pheno.prs,family = "binomial")
    model1 <- glm(y~prs, data=pheno.prs,family = "binomial")
    predicted1 <- predict(model1,pheno.prs, type="response")
    alltunauc[p_ind]= as.numeric(pROC::auc(pheno.prs$y,predicted1,quiet=T))
  }
  
  prs_tun_order = order(alltunauc)
  #order prs in terms of AUC smallest to largest
  super_tune1=super_tune[,prs_tun_order]
  removehighcorr=function(dat=0,corcutoff=0.98)
  {
    tmp <- cor(dat)
    tmp[upper.tri(tmp)] <- 0
    diag(tmp) <- 0
    #first remove prs (with smaller AUC) and have high corr with others
    datnew <- dat[,!apply(tmp,2,function(x) any(abs(x) > corcutoff))]
    #tmp2 =tmp[,!apply(tmp,2,function(x) any(abs(x) > corcutoff))]
    return(datnew)
  }
  
  super_tune_clean=removehighcorr(dat=super_tune1,corcutoff = corcutoff)
  print(paste0(ncol(super_tune_clean), ' independent PRS'))
  super_validate_clean=super_validate[,match(colnames(super_tune_clean),colnames(super_validate))]
  all(colnames(super_tune_clean)==colnames(super_validate_clean))
  
  return_list <- list("tune" = super_tune_clean,
                      "validate" = super_validate_clean)
  return(return_list)
}
train_val_list <- PRSTrainValSplit(x = prs_mat_eb,famtrain = famtrain,famtest=famtest)
prs_tune_sl <- train_val_list[[1]]
prs_valid_sl <- train_val_list[[2]]

#get the used snps 

colnames(prs_tune_sl)=gsub(".","_",colnames(prs_tune_sl),fixed=T)
colnames(prs_tune_sl)=gsub("-","_",colnames(prs_tune_sl),fixed=T)
colnames(prs_valid_sl)=gsub(".","_",colnames(prs_valid_sl),fixed=T)
colnames(prs_valid_sl)=gsub("-","_",colnames(prs_valid_sl),fixed=T)
rownames(prs_tune_sl)=famtrain$V1

#use four linear learners

# phenotype1=phenotype[match(famtrain[,1],phenotype$GWAS_ID),]
# traindat=cbind.data.frame(prs_tune_sl)
# phenotype2=phenotype[match(famtest[,1],phenotype$GWAS_ID),]
idx=match(famtrain[,1],phenoonco$ID)
phenotype1=phenoonco[idx,]
traindat=cbind.data.frame(prs_tune_sl)

get_pheno_test=function(famtest)
{
  phenotype=phenoicogs[,c("ID",paste0("pc",1:10),"age","Behaviour1","y")]
  if (!all(famtest$V2 %in% phenoicogs$ID))
  {
    tmp=famtest$V2[!famtest$V2 %in% phenoicogs$ID]
    idx=match(tmp,phenoonco$ID)
    phenotype0=phenoonco[idx,c("ID",paste0("pc",1:10),"age","Behaviour1","y")]
    phenotype=rbind(phenotype,phenotype0)
  }
  phenotype=phenotype[match(famtest$V2,phenotype$ID),]
  phenotype$y=famtest$V6-1
  phenotype=phenotype[,c("ID",paste0("pc",1:10),"age","Behaviour1","y")]
  return(phenotype)
}

phenotest=get_pheno_test(famtest)
idx=match(famtest[,1],phenotest$ID)
phenotype2=phenotest[idx,]

valdat=cbind.data.frame(prs_valid_sl)
set.seed(1000)
# sl <- SuperLearner(Y = phenotype1$y, 
#                    X = traindat, 
#                    family = binomial(),
#                    method="method.AUC",
#                    SL.library = c("SL.glm","SL.step.forward","SL.glmnet","SL.rpart"))
# sl_linear1 <- SuperLearner(Y = phenotype1$y, 
#                    X = traindat, 
#                    family = binomial(),
#                    method="method.AUC",
#                    SL.library = c("SL.glmnet"))
# sl_linear2 <- SuperLearner(Y = phenotype1$y, 
#                            X = traindat, 
#                            family = binomial(),
#                            method="method.AUC",
#                            SL.library = c("SL.bayesglm"))
# sl_linear3 <- SuperLearner(Y = phenotype1$y, 
#                            X = traindat, 
#                            family = binomial(),
#                            method="method.AUC",
#                            SL.library = c("SL.biglasso"))
# sl_linear4 <- SuperLearner(Y = phenotype1$y, 
#                            X = traindat, 
#                            family = binomial(),
#                            method="method.AUC",
#                            SL.library = c("SL.lda"))
# sl_linear5 <- SuperLearner(Y = phenotype1$y, 
#                            X = traindat, 
#                            family = binomial(),
#                            method="method.AUC",
#                            SL.library = c("SL.glm"))
# sl_linear6 <- SuperLearner(Y = phenotype1$y, 
#                            X = traindat, 
#                            family = binomial(),
#                            method="method.AUC",
#                            SL.library = c("SL.lm"))

sl_linear1 <- SuperLearner(Y = phenotype1$y, 
                           X = traindat, 
                           family = binomial(),
                           method="method.AUC",
                           SL.library = c("SL.glmnet"))
sl_linear2 <- SuperLearner(Y = phenotype1$y, 
                           X = traindat, 
                           family = binomial(),
                           method="method.AUC",
                           SL.library = c("SL.glmnet","SL.glm"))
sl_linear3 <- SuperLearner(Y = phenotype1$y, 
                           X = traindat, 
                           family = binomial(),
                           method="method.AUC",
                           SL.library = c("SL.glmnet","SL.lm"))
sl_linear4 <- SuperLearner(Y = phenotype1$y, 
                           X = traindat, 
                           family = binomial(),
                           method="method.AUC",
                           SL.library = c("SL.glmnet","SL.ridge"))
sl_linear5 <- SuperLearner(Y = phenotype1$y, 
                           X = traindat, 
                           family = binomial(),
                           method="method.AUC",
                           SL.library = c("SL.glmnet","SL.bayesglm"))
sl_linear6 <- SuperLearner(Y = phenotype1$y, 
                           X = traindat, 
                           family = binomial(),
                           method="method.AUC",
                           SL.library = c("SL.biglasso"))

#get auc for tuning/validaiton/all
getauc=function(mymodel=sl_nonl)
{
  auc=data.frame(tun=0,val=0,all=0)
  y_pred_tun <- predict(mymodel,traindat, onlySL = TRUE)[[1]]
  
  y_pred_val <- predict(mymodel,valdat, onlySL = TRUE)[[1]]
  pheno.prs=cbind.data.frame(phenotype2,prs=y_pred_val[[1]])
  allprs=data.frame(GWAS_ID=c(famtrain$V1,famtest$V1),prs=c(y_pred_tun,y_pred_val))
  pheno1.prs=merge(phenotype,allprs,by="GWAS_ID")
  #model1 <- glm(y~prs+AGE+EV1+EV2+EV3+EV4+EV5+EV6+EV7+EV8+EV9+EV10, data=pheno1.prs,family = "binomial")
  model1 <- glm(y~prs, data=pheno1.prs,family = "binomial")
  predicted1 <- predict(model1,pheno1.prs, type="response")
  auc$all=as.numeric(auc(pheno1.prs$y,predicted1,quiet=T))
  
  pheno2.prs=pheno1.prs[pheno1.prs$GWAS_ID %in% famtrain[,1],]
  #model2 <- glm(y~prs+AGE+EV1+EV2+EV3+EV4+EV5+EV6+EV7+EV8+EV9+EV10, data=pheno2.prs,family = "binomial")
  model2 <- glm(y~prs, data=pheno2.prs,family = "binomial")
  predicted2 <- predict(model2,pheno2.prs, type="response")
  auc$tun=as.numeric(auc(pheno2.prs$y,predicted2,quiet=T))
  
  pheno3.prs=pheno1.prs[pheno1.prs$GWAS_ID %in% famtest[,1],]
  #model3 <- glm(y~prs+AGE+EV1+EV2+EV3+EV4+EV5+EV6+EV7+EV8+EV9+EV10, data=pheno3.prs,family = "binomial")
  model3 <- glm(y~prs, data=pheno3.prs,family = "binomial")
  predicted3 <- predict(model3,pheno3.prs, type="response")
  auc$val=as.numeric(auc(pheno3.prs$y,predicted3,quiet=T))
  print(auc)
  auc_val=data.frame(auc=auc$val,auc_low=0,auc_high=0)
  library(boot)
  AUCBoot = function(data,indices){
    boot_data = data[indices, ]
    #model4 <- glm(y~prs+AGE+EV1+EV2+EV3+EV4+EV5+EV6+EV7+EV8+EV9+EV10, data=boot_data,family = "binomial")
    model4 <- glm(y~prs, data=boot_data,family = "binomial")
    predicted4 <- predict(model4,boot_data, type="response")
    auc=as.numeric(auc(boot_data$y,predicted4,quiet=T))
    return(c(auc))
  }
  # boot_auc = boot(data =pheno3.prs, statistic = AUCBoot, R = 10000)
  # tmp=boot.ci(boot_auc,type="bca")
  # auc_val$auc_low=tmp$bca[4]
  # auc_val$auc_high=tmp$bca[5]
  # print(auc_val)
  return(list(auc=auc,auc_val=auc_val,allprs=allprs))
  
}

#aucres=getauc(mymodel = sl)
aucres_linear1=getauc(mymodel = sl_linear1)
aucres_linear2=getauc(mymodel = sl_linear2)
aucres_linear3=getauc(mymodel = sl_linear3)
aucres_linear4=getauc(mymodel = sl_linear4)
aucres_linear5=getauc(mymodel = sl_linear5)
aucres_linear6=getauc(mymodel = sl_linear6)


# res=list(prs_mat_eb=prs_mat_eb,best_snps=best_snps,sl=sl,sl_linear1=sl_linear1,sl_linear2=sl_linear2,sl_linear3=sl_linear3,
#          sl_linear4=sl_linear4,sl_linear5=sl_linear5,sl_linear6=sl_linear6,
#          aucres=aucres,aucres_linear1=aucres_linear1,aucres_linear2=aucres_linear2,aucres_linear3=aucres_linear3,aucres_linear4=aucres_linear4,
#          aucres_linear5=aucres_linear5,aucres_linear6=aucres_linear6,
#          train_val_list=train_val_list)
# save(res,file=paste0(out_dir,outprefix,"_runsuper_linear.RData"))
res=list(prs_mat_eb=prs_mat_eb,best_snps=best_snps,sl_linear1=sl_linear1,sl_linear2=sl_linear2,sl_linear3=sl_linear3,
         sl_linear4=sl_linear4,sl_linear5=sl_linear5,sl_linear6=sl_linear6,
         aucres_linear1=aucres_linear1,aucres_linear2=aucres_linear2,aucres_linear3=aucres_linear3,aucres_linear4=aucres_linear4,
         aucres_linear5=aucres_linear5,aucres_linear6=aucres_linear6,
         train_val_list=train_val_list)
save(res,file=paste0(out_dir,outprefix,"_runsuper_linear_comb.RData"))
#runsuper_new1.RData: AUC adjust for all covariates
#runsuper_prs.RData: AUC adjust for nothing
#runsuper_linear.RData: linear learner
