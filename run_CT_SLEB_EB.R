#!/usr/bin/env Rscript
#step 3: compute EB PRS, no need. It was done in run_CT_SLEB_super.R
.libPaths(c("/data/BB_Bioinformatics/Kevin/tools/Rpackages",.libPaths()))
args = commandArgs(trailingOnly=TRUE)
taski=as.integer(args[1]) #taski=1:(length(PRS_farm$pthres)
setwd("/data/BB_Bioinformatics/Kevin/BCAC/code")
library(CTSLEB)
library(data.table)
library(dplyr)

set.seed(1000)
print("Empical Bayes estimation")
print(paste0("taski=",taski))

PreparePlinkFile=function (snp_list = snp_list, sum_com, pthres = c(5e-08, 5e-07, 
                                                   5e-06, 5e-05, 5e-04, 0.005, 0.05, 0.5, 1), results_dir, params_farm = as.null()) 
{
  print("executing PreparePlinkFile()... ")
  if (is.null(params_farm)) {
    print("no params_farm")
  }
  else {
    print("params_farm list will be used")
    pthres <- as.numeric(unlist(params_farm["pthres"]))
  }
  unique_id <- unique(rbindlist(snp_list, use.name = FALSE))
  names(unique_id) <- "SNP"
  unique_infor <- left_join(unique_id, sum_com, by = "SNP")
  print("unique_infor dataframe complete")
  n_col <- length(snp_list)
  n_row <- nrow(unique_infor)
  beta_mat <- matrix(unique_infor$BETA, nrow = n_row, ncol = n_col)
  names <- rep("c", n_col)
  temp <- 1
  for (ldx in seq(n_col)) {
    LD <- snp_list[[ldx]]
    names(LD) <- "SNP"
    idx <- which(unique_infor$SNP %in% LD$SNP == F)
    beta_mat[idx, ldx] <- 0
    names[ldx] <- names(snp_list[[ldx]])
  }
  colnames(beta_mat) <- names
  scores <- data.frame(SNP = unique_id, A1 = unique_infor$A1, 
                       beta_mat)
  print("scores dataframe complete")
  p_values <- data.frame(SNP = unique_id, P = unique_infor$P)
  print("p_values dataframe complete")
  q_range <- CreateQRange(pthres)
  names <- c("scores", "p_values", "unique_infor", "q_range")
  values <- list(scores, p_values, unique_infor, q_range)
  list <- setNames(values, names)
  print("PreparePlinkFile() complete ...")
  return(list)
}

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

plink2score <- function(plink2_exec = "plink2 ",
                        bfile,
                        q_range_file,
                        p_value_file,
                        score_file,
                        score_col_nums,
                        results_dir,
                        pthres_idx,
                        threads = 2,
                        memory = 8000,
                        out,
                        params_farm=as.null())
{
  
  if (is.null(params_farm)) {
    #print("no params_farm")
  } else {
    #print("params_farm list will be used")
    plink2_exec <- as.character(unlist(params_farm["plink2_exec"]))
    mem <- as.character(unlist(params_farm["mem"]))
    threads <- as.character(unlist(params_farm["threads"]))
  }
  
  prs_out <- paste0(prs_p_other_,pthres_idx)
  cmd=paste0(plink2_exec, " ",
             "--bfile ", bfile, " ",
             "--q-score-range ", q_range_file, " ", p_value_file, " ",
             "--score-col-nums 3-", score_col_nums, " ",
             "--score ", score_file, " cols=+scoresums,-scoreavgs ",
             "--threads ", threads, " ",
             "--memory ", memory, " ",
             "--out ", prs_out)
  print(cmd)
  system(paste0(plink2_exec, " ",
                "--bfile ", bfile, " ",
                "--q-score-range ", q_range_file, " ", p_value_file, " ",
                "--score-col-nums 3-", score_col_nums, " ",
                "--score ", score_file, " cols=+scoresums,-scoreavgs ",
                "--threads ", threads, " ",
                "--memory ", memory, " ",
                "--out ", prs_out)
  )
}

helper_ebscore_loop_swarm <- function(taski=1,
                                      plink2_exec,
                                      bfile,
                                      eb_plink_list,
                                      plink_list,
                                      pthres,
                                      threads,
                                      memory,
                                      results_dir,
                                      out_prefix = as.null(),
                                      params_farm = as.null()){
  
  temp_dir <- paste0(results_dir,"temp/")
  
  outfile_prefix <- helper_prefix(out_prefix = out_prefix,
                                  ebayes = TRUE)
  
  scores <- eb_plink_list[[1]]
  p_values_eb <- eb_plink_list[[2]]
  clump_info <- plink_list[[3]]
  p_values_eb_file <- paste0(as.character(unlist(eb_plink_list["p_values_eb_file"])))
  score_eb_file <- paste0(as.character(unlist(eb_plink_list["score_eb_file"])))
  q_range_file <- as.character(unlist(plink_list["q_range_file"]))
  prs_p_other_ <- paste0(temp_dir, outfile_prefix, "prs_p_other_")
  assign("prs_p_other_", prs_p_other_, envir = .GlobalEnv)
  p_values_temp <- p_values_eb
  k1=taski
  idx <- which(clump_info$P_ref <= pthres[k1])
  # print(paste0("pthres: ", pthres[k1]))
  #print(class(pthres[k1]))
  p_values_temp$P[idx] <- 0
  # print(paste0("Number of variants less than pthres :",
  #              sum(p_values_temp$P == 0)))
  # print(paste0("writing ", p_values_eb_file))
  write.table(p_values_temp,
              file = p_values_eb_file,
              col.names = F,
              row.names = F,
              quote=F)
  score_col_nums <- ncol(scores)
  plink2score(params_farm = params_farm,
              plink2_exec = plink2_exec,
              bfile = bfile,
              q_range_file = q_range_file,
              p_value_file = p_values_eb_file,
              score_file = score_eb_file,
              score_col_nums = score_col_nums,
              results_dir = results_dir,
              pthres_idx = k1,
              threads = threads,
              out = prs_p_other_,
              memory = memory)
  
  return(prs_p_other_)
}

PRSscoreEBayes_swarm <- function(taski=1,
                                 plink2_exec = "plink2 ",
                                 bfile,
                                 eb_plink_list,
                                 plink_list,
                                 pthres = c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0),
                                 threads = 2,
                                 memory = 8000,
                                 results_dir,
                                 out_prefix = as.null(),
                                 params_farm=as.null()){
  
  print("executing PRSscoreEbayes()... ")
  if (is.null(params_farm)) {
    #print("no params_farm")
  } else {
    #print("params_farm list will be used")
    plink2_exec <- as.character(unlist(params_farm["plink2_exec"]))
    memory <- as.integer(unlist(params_farm["mem"]))
    threads <- as.integer(unlist(params_farm["threads"]))
    pthres <- as.numeric(unlist(params_farm["pthres"]))
  }
  
  this_prs_p_other_ <- helper_ebscore_loop_swarm(taski=taski,
                                                 params_farm = params_farm,
                                                 plink2_exec = plink2_exec,
                                                 bfile = bfile,
                                                 eb_plink_list = eb_plink_list,
                                                 plink_list = plink_list,
                                                 pthres = pthres,
                                                 threads = threads,
                                                 memory = memory,
                                                 results_dir = results_dir,
                                                 out_prefix = out_prefix)
  # scores_eb <- eb_plink_list[[1]]
  # this_prs_mat <- helper_CombinePRS(scores = scores_eb,
  #                                   pthres = pthres,
  #                                   prs_p_other_ = this_prs_p_other_)
  print("prs_mat_eb object created")
  #return(this_prs_mat)
}

PreparePlinkFileEBayes <- function(taski=taski,
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
  score_eb_file <- paste0(temp_dir,"score_eb_file",taski)
  this_p_values <- data.frame(SNP = unique_id,
                              P = post_clump_info$P)
  p_values_eb_file <-paste0(temp_dir,"p_values_eb_file",taski)
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

CalculateEBEffectSize_swarm <- function(taski=1,
                                        bfile,
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
  
  plinklist_eb <- PreparePlinkFileEBayes(taski=taski,
                                         snp_list = snp_list,
                                         clump_info = clump_info,
                                         post_clump_info = clump_info_post,
                                         post_beta = post_coef_mat,
                                         results_dir = results_dir)
  
  assign("best_snps_set", best_snp_set, envir = .GlobalEnv)
  assign("unique_infor_post", clump_info_post, envir = .GlobalEnv)
  assign("plink_list_eb", plinklist_eb, envir = .GlobalEnv)
  
  scores_eb <- plinklist_eb[[1]]
  score_eb_file <- paste0(as.character(unlist(plinklist_eb["score_eb_file"])))
  write.table(scores_eb,
              file = score_eb_file,
              row.names = F,
              col.names = F,
              quote=F)
  p_values_eb <- plinklist_eb[[2]]
  
  ebayes_prs <- PRSscoreEBayes_swarm(taski=taski,
                                     bfile = bfile,
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
CreateQRange <- function(x)
{
  print("executing CreateQRange()... ")
  if (is.null(x)) {
    pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,
                5E-02,5E-01,1.0)
    print(paste0("pthres is NULL...using default values "))
  } else {
    pthres <- x
  }
  
  n_pthres <- length(pthres)
  q_range <- data.frame(filename = rep("p_tar",n_pthres),
                        small_P = rep(0,n_pthres),
                        max_P = rep(1,n_pthres))
  
  temp <- 1
  
  for(k2 in 1:length(pthres)){
    
    q_range[temp,1] = paste0("p_tar_",k2)
    q_range[temp,3] = pthres[k2]
    temp = temp+1
  }
  q_range = q_range[1:(temp-1),]
  print("CreateQRange() complete... ")
  return(q_range)
}

config=read.table(args[2]) #config input file
#config=read.table("CT_SLEB_EAS.config")
if (!file.exists(args[2])) stop("no input file.")
#if (nrow(config)!=9 || ncol(config)!=2) stop("input file format is not right")

#Specify the directory for the summary statistics
EUR_sumstats_file <- config[1,2] #"/data/BB_Bioinformatics/Kevin/BCAC/result/PRS1/ctsleb/EUR/EUR_CTSLEB.sumdat" #  reference population
Target_sumstats_file <- config[2,2]  #"/data/BB_Bioinformatics/Kevin/BCAC/result/PRS1/ctsleb/EAS/EAS_CTSLEB.sumdat" #  target population

#Specify the directory to the reference data for clumping purpose
EUR_ref_plinkfile <- config[3,2] #"../result/PRS1/euro_onco_tuning"
Target_ref_plinkfile <- config[4,2] #"../result/PRS1/asian_onco_tuning"
#Specify the tuning and validation data directory
Target_test_plinkfile <- config[5,2] #"../result/PRS1/tuning_validation"
#output folder
out_dir = config[6,2] # "/data/BB_Bioinformatics/Kevin/BCAC/result/PRS1/ctsleb/EAS/"
outprefix <- config[7,2] #"EAS_CTSLEB"
plink19_exec=config[8,2] #"/usr/local/apps/plink/1.9/plink"
plink2_exec=config[9,2] #"/usr/local/apps/plink/2.3-alpha/plink2"
#famtrainfile=config[10,2] # /data/BB_Bioinformatics/Kevin/PRS_EASLC/result/EAS_never_train1.fam

#used to select tuning
#famtrain=read.table(famtrainfile)

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
                          plink2_exec = plink2_exec,threads = 2,mem=64000)
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
scores <- plink_list[[1]]
#prs_p_other_="../result/CTSLEB_Never/temp/EASnever_CTSLEB_train_prs_p_other_"
prs_p_other_ <- paste0(temp.dir, out.prefix, "prs_p_other_")
prs_mat <- helper_CombinePRS(scores = scores,
                             pthres = pthres,
                             prs_p_other_ = prs_p_other_)

famtrain=read.table(paste0(Target_ref_plinkfile,".fam"))
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
#calculate eb effect using EB coefficients
prs_mat_eb <- CalculateEBEffectSize_swarm(taski=taski,
                                    bfile = Target_ref_plinkfile,
                                    snp_ind = best_snps,
                                    plink_list = plink_list,
                                    out_prefix = outprefix,
                                    results_dir = out_dir,
                                    params_farm = PRS_farm)
save(best_snps,file=paste0(out_dir,"best_snps.RData"))
#cmd=data.frame(cmd="/data/BB_Bioinformatics/Kevin/BCAC/code/run_CT_SLEB_EB.R ",taski=1:(length(PRS_farm$pthres)),config="/data/BB_Bioinformatics/Kevin/BCAC/code/CT_SLEB_EAS.config")
#target="EAS"
#write.table(cmd,file=paste0("CTSLEB_",target,"_runEB.swarm"),row.names = F,col.names = F,quote=F,sep="\t")
#swarm  -f  /data/BB_Bioinformatics/Kevin/BCAC/code/CTSLEB_EAS_runEB.swarm --module R -g 72 --time=3-00:00:00 --gres=lscratch:72
#cmd=data.frame(cmd="/data/BB_Bioinformatics/Kevin/BCAC/code/run_CT_SLEB_EB.R ",taski=1:(length(PRS_farm$pthres)),config="/data/BB_Bioinformatics/Kevin/BCAC/code/CT_SLEB_EUR.config")
#target="EUR"
#write.table(cmd,file=paste0("CTSLEB_",target,"_runEB.swarm"),row.names = F,col.names = F,quote=F,sep="\t")
#swarm  -f  /data/BB_Bioinformatics/Kevin/BCAC/code/CTSLEB_EUR_runEB.swarm --module R -g 72 --time=3-00:00:00 --gres=lscratch:72

