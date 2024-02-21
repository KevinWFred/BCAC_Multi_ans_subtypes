#!/usr/bin/env Rscript
#second step,compute 2D CT PRS 
.libPaths(c("/data/BB_Bioinformatics/Kevin/tools/Rpackages",.libPaths()))

args = commandArgs(trailingOnly=TRUE)
taski=as.integer(args[1])
#taski is 1:length(PRS_farm$pthres)

setwd("/data/BB_Bioinformatics/Kevin/BCAC/code")
library(CTSLEB)
library(data.table)
library(dplyr)
print("scoring")
print(paste0("taski=",taski))


helper_PreparePlinkFile <- function(taski=1,
                                    plink_list,
                                    results_dir)
  
{
  
  temp.dir <- paste0(results_dir,"temp/")
  scores <- plink_list[[1]]
  score_file <- paste0(temp.dir,"score_file")
  write.table(scores,
              file = paste0(score_file,taski) ,
              row.names = F,
              col.names = F,
              quote=F)
  
  p_values <- plink_list[[2]]
  p_value_file <- paste0(temp.dir,"p_value_file")
  
  unique_infor <- plink_list[[3]]
  
  q_range <- plink_list[[4]]
  q_range_file <- paste0(temp.dir,"q_range_file")
  if (!file.exists(q_range_file))
    write.table(q_range,
              file = q_range_file,
              row.names = F,
              col.names = F,
              quote=F)
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

helper_score_loop_swarm <- function(taski=1,
                              plink2_exec,
                              bfile,
                              plink_list,
                              pthres,
                              threads,
                              memory,
                              results_dir,
                              params_farm = as.null(),
                              out_prefix = as.null()){
  
  if (is.null(params_farm)) {
    #print("no params_farm")
  } else {
    # print("params_farm list will be used")
    memory <- as.integer(unlist(params_farm["mem"]))
    threads <- as.integer(unlist(params_farm["threads"]))
    pthres <- as.numeric(unlist(params_farm["pthres"]))
  }
  
  temp.dir <- paste0(results_dir,"temp/")
  out.prefix <- helper_prefix(out_prefix = out_prefix,
                              ebayes = FALSE)
  
  scores <- plink_list[[1]]
  p_values <- plink_list[[2]]
  unique_infor <- plink_list[[3]]
  p_value_file <- paste0(as.character(unlist(plink_list["p_value_file"])),taski)
  score_file <- paste0(as.character(unlist(plink_list["score_file"])),taski)
  q_range_file <- as.character(unlist(plink_list["q_range_file"]))
  prs_p_other_ <- paste0(temp.dir, out.prefix, "prs_p_other_")
  pthres <- pthres
  assign("prs_p_other_", prs_p_other_, envir = .GlobalEnv)
  p_values_temp <- p_values
  
  k1=taski
  idx <- which(unique_infor$P_ref <= pthres[k1])
  print(paste0("pthres: ", pthres[k1]))
  # print(class(pthres[k1]))
  p_values_temp$P[idx] <- 0
  # print(paste0("Number of variants less than pthres :",
  #              sum(p_values_temp$P == 0)))
  # print(paste0("writing ", p_value_file))
  write.table(p_values_temp,
              file = p_value_file,
              col.names = F,
              row.names = F,
              quote=F)
  score.col_nums <- ncol(scores)
  print(paste0("pthres_idx=",k1))
  plink2score(plink2_exec = plink2_exec,
              bfile = bfile,
              q_range_file = q_range_file,
              p_value_file = p_value_file,
              score_file = score_file,
              score_col_nums = score.col_nums,
              results_dir = results_dir,
              pthres_idx = k1,
              threads = threads,
              memory = memory,
              out = prs_p_other_,
              params_farm=params_farm)
  
  
  return(prs_p_other_)
}
PRSscore_swarm <- function(taski=1,
                     plink2_exec = "plink2 ",
                     bfile,
                     plink_list,
                     pthres = c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0),
                     threads = 4,
                     memory = 8000,
                     results_dir,
                     out_prefix = as.null(),
                     params_farm=as.null()){
  
  print("executing PRSscore()... ")
  print(paste0("pthres=",pthres[taski]))
  if (is.null(params_farm)) {
    #print("no params_farm")
  } else {
    # print("params_farm list will be used")
    plink2_exec <- as.character(unlist(params_farm["plink2_exec"]))
    memory <- as.integer(unlist(params_farm["mem"]))
    threads <- as.integer(unlist(params_farm["threads"]))
    pthres <- as.numeric(unlist(params_farm["pthres"]))
  }
  
  plink2_exec <- plink2_exec
  bfile <- bfile
  pthres <- pthres
  threads <- threads
  memory <- memory
  results_dir <- results_dir
  out_prefix <- out_prefix
  params_farm <- params_farm
  plink_list <- plink_list
  
  prs_p_other_ <- helper_score_loop_swarm(taski=taski,
                                    plink2_exec = plink2_exec,
                                    bfile = bfile,
                                    pthres = pthres,
                                    threads = threads,
                                    memory = memory,
                                    plink_list = plink_list,
                                    results_dir = results_dir,
                                    out_prefix = out_prefix)
  
}



#Read configure file
# An input file looks like this:
# reference_sumstat /data/BB_Bioinformatics/Kevin/PRS_EASLC/result/EUnever_CTSLEB.sumdat
# target_sumstat /data/BB_Bioinformatics/Kevin/PRS_EASLC/result/EASnever_CTSLEB.sumdat
# reference_clumping_plink /data/BB_Bioinformatics/Kevin/tools/CTSLEB/EUR/chr_all_FLCCA
# target_clumping_plink /data/BB_Bioinformatics/Kevin/PRS_EASLC/result/FLCCA_1kg
# target_test_plink /data/BB_Bioinformatics/Kevin/PRS_EASLC/result/FLCCA_1kg
# output_dir /data/BB_Bioinformatics/Kevin/PRS_EASLC/result/CTSLEB_Never_swarm/
#   outprefix EASnever_CTSLEB
# plink19_exec /usr/local/apps/plink/1.9/plink
# plink2_exec /usr/local/apps/plink/2.3-alpha/plink2

config=read.table(args[2]) #config input file
#config=read.table("CT_SLEB_EAS.config")
if (!file.exists(args[2])) stop("no input file.")
if (nrow(config)!=9 || ncol(config)!=2) stop("input file format is not right")

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

if (! dir.exists(out_dir))
{
  dir.create(out_dir)
  dir.create(paste0(out_dir,"temp"))
}

sum_EUR <- fread(EUR_sumstats_file,header=T)
sum_Target <- fread(Target_sumstats_file,header=T)
sum_com <- AlignSum(sum_target = sum_Target,
                    sum_ref = sum_EUR)

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
# names(plink_list):data used to compute 2D CT PRS
# [1] "scores"       "p_values"     "unique_infor" "q_range" 
print(dim(plink_list$scores)) #1799881      14

file_list <- helper_PreparePlinkFile(taski=taski,plink_list = plink_list,
                                     results_dir = out_dir)
# $score_file
# [1] "/data/BB_Bioinformatics/Kevin/BCAC/result/PRS1/ctsleb/EAS/temp/score_file"
# 
# $p_value_file
# [1] "/data/BB_Bioinformatics/Kevin/BCAC/result/PRS1/ctsleb/EAS/temp/p_value_file"
# 
# $q_range_file
# [1] "/data/BB_Bioinformatics/Kevin/BCAC/result/PRS1/ctsleb/EAS/temp/q_range_file"
plink_list <- c(plink_list,file_list)
# names(plink_list)
# [1] "scores"       "p_values"     "unique_infor" "q_range"      "score_file"  
# [6] "p_value_file" "q_range_file"
prs_mat <- PRSscore_swarm(taski=taski,
                    params_farm = PRS_farm,
                    plink2_exec = plink2_exec,
                    bfile = Target_test_plinkfile,
                    plink_list = plink_list,
                    threads = PRS_farm$threads,
                    memory = PRS_farm$mem,
                    out_prefix = outprefix,
                    results_dir = out_dir)


# cmd=data.frame(cmd="/data/BB_Bioinformatics/Kevin/BCAC/code/run_CT_SLEB_score.R ",taski=1:(length(PRS_farm$pthres)),config="/data/BB_Bioinformatics/Kevin/BCAC/code/CT_SLEB_EAS.config")
#target="EAS"
#write.table(cmd,file=paste0("CTSLEB_",target,"_runscore.swarm"),row.names = F,col.names = F,quote=F,sep="\t")
#swarm  -f  /data/BB_Bioinformatics/Kevin/BCAC/code/CTSLEB_EAS_runscore.swarm --module R -g 72 --time=3-00:00:00 --gres=lscratch:72

# cmd=data.frame(cmd="/data/BB_Bioinformatics/Kevin/BCAC/code/run_CT_SLEB_score.R ",taski=1:(length(PRS_farm$pthres)),config="/data/BB_Bioinformatics/Kevin/BCAC/code/CT_SLEB_EUR.config")
#target="EUR"
#write.table(cmd,file=paste0("CTSLEB_",target,"_runscore.swarm"),row.names = F,col.names = F,quote=F,sep="\t")
#swarm  -f  /data/BB_Bioinformatics/Kevin/BCAC/code/CTSLEB_EUR_runscore.swarm --module R -g 72 --time=3-00:00:00 --gres=lscratch:72
