#!/usr/bin/env Rscript
#first step, 2D clumping
#run_CT_SLEB_clump: 2D clumping save results in temp folder
#run_CT_SLEB_score: compute PRS based on 2D clumping results
#run_CT_SLEB_EB_score: use the best snps based on tuning data, get EB weights, compute PRS on all samples (tuning+valiation)
#run_CT_SLEB_super: use tuning data to find the best combination of PRS, and compute PRS on validation

.libPaths(c("/data/BB_Bioinformatics/Kevin/tools/Rpackages",.libPaths()))
setwd("/data/BB_Bioinformatics/Kevin/BCAC/code")

# library(devtools)
# install_github("andrewhaoyu/CTSLEB")
args = commandArgs(trailingOnly=TRUE)
#taski is 1:(length(r2_vec)*length(wc_base_vec))

taski=as.integer(args[1])
library(CTSLEB)
library(data.table)
library(dplyr)
print("clumping")
print(paste0("taski=",taski))

RunClump_swarm <- function(taski=1,
                     plink19_exec,
                     ref_plink,
                     target_plink,
                     out_prefix = as.null(),
                     results_dir,
                     ref_splitfile,
                     target_splitfile,
                     r2_vec = c(0.01,0.05,0.1,0.2,0.5,0.8),
                     wc_base_vec = c(50,100),
                     mem = 8000,
                     threads = 2,
                     params_farm = as.null()) {
  print("executing RunClump()... ")
  if (is.null(params_farm)) {
    print("RunClump() no params_farm")
  } else {
    print("RunClump() params_farm list will be used")
    plink19_exec <- as.character(unlist(params_farm["plink19_exec"]))
    r2_vec <- as.numeric(unlist(params_farm["r2_vec"]))
    wc_base_vec <- as.integer(unlist(params_farm["wc_base_vec"]))
    mem <- as.integer(unlist(params_farm["mem"]))
    threads <- as.integer(unlist(params_farm["threads"]))
  }
  
  temp.dir <- paste0(results_dir,"temp/")
  if (!dir.exists(results_dir)) dir.create(results_dir)
  if (!dir.exists(temp.dir)) dir.create(temp.dir)
  if (is.null(out_prefix)) {
    print("not out_prefix")
    out.prefix <- ""
  } else {
    print(paste0("out_prefix: ", out_prefix))
    out.prefix <- paste0(out_prefix, "_")
  }

  r_ind=ceiling(taski /length(wc_base_vec))
  print(paste0("r_ind:",r_ind))
  #print(r2_vec[r_ind])
  wc_vec <- wc_base_vec/r2_vec[r_ind]
  w_ind=taski %% length(wc_base_vec)+1
  print(w_ind)
  
  pthr <-1
  r2thr <- r2_vec[r_ind]
  kbpthr <- wc_vec[w_ind]
  ref_outfile <- paste0(temp.dir, out.prefix, "ref_CT_rind_",r_ind,"_wcind_",w_ind)
  target_outfile <- paste0(temp.dir, out.prefix,"target_CT_rind_",r_ind,"_wcind_",w_ind)
  Plink19Clump(plink19_exec = plink19_exec,
               bfile = ref_plink,
               clump = ref_splitfile,
               clump_p1 = pthr,
               clump_r2 = r2thr,
               clump_kb = kbpthr,
               threads = threads,
               memory = mem,
               out = ref_outfile,
               params_farm = params_farm
  )
  Plink19Clump(plink19_exec = plink19_exec,
               bfile = target_plink,
               clump = target_splitfile,
               clump_p1 = pthr,
               clump_r2 = r2thr,
               clump_kb = kbpthr,
               threads = threads,
               memory = mem,
               out = target_outfile,
               params_farm = params_farm
  )
  print("RunClump() complete ...")
  
}

WriteSplitTables <- function(x,
                             results_dir,
                             ref_split_file = "sum_ref.txt",
                             target_split_file = "sum_target.txt")
{
  print("executing WriteSplitTables()... ")
  ref_split_file <- paste0(results_dir, ref_split_file)
  target_split_file <- paste0(results_dir,target_split_file)
  print("creating sum_ref object")
  sum_ref <- x[[1]]
  print("creating sum_target object")
  sum_target <- x[[2]]
  
  assign("ref_split_file", ref_split_file, envir = .GlobalEnv)
  assign("target_split_file", target_split_file, envir = .GlobalEnv)
  
  write.table(sum_ref,ref_split_file,col.names = T,row.names = F,quote=F)
  write.table(sum_target,target_split_file,col.names = T,row.names = F,quote=F)
  names <- c("ref_split_file",
             "target_split_file")
  values <- list(ref_split_file,
                 target_split_file)
  write_list <- setNames(values, names)
  print("WriteSplitTables() complete ... ")
  return(write_list)
}

SplitSum_swarm <- function(taski=1,
                           x,
                     results_dir,
                     write_tables = TRUE){
  print("executing SplitSum() ...")
  sum_com <- x
  sum_com_select <- sum_com %>%
    mutate(split_ind =
             ifelse(
               (P<P_ref)|is.na(P_ref),T,F)
    )%>%
    select(SNP,P,P_ref,split_ind)
  
  sum_com_select_other_ref <- sum_com_select %>%
    filter(split_ind==F) %>%
    select(SNP,P_ref) %>%
    rename(P = P_ref)
  
  sum_com_select_tar_ref <- sum_com_select %>%
    filter(split_ind==T) %>%
    select(SNP,P)
  
  sum_list <- list(sum_com_select_other_ref,
                   sum_com_select_tar_ref)
  
  if (write_tables) {
    write_list <- WriteSplitTables(x = sum_list,
                                   results_dir = results_dir,
                                   ref_split_file = paste0("sum_ref",taski,".txt"),
                                   target_split_file = paste0("sum_target",taski,".txt"))
    
    return(write_list)
  } else {
    print(paste0("WriteSplitTables() was not performed"))
    return(sum_list)
  }
  print("SplitSum() complete ... ")
  return(sum_list)
  
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
EUR_ref_plinkfile <- config[3,2] #"../result/PRS1/euro_onco_tuning",used for clumping
Target_ref_plinkfile <- config[4,2] #"../result/PRS1/asian_onco_tuning",used for clumping and tuning
#Specify the tuning and validation data directory
Target_test_plinkfile <- config[5,2] #"../result/PRS1/tuning_validation",used for computing PRS,including tuning and validation.
#output folder
out_dir = config[6,2] # "/data/BB_Bioinformatics/Kevin/BCAC/result/PRS1/ctsleb/EAS/"
outprefix <- config[7,2] #"EAS_CTSLEB"
plink19_exec=config[8,2] #"/usr/local/apps/plink/1.9.0-beta4.4/plink"
plink2_exec=config[9,2] #"/usr/local/apps/plink/2.3-alpha/plink2"
if (! dir.exists(out_dir))
{
  dir.create(out_dir)
  dir.create(paste0(out_dir,"temp")) #the clumped results saved here
}

if (!file.exists(EUR_sumstats_file) || !file.exists(Target_sumstats_file) || !file.exists(paste0(EUR_ref_plinkfile,".bim")) ||
    !file.exists(paste0(Target_ref_plinkfile,".bim")) || !file.exists(paste0(Target_test_plinkfile,".bim"))||
    !file.exists(plink19_exec) || !file.exists(plink2_exec))
  stop("some inputs files are missing")

#Read sumstats
sum_EUR <- fread(EUR_sumstats_file,header=T)
sum_Target <- fread(Target_sumstats_file,header=T)
sum_com <- AlignSum(sum_target = sum_Target,
                    sum_ref = sum_EUR)


PRS_farm <- SetParamsFarm(plink19_exec = plink19_exec,
                          plink2_exec = plink2_exec,threads = 2,mem=64000)

write_list <- SplitSum_swarm(taski=taski,
                      x = sum_com,
                       results_dir = out_dir,
                       write_tables = TRUE)

ref_splitfile <- unlist(write_list["ref_split_file"])
target_splitfile <- unlist(write_list["target_split_file"])

snp_list <- RunClump_swarm(taski=taski,
                     params_farm = PRS_farm,
                     plink19_exec = plink19_exec,
                     ref_plink = EUR_ref_plinkfile, #ref tuning
                     target_plink = Target_ref_plinkfile, #target tuning
                     ref_splitfile = ref_splitfile,
                     target_splitfile = target_splitfile,
                     out_prefix = outprefix,
                     results_dir = out_dir)

#cmd=data.frame(cmd="/data/BB_Bioinformatics/Kevin/BCAC/code/run_CT_SLEB_clump.R ",taski=1:(length(PRS_farm$r2_vec)*length(PRS_farm$wc_base_vec)),config="/data/BB_Bioinformatics/Kevin/BCAC/code/CT_SLEB_EAS.config")
#target="EAS"
# write.table(cmd,file=paste0("CTSLEB_",target,"_runclump.swarm"),row.names = F,col.names = F,quote=F,sep="\t")

#swarm  -f  /data/BB_Bioinformatics/Kevin/BCAC/code/CTSLEB_EAS_runclump.swarm --module R -g 72 --time=3-00:00:00 --gres=lscratch:72
#replace EAS with EUR to get the new swarmfile
#swarm  -f  /data/BB_Bioinformatics/Kevin/BCAC/code/CTSLEB_EUR_runclump.swarm --module R -g 72 --time=3-00:00:00 --gres=lscratch:72
