#!/usr/bin/Rscript

#### test correlation between parallel chains of peaky_fs
#### update to run all LIBs

# load
library(peaky)
library(data.table)

# functions
MPPC_cor = function(result_A, result_B){
  cor(result_A$rjmcmc_pos, result_B$rjmcmc_pos)
}

MPPC_cor_by_dir = function(baits_rjmcmc_dir_A, baits_rjmcmc_dir_B, baits_dir, cor_threshold=.75){
  files_A = grep("bait_rjmcmc_[0-9]+.rds$",list.files(baits_rjmcmc_dir_A),value=TRUE)
  files_B = grep("bait_rjmcmc_[0-9]+.rds$",list.files(baits_rjmcmc_dir_B),value=TRUE)
  
  if (!dir.exists(baits_dir)) {
    stop("Baits directory not found.")
  }
  
  L = paste0(baits_dir, "/MPPCcor_log.txt")
  
  peaky:::note(L,TRUE,"Calculating MPPC cor for interpreted peaky results in:\nA: ",
               baits_rjmcmc_dir_A,"\nB: ",baits_rjmcmc_dir_B)
  
  files_all = union(files_A,files_B)
  files_both = intersect(files_A,files_B)
  
  paths_A = paste0(baits_rjmcmc_dir_A,"/",files_both)
  paths_B = paste0(baits_rjmcmc_dir_B,"/",files_both)
  peaky:::note(L,FALSE,"Unique bait results across A and B: ",length(files_all),
               "\nBait results in both A and B: ", length(files_both))
  
  peaky:::note(L,TRUE,"Loading common results and calculating MPPC...")
  
  results_A = lapply(paths_A,readRDS)
  results_B = lapply(paths_B,readRDS)
  
  cors = mapply(MPPC_cor, results_A, results_B)
  
  cor_table = data.table(results_file=files_both, MPPC_cor=cors)
  file_table = data.table(results_file=files_all, in_A=files_all%in%files_A, in_B=files_all%in%files_B,
                          bait_path=paste0(baits_dir,"/",gsub("rjmcmc_","",files_all)))
  
  file_cor_table = merge(file_table, cor_table, by="results_file", all.x = TRUE)
  file_cor_table[,redo:=!is.na(MPPC_cor) & MPPC_cor<cor_threshold]
  
  table_path = paste0(baits_dir,"/","MPPCcor.csv")
  baitlist_path = paste0(baits_dir,"/","baitlist_MPPCcor_sub_",cor_threshold,".txt")
  
  peaky:::note(L,TRUE,"Saving MPPC correlation overview:\n",table_path,
               "\nSaving baitlist for baits with MPPC correlation below ",cor_threshold,":\n",baitlist_path)
  fwrite(file_cor_table, table_path)
  write(file_cor_table[redo==TRUE,bait_path], baitlist_path)
  peaky:::note(L,FALSE,"Done.")
  
  return(file_cor_table)
}

# CHIC samples
sample_list = c("B80T5", "Hs578T", "MCF10A", "MCF7", "MDAMB231", "T47D")

for(k in c("SCHIC")){
  
  run_dir  = paste0("/working/lab_georgiat/jonathB/CaptureHiC/peaky/peaky_fs/peaky_chain_runs/",k,"_split_baits/")
  
  for(i in sample_list){
    
    # use full paths
    baits_rjmcmc_dir_A = paste0(run_dir,i,"/baits_rjmcmc1/") 
    baits_rjmcmc_dir_B = paste0(run_dir,i,"/baits_rjmcmc2/")
    baits_dir          = paste0(run_dir,i,"/baits/")
    cor_threshold      = 0.75
    
    # run
    MPPC_cor_by_dir(baits_rjmcmc_dir_A, baits_rjmcmc_dir_B, baits_dir, 0.75)
    
  }
}