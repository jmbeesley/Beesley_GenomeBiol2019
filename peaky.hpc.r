
#!/usr/bin/Rscript

## Filesystem version
# To make the software more convenient for use on a cluster, Peaky has wrapper functions that interact directly with the filesystem. 
# These wrappers generally take integer arguments specifying which inputs (out of some directory) to work with, and save results directly to disk.
# The workflow is similar to (and arguably simpler than) that of the interactive version, and allows for parallelization at the following computationally intensive steps:
### Removing initial noise (going from raw readcounts to adjusted readcounts) with model_bin_fs
### Separating direct from collateral contacts with peaky_fs()
#These two functions are demonstrated within for-loops in the example below, but would normally be run in parallel (e.g. via a scheduler).


library(peaky)

args     = commandArgs(TRUE)
args[1] -> lib

# read in data
Ffile               = "/working/lab_georgiat/jonathB/CaptureHiC/peaky/peaky_fs/raw/fragments.bed"
#colnames(fragments) = c("chrom","chromStart","chromEnd","ID")
Ifile               = paste0("/working/lab_georgiat/jonathB/CaptureHiC/peaky/peaky_fs/raw/interactions_", lib, ".tsv")

# set directories
results_dir         = paste0("/working/lab_georgiat/jonathB/CaptureHiC/peaky/peaky_fs/peakyRuns_specifyBaits/results")
lib_dir             = paste0("/working/lab_georgiat/jonathB/CaptureHiC/peaky/peaky_fs/peakyRuns_specifyBaits/", lib)
bins_dir            = paste0(lib_dir, "/bins")
fits_dir            = paste0(lib_dir, "/fits")
baits_dir           = paste0(lib_dir, "/baits")
rjmcmc_dir          = paste0(lib_dir, "/rjmcmc")
baits_of_interest   = "/mnt/lustre/working/lab_georgiat/jonathB/CaptureHiC/peaky/peaky_fs/peakyRuns_specifyBaits/baitdir/baits.of.interest.txt"
baits_filename      = paste0(baits_dir,"/bait_", readLines(baits_of_interest), ".rds")
baitlist            = paste0(baits_dir,"/baitlist.txt")
omega_power         = -5

# bin interactions based on distance
bin_interactions_fs(Ifile, Ffile, output_dir = bins_dir, bins = 5)


# null models
for(bin_index in 1:5){
  model_bin_fs(bins_dir, bin_index, output_dir = fits_dir, subsample_size = 1000)
}

# 
split_baits_fs(bins_dir, residuals_dir = fits_dir, indices=1:5, output_dir = baits_dir)




bait_indexes = which(readLines(baitlist) %in% c(paste0("/working/lab_georgiat/jonathB/CaptureHiC/peaky/peaky_fs/peakyRuns_specifyBaits/", lib, "/baits/bait_486406.rds"),
                                                paste0("/working/lab_georgiat/jonathB/CaptureHiC/peaky/peaky_fs/peakyRuns_specifyBaits/", lib, "/baits/bait_596030.rds"),
                                                paste0("/working/lab_georgiat/jonathB/CaptureHiC/peaky/peaky_fs/peakyRuns_specifyBaits/", lib, "/baits/bait_596031.rds")))


###
for(i in bait_indexes){
  peaky_fs(baitlist, i, output_dir = rjmcmc_dir, omega_power = omega_power, iterations = 1e+07)
}




rjmcmc_list(rjmcmc_dir)

rjmcmclist       = paste0(rjmcmc_dir,"/rjmcmclist.txt")
baits_rjmcmc_dir = paste0(lib_dir,"/baits_rjmcmc")

for(i in 1:3){

  P = interpret_peaky_fs(rjmcmclist, i, baits_dir ,baits_rjmcmc_dir, omega_power)
  
  P2 = P$output

  pdf(paste0(results_dir, "/", lib, "_bait_", P2$baitID[1], ".pdf"))  
  
  par(mfrow=c(3,1))
  zoom = P2[abs(P2$dist)<2e6]
  plot(x=zoom$dist, xlab="Distance from bait (bp)",
       y=zoom$residual, ylab="Adjusted readcount",
       main = paste0("Bait_", P2$baitID[1]))
  
  plot(x=zoom$dist, xlab="Distance from bait (bp)",
       y=zoom$beta_mean, ylab="Mean contact strength",
       col="green")
  
  plot(x=zoom$dist, xlab="Distance from bait (bp)",
       y=zoom$rjmcmc_pos, ylab="MPPC",
       col="blue")
  
  dev.off()

  P3 = P2[P2$rjmcmc_pos>0.2]
  write.table(P2, paste0(results_dir, "/", lib, ".rjmcmc_all", i, ".tsv"), quote=F, row.names=F, sep='\t')
  write.table(P3, paste0(results_dir, "/", lib, ".rjmcmc_high", i, ".tsv"), quote=F, row.names=F, sep='\t')
  
}





