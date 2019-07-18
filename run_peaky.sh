#!/bin/bash

module load R/3.4.1

while read -r LIB ; do

  Rscript /working/lab_georgiat/jonathB/CaptureHiC/peaky/peaky_fs/scripts/peaky.hpc.r ${LIB}

done < /working/lab_georgiat/jonathB/CaptureHiC/peaky/peaky_fs/library.list.txt
