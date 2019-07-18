#/bin/bash

module load R/3.4.1

Rscript /mnt/lustre/working/lab_georgiat/jonathB/CaptureHiC/peaky/peaky_fs/scripts/peaky_correlation/MPPC_cor.r

### collate resutls

PEAKY_DIR=/working/lab_georgiat/jonathB/CaptureHiC/peaky/peaky_fs/peaky_chain_runs/
COR_DIR=/working/lab_georgiat/jonathB/CaptureHiC/peaky/peaky_fs/peaky_chain_runs/peaky_correlation/

for CAP in PCHIC SCHIC ; do
  
  while read -r LIB ; do

    # copy results
    cp ${PEAKY_DIR}/${CAP}_split_baits/${LIB}/baits/MPPCcor* ${COR_DIR}/${CAP}/${LIB}/
    cp ${PEAKY_DIR}/${CAP}_split_baits/${LIB}/baits/baitlist_MPPCcor_sub_0.75.txt ${COR_DIR}/${CAP}/${LIB}/

  done < /working/lab_georgiat/jonathB/CaptureHiC/peaky/peaky_fs/peaky_chain_runs/sample.list.txt

done

## make single combined table for each capture
for CAP in PCHIC SCHIC ; do
  while read -r LIB ; do
    awk -F'\t' -v LIB=${LIB} 'BEGIN{ print "results_file"",""in_A"",""in_B"",""bait_path"",""MPPC_cor"",""redo"",""LIB" }; NR>1{ print $0","LIB }' ${COR_DIR}/${CAP}/${LIB}/MPPCcor.csv 
  done < /working/lab_georgiat/jonathB/CaptureHiC/peaky/peaky_fs/peaky_chain_runs/sample.list.txt > ${COR_DIR}/${CAP}.combined.csv
done 

# final tidy
for CAP in PCHIC SCHIC ; do
awk 'NR == 1 { print } ; !/results/{ print }' ${COR_DIR}/${CAP}.combined.csv > ${COR_DIR}/xx ; mv ${COR_DIR}/xx ${COR_DIR}/${CAP}.combined.csv
done