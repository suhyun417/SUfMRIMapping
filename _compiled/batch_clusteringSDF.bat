#!/bin/bash
# This file is batch_clusteringSDF.bat.
# to submit the compiled executable shell script as a batch job 
# on the Biowulf computational nodes
 
cd /data/parks20/analysis/NeuroMRI/_compiled
./run_doClusteringSDF_multipleSubj_prob.sh /usr/local/matlab-compiler/v91 1 1
