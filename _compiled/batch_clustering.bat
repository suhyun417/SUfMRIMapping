#!/bin/bash
# This file is batch_Clustering.bat.
# to submit the compiled executable shell script as a batch job 
# on the Biowulf computational nodes
 
cd /data/parks20/analysis/NeuroMRI/_compiled
./run_doClusteringCorrMap_multipleSubj_prob_1000_corrcrit.sh /usr/local/matlab-compiler/v91 1 1
