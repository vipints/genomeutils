#!/bin/bash
#
# ORG_NAMES text file contains the list of organism names included in the experiment  
#
for org_name in `cat ORG_NAMES`
do 
    echo $org_name
    cd $org_name/trans_pred
    echo "bash run_TranSkimmer.sh" | qsub -l nodes=1:ppn=1 -l mem=18gb -l vmem=18gb -l pmem=18gb -l walltime=05:00:00 -d `pwd` -N ${org_name}_TSkm 
    cd ../..
done 
