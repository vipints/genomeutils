#!/bin/bash
#
# ORG_NAMES text file contains the list of organism names included in the experiment  
#
for org_name in `cat ORG_NAMES`
do 
    ## expecting the read alignments in read_mapping sub folder 
    cd $org_name/read_mapping/
    ## submitting the job to the cluster computing resource 
    echo "samtools view -h Aligned.sort.out.bam | grep -e "^@" -e "NH:i:1[[:space:]]" | samtools view -bS - > unique.bam; samtools index unique.bam" | qsub -l nodes=1:ppn=1 -l mem=4gb,vmem=4gb,pmem=4gb,walltime=05:00:00 -d `pwd` -N ${org_name}_UNQ
    cd ../..
done 
