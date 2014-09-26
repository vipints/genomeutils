#!/bin/bash

for org_name in H_sapiens P_troglodytes M_mulatta R_norvegicus M_musculus O_cuniculus E_caballus B_taurus C_familiaris S_scrofa G_gallus M_gallopavo A_carolinensis T_nigroviridis O_latipes D_rerio X_tropicalis M_domestica A_gambiae A_mellifera D_melanogaster D_simulans C_elegans C_briggsae O_sativa Z_mays G_max B_rapa A_thaliana V_vinifera
do 
    ## expecting the read alignments in read_mapping sub folder 
    cd $org_name/read_mapping/
    ## submitting the job to the cluster computing resource 
    echo "samtools view -h Aligned.sort.out.bam | grep -e "^@" -e "NH:i:1[[:space:]]" | samtools view -bS - > unique.bam; samtools index unique.bam" | qsub -l nodes=1:ppn=1 -l mem=4gb,vmem=4gb,pmem=4gb,walltime=05:00:00 -d `pwd` -N ${org_name}_UNQ
    cd ../..
done 
