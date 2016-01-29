#from signal_labels import fetch_signal_seq_examples as fetch_seq 
#from data_quality_check import split_train_test_examples

#fetch_seq.get_feature_seq("H_sapiens/ensembl_release_79/STARgenome/hg19_chrOnly.fa", "../H_sapiens_stringtie_genes.gff")
#split_train_test_examples.split_data_random_non_overlap() 

#split_train_test_examples.get_matching_neg_example("tss_sig_pos_example_1.fa", "tss_sig_neg_example.fa", "tss_sig_neg_example_1.fa") 
#split_train_test_examples.get_matching_neg_example("tss_sig_pos_example_2.fa", "tss_sig_neg_example.fa", "tss_sig_neg_example_2.fa") 
#split_train_test_examples.get_matching_neg_example("tss_sig_pos_example_3.fa", "tss_sig_neg_example.fa", "tss_sig_neg_example_3.fa") 


import os 
import sys 
import yaml

from collections import defaultdict
import libpyjobrunner as pg


def distribute_model_train(args_list):
    """
    wrapper function to distribure model training for number of organisms 
    """
   
    import train_signal_model as tsm

    organism, work_path, data_path = args_list
    
    os.chdir(work_path) 
    
    tsm.train_wdspeck_svm(organism, "tss", data_path) #FIXME the signal type 
    return "done"


def main(yaml_config):
    """
    """

    config_map = yaml.safe_load(open(yaml_config, "rU"))
    exp_path_pfx = config_map['experiment_data_path']['dir']
    
    org_db = defaultdict()
    for ent in config_map['experiment']:
        species_name = ent['organism_name']

        genus, species = species_name.strip().split("_")
        short_name = "%s_%s" % (genus[0].upper(), species)

        org_db[short_name] = dict(short_name = short_name)  
        org_db[short_name]['work_dir'] = "%s/%s/set_union_refix" % (exp_path_pfx, short_name) 
        org_db[short_name]['data_dir'] = "%s/%s/set_2" % (exp_path_pfx, short_name) 
    
    ## prepare jobs 
    Jobs = [] 
    for org_name, det in org_db.items():

        arg = [[org_name, det['work_dir'], det['data_dir']]]
        
        job = pg.cBioJob(distribute_model_train, arg) 
        job.mem="4gb"
        job.vmem="4gb"
        job.pmem="4gb"
        job.pvmem="4gb"
        job.nodes = 1
        job.ppn = 1
        job.walltime = "2:00:00"

        Jobs.append(job)
    
    compute_local = True
    print "sending jobs to worker" 
    processedJobs = pg.process_jobs(Jobs, local=compute_local) 


if __name__ == "__main__":
    
    ## run for multiple organisms 

    ## read the config file 
    try:
        yaml_config = sys.argv[1]
    except:
        print __doc__

    main(yaml_config)
