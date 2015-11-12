from signal_labels import fetch_signal_seq_examples as fetch_seq 
from data_quality_check import split_train_test_examples

fetch_seq.get_feature_seq("H_sapiens/ensembl_release_79/STARgenome/hg19_chrOnly.fa", "../H_sapiens_stringtie_genes.gff")
split_train_test_examples.split_data_random_non_overlap() 

split_train_test_examples.get_matching_neg_example("tss_sig_pos_example_1.fa", "tss_sig_neg_example.fa", "tss_sig_neg_example_1.fa") 
split_train_test_examples.get_matching_neg_example("tss_sig_pos_example_2.fa", "tss_sig_neg_example.fa", "tss_sig_neg_example_2.fa") 
split_train_test_examples.get_matching_neg_example("tss_sig_pos_example_3.fa", "tss_sig_neg_example.fa", "tss_sig_neg_example_3.fa") 
