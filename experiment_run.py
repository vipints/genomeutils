#!/usr/bin/env python 
"""
master script to execute the pipeline
"""

import os 
import sys 


if __name__=="__main":
    
    infile = ""
    data_path = ""
    exp_path = ""
   
    from signal_labels import org_details_db as odb

    org_details = odb.make_org_db(infile, data_path, exp_path) 

    ## iterating through each organism 
