#!/usr/bin/env python
"""
Extract genome annotation from a GFF (a tab delimited format for storing sequence features and annotations) file.

Usage: GFFParser.py in.gff3  

Requirements: 
    Numpy :- http://numpy.org/ 
    Scipy :- http://scipy.org/ 

Copyright (C)	

2009-2012 Friedrich Miescher Laboratory of the Max Planck Society, Tubingen, Germany. 
2012-2013 Memorial Sloan-Kettering Cancer Center, New York City, USA.
"""

import re
import os
import sys
import urllib
import numpy as np
import scipy.io as sio
from collections import defaultdict
import helper as utils 

def _attribute_tags(col9):
    """ 
    Split the key-value terms from the attribute column
    """
    
    info = defaultdict(list)
    is_gff = False
    
    if not col9:
        return is_gff, info
        
    # trim the line ending semi-colon  ucsc may have some white-space  
    col9 = col9.rstrip(';| ')
    # attributes from 9th column 
    atbs = col9.split(" ; ")
    if len(atbs) == 1:
        atbs = col9.split("; ")
        if len(atbs) == 1:
            atbs = col9.split(";")
    # check the GFF3 pattern which has key value pairs like:
    gff3_pat = re.compile("\w+=")
    # sometime GTF have: gene_id uc002zkg.1;
    gtf_pat = re.compile("\s?\w+\s")

    key_vals = []

    if gff3_pat.match(atbs[0]): # gff3 pattern 
        is_gff = True
        key_vals = [at.split('=') for at in atbs]
    elif gtf_pat.match(atbs[0]): # gtf pattern
        for at in atbs:
            key_vals.append(at.strip().split(" "))
    else:
        # to handle attribute column has only single value 
        key_vals.append(['ID', atbs[0]])
    # get key, val items 
    for item in key_vals:
        key, val = item
        if val[0] == '"' and val[-1] == '"':
            val = val[1:-1] 
        # replace the web formating place holders to plain text format 
        info[key].extend([urllib.unquote(v) for v in val.split(',') if v])

    return is_gff, info
                
def _spec_features_keywd(gff_parts):
    """
    Specify the feature key word according to the GFF specifications
    """

    for t_id in ["transcript_id", "transcriptId", "proteinId"]:
        try:
            gff_parts["info"]["Parent"] = gff_parts["info"][t_id]
            break
        except KeyError:
            pass
    ## TODO key words
    for flat_name in ["Transcript", "CDS"]:
        if gff_parts["info"].has_key(flat_name):
            # parents
            if gff_parts['type'] in [flat_name] or re.search(r'transcript', gff_parts['type'], re.IGNORECASE):
                if not gff_parts['id']:
                    gff_parts['id'] = gff_parts['info'][flat_name][0]
                    #gff_parts["info"]["ID"] = [gff_parts["id"]]
            # children 
            elif gff_parts["type"] in ["intron", "exon", "three_prime_UTR",
                        "coding_exon", "five_prime_UTR", "CDS", "stop_codon",
                        "start_codon"]:
                gff_parts["info"]["Parent"] = gff_parts["info"][flat_name]
            break
    return gff_parts

def GFFParse(ga_file):
    """
    Parsing GFF/GTF file based on feature relationship.
    """

    ga_handle = utils._open_file(ga_file)

    for rec in ga_handle:
        rec = rec.strip('\n\r')

        # skip empty line fasta identifier and commented line
        if not rec or rec[0] in  ['#', '>']:
            continue
        # skip the genome sequence 
        if not re.search('\t', rec):
            continue

        parts = rec.split('\t')
        assert len(parts) >= 8, rec

        # process the attribute column (9th column)
        ftype, tags = _attribute_tags(parts[-1])
        if not tags: # skip the line if no attribute column.
	        continue 
        
        # extract fields  
        gff_info = dict()
        gff_info["is_gff3"] = ftype
        if parts[1]:
            tags["source"].append(parts[1])
        if parts[5]:
            tags["score"].append(parts[5])
        if parts[7]:
            tags["phase"].append(parts[7])
        gff_info['info'] = dict(tags)
        gff_info['chr'] = parts[0]
        if parts[3] and parts[4]:
            gff_info['location'] = [int(parts[3]) ,
                        int(parts[4])]
            gff_info['type'] = parts[2]
            gff_info['id'] = tags.get('ID', [''])[0]
            if parts[6] in ['?', '.']:
                parts[6] = None 
            gff_info['strand'] = parts[6]

            # key word according to the GFF spec.
            if not ftype:
                gff_info = _spec_features_keywd(gff_info)
        
            # link the feature relationships
            if gff_info['info'].has_key('Parent'): 
                for p in gff_info['info']['Parent']:
                    if p == gff_info['id']:
                        gff_info['id'] = ''
                        break
                rec_category = 'child'
            elif gff_info['id']:
                rec_category = 'Parent'
            else:
                rec_category = 'record'

        #TODO infer the parent child relationship
    ga_handle.close()
    
    #return 

def __main__():
    """
    extract genome feature information main factory  
    """
    try:
        gff_file = sys.argv[1]
    except:
        print __doc__
        sys.exit(-1)
    
    GFFParse(gff_file)
    
if __name__=='__main__':
    __main__()
