#!/usr/bin/env python 
"""
Program to filter BLAT alignment transcript file 

Usage : GFF file created by psl_to_gff.pl program 
Here the upper bound mismatch as 4 and minimum matching length of nucleotide as 25 

based on the criteria, program filters the fragments and writes a GFF file 
"""
import sys 

def gff_clean():
    """
    clean the gff file based on the filter
    """

    NM=0
    seqLen=25
    featdb=dict()

    gfh=open(gffname, 'rU')
    for rec in gfh:
        rec = rec.strip('\n\r').split('\t')
        if rec[2]=='transcript':
            if not int(rec[7])<=NM:
                continue
            mLen,fid=0,None
            for atb in rec[-1].split(';'):
                key,val=atb.split('=')
                if key=='ID':
                    fid=val
                if key=='Match':
                    mLen=int(val)
            if mLen >= seqLen: 
                featdb[fid]=1
        else:
            parent,codLen=None,0
            for atb in rec[-1].split(';'):
                key,val=atb.split('=')
                if key=='Parent':
                    parent=val
                if key=='Len':
                    codLen=int(val)
            if parent in featdb:
                if codLen >= seqLen:
                    print '\t'.join(rec)
            #break
    gfh.close()

if __name__=="__main__":

    try:
        gffname=sys.argv[1]
    except:
        print __doc__
        sys.exit(-1)
    
    gff_clean(gffname)
