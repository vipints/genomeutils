#!/usr/bin/env python
"""
Transform contig and feature co-ordinate based on the mapping file. 

Usage: python swap_contigs.py map_file switch_type 

map_file can be created from compact_chrom.py
switch_type annoGFF, predGFF, predSPF, asciiSPF, wigSPF


"""

import re, sys, os, time
import struct

def annotationInGFF3(occ, c_info):
    """ 
    Transform contigs in genome annotation (in GFF3) based on a mapping file created by MergeContigs tool
    """
    try:
        annof = open(sys.argv[3], "rU")
        info = open(sys.argv[4], "w")
    except:
        sys.stderr.write("Check input parameters\n\tGenome Annotation in GFF3 file\n\tResult file name in GFF3 format\n")
        sys.exit(-1)
   
    ## checking consistency of contig numbers in both files.
    gff_contig = os.popen('cut -f 1 ' + sys.argv[3] + ' | sort | uniq').read()
    gff_contig = gff_contig.strip()
    
    num = 0
    for chr in gff_contig.split('\n'):
        if re.match(r'#', chr): continue
        num += 1
    
    if len(c_info) != num:
        sys.stdout.write("Contig numbers in Mapping file: " + str(len(c_info)) + "\tand GFF3 file: " + str(num) + "\n")
        sys.stdout.write("Warning: Found mismatch in contig numbers from provided Mapping and GFF3 file\n")

    print 
    print '2. GFF file validation over.'
    print

    for l in annof:
        l = l.strip()
        if re.match(r'#', l): info.write(l + "\n");continue
        
        l = l.split('\t')
        if l[2] == 'region': continue

        ## validating each line in the file.
        if len(l) != 9:
            sys.stdout.write("Skipping invalid GFF line: ")
            l = '\t'.join(l)
            sys.stdout.write(l + "\n")
            continue
        
        ## checking merged contig contain only one chromosome.
        if occ[c_info[l[0]]['lmark']] == 1: 
            l[0] = c_info[l[0]]['lmark']
            l = '\t'.join(l)
            info.write(l + "\n")
            continue
        
        ## mapping the old contig postions to new contig positions.
        f_len = int(l[4]) - int(l[3])
        l[3] = str(c_info[l[0]]['start'] + int(l[3]) - 1) 
        l[4] = str(int(l[3]) + f_len)
        l[0] = c_info[l[0]]['lmark']
        l = '\t'.join(l)
        info.write(l + "\n")
    
    annof.close()
    info.close()

    print
    print '3. Transformed all feature co-ordinate based on Contig Mapping file.'
    print

    print
    print '4. Created new Annotation file in GFF format.'
    print

def contigPOS(file):
    """ 
    Reading position information from SPF binary file 
    """
    loc = []
    try:
        inh = open(sys.argv[3] + '/pred/' + file, "rb")
    except:
        sys.stderr.write("SPF Co-ordinate file not found. Program terminating.\n")
        sys.exit(-1)
    while True:
        rec = inh.read(4)
        if len(rec) != 4:break
        (p,) = struct.unpack('i', rec)
        loc.append(p)
    inh.close()
    return loc

def posSCORE(file):
    """ 
    Reading score values from SPF binary format 
    """
    score = []
    sfile = re.search(r'(.+)\.pos', file).group(1)
    try:
        inh = open(sys.argv[3] + '/pred/' + sfile + '.Conf_cum', "rb")
    except:
        sys.stderr.write("SPF Score file not found, Program terminating.\n")
        sys.exit(-1)
    while True:
        rec = inh.read(4)
        if len(rec) != 4: break
        (s,) = struct.unpack('f', rec)
        s = '%.3f' % s
        score.append(s)
    inh.close()
    return score 

def posBINformat(orient, contig_pos):
    """ 
    Writing Position information in SPF binary format 
    """
    for chr in contig_pos:
        num = re.search(r'(\d+)$', chr).group(1)
        inh = open(sys.argv[4] + '/pred/contig_' + str(num) + orient + '.pos', "wb")
        for l in contig_pos[chr]: inh.write(struct.pack("i", l))
        inh.close()
    contig_pos.clear()        

def scoreBINformat(orient, contig_score):
    """ 
    Writing score information in SPF binary format 
    """
    for chr in contig_score:
        num = re.search(r'(\d+)$', chr).group(1)
        inh = open(sys.argv[4] + '/pred/contig_' + str(num) + orient + '.Conf_cum', "wb")
        for l in contig_score[chr]: inh.write(struct.pack("f", float(l)))
        inh.close()
    contig_score.clear()        
        
def predictionInSPF(c_info, occp_con):
    """
    Transform contigs in Signal Prediction Format file based on mapping file created by MergeContigs tool
    """
    try:
        files = os.listdir(sys.argv[3] + '/pred')
        os.system('mkdir -p ' + sys.argv[4] + '/pred/')
    except:
        sys.stderr.write("Check input parameters\n\tSPF prediction file\n\tPath name to store transformed SPF files\n")
        sys.exit(-1)

    if len(files) == 0:
        sys.stderr.write("Empty SPF file. Cannot continue, program terminating....\n")
        sys.exit(-1)

    print 
    print '2. SPF file validation over.'
    print

    pos_file = []
    for file in files:
        if re.search(r'.+\.pos$', file):pos_file.append(file)

    for file in pos_file:
        ## reading position and score value 
        loc = contigPOS(file)
        score = posSCORE(file)

        orient = re.search(r'.*(.+)\.pos', file).group(1)
        conid = 'Contig' + re.search(r'.*\_(.+)[+|-]\.pos', file).group(1)
        j = 0
        (contig_pos, contig_score) = ({}, {})

        ## remapping the new contig to old genome contigs.
        for pos in loc: 
            for cid in occp_con[conid]:
                if pos >= c_info[cid]['start'] and pos <= c_info[cid]['stop']:
                    point = pos - (c_info[cid]['start'] -1)
                    if cid in contig_pos:                    
                        contig_pos[cid].append(point)
                    else:
                        contig_pos[cid] = [point]
                    if cid in contig_score:
                        contig_score[cid].append(score[j])
                    else:
                        contig_score[cid] = [score[j]]
                    break
            j += 1

        ## writing new SPF binary files for position and score value.
        posBINformat(orient, contig_pos) 
        scoreBINformat(orient, contig_score)

    print 
    print '3. Created new Signal Prediction file in SPF format.'
    print

    try:
        trunkinh = open(sys.argv[5], 'rU')
        trunkouth = open(sys.argv[6], 'w')
    except:
        sys.stderr.write("Truncated SPF file missing\n")
        sys.exit(-1)

    for l in trunkinh:
        l = l.strip()
        if re.match(r'#', l): trunkouth.write(l + "\n");continue
        l = l.split('\t')

        if len(l) != 6:
            sys.stdout.write("Skipping invalid ASCII-SPF file line: ")
            l = '\t'.join(l)
            sys.stdout.write(l + "\n")
            continue

        ## Merged contig contains only one chromosome
        if len(occp_con[l[0]]) == 1: 
            l[0] = occp_con[l[0]][0]
            l = '\t'.join(l)
            trunkouth.write(l + "\n")
            continue

        ## remaining contigs appear here and mapping back the coordinates.
        for chr in occp_con[l[0]]:
            if int(l[3]) >= c_info[chr]['start'] and int(l[3]) <= c_info[chr]['stop']:
                pos = int(l[3]) - (c_info[chr]['start'] - 1)
                l[3] = str(pos)
                l[0] = chr
                line = '\t'.join(l)
                trunkouth.write(line + "\n")
                break
            
    trunkinh.close()
    trunkouth.close()

def predictionInASCIISPF(c_info, occp_con):
    """
    Transform contigs in ASCIISPF prediction file based on Mapping file created by MergeContigs tool
    """
    try:
        spfh = open(sys.argv[3], "rU")
        aspfh = open(sys.argv[4], "w")
    except:
        sys.stderr.write("Check input parameters\n\tASCII-SPF prediction file\n\tResult file name\n")
        sys.exit(-1)
    
    for l in spfh:
        l = l.strip()
        if re.match(r'#', l): aspfh.write(l + "\n");continue
        l = l.split('\t')

        if len(l) != 6:
            sys.stdout.write("Skipping invalid ASCII-SPF file line: ")
            l = '\t'.join(l)
            sys.stdout.write(l + "\n")
            continue

        ## Merged contig contains only one chromosome
        if len(occp_con[l[0]]) == 1: 
            l[0] = occp_con[l[0]][0]
            l = '\t'.join(l)
            aspfh.write(l + "\n")
            continue

        ## remaining contigs appear here and mapping back the coordinates.
        for chr in occp_con[l[0]]:
            if int(l[3]) >= c_info[chr]['start'] and int(l[3]) <= c_info[chr]['stop']:
                pos = int(l[3]) - (c_info[chr]['start'] - 1)
                l[3] = str(pos)
                l[0] = chr
                line = '\t'.join(l)
                aspfh.write(line + "\n")
                break

    spfh.close()
    aspfh.close()
    
    print
    print '2. Transformed all feature co-ordinate based on Contig Mapping file.'
    print 

    print 
    print '3. Created new Signal Prediction file in ASCII-SPF.'
    print 


def predictionInWigSPF(c_info, occp_con):
    """
    Transform contigs in Wiggle SPF prediction file based on mapping file created by MergeContigs tool
    """
    try:
        wigh = open(sys.argv[3], "rU")
        wspfh = open(sys.argv[4], "w")
    except:
        sys.stderr.write("Check input parameters\n\tWIG-SPF prediction file\n\tResult file name.\n")
        sys.exit(-1)

    wigpos = dict()
    (track_name, pchr) = ('', '')
    
    for l in wigh:
        l = l.strip()
        if re.match(r'#', l): continue
        if re.search(r'^track', l): track_name = l;continue
        if re.search(r'^variableStep', l): 
            chr = re.search(r'chrom=(.+)\s.+', l).group(1)
            if chr == pchr:
                orient = '-'
            else:
                orient = '+'
            pchr = chr    
            continue
        if len(occp_con[chr]) == 1:
            ke = occp_con[chr][0] + orient
            if ke in wigpos:
                wigpos[ke].append(l)
            else:
                wigpos[ke] = [l]
            continue
        l = l.split(' ')

        if len(l) != 2:
            sys.stdout.write("Skipping invalid Wiggle SPF file line: ")
            l = ' '.join(l)
            sys.stdout.write(l + "\n")
            sys.exit(-1)
        
        for cid in occp_con[chr]:
            if int(l[0]) >= c_info[cid]['start'] and int(l[0]) <= c_info[cid]['stop']:
                pos = int(l[0]) - (c_info[cid]['start'] - 1)
                l[0] = str(pos)
                ke = cid + orient
                l = ' '.join(l)
                if ke in wigpos:
                    wigpos[ke].append(l)
                else:
                    wigpos[ke] = [l]
                break                    
    
    print
    print '2. Transformed all feature co-ordinate based on Contig Mapping file.'
    print

    track_name = re.search(r'(.*)\s[+|-]$', track_name).group(1)
    wspfh.write(track_name + "\n")
    for contig in wigpos:
        wspfh.write("variableStep chrom=" + contig + " span=1\n")
        for line in wigpos[contig]:wspfh.write(line + "\n")
    wigpos.clear()
    wspfh.close()
    wigh.close()        

    print 
    print '3. Created new Signal Prediction file in ASCII-SPF.'
    print 
        
def predictionInGFF(c_info, occp_con):
    """
    Transform contig location of a GFF3 prediction file based on the mapping file created using MergeContig tool 
    """
    try:
        pred_gff_hand = open(sys.argv[3], "rU") 
        out_gff_hand = open(sys.argv[4], "w")
    except:
        sys.stderr.write("Check Input parameters:\n\tPrediction file in GFF3 format\n\tResult file name\n")
        sys.exit(-1)

    for l in pred_gff_hand:
        l = l.strip()
        if re.match(r'#', l): out_gff_hand.write(l + "\n"); continue
        
        l = l.split('\t')
        if len(l) != 9:
            sys.stdout.write("Skipping invalid GFF format line: ")
            l = '\t'.join(l)
            sys.stdout.write(l + "\n")
            continue 
        
        ## checking merged contig contains only one chromosome. 
        if len(occp_con[l[0]]) == 1:
            l[0] = occp_con[l[0]][0]
            l = '\t'.join(l)
            out_gff_hand.write(l + "\n")
            continue

        ## several small contigs merged on one big contig will appear here for transforming their coordinates.
        for chr in occp_con[l[0]]:
            if int(l[3]) >= c_info[chr]['start'] and int(l[3]) <= c_info[chr]['stop']:
                flen = int(l[4]) - int(l[3])
                start = int(l[3]) - (c_info[chr]['start'] - 1)
                stop = start + flen
                l[0] = chr
                l[3] = str(start)
                l[4] = str(stop)
                l = '\t'.join(l)
                out_gff_hand.write(l + "\n")
                break

    pred_gff_hand.close()
    out_gff_hand.close()

    print 
    print '2. Transformed all feature co-ordinate based on Contig Mapping file.' 
    print 

    print 
    print '3. Created new Prediction file in GFF format.'
    print 

def __main__():
   
    stime = time.asctime(time.localtime(time.time()))
    print '-------------------------------------------------'
    print 'Switch Contigs started on ' + stime 
    print '-------------------------------------------------'
    
    ## retrieving mapping file information.
    try:
        maph = open(sys.argv[1], "rU")
    except:
        print __doc__
        sys.exit(-1)

    occ, c_info, occp_con = dict(), dict(), dict() ## Occurence number, Information, Occupied short contigs
    for l in maph:
        l = l.strip()
        if re.match(r'#', l): continue
        if re.search(r'NSPACER', l): continue
        l = l.split('\t')

        if len(l) != 5:
            l = '\t'.join(l)
            sys.stderr.write("Skipping invalid Contig Mapping line: " + l + "\n")
            continue
            
        info = {}
        info['lmark'] = l[0]
        info['start'] = int(l[2])
        info['stop'] = int(l[3])
        c_info[l[1]] = info
        if l[0] in occ:
            occ[l[0]] += 1
            occp_con[l[0]].append(l[1])
        else:
            occ[l[0]] = 1
            occp_con[l[0]] = [l[1]]
    maph.close()
    ## information retrieved from mapping file.

    print 
    print '1. Parsed information from Contig Mapping file.'
    print 

    try:
        con_type = sys.argv[2]
    except:
        print __doc__
        sys.exit(-1)
    
    if con_type == 'annoGFF':
        print '- Contig transformation type :  Genome annotation in GFF'
        annotationInGFF3(occ, c_info)
    elif con_type == 'predSPF':
        print '- Contig transformation type :  Signal Prediction in SPF'
        predictionInSPF(c_info, occp_con)
    elif con_type == 'asciiSPF':
        print '- Contig transformation type :  Signal Prediction in ASCII-SPF'
        predictionInASCIISPF(c_info, occp_con)
    elif con_type == 'wigSPF':
        print '- Contig transformation type :  Signal Prediction in WIG-SPF'
        predictionInWigSPF(c_info, occp_con)
    elif con_type =='predGFF':
        print '- Contig transformation type :  Gene Prediction in GFF'
        predictionInGFF(c_info, occp_con)
    else:
        sys.stderr.write("Contig transformation type is invalid\n")
        sys.exit(-1)
    
    occ.clear()
    c_info.clear()
    occp_con.clear()
    
    print 
    print 'Done.'
    etime = time.asctime( time.localtime(time.time()) )
    print '--------------------------------------------------'
    print 'Switch Contigs finished on ' + etime
    print '--------------------------------------------------'

if __name__=="__main__":__main__()
