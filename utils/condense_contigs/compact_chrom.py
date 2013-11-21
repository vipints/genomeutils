#!/usr/bin/env python
"""
Genomes that come in many small contigs (> several hundreds) this program merge 
into few larger ones.

Usage: python compact_chrom.py gio_file_path no_chrom merged_gio_map_file merged_gio_file_path
"""

import sys, os, re, time
from operator import itemgetter

def write_map_file(new_path, flat_seq, chr_id, cname, map, chr_no):
    """
    Coordinate Mapping file from OLD to NEW Contigs 
    """
    f_flat = open(new_path + '/genome/Contig' + str(chr_id[0]+1) + '.flat', 'w') ## create new contig file in .flat format in target directory
    f_flat.write(flat_seq)
    f_flat.close()

    i = 0
    tc = len(cname)
    for ele in cname: ## writing a mapping file old contig to new contig information.
        if i == 0:
            start = 1
            stop = map[i]
        else:
            start = p_stop
            stop = start + map[i] - 1
        print 'Contig' + str(chr_id[0]+1) + '\t' + ele + '\t' + str(start) + '\t' + str(stop) + '\t' + str(map[i])
        if i==(tc-1): break
        print 'Contig' + str(chr_id[0]+1) + '\t' + 'NSPACER\t' + str(stop+1) + '\t' + str(stop+25000) + '\t' + str(25000) 
        p_stop = stop + 25001 # default spacer nts 
        i += 1    
    
    (flat_seq, cname, map) = ('', [], [])
    chr_no.append(chr_id[0]+1)
    chr_id = chr_id[1:]

    return (flat_seq, cname, map, chr_no, chr_id)

def __main__():
    
    stime = time.asctime(time.localtime(time.time()))
    try:
        old_path = sys.argv[1] ## GIO file path with small several hundred contigs
        no_chr = sys.argv[2] ## Target number of new chromosomes
        giolog = open(sys.argv[3], 'w') ## New GIO file
        new_path = sys.argv[4] ## New GIO path
    except:
        print __doc__
        sys.exit(-1)

    giolog.write('------------------------------------------------\n')
    giolog.write('Merge Contigs started on ' + stime + '\n')
    giolog.write('------------------------------------------------\n\n')
    
    ## calculate the total sequence length from all contigs
    (seq_len, old_c, contigs) = (0, 0, {})
    contig_files = os.listdir(old_path + '/genome/')
    for file in contig_files:
        old_c += 1 ## getting number of contigs
        fsize = os.path.getsize(old_path + '/genome/' + file) ## getting file size
        contigs[file] = fsize
        inh = open(old_path + '/genome/' + file, "rU")
        for line in inh:
            line = line.strip()
            seq_len += len(line)
        inh.close()
    if int(no_chr) >= old_c:
        sys.stderr.write("Resulting merge contigs are greater than or equal to original contigs\n")
        sys.exit(-1)

    ## calculate sequence quantity for each new contigs
    chr_quan = seq_len/int(no_chr)
    os.system('mkdir -p ' + new_path + '/genome/')
    chr_id = [x for x in range(int(no_chr))]

    print '##Contig merging information file'
    print '##New Seq ID\tMap element\tStart\tStop\tLength'

    sorted_contigs = sorted(contigs.iteritems(), key=itemgetter(1), reverse=True) ## sort contigs based on size in reverse order.
    (flat_seq, sl, cname, map, chr_no, seq, capture) = ('', 0, [], [], [], '', 0)

    for file in sorted_contigs:
        cn = re.search(r'(.+)\.flat', file[0]).group(1)
        cname.append(cn)
 
        seq = ''
        inh = open(old_path + '/genome/' + file[0], "rU")
        for line in inh:
            line = line.strip()
            seq += line 
        inh.close()
        map.append(len(seq))  
        sl += len(seq)

        if sl < chr_quan:
            n_spacer = 'N'*25000
            flat_seq += seq + n_spacer
            capture = 0
            continue

        flat_seq += seq
        sl = 0
        (flat_seq, cname, map, chr_no, chr_id) = write_map_file(new_path, flat_seq, chr_id, cname, map, chr_no) ## writing new GIO file and corresponding mapping file.
        capture = 1
       
    if capture == 0:
        flat_seq += seq
        seq = ''
        (flat_seq, cname, map, chr_no, chr_id) = write_map_file(new_path, flat_seq, chr_id, cname, map, chr_no) ## final catching.

    ## GIO file
    giolog.write('Genome properties:\n')
    giolog.write(' * ' + str(len(chr_no)) + ' contigs\n')
    giolog.write(' * ' + str(seq_len/1000) + 'kbs total length\n\n') 
    giolog.write('Contig list:\n')
    for chr in chr_no:
        giolog.write('    Contig' + str(chr) + '\n')
    
    gconfig = open(new_path + '/genome.config', 'w')
    gconfig.write('BASEDIR ' + new_path + '\n\n')
    gconfig.write('CONTIGS ' + str(len(chr_no)) + '\n')
    for chr in chr_no:
        gconfig.write('Contig' + str(chr) + '\t' + 'genome/Contig' + str(chr) + '.flat\t' + 'genome/Contig' + str(chr) + '.dna\n')
    gconfig.write('\n')
    gconfig.write('ALPHABET acgt\n\n')
    gconfig.write('ESTFILES 0\n\n')
    gconfig.write('CDNAFILES 0\n\n')
    gconfig.write('ANNOTATIONFILES 0\n')
    gconfig.close()
    etime = time.asctime( time.localtime(time.time()) )
    giolog.write('\n-------------------------------------------------\n')
    giolog.write('Merge Contigs finished on ' + etime + '\n')
    giolog.write('-------------------------------------------------\n')	

if __name__=="__main__":__main__()
