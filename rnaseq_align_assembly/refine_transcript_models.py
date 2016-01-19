#!/usr/bin/env python 
"""
filter out gene models from an genome annotation file based on 
the splice-site sequence consensus, length of the ORF and 
sequencing read coverage to the trancript. 

program requires genome annotation in gtf/gff and genome sequence 
file in fasta format. 

usage:
    python refine_transcript_models.py in.gff in.fasta  

Requirement:
    numpy       :- http://numpy.org 
    biopython   :- http://biopython.org
    gfftools    :- http://github.com/vipints/gfftools 
"""

from __future__ import division
import sys 
import numpy 
from Bio import SeqIO 
from collections import defaultdict
from gfftools import GFFParser, helper 


def filter_gene_models(gff_name, fas_file, outFile):
    """
    check the sequence consistency/quality of predicted fragment

    @args gff_name: result file gff format from TranscriptSkimmer
    @type gff_name: str
    @args fas_file: genome sequence in fasta format
    @type fas_file: str 
    @args outFile: filtered gene output file 
    @type outFile: str 
    """

    sys.stdout.write( 'using genome sequence file %s\n' % fas_file)
    sys.stdout.write( 'using genome annotation file %s\n' % gff_name)

    sys.stdout.write( "parsing genome annotation file...\n")
    gff_content = GFFParser.Parse(gff_name) ## getting the genome annotation from GFF file 
    sys.stdout.write( " ...done\n") 

    sys.stdout.write( "screening for spliced transcripts...\n")
    orf_short = 0 
    spliced_cand = 0 
    sing_exon_gen = 0
    transcript_cov = 0 
    min_orf_length = 400 

    transcripts_region = defaultdict(list)
    for gene_recd in gff_content: ## screening the spliced transcripts
        spliced_transcript = defaultdict(list)

        for idx, sub_rec in enumerate(gene_recd['transcripts']):
            try:
                exon_cnt = len(gene_recd['exons'][idx])
            except:
                continue

            if exon_cnt > 1: ## skipping the single-exon transcripts 
                if gene_recd['transcript_info'][idx]: ## discarding the transcript based on the read coverage value
                    if float(numpy.atleast_1d(gene_recd['transcript_info'][idx])[0]) < 10: ## read coverage value to consider  
                        transcript_cov += 1 
                        continue

                orf_length = 0 
                for idk, ex in enumerate(gene_recd['exons'][idx]):
                    orf_length += ex[1]-(ex[0]-1)

                    if idk == 0:
                        #ex[0] = None 
                        spliced_transcript[(gene_recd['name'], sub_rec[0], gene_recd['strand'])].append((None, ex[1]))
                    elif exon_cnt-1 == idk:
                        #ex[1] = None
                        spliced_transcript[(gene_recd['name'], sub_rec[0], gene_recd['strand'])].append((ex[0], None))
                    else:
                        spliced_transcript[(gene_recd['name'], sub_rec[0], gene_recd['strand'])].append((ex[0], ex[1]))

                if orf_length < min_orf_length: ## min orf length for the transcripts 
                    
                    del spliced_transcript[(gene_recd['name'], sub_rec[0], gene_recd['strand'])] ## clearing that transcript details
                    orf_short += 1 
                    continue
                    
                spliced_cand +=1
            else:
                sing_exon_gen +=1 
                #TODO orf length of the single exon gene will be good 
                # to look, some histone genes are long enough to have 
                # strong TSS region 
                """
                single_exon_len = 0 
                for idk, ex in enumerate(gene_recd['exons'][idx]):
                    single_exon_len = ex[1]-(ex[0]-1)

                if single_exon_len > 1600:
                    spliced_transcript[(gene_recd['name'], sub_rec[0], gene_recd['strand'])].append((ex[0], ex[1]))
                """

        if spliced_transcript: 
            transcripts_region[gene_recd['chr']].append(spliced_transcript)
    
    sys.stdout.write( "...considering %d spliced transcripts\n" % spliced_cand) 
    sys.stdout.write( "discarding transcripts...\n\t%d transcripts with single exon\n" % sing_exon_gen)
    sys.stdout.write( "\t%d transcripts with read coverage value less than 10\n" % transcript_cov)
    sys.stdout.write( "\t%d transcripts with orf region less than 400 nucleotides\n" % orf_short)

    genemodels = check_splice_site_consensus(fas_file, transcripts_region)

    write_filter_gene_models(gff_content, genemodels, outFile)


def check_splice_site_consensus(fas_file, splice_region):
    """
    splice site consensus check
    """

    sys.stdout.write( "splice site sequence consensus check started...\n")
    get_gene_models = defaultdict()
    splice_site_con = 0 
    fas_fh = helper.open_file(fas_file)
    for fas_rec in SeqIO.parse(fas_fh, "fasta"):
        if fas_rec.id in splice_region:
            for details in splice_region[fas_rec.id]:
                for genes, regions in details.items():
                    acc_cons_cnt = 0 
                    don_cons_cnt = 0 

                    if len(regions) == 1:## single exon long transcripts no checks 
                        get_gene_models[(fas_rec.id, genes[0], genes[1], genes[2])] = 1   
                        continue

                    for region in regions:
                        if genes[-1] == '+':
                            #if not numpy.isnan(region[0]):## acceptor splice site 
                            if region[0]:## acceptor splice site 
                                acc_seq = fas_rec.seq[int(region[0])-3:int(region[0])-1]
                                if str(acc_seq).upper() == "AG":
                                    acc_cons_cnt += 1 

                            if region[1]:
                                don_seq = fas_rec.seq[int(region[1]):int(region[1])+2]
                                if str(don_seq).upper() == "GT":
                                    don_cons_cnt +=1 

                        elif genes[-1] == '-':
                            if region[0]: ## donor splice site 
                                don_seq = fas_rec.seq[int(region[0])-3:int(region[0])-1]
                                don_seq = don_seq.reverse_complement()
                                if str(don_seq).upper() == "GT":
                                    don_cons_cnt +=1 
                            
                            if region[1]:
                                acc_seq = fas_rec.seq[int(region[1]):int(region[1])+2]
                                acc_seq = acc_seq.reverse_complement()
                                if str(acc_seq).upper() == "AG":
                                    acc_cons_cnt += 1 
                    ## check for half of the consensus sites 
                    if acc_cons_cnt > (len(regions)/2) and don_cons_cnt > (len(regions)/2):
                        get_gene_models[(fas_rec.id, genes[0], genes[1], genes[2])] = 1   
                    else:
                        splice_site_con +=1 
    fas_fh.close()
    sys.stdout.write( "...considering %d best transcripts\n" % len(get_gene_models))
    sys.stdout.write( "discarding transcripts...\n")
    sys.stdout.write( "\t%d splice-site consensus sequence missing\n" % splice_site_con)

    return get_gene_models


def write_filter_gene_models(gff_cont, gene_models, outFile):
    """
    writing the filtered gene models to the result file
    """

    sys.stdout.write( "writing filtered gene models to %s ...\n" % outFile)
    true_genes = 0 
    true_transcripts = 0 
    out_fh = open(outFile, "w")
    for recd in gff_cont:
        trans_indices = [] 

        for idx, sub_rec in enumerate(recd['transcripts']):
            if (recd['chr'], recd['name'], sub_rec[0], recd['strand']) in gene_models:
                trans_indices.append(idx)

        if trans_indices:
            true_genes += 1 
            chr_name = recd['chr']
            strand = recd['strand']
            start = recd['start']
            stop = recd['stop']
            source = recd['source']
            ID = recd['name']
            Name = recd['gene_info']['Name']
            Name = ID if Name != None else Name  
            out_fh.write('%s\t%s\tgene\t%d\t%d\t.\t%s\t.\tID=%s;Name=%s\n' \
                % (chr_name, source, start, stop, strand, ID, Name))
                
            for idz, tid in enumerate(recd['transcripts']):
                if idz in trans_indices:

                    true_transcripts += 1 
                    t_start = recd['exons'][idz][0][0]
                    t_stop = recd['exons'][idz][-1][-1]
                    t_type = recd['transcript_type'][idz] 

                    out_fh.write('%s\t%s\t%s\t%d\t%d\t.\t%s\t.\tID=%s;Parent=%s\n' \
                        % (chr_name, source, t_type, t_start, t_stop, strand, tid[0], ID))
                    
                    for ex_cod in recd['utr5_exons'][idz]:
                        out_fh.write('%s\t%s\tfive_prime_UTR\t%d\t%d\t.\t%s\t.\tParent=%s\n' \
                            % (chr_name, source, ex_cod[0], ex_cod[1], strand, tid[0])) 
                    for ex_cod in recd['cds_exons'][idz]:
                        out_fh.write('%s\t%s\tCDS\t%d\t%d\t.\t%s\t%d\tParent=%s\n' \
                            % (chr_name, source, ex_cod[0], ex_cod[1], strand, ex_cod[2], tid[0])) 
                    for ex_cod in recd['utr3_exons'][idz]:
                        out_fh.write('%s\t%s\tthree_prime_UTR\t%d\t%d\t.\t%s\t.\tParent=%s\n' \
                            % (chr_name, source, ex_cod[0], ex_cod[1], strand, tid[0]))
                    for ex_cod in recd['exons'][idz]:
                        out_fh.write('%s\t%s\texon\t%d\t%d\t.\t%s\t.\tParent=%s\n' \
                            % (chr_name, source, ex_cod[0], ex_cod[1], strand, tid[0])) 
    out_fh.close()
    sys.stdout.write( "...done\n")
    sys.stdout.write( "number of genes considered %d\n" % true_genes)
    sys.stdout.write( "number of transcripts considered %d\n" % true_transcripts)


if __name__ == "__main__":
    try:
        gff_name = sys.argv[1]
        fas_file = sys.argv[2] 
    except:
        print __doc__
        sys.exit(-1) 

    outfile = "filter_genes.gff"
    filter_gene_models(gff_name, fas_file, outfile)
