#!/usr/bin/env python
"""
Merge contig locations in a BED file to a bin of YY nts location on a single chromosome. 
The score of that bin will be the maximum score obtained in the 500 nts bin.

Usage:
cat in.bed | python CondenseContigRegions.py bin_size > bin-out.bed 

in.bed looks like:
track name=junctions description="TSS score"
chrC    70043   70044   44      -
chrC    105671  105913  1       +
chrC    132687  132980  2       -
chrC    132735  132977  1       -
chrM    11450   11722   2       -
chr2    6149    8097    2       +
chr2    6925    8079    3       +
chr2    77246   77428   1       +
chr2    87142   87322   1       +

"""
import sys 

if __name__=="__main__":

    try:
        binsize = int(sys.argv[1])
    except:
        print 'Incorrect argument supplied'
        print __doc__
        sys.exit()

    win_score = []
    wind_cnt, condenseStart, LastPos = 0, 0, 0
    chr_change = None 

    for line in sys.stdin
        line = line.strip('\n\r').split('\t')
        line[3] = round(float(line[3]), 4)
        win_score.append(line[3])

        if wind_cnt == (binsize-1): 
            win_score.sort()
            dense_line = [line[0],
                    str(condenseStart),
                    str(condenseStart+1),
                    str(win_score[-1]),
                    line[4]
                    ]
            print '\t'.join(dense_line)
            condenseStart += binsize
            wind_cnt = 0
            win_score = []
            continue
        wind_cnt += 1
        LastPos = int(line[1])

        if chr_change != line[0]:
            print_last_bin()
            chr_change=None
            continue
        chr_change=line[0]   

def print_last_bin():
    """
    print the final bin to the merged list in the file
    """
    if wind_cnt != 0:
        win_score.sort()
        dense_line = [line[0],
                    str(LastPos),
                    str(LastPos+1),
                    str(win_score[-1]),
                    line[4]
                    ]
        print '\t'.join(dense_line)
