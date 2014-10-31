#!/usr/bin/perl
#
# convert a blat output file (tab format) for est alignments
# to a gbrowse gff file and the multiple fasta file with the matching ests

use strict;
use Getopt::Long;

my $usage = "$0 -- convert blat file to gbrowse file\n";
$usage .= "\n";
$usage .= "Usage: $0 blat.psl gbrowse.gff\n";
$usage .= "Options:\n";
$usage .= "    --estnames=file    output file with the names of the ESTs\n";
$usage .= "    --source=name      identifyier in the source column\n";
$usage .= "\n";

my $coloffset=0;
my $source;
my $estfilename;

if ($#ARGV <1) {
    die "Unknown option \n\n$usage";
}

GetOptions('estnames:s'=>\$estfilename,
	   'source:s'=>\$source);

$source = "blat" unless (defined($source));

my $blatfilename = $ARGV[0];
my $hintsfilename = $ARGV[1];


open(BLAT, "<$blatfilename") || die "Couldn't open $blatfilename\n";
open(GFF, ">$hintsfilename") || die "Could not open $hintsfilename";
if (defined $estfilename){
    open(EST, ">$estfilename") || die "Could not open $estfilename";
}
my ($i, $j, $mstart, $mend);
my (@dsshints, @asshints, @exonhints, @exonparthints, @intronhints);
my ($match,$TgapCount,$strand,$qname,$qsize,$blockSizes,$tStarts, $qStarts);
my (@f,@b,@t,@q);
my (@blockbegins, @blockends);
my $numBlocks;

# hint lists are sorted by by increasing begin position
my @hint; # (begin, end, strand, tname, qname)
my $hintref;
my $targetname;
my $skiplines=0;
my %estnames;
my $counter=0;
my $prev_tid;
while (<BLAT>) {
    if (/psLayout/){
	$skiplines=5;
    }
    if ($skiplines>0) {
	$skiplines--;
	next;
    }
    s/#.*//;
    next unless /\S/;

    @f = split /\t/, $_, $coloffset+21;
    if (@f < $coloffset+20) { warn "Not BLAT format"; next } # blat format from the GenomeBrowser has an additional first column
    
    $match = $f[$coloffset+0];
    $TgapCount = $f[$coloffset+6];
    $strand = $f[$coloffset+8];
    $qname = $f[$coloffset+9];
    $qsize = $f[$coloffset+10];
    $targetname = $f[$coloffset+13];
    $blockSizes = $f[$coloffset+18];
    $qStarts = $f[$coloffset+19];
    $tStarts = $f[$coloffset+20];
    
    my $match_score=sprintf "%.2f",  (100 * ( $match ) / $qsize);
    my $NM=$f[$coloffset+1];
    #print " TgapCount=", $TgapCount;
    #print " strand=", $strand;
    #print " qname=", $qname;
    #print " blockSizes=", $blockSizes;
    #print " tStarts=", $tStarts, "\n";

    $blockSizes =~ s/[, ]$//;
    $tStarts =~ s/[, ]$//;
    @b = split /,/, $blockSizes;
    @t = split /,/, $tStarts;
    @q = split /,/, $qStarts;
    
    #print "blocksizes ", (join ", ", @b), " blockbegins ", (join ", ", @t) , "\n"; 
    
    $numBlocks = scalar @t;
    # Go throught the line
    #
    @blockbegins=();
    @blockends=();
    for ($i=0; $i<$numBlocks; $i++) {
	$mstart = $t[$i]+1; # blat is 0-based
	$mend = $mstart + $b[$i] - 1;

	push @blockbegins, $mstart;
	push @blockends, $mend;
    }
    $numBlocks = scalar @blockbegins;
    
    my $a, $b;
    if ($prev_tid eq $qname){ # naming transcripts in a meaning full way.
        $counter++;
    }
    else {
        $counter=1;
    }
    #if ($strand ne "-") {
	#print GFF "$targetname\t$source\ttranscript\t",($t[0]+1),"\t", ($t[$numBlocks-1]+$b[$numBlocks-1]), "\t",$match_score,"\t$strand\t.\tID=$qname","_","$counter;QS=$qsize;M=$match;MS=$match_score;NM=$NM","\n";
    #} else {
	print GFF "$targetname\t$source\ttranscript\t",($t[0]+1),"\t", ($t[$numBlocks-1]+$b[$numBlocks-1]), "\t",$match_score,"\t$strand\t",$NM,"\tID=$qname", "_", $counter, ";Match=$match\n";
    #}
	
    for ($i=0; $i<$numBlocks; $i++) {
	#if ($strand ne "-") {
	#    print GFF "$targetname\t$source\texon\t",($t[$i]+1),"\t", ($t[$i]+$b[$i]), "\t",$match_score,"\t$strand\t.\tTarget $source:$qname ", ($q[$i]+1), " ", ($q[$i]+$b[$i]),"\n";
	#} else {
	    print GFF "$targetname\t$source\texon\t",($t[$i]+1),"\t", ($t[$i]+$b[$i]), "\t", $match_score, "\t$strand\t.\tParent=$qname","_","$counter",";Len=",$b[$i],"\n";
	#}
    }
    $estnames{$qname}++;
    $prev_tid=$qname;
}

if (defined $estfilename) {
    foreach (keys %estnames) {
	print EST;
	print EST "\n";
    }
}
