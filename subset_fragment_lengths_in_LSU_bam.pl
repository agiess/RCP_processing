#!/usr/bin/perl -w
use strict;

#19/06/18
#script to select fragment lengths corrosponding to translating peaks in LSU TCP-seq libraries
#there wll be a second script to remove the translating lengths from the SSU bam

my $bam_TCPseq=$ARGV[0];
my $MIN_LENGTH=$ARGV[1];
my $MAX_LENGTH=$ARGV[2];

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
open BAM,"samtools view -h $bam_TCPseq |";

while(<BAM>){

    if(/^(\@)/){  #header lines (if you used -h in the samools command)

        print $_;

    }else{

        s/\n//;  s/\r//;  ## removing new line
        my @sam = split(/\t+/);  ## splitting SAM line into array

        my $leftMost=$sam[3]; #leftmost position of match 5' for fwd, 3' for rev
        my $flag=$sam[1];
        my $chr=$sam[2];
        my $mapq=$sam[4];
        my $cigar=$sam[5];
        my $seq=$sam[9];

        my $fivePrime;
        my $threePrime;

        if ($chr=~/chr(.*)/){ #ensembl chromosome names do not contain the "chr" prefix
            $chr=$1;
        }

        unless ($flag & 0x4){   #if aligned
            my $length=length($seq);
            if ($length >= $MIN_LENGTH && $length <= $MAX_LENGTH){
                print "$_\n";
            }
        }
    }
}
close (BAM);

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

exit;
