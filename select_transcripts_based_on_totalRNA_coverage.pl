#!/usr/bin/perl -w
use strict;

#27/09/18
#script to find the migh highel expressed transcript of each gene

my $inGtf=$ARGV[0]; 
my $bam_totalRNA=$ARGV[1];
my $outfile=$ARGV[2];

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#pass through the genome, find annotated start codons and setup transcript models for all isoforms of each gene
my %transcript_start_codon_fwd;
my %transcript_stop_codon_fwd;
my %transcript_exons_fwd;
my %transcript_start_codon_rev;
my %transcript_stop_codon_rev;
my %transcript_exons_rev;
my %transcript_2_chr; #key = gene_id; value = chr

my %gene_to_transcript;

open(GENES2,$inGtf) || die "can't open $inGtf";      #gft is 1 based
while (<GENES2>){
    unless(/^#/){
        my @b=split("\t");
        my $chr=$b[0];
        my $class=$b[2];
        my $start=$b[3];
        my $end=$b[4];
        my $dir=$b[6];
        my ($gene_id) = $b[8] =~ /gene_id\s"([^\"]+)";/;
        my ($transcript_id) = $b[8] =~ /transcript_id\s"([^\"]+)";/;
        my $gene_biotype="NA";
        if ($b[8] =~ /gene_biotype\s"([^\"]+)";/){
            $gene_biotype=$1;
        }

        if ($gene_id && $transcript_id){

            if ($gene_biotype eq "protein_coding"){ #restrict to protien coding genes (control for ncRNAs)

                $gene_to_transcript{$gene_id}{$transcript_id}=1;
                $transcript_2_chr{$transcript_id}=$chr;

                if ($dir eq "+"){ #fwd cases. Use start positions as 5'

                    if ($class eq "start_codon"){
                        if (exists ($transcript_start_codon_fwd{$transcript_id})){ #if multiple start codon line take the lower
                            if ($start < $transcript_start_codon_fwd{$transcript_id}){
                                $transcript_start_codon_fwd{$transcript_id}=$start;
                            }
                        }else{
                            $transcript_start_codon_fwd{$transcript_id}=$start;
                        }
                    }

                    if ($class eq "stop_codon"){
                        if (exists ($transcript_stop_codon_fwd{$transcript_id})){ #if multiple stop codon line take the lower
                            if ($start < $transcript_stop_codon_fwd{$transcript_id}){
                                $transcript_stop_codon_fwd{$transcript_id}=$start;
                            }
                        }else{
                            $transcript_stop_codon_fwd{$transcript_id}=$start;
                        }
                    }

                    if ($class eq "exon"){
                        $transcript_exons_fwd{$transcript_id}{$start}=$end;
                    }

                }else{ #revese cases use end as 5'

                    if ($class eq "start_codon"){
                        if (exists ($transcript_start_codon_rev{$transcript_id})){ #if multiple start codon line take the higher
                            if ($end > $transcript_start_codon_rev{$transcript_id}){
                                $transcript_start_codon_rev{$transcript_id}=$end;
                            }
                        }else{
                            $transcript_start_codon_rev{$transcript_id}=$end;
                        }
                    }

                    if ($class eq "stop_codon"){
                        if (exists ($transcript_stop_codon_rev{$transcript_id})){ #if multiple stop codon line take the higher
                            if ($end > $transcript_stop_codon_rev{$transcript_id}){
                                $transcript_stop_codon_rev{$transcript_id}=$end;
                            }
                        }else{
                            $transcript_stop_codon_rev{$transcript_id}=$end;
                        }
                    }

                    if ($class eq "exon"){
                        $transcript_exons_rev{$transcript_id}{$start}=$end;
                    }
                }
            }
        }
    }
}
close(GENES2);

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#setup trascript models

my %transcript_model_fwd;
my %start_coord_fwd;
my %stop_coord_fwd;

my %transcript_exon_exon_in_leader;

#5' is the #transcript_model{$transcript}{0}
#3' is the last coord
my %three_prime_most_coord_fwd;
my %five_prime_most_coord_fwd;

#transcript_model{$transcript}{0}=12345   #1st nt of start codon
#transcript_model{$transcript}{1}=12346
#transcript_model{$transcript}{2}=12347
#transcript_model{$transcript}{3}=12348
#transcript_model{$transcript}{4}=12349
#...
#to end of exons             #last nt of stop codon

for my $transcript (keys %transcript_exons_fwd){
    if ( (exists ($transcript_start_codon_fwd{$transcript})) && (exists ($transcript_stop_codon_fwd{$transcript})) ) { #restrict to transcripts with annotated start + stop codon

        $five_prime_most_coord_fwd{$transcript}=0;
        my $model_pos=0;
        my $exon_count=0; 

        for my $exon_start (sort {$a <=> $b} keys %{ $transcript_exons_fwd{$transcript} } ){
             my $exon_end=$transcript_exons_fwd{$transcript}{$exon_start};

            #fwd exons are in ascending order
            # start(-1)-> 100958 100975
            #             101077 101715 <-end(+1)

            for ($exon_start .. $exon_end){
                $transcript_model_fwd{$transcript}{$model_pos}=$_;

                if ($_ == $transcript_stop_codon_fwd{$transcript}){
                    $stop_coord_fwd{$transcript}=$model_pos;    #find the index of the stop codon per transcript
                }

                if ($_ == $transcript_start_codon_fwd{$transcript}){
                    $start_coord_fwd{$transcript}=$model_pos;    #find the index of the start codon per transcript
                    $transcript_exon_exon_in_leader{$transcript}=$exon_count;
                }
                $model_pos++;
            }
            $exon_count++;
        }
        $three_prime_most_coord_fwd{$transcript}=$model_pos-1; #store the 3 prime most position of each transcript
    }
}

my %transcript_model_rev;
my %start_coord_rev;
my %stop_coord_rev;

#5' is the #transcript_model{$transcript}{0}
#3' is the last coord
my %three_prime_most_coord_rev;
my %five_prime_most_coord_rev;

for my $transcript (keys %transcript_exons_rev){
    if ( (exists ($transcript_start_codon_rev{$transcript})) && (exists ($transcript_stop_codon_rev{$transcript})) ) { #restrict to transcripts with annotated start + stop codon

        $five_prime_most_coord_rev{$transcript}=0;
        my $model_pos=0;
        my $exon_count=0;

        for my $exon_end (reverse (sort {$a <=> $b} keys %{ $transcript_exons_rev{$transcript} } )){
            my $exon_start=$transcript_exons_rev{$transcript}{$exon_end};

            #rev exons are sorted in decending order  
            #           447087 447794 <-start(+1)
            # end(-1)-> 446060 446254

            while ($exon_start >= $exon_end){
                $transcript_model_rev{$transcript}{$model_pos}=$exon_start;

                if ($exon_start == $transcript_stop_codon_rev{$transcript}){
                    $stop_coord_rev{$transcript}=$model_pos;    #find the index of the stop codon per transcript
                }
                if ($exon_start == $transcript_start_codon_rev{$transcript}){
                    $start_coord_rev{$transcript}=$model_pos;    #find the index of the start codon per transcript
                    $transcript_exon_exon_in_leader{$transcript}=$exon_count;
                }
                $model_pos++;
                $exon_start--;
            }
            $exon_count++;

        }
        $three_prime_most_coord_rev{$transcript}=$model_pos-1; #store the 3 prime most position of each transcript
    }
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#loop though transcripts and assign exonic regions to hashes for quick searching.
my %exon_search_fwd; #key1=chr, key2=pos, value=transcript
my %exon_search_rev; #key1=chr, key2=pos, value=transcript
my %exon_counts_totalRNA;
my %transcript_lengths;

for my $transcript (keys %transcript_model_fwd){
    my $chr=$transcript_2_chr{$transcript};
    my $five_prime_coord=$five_prime_most_coord_fwd{$transcript};
    my $three_prime_coord=$three_prime_most_coord_fwd{$transcript};

    $transcript_lengths{$transcript}=$three_prime_coord-$five_prime_coord; 
    $exon_counts_totalRNA{$transcript}=0;

    for my $coord (sort {$a <=> $b} keys %{ $transcript_model_fwd{$transcript} } ){
        my $pos=$transcript_model_fwd{$transcript}{$coord};

        if (exists ($exon_search_fwd{$chr}{$pos})){
            $exon_search_fwd{$chr}{$pos}.=";".$transcript; #concaternate overlapping transcripts
        }else{
            $exon_search_fwd{$chr}{$pos}=$transcript; 
        }    
    }
}

for my $transcript (keys %transcript_model_rev){
    my $chr=$transcript_2_chr{$transcript};
    my $five_prime_coord=$five_prime_most_coord_rev{$transcript};
    my $three_prime_coord=$three_prime_most_coord_rev{$transcript};

    $transcript_lengths{$transcript}=$three_prime_coord-$five_prime_coord;
    $exon_counts_totalRNA{$transcript}=0;

    for my $coord (sort {$a <=> $b} keys %{ $transcript_model_rev{$transcript} } ){
        my $pos=$transcript_model_rev{$transcript}{$coord};

        if (exists ($exon_search_rev{$chr}{$pos})){
            $exon_search_rev{$chr}{$pos}.=";".$transcript; #concaternate overlapping transcripts
        }else{
            $exon_search_rev{$chr}{$pos}=$transcript; 
        }           
    }
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤t¤
#open and store totalRNA counts

open BAM,"samtools view $bam_totalRNA |";

while(<BAM>){

    next if(/^(\@)/);  ## skipping the header lines (if you used -h in the samools command)
    s/\n//;  s/\r//;  ## removing new line
    my @sam = split(/\t+/);  ## splitting SAM line into array

    my $leftMost=$sam[3]; #leftmost position of match 5' for fwd, 3' for rev
    my $flag=$sam[1];
    my $chr=$sam[2];
    my $mapq=$sam[4];
    my $cigar=$sam[5];

    if ($chr=~/chr(.*)/){ #ensembl chromosome names do not contain the "chr" prefix
        $chr=$1;
    }

    unless ($flag & 0x4){   #if aligned

        if ($mapq >= 10){     #mapping uniqueness filter

            if ($flag & 0x10){  #Reverse reads. Starting from the leftmost position parse the cigar and check if matching positions overlap leaders or cds's

                while ($cigar !~ /^$/){
                    if ($cigar =~ /^([0-9]+[MIDN])/){
                        my $cigar_part = $1;
                        if ($cigar_part =~ /(\d+)M/){   #alignment matching 
                            for my $pos ($leftMost .. ($leftMost+$1-1)){ #search though this position

                                if( exists($exon_search_rev{$chr}{$pos})){
                                    my @transcript=split(";", $exon_search_rev{$chr}{$pos});
                                    for my $transcript (@transcript){
                                        $exon_counts_totalRNA{$transcript}++;
                                    }
                                }

                                if( exists($exon_search_fwd{$chr}{$pos})){
                                    my @transcript=split(";", $exon_search_fwd{$chr}{$pos});
                                    for my $transcript (@transcript){
                                        $exon_counts_totalRNA{$transcript}++;
                                    }
                                }

                            }
                            $leftMost+=$1;
                        } elsif ($cigar_part =~ /(\d+)I/){  #insertion (to the reference) #do nothing this region is not in the reference
                        } elsif ($cigar_part =~ /(\d+)D/){  #deletion (from the reference)
                            $leftMost+=$1; #skip this position. Add to position count but do not search
                        } elsif ($cigar_part =~ /(\d+)N/){  #Skipped region from the reference
                            $leftMost+=$1; #skip this position. Add to position count but do not search  
                        }
                        $cigar =~ s/$cigar_part//;
                    }
                }

            }else{ #Forward reads. Starting from the leftmost position parse the cigar and check if matching positions overlap leaders or cds's
                 while ($cigar !~ /^$/){
                    if ($cigar =~ /^([0-9]+[MIDN])/){
                        my $cigar_part = $1;
                        if ($cigar_part =~ /(\d+)M/){   #alignment matching 
                            for my $pos ($leftMost .. ($leftMost+$1-1)){ #search though this position

                                if( exists ($exon_search_rev{$chr}{$pos})){
                                    my @transcript=split(";", $exon_search_rev{$chr}{$pos});
                                    for my $transcript (@transcript){
                                        $exon_counts_totalRNA{$transcript}++;
                                    }
                                }

                                if( exists ($exon_search_fwd{$chr}{$pos})){
                                    my @transcript=split(";", $exon_search_fwd{$chr}{$pos});
                                    for my $transcript (@transcript){
                                        $exon_counts_totalRNA{$transcript}++;
                                    }
                                }
                            }
                            $leftMost+=$1;
                        } elsif ($cigar_part =~ /(\d+)I/){  #insertion (to the reference) #do nothing this region is not in the reference
                        } elsif ($cigar_part =~ /(\d+)D/){  #deletion (from the reference)
                            $leftMost+=$1; #skip this position. Add to position count but do not search
                        } elsif ($cigar_part =~ /(\d+)N/){  #Skipped region from the reference
                            $leftMost+=$1; #skip this position. Add to position count but do not search  
                        }
                        $cigar =~ s/$cigar_part//;
                    }
                }
            }
        }
    }
}
close (BAM);


#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
#For each gene find the most highly expressed transcript

my %highest_coverage;
my %most_highly_translated;

#header
#print "#gene_id,transcript_id,chr,transcript_length,total_RNA_sum,total_RNA_coverage\n";

for my $gene (keys %gene_to_transcript){

   for my $transcript (keys %{ $gene_to_transcript{$gene} }){

        if (exists ($transcript_lengths{$transcript})){ #these will have annottated start and stop codons
            my $transcript_length=$transcript_lengths{$transcript};
            my $chr=$transcript_2_chr{$transcript};
            my $transcript_total_RNA_sum=$exon_counts_totalRNA{$transcript};
            my $transcript_coverage=eval { $transcript_total_RNA_sum/$transcript_length; } || 0;

#            print "$gene,$transcript,$chr,$transcript_length,$transcript_total_RNA_sum,$transcript_coverage\n";

            if (exists ( $most_highly_translated{$gene} ) ){
                if ($transcript_coverage > $highest_coverage{$gene} ){
                    $highest_coverage{$gene}=$transcript_coverage;
                    $most_highly_translated{$gene}=$transcript;
                }
            }else{
                $highest_coverage{$gene}=$transcript_coverage;
                $most_highly_translated{$gene}=$transcript;
            }
        }
    }
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#print highest coverage transcript per gene
open (OUT, ">$outfile") || die;

for my $gene (keys %most_highly_translated){
    print OUT "$gene,$most_highly_translated{$gene}\n";
}
close(OUT);

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

exit;

