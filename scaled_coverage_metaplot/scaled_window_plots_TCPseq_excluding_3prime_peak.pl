#!/usr/bin/perl -w
use strict;

#11/02/19
#script to calculate the coverage, 5' counts and 3' counts across transcripts regions. There regions are then scaled to 100nt.

my $inGtf=$ARGV[0]; 
my $leaders=$ARGV[1];
my $bam_SSU=$ARGV[2];  
my $bam_LSU=$ARGV[3];  
my $outDir=$ARGV[4];
my $most_highly_expressed=$ARGV[5];
my $use_cage=$ARGV[6];  # unless($use_cage eq "exclude")
my ($prefix)=$bam_SSU=~/([^\/]+).bam$/;


#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#open list of most highly expressed transcripts
my %most_expressed_transcript; #key = gene_id, transcript_id, #value = sum_exon_lengths;

open(EXP,$most_highly_expressed) || die "can't open $most_highly_expressed";
while (<EXP>){

    chomp();

    unless(/^#/){
        my @b=split(",");
        my $gene=$b[0];
        my $transcript=$b[1];
        $most_expressed_transcript{$gene}=$transcript;
    }
}
close (EXP);

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#pass through the GTF, find annotated start codons and setup transcript models for the most highly expressed transcript of each gene
my %gene_start_codon_fwd;
my %gene_stop_codon_fwd;
my %gene_exons_fwd;
my %gene_start_codon_rev;
my %gene_stop_codon_rev;
my %gene_exons_rev;
my %gene_2_chr; #key = gene_id; value = chr

my %overlaps_inframe_gene;
my %leader_overlaps_upstream;
my %gene_overlaps_downstream_leader;
my %gene_overlaps_nc;
my %has_cage_defined_leader;

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

        if ($gene_id && $transcript_id){

            if (exists ( $most_expressed_transcript{$gene_id} )){    #if the transcript is in the list of longest transcripts

                if ($transcript_id eq $most_expressed_transcript{$gene_id}){

                    $gene_2_chr{$gene_id}=$chr;
                    $overlaps_inframe_gene{$gene_id}=0;
                    $leader_overlaps_upstream{$gene_id}=0;
                    $gene_overlaps_downstream_leader{$gene_id}=0;
                    $gene_overlaps_nc{$gene_id}=0;
                    $has_cage_defined_leader{$gene_id}=0;

                    if ($dir eq "+"){ #fwd cases. Use start positions as 5'

                        if ($class eq "start_codon"){
                            if (exists ($gene_start_codon_fwd{$gene_id})){ #if multiple start codon line take the lower
                                if ($start < $gene_start_codon_fwd{$gene_id}){
                                    $gene_start_codon_fwd{$gene_id}=$start;
                                }
                            }else{
                                $gene_start_codon_fwd{$gene_id}=$start;
                            }
                        }

                        if ($class eq "stop_codon"){
                            if (exists ($gene_stop_codon_fwd{$gene_id})){ #if multiple stop codon line take the lower
                                if ($start < $gene_stop_codon_fwd{$gene_id}){
                                    $gene_stop_codon_fwd{$gene_id}=$start;
                                }
                            }else{
                                $gene_stop_codon_fwd{$gene_id}=$start;
                            }
                        }

                        if ($class eq "exon"){
                            $gene_exons_fwd{$gene_id}{$start}=$end;
                        }

                    }else{ #revese cases use end as 5'

                        if ($class eq "start_codon"){
                            if (exists ($gene_start_codon_rev{$gene_id})){ #if multiple start codon line take the higher
                                if ($end > $gene_start_codon_rev{$gene_id}){
                                    $gene_start_codon_rev{$gene_id}=$end;
                                }
                            }else{
                                $gene_start_codon_rev{$gene_id}=$end;
                            }
                        }

                        if ($class eq "stop_codon"){
                            if (exists ($gene_stop_codon_rev{$gene_id})){ #if multiple stop codon line take the higher
                                if ($end > $gene_stop_codon_rev{$gene_id}){
                                    $gene_stop_codon_rev{$gene_id}=$end;
                                }
                            }else{
                                $gene_stop_codon_rev{$gene_id}=$end;
                            }
                        }

                        if ($class eq "exon"){
                            $gene_exons_rev{$gene_id}{$start}=$end;
                        }
                    }
                }
            }
        }
    }
}
close(GENES2);

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#setup transcript models

my %gene_model_fwd;
my %start_coord_fwd;
my %stop_coord_fwd;

#5' is the #gene_model{$gene}{0}
#3' is the last coord
my %three_prime_most_coord_fwd;
my %five_prime_most_coord_fwd;

#gene_model{$gene}{0}=12345   #1st nt of start codon
#gene_model{$gene}{1}=12346
#gene_model{$gene}{2}=12347
#gene_model{$gene}{3}=12348
#gene_model{$gene}{4}=12349
#...
#to end of exons             #last nt of stop codon

for my $gene (keys %gene_exons_fwd){
    if ( (exists ($gene_start_codon_fwd{$gene})) && (exists ($gene_stop_codon_fwd{$gene})) ) { #restrict to genes with annotated start + stop codon

        my $model_pos=0;
        $five_prime_most_coord_fwd{$gene}=$model_pos;  #initalise the 5' to the first coord coord

        for my $exon_start (sort {$a <=> $b} keys %{ $gene_exons_fwd{$gene} } ){
            my $exon_end=$gene_exons_fwd{$gene}{$exon_start};

            #fwd exons are in ascending order
            # start(-1)-> 100958 100975
            #             101077 101715 <-end(+1)

            for ($exon_start .. $exon_end){
                $gene_model_fwd{$gene}{$model_pos}=$_;

                if ($_ == $gene_stop_codon_fwd{$gene}){
                    $stop_coord_fwd{$gene}=$model_pos;    #find the index of the stop codon per gene
                }

                if ($_ == $gene_start_codon_fwd{$gene}){
                    $start_coord_fwd{$gene}=$model_pos;    #find the index of the start codon per gene
                }
                $model_pos++;
            }
        }
        $three_prime_most_coord_fwd{$gene}=$model_pos-1; #store the 3 prime most position of each gene
    }
}

my %gene_model_rev;
my %start_coord_rev;
my %stop_coord_rev;

#5' is the #gene_model{$gene}{0}
#3' is the last coord
my %three_prime_most_coord_rev;
my %five_prime_most_coord_rev;

for my $gene (keys %gene_exons_rev){
    if ( (exists ($gene_start_codon_rev{$gene})) && (exists ($gene_stop_codon_rev{$gene})) ) { #restrict to genes with annotated start + stop codon

        my $model_pos=0;

        $five_prime_most_coord_rev{$gene}=$model_pos;  #initalise the 5' to the first coord coord
 
        for my $exon_end (reverse (sort {$a <=> $b} keys %{ $gene_exons_rev{$gene} } )){
            my $exon_start=$gene_exons_rev{$gene}{$exon_end};

            #rev exons are sorted in decending order  
            #           447087 447794 <-start(+1)
            # end(-1)-> 446060 446254

            while ($exon_start >= $exon_end){
                $gene_model_rev{$gene}{$model_pos}=$exon_start;

                if ($exon_start == $gene_stop_codon_rev{$gene}){
                    $stop_coord_rev{$gene}=$model_pos;    #find the index of the stop codon per gene
                }
                if ($exon_start == $gene_start_codon_rev{$gene}){
                    $start_coord_rev{$gene}=$model_pos;    #find the index of the start codon per gene
                }
                $model_pos++;
                $exon_start--;
            }
        }
        $three_prime_most_coord_rev{$gene}=$model_pos-1; #store the 3 prime most position of each gene
    }
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#parse leaders
my %leader_positions_fwd;
my %leader_positions_rev;

#filters:
my %leader_length;
my %cage_peak_value;

#gene_id,       0
#transcript_id, 1
#chr,           2
#dir,           3
#overlaps_inframe_gene,           4
#leader_overlaps_upstream,        5
#gene_overlaps_downstream_leader, 6
#highest_cage_peak,               7
#count_at_highest_cage_peak,      8
#leader_length                    9

#ENSDARG00000037917,ENSDART00000161963,3,fwd,FALSE,FALSE,FALSE,34716685,2.30440468389324,143
#ENSDARG00000104069,ENSDART00000167982,5,fwd,FALSE,FALSE,FALSE,337237,9.98331428397882,122
#ENSDARG00000037925,ENSDART00000130591,3,fwd,FALSE,FALSE,FALSE,36250098,0.346429817638395,-679
#ENSDARG00000029263,ENSDART00000078466,3,fwd,FALSE,FALSE,FALSE,NaN,0,NaN

open(LEAD, $leaders) || die "can't open $leaders";
while (<LEAD>){
    unless(/^#/){
 
        chomp;
        my @b=split(",");

        my $gene=$b[0];
        my $transcript=$b[1];
        my $chr=$b[2];
        my $dir=$b[3];
        my $overlaps_inframe_gene=$b[4];
        my $leader_overlaps_upstream=$b[5];
        my $gene_overlaps_downstream_leader=$b[6];
        my $gene_overlaps_nc=$b[7];
        my $highest_cage_peak=$b[8];
        my $count_at_highest_cage_peak=$b[9];
        my $leader_length=$b[10];

        if ($overlaps_inframe_gene eq "TRUE"){ $overlaps_inframe_gene{$gene}=1; }
        if ($leader_overlaps_upstream eq "TRUE"){ $leader_overlaps_upstream{$gene}=1; }
        if ($gene_overlaps_downstream_leader eq "TRUE"){ $gene_overlaps_downstream_leader{$gene}=1; }
        if ($gene_overlaps_nc eq "TRUE"){ $gene_overlaps_nc{$gene}=1; }

        unless ($leader_length eq "NaN"){  #only take genes that have a detectable cage peak

            unless ($leader_length < 0 ){  #exlude genes with negative leader sizes, as they cuase problems with FPKM

                $has_cage_defined_leader{$gene}=1;

                if ($dir eq "fwd"){ 
                    if (exists ($start_coord_fwd{$gene})){
                        unless ($highest_cage_peak >=  $gene_model_fwd{$gene}{$start_coord_fwd{$gene}}){  #exclude genes where the TSS is downstream of the start codon
                            $leader_positions_fwd{$gene}=$highest_cage_peak;
                            $cage_peak_value{$gene}=$count_at_highest_cage_peak;
                            $leader_length{$gene}=$leader_length;
                            $has_cage_defined_leader{$gene}=1;
                        } 
                    }
                }else{
                    if (exists ($start_coord_rev{$gene})){
                        unless ($highest_cage_peak <=  $gene_model_rev{$gene}{$start_coord_rev{$gene}}){  #exclude genes where the TSS is downstream of the start codon
                            $leader_positions_rev{$gene}=$highest_cage_peak;
                            $cage_peak_value{$gene}=$count_at_highest_cage_peak;
                            $leader_length{$gene}=$leader_length;
                            $has_cage_defined_leader{$gene}=1;
                        }
                    }
                }     
            }
        }
    }
}
close(LEAD);

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#extend transcript models to incorportate cage derived leaders
unless($use_cage eq "exclude"){

    for my $gene (keys %leader_positions_fwd){

        if (exists ($gene_2_chr{$gene})){  #to restict to protien_coding genes
            my $leader_start=$leader_positions_fwd{$gene};
            my $three_prime_coord=$three_prime_most_coord_fwd{$gene};
            my $three_prime_pos=$gene_model_fwd{$gene}{$three_prime_coord};
            my $five_prime_coord=0;
            my $five_prime_pos=$gene_model_fwd{$gene}{$five_prime_coord};
 
            if ($leader_start >= $five_prime_pos){  #great, just find and save the coordinate
 
                for my $coord ($five_prime_coord .. $three_prime_coord){
                    my $pos=$gene_model_fwd{$gene}{$coord};       
                    if ($pos == $leader_start){
                        $five_prime_most_coord_fwd{$gene}=$coord;
#                       if ($gene eq "ENSDARG00000005026"){ print "new_5 shorter: $coord\n"; }
                        last; 
                    }
                }

            }else{  #extend the coords

                my $extended_coord=0; 
                while ($five_prime_pos > $leader_start){        
                    $extended_coord--;
                    $five_prime_pos--;
                    $gene_model_fwd{$gene}{$extended_coord}=$five_prime_pos;
                }
                $five_prime_most_coord_fwd{$gene}=$extended_coord;
   #            if ($gene eq "ENSDARG00000005026"){ print "new_5 extended: $extended_coord\n"; }
            }
        }
    } 

    for my $gene (keys %leader_positions_rev){

        if (exists ($gene_2_chr{$gene})){  #to restict to protien_coding genes
            my $leader_start=$leader_positions_rev{$gene};
            my $three_prime_coord=$three_prime_most_coord_rev{$gene};
            my $three_prime_pos=$gene_model_rev{$gene}{$three_prime_coord};
            my $five_prime_coord=0;
            my $five_prime_pos=$gene_model_rev{$gene}{$five_prime_coord};
 
            if ($leader_start <= $five_prime_pos){  #great, just find and save the coordinate

                for my $coord ($five_prime_coord .. $three_prime_coord){
                    my $pos=$gene_model_rev{$gene}{$coord};
                    if ($pos == $leader_start){
                        $five_prime_most_coord_rev{$gene}=$coord;
                        last;
                    }
                }

            }else{   #extend the coords

                my $extended_coord=0;
                while ($five_prime_pos < $leader_start){
                    $extended_coord--;
                    $five_prime_pos++;
                    $gene_model_rev{$gene}{$extended_coord}=$five_prime_pos;
                }
                $five_prime_most_coord_rev{$gene}=$extended_coord;
            }
        }
    }

}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#Transcript coord key
#Transcript 5' coord   $five_prime_most_coord_fwd{$gene}  #1st nt of annotated transcript
#Transcript 3' coord   $three_prime_most_coord_???{$gene} #last nt of transcript
#Start codon coord     $start_coord_???{$gene}            #1st nt in start codon
#Stop codon coord      $stop_coord_???{$gene}             #1st nt in stop codon
#$gene_model_fwd{$gene}{$coord}==genomic position

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#loop though genes and assign leader, CDS and trailer regions to hashes for quick searching. (added start + stop codons).
my %whole_transcripts_fwd;
my %whole_transcripts_rev;

my %leader_search_fwd; #key1=chr, key2=pos, value=gene
my %leader_search_rev; #key1=chr, key2=pos, value=gene
my %CDS_search_fwd; #key1=chr, key2=pos, value=gene
my %CDS_search_rev; #key1=chr, key2=pos, value=gene
my %trailer_search_fwd; #key1=chr, key2=pos, value=gene
my %trailer_search_rev; #key1=chr, key2=pos, value=gene

my %TTS_ending_search_fwd;
my %TTS_ending_search_rev;

my %leader_counts_five_SSU;
my %leader_counts_three_SSU;
my %leader_counts_coverage_SSU;
my %CDS_counts_five_SSU;
my %CDS_counts_three_SSU;
my %CDS_counts_coverage_SSU;
my %trailer_counts_five_SSU;
my %trailer_counts_three_SSU;
my %trailer_counts_coverage_SSU;

my %leader_counts_five_LSU;
my %leader_counts_three_LSU;
my %leader_counts_coverage_LSU;
my %CDS_counts_five_LSU;
my %CDS_counts_three_LSU;
my %CDS_counts_coverage_LSU;
my %trailer_counts_five_LSU;
my %trailer_counts_three_LSU;
my %trailer_counts_coverage_LSU;

my %length_leader;
my %length_CDS;
my %length_trailer;

for my $gene (keys %gene_model_fwd){
    my $chr=$gene_2_chr{$gene};
    my $start_coord=$start_coord_fwd{$gene};
    my $stop_coord=$stop_coord_fwd{$gene};
    my $relational_position_leader=0;
    my $relational_position_CDS=0;
    my $relational_position_trailer=0;
    my $three_prime_coord=$three_prime_most_coord_fwd{$gene};

    for my $coord (sort {$a <=> $b} keys %{ $gene_model_fwd{$gene} } ){

        if ($coord >= $five_prime_most_coord_fwd{$gene}){

            my $pos=$gene_model_fwd{$gene}{$coord};

            if (exists ($whole_transcripts_fwd{$chr}{$pos})){ #mark overlapping genes
                my @gene = split("_", $whole_transcripts_fwd{$chr}{$pos});
                for (@gene){ $overlaps_inframe_gene{$_}=1; }
                $overlaps_inframe_gene{$gene}=1;
                $whole_transcripts_fwd{$chr}{$pos}.="_".$gene;
            }else{
                $whole_transcripts_fwd{$chr}{$pos}=$gene;
            }

            if ($coord < $start_coord){
                $leader_search_fwd{$chr}{$pos}=$gene.";".$relational_position_leader;
                $leader_counts_five_SSU{$gene}{$relational_position_leader}=0;  #initalise all positions
                $leader_counts_five_LSU{$gene}{$relational_position_leader}=0;
                $leader_counts_three_SSU{$gene}{$relational_position_leader}=0;
                $leader_counts_three_LSU{$gene}{$relational_position_leader}=0;
                $leader_counts_coverage_SSU{$gene}{$relational_position_leader}=0;
                $leader_counts_coverage_LSU{$gene}{$relational_position_leader}=0;
#                print "$gene,fwd,LEAD,$chr,$pos,$relational_position_leader\n";
                $relational_position_leader++;

            }elsif($coord <= ($stop_coord+2)){
                $CDS_search_fwd{$chr}{$pos}=$gene.";".$relational_position_CDS;
                $CDS_counts_five_SSU{$gene}{$relational_position_CDS}=0;  #initalise all positions
                $CDS_counts_five_LSU{$gene}{$relational_position_CDS}=0;
                $CDS_counts_three_SSU{$gene}{$relational_position_CDS}=0;
                $CDS_counts_three_LSU{$gene}{$relational_position_CDS}=0;
                $CDS_counts_coverage_SSU{$gene}{$relational_position_CDS}=0;
                $CDS_counts_coverage_LSU{$gene}{$relational_position_CDS}=0;
#                print "$gene,fwd,CDS,$chr,$pos,$relational_position_CDS\n";
                $relational_position_CDS++;

            }elsif( $coord >= $three_prime_coord-9){   #exlude read overlapping with the last 10nt of the transcript
                $TTS_ending_search_fwd{$chr}{$pos}=$gene;  
                          
            }elsif( $coord > ($stop_coord+2)){
                $trailer_search_fwd{$chr}{$pos}=$gene.";".$relational_position_trailer;
                $trailer_counts_five_SSU{$gene}{$relational_position_trailer}=0;  #initalise all positions
                $trailer_counts_five_LSU{$gene}{$relational_position_trailer}=0;
                $trailer_counts_three_SSU{$gene}{$relational_position_trailer}=0;
                $trailer_counts_three_LSU{$gene}{$relational_position_trailer}=0;
                $trailer_counts_coverage_SSU{$gene}{$relational_position_trailer}=0;
                $trailer_counts_coverage_LSU{$gene}{$relational_position_trailer}=0;
#                print "$gene,fwd,TRAIL,$chr,$pos,$relational_position_trailer\n";
                $relational_position_trailer++;
            }
        }      
    }
    $length_leader{$gene}=$relational_position_leader;
    $length_CDS{$gene}=$relational_position_CDS;
    $length_trailer{$gene}=$relational_position_trailer;
}

for my $gene (keys %gene_model_rev){
    my $chr=$gene_2_chr{$gene};
    my $start_coord=$start_coord_rev{$gene};
    my $stop_coord=$stop_coord_rev{$gene};
    my $relational_position_leader=0;
    my $relational_position_CDS=0;
    my $relational_position_trailer=0;
    my $three_prime_coord=$three_prime_most_coord_rev{$gene};

    for my $coord (sort {$a <=> $b} keys %{ $gene_model_rev{$gene} } ){

        if ($coord >= $five_prime_most_coord_rev{$gene}){
 
            my $pos=$gene_model_rev{$gene}{$coord};

            if (exists ($whole_transcripts_rev{$chr}{$pos})){ #mark overlapping genes
                my @gene = split("_", $whole_transcripts_rev{$chr}{$pos});
                for (@gene){ $overlaps_inframe_gene{$_}=1; }
                $overlaps_inframe_gene{$gene}=1;
                $whole_transcripts_rev{$chr}{$pos}.="_".$gene;
            }else{
                $whole_transcripts_rev{$chr}{$pos}=$gene;
            }

            if ($coord < $start_coord){
                $leader_search_rev{$chr}{$pos}=$gene.";".$relational_position_leader;
                $leader_counts_five_SSU{$gene}{$relational_position_leader}=0;  #initalise all positions
                $leader_counts_five_LSU{$gene}{$relational_position_leader}=0;
                $leader_counts_three_SSU{$gene}{$relational_position_leader}=0;
                $leader_counts_three_LSU{$gene}{$relational_position_leader}=0;
                $leader_counts_coverage_SSU{$gene}{$relational_position_leader}=0;
                $leader_counts_coverage_LSU{$gene}{$relational_position_leader}=0;
#                print "$gene,rev,LEAD,$chr,$pos,$relational_position_leader\n";                   
                $relational_position_leader++;

            }elsif($coord <= ($stop_coord+2)){
                $CDS_search_rev{$chr}{$pos}=$gene.";".$relational_position_CDS;
                $CDS_counts_five_SSU{$gene}{$relational_position_CDS}=0;  #initalise all positions
                $CDS_counts_five_LSU{$gene}{$relational_position_CDS}=0;
                $CDS_counts_three_SSU{$gene}{$relational_position_CDS}=0;
                $CDS_counts_three_LSU{$gene}{$relational_position_CDS}=0;
                $CDS_counts_coverage_SSU{$gene}{$relational_position_CDS}=0;
                $CDS_counts_coverage_LSU{$gene}{$relational_position_CDS}=0;
#                print "$gene,rev,CDS,$chr,$pos,$relational_position_CDS\n";
                $relational_position_CDS++;

            }elsif ($coord >= $three_prime_coord-9){   #exlude read overlapping with the last 10nt of the transcript
                $TTS_ending_search_rev{$chr}{$pos}=$gene;    

            }elsif( $coord > ($stop_coord+2)){
                $trailer_search_rev{$chr}{$pos}=$gene.";".$relational_position_trailer;
                $trailer_counts_five_SSU{$gene}{$relational_position_trailer}=0;  #initalise all positions
                $trailer_counts_five_LSU{$gene}{$relational_position_trailer}=0;
                $trailer_counts_three_SSU{$gene}{$relational_position_trailer}=0;
                $trailer_counts_three_LSU{$gene}{$relational_position_trailer}=0;
                $trailer_counts_coverage_SSU{$gene}{$relational_position_trailer}=0;
                $trailer_counts_coverage_LSU{$gene}{$relational_position_trailer}=0;
#                print "$gene,rev,TRAIL,$chr,$pos,$relational_position_trailer\n";
                $relational_position_trailer++;
            }
        }   
    }
    $length_leader{$gene}=$relational_position_leader;
    $length_CDS{$gene}=$relational_position_CDS;
    $length_trailer{$gene}=$relational_position_trailer;
}

#for my $g (keys %length_leader){
#    print "$g,$length_leader{$g},$length_CDS{$g},$length_trailer{$g}\n";
#}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#open SSU and assign
open SSU,"samtools view $bam_SSU |";
while(<SSU>){

    next if(/^(\@)/);  ## skipping the header lines (if you used -h in the samools command)
    s/\n//;  s/\r//;  ## removing new line
    my @sam = split(/\t+/);  ## splitting SAM line into array

    my $leftMost=$sam[3]; #leftmost position of match 5' for fwd, 3' for rev
    my $flag=$sam[1];
    my $chr=$sam[2];
    my $mapq=$sam[4];
    my $cigar=$sam[5];
    my $seq=$sam[9];
    my $threePrime;
    my $fivePrime;

    if ($chr=~/chr(.*)/){ #ensembl chromosome names do not contain the "chr" prefix
        $chr=$1;
    }

    unless ($flag & 0x4){   #if aligned

        if ($mapq >= 10){     #mapping uniqnes filter

            #both chew_2013 and subtelney_2014 riboseq are directional
            if ($flag & 0x10){  #Reverse reads. Starting from the leftmost position parse the cigar and check if it matches a uORF

                $threePrime=$leftMost;               #assign 5' amd 3' to positions
                my $length=length($seq);             #parse cigar for indels and adjust the length of the alignment     
                while ($cigar =~/(\d+)I/g){          #add to length for insertions
                    $length+=$1;
                }
                while ($cigar =~/(\d+)D/g){          #substact from length for deletions
                    $length-=$1;
                }
                $fivePrime=$leftMost+($length-1);    #SAM is 1 based

                unless ( (exists ($TTS_ending_search_rev{$chr}{$fivePrime})) || (exists ($TTS_ending_search_rev{$chr}{$threePrime})) ){ #exlude reads overlapping the 3' end of the transcript
              
                    if (exists ($leader_search_rev{$chr}{$fivePrime})){
                        my ($gene,$meta_pos)=$leader_search_rev{$chr}{$fivePrime}=~/^([^\;]+)\;(\d+)$/;
                        $leader_counts_five_SSU{$gene}{$meta_pos}++;
                    }

                    if (exists ($leader_search_rev{$chr}{$threePrime})){
                        my ($gene,$meta_pos)=$leader_search_rev{$chr}{$threePrime}=~/^([^\;]+)\;(\d+)$/;
                        $leader_counts_three_SSU{$gene}{$meta_pos}++;
                    }

                    if (exists ($CDS_search_rev{$chr}{$fivePrime})){
                         my ($gene,$meta_pos)=$CDS_search_rev{$chr}{$fivePrime}=~/^([^\;]+)\;(\d+)$/;
                        $CDS_counts_five_SSU{$gene}{$meta_pos}++;
                    }

                    if (exists ($CDS_search_rev{$chr}{$threePrime})){
                        my ($gene,$meta_pos)=$CDS_search_rev{$chr}{$threePrime}=~/^([^\;]+)\;(\d+)$/;
                        $CDS_counts_three_SSU{$gene}{$meta_pos}++;
                    }   
            
                    if (exists ($trailer_search_rev{$chr}{$fivePrime})){
                        my ($gene,$meta_pos)=$trailer_search_rev{$chr}{$fivePrime}=~/^([^\;]+)\;(\d+)$/;
                        $trailer_counts_five_SSU{$gene}{$meta_pos}++;
                    }

                    if (exists ($trailer_search_rev{$chr}{$threePrime})){
                        my ($gene,$meta_pos)=$trailer_search_rev{$chr}{$threePrime}=~/^([^\;]+)\;(\d+)$/;
                        $trailer_counts_three_SSU{$gene}{$meta_pos}++;
                    }   
                
                    while ($cigar !~ /^$/){
                        if ($cigar =~ /^([0-9]+[MIDN])/){
                            my $cigar_part = $1;
                            if ($cigar_part =~ /(\d+)M/){   #alignment matching 
                                for my $pos ($leftMost .. ($leftMost+$1-1)){ #search though this position

                                    if (exists ($leader_search_rev{$chr}{$pos})){
                                        my ($gene,$meta_pos)=$leader_search_rev{$chr}{$pos}=~/^([^\;]+)\;(\d+)$/;
                                        $leader_counts_coverage_SSU{$gene}{$meta_pos}++;
                                    }

                                    if (exists ($CDS_search_rev{$chr}{$pos})){
                                        my ($gene,$meta_pos)=$CDS_search_rev{$chr}{$pos}=~/^([^\;]+)\;(\d+)$/;
                                        $CDS_counts_coverage_SSU{$gene}{$meta_pos}++;
                                    }
   
                                    if (exists ($trailer_search_rev{$chr}{$pos})){
                                        my ($gene,$meta_pos)=$trailer_search_rev{$chr}{$pos}=~/^([^\;]+)\;(\d+)$/;
                                        $trailer_counts_coverage_SSU{$gene}{$meta_pos}++;
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
            }else{ #Forward reads. Starting from the leftmost position parse the cigar and check if it matches a uORF

                $fivePrime=$leftMost;
                my $length=length($seq);
                while ($cigar =~/(\d+)I/g){           #add to length for insertions
                    $length+=$1;
                }
                while ($cigar =~/(\d+)D/g){           #substact from length for deletions
                    $length-=$1;
                }
                $threePrime=$leftMost+($length-1);    #SAM is 1 based

                unless ( (exists ($TTS_ending_search_fwd{$chr}{$fivePrime})) || (exists ($TTS_ending_search_fwd{$chr}{$threePrime})) ){ #exlude reads overlapping the 3' end of the transcript

                    if (exists ($leader_search_fwd{$chr}{$fivePrime})){
                        my ($gene,$meta_pos)=$leader_search_fwd{$chr}{$fivePrime}=~/^([^\;]+)\;(\d+)$/;
                        $leader_counts_five_SSU{$gene}{$meta_pos}++;
                    }
    
                    if (exists ($leader_search_fwd{$chr}{$threePrime})){
                        my ($gene,$meta_pos)=$leader_search_fwd{$chr}{$threePrime}=~/^([^\;]+)\;(\d+)$/;
                        $leader_counts_three_SSU{$gene}{$meta_pos}++;
                    }   
    
                    if (exists ($CDS_search_fwd{$chr}{$fivePrime})){
                        my ($gene,$meta_pos)=$CDS_search_fwd{$chr}{$fivePrime}=~/^([^\;]+)\;(\d+)$/;
                        $CDS_counts_five_SSU{$gene}{$meta_pos}++;
                    }
    
                    if (exists ($CDS_search_fwd{$chr}{$threePrime})){
                        my ($gene,$meta_pos)=$CDS_search_fwd{$chr}{$threePrime}=~/^([^\;]+)\;(\d+)$/;
                        $CDS_counts_three_SSU{$gene}{$meta_pos}++;
                    }
    
                    if (exists ($trailer_search_fwd{$chr}{$fivePrime})){
                        my ($gene,$meta_pos)=$trailer_search_fwd{$chr}{$fivePrime}=~/^([^\;]+)\;(\d+)$/;
                        $trailer_counts_five_SSU{$gene}{$meta_pos}++;
                    }
    
                    if (exists ($trailer_search_fwd{$chr}{$threePrime})){
                        my ($gene,$meta_pos)=$trailer_search_fwd{$chr}{$threePrime}=~/^([^\;]+)\;(\d+)$/;
                        $trailer_counts_three_SSU{$gene}{$meta_pos}++;
                    }
    
                    while ($cigar !~ /^$/){
                        if ($cigar =~ /^([0-9]+[MIDN])/){
                            my $cigar_part = $1;
                            if ($cigar_part =~ /(\d+)M/){   #alignment matching 
                                for my $pos ($leftMost .. ($leftMost+$1-1)){ #search though this position
            
                                    if (exists ($leader_search_fwd{$chr}{$pos})){
                                        my ($gene,$meta_pos)=$leader_search_fwd{$chr}{$pos}=~/^([^\;]+)\;(\d+)$/;
                                        $leader_counts_coverage_SSU{$gene}{$meta_pos}++;
                                    }
    
                                    if (exists ($CDS_search_fwd{$chr}{$pos})){
                                        my ($gene,$meta_pos)=$CDS_search_fwd{$chr}{$pos}=~/^([^\;]+)\;(\d+)$/;
                                        $CDS_counts_coverage_SSU{$gene}{$meta_pos}++;
                                    }
    
                                    if (exists ($trailer_search_fwd{$chr}{$pos})){
                                        my ($gene,$meta_pos)=$trailer_search_fwd{$chr}{$pos}=~/^([^\;]+)\;(\d+)$/;
                                        $trailer_counts_coverage_SSU{$gene}{$meta_pos}++;
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
}
close(SSU);

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#again for LSU

open LSU,"samtools view $bam_LSU |";
while(<LSU>){

    next if(/^(\@)/);  ## skipping the header lines (if you used -h in the samools command)
    s/\n//;  s/\r//;  ## removing new line
    my @sam = split(/\t+/);  ## splitting SAM line into array

    my $leftMost=$sam[3]; #leftmost position of match 5' for fwd, 3' for rev
    my $flag=$sam[1];
    my $chr=$sam[2];
    my $mapq=$sam[4];
    my $cigar=$sam[5];
    my $seq=$sam[9];
    my $threePrime;
    my $fivePrime;

    if ($chr=~/chr(.*)/){ #ensembl chromosome names do not contain the "chr" prefix
        $chr=$1;
    }

    unless ($flag & 0x4){   #if aligned

        if ($mapq >= 10){     #mapping uniqnes filter

            if ($flag & 0x10){  #Reverse reads. Starting from the leftmost position parse the cigar and check if it matches a uORF

                $threePrime=$leftMost;               #assign 5' amd 3' to positions
                my $length=length($seq);             #parse cigar for indels and adjust the length of the alignment     
                while ($cigar =~/(\d+)I/g){          #add to length for insertions
                    $length+=$1;
                }
                while ($cigar =~/(\d+)D/g){          #substact from length for deletions
                    $length-=$1;
                }
                $fivePrime=$leftMost+($length-1);    #SAM is 1 based

                unless ( (exists ($TTS_ending_search_rev{$chr}{$fivePrime})) || (exists ($TTS_ending_search_rev{$chr}{$threePrime})) ){ #exlude reads overlapping the 3' end of the transcript

                    if (exists ($leader_search_rev{$chr}{$fivePrime})){
                        my ($gene,$meta_pos)=$leader_search_rev{$chr}{$fivePrime}=~/^([^\;]+)\;(\d+)$/;
                        $leader_counts_five_LSU{$gene}{$meta_pos}++;
                    }
    
                    if (exists ($leader_search_rev{$chr}{$threePrime})){
                        my ($gene,$meta_pos)=$leader_search_rev{$chr}{$threePrime}=~/^([^\;]+)\;(\d+)$/;
                        $leader_counts_three_LSU{$gene}{$meta_pos}++;
                    }   
    
                    if (exists ($CDS_search_rev{$chr}{$fivePrime})){
                        my ($gene,$meta_pos)=$CDS_search_rev{$chr}{$fivePrime}=~/^([^\;]+)\;(\d+)$/;
                        $CDS_counts_five_LSU{$gene}{$meta_pos}++;
                    }
    
                    if (exists ($CDS_search_rev{$chr}{$threePrime})){
                        my ($gene,$meta_pos)=$CDS_search_rev{$chr}{$threePrime}=~/^([^\;]+)\;(\d+)$/;
                        $CDS_counts_three_LSU{$gene}{$meta_pos}++;
                    }
    
                    if (exists ($trailer_search_rev{$chr}{$fivePrime})){
                        my ($gene,$meta_pos)=$trailer_search_rev{$chr}{$fivePrime}=~/^([^\;]+)\;(\d+)$/;
                        $trailer_counts_five_LSU{$gene}{$meta_pos}++;
                    }
    
                    if (exists ($trailer_search_rev{$chr}{$threePrime})){
                        my ($gene,$meta_pos)=$trailer_search_rev{$chr}{$threePrime}=~/^([^\;]+)\;(\d+)$/;
                        $trailer_counts_three_LSU{$gene}{$meta_pos}++;
                    }
    
                    while ($cigar !~ /^$/){
                        if ($cigar =~ /^([0-9]+[MIDN])/){
                            my $cigar_part = $1;
                            if ($cigar_part =~ /(\d+)M/){   #alignment matching 
                                for my $pos ($leftMost .. ($leftMost+$1-1)){ #search though this position
    
                                    if (exists ($leader_search_rev{$chr}{$pos})){
                                        my ($gene,$meta_pos)=$leader_search_rev{$chr}{$pos}=~/^([^\;]+)\;(\d+)$/;
                                        $leader_counts_coverage_LSU{$gene}{$meta_pos}++;
                                    }
    
                                    if (exists ($CDS_search_rev{$chr}{$pos})){
                                        my ($gene,$meta_pos)=$CDS_search_rev{$chr}{$pos}=~/^([^\;]+)\;(\d+)$/;
                                        $CDS_counts_coverage_LSU{$gene}{$meta_pos}++;
                                    }
    
                                    if (exists ($trailer_search_rev{$chr}{$pos})){
                                        my ($gene,$meta_pos)=$trailer_search_rev{$chr}{$pos}=~/^([^\;]+)\;(\d+)$/;
                                        $trailer_counts_coverage_LSU{$gene}{$meta_pos}++;
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
            }else{ #Forward reads. Starting from the leftmost position parse the cigar and check if it matches a uORF

                $fivePrime=$leftMost;
                my $length=length($seq);
                while ($cigar =~/(\d+)I/g){           #add to length for insertions
                    $length+=$1;
                }
                while ($cigar =~/(\d+)D/g){           #substact from length for deletions
                    $length-=$1;
                }
                $threePrime=$leftMost+($length-1);    #SAM is 1 based

                unless ( (exists ($TTS_ending_search_fwd{$chr}{$fivePrime})) || (exists ($TTS_ending_search_fwd{$chr}{$threePrime})) ){ #exlude reads overlapping the 3' end of the transcript

                    if (exists ($leader_search_fwd{$chr}{$fivePrime})){
                        my ($gene,$meta_pos)=$leader_search_fwd{$chr}{$fivePrime}=~/^([^\;]+)\;(\d+)$/;
                        $leader_counts_five_LSU{$gene}{$meta_pos}++;
                    }
    
                    if (exists ($leader_search_fwd{$chr}{$threePrime})){
                        my ($gene,$meta_pos)=$leader_search_fwd{$chr}{$threePrime}=~/^([^\;]+)\;(\d+)$/;
                        $leader_counts_three_LSU{$gene}{$meta_pos}++;
                    }   
    
                    if (exists ($CDS_search_fwd{$chr}{$fivePrime})){
                        my ($gene,$meta_pos)=$CDS_search_fwd{$chr}{$fivePrime}=~/^([^\;]+)\;(\d+)$/;
                        $CDS_counts_five_LSU{$gene}{$meta_pos}++;
                    }
    
                    if (exists ($CDS_search_fwd{$chr}{$threePrime})){
                        my ($gene,$meta_pos)=$CDS_search_fwd{$chr}{$threePrime}=~/^([^\;]+)\;(\d+)$/;
                        $CDS_counts_three_LSU{$gene}{$meta_pos}++;
                    }
    
                    if (exists ($trailer_search_fwd{$chr}{$fivePrime})){
                        my ($gene,$meta_pos)=$trailer_search_fwd{$chr}{$fivePrime}=~/^([^\;]+)\;(\d+)$/;
                        $trailer_counts_five_LSU{$gene}{$meta_pos}++;
                    }
    
                    if (exists ($trailer_search_fwd{$chr}{$threePrime})){
                        my ($gene,$meta_pos)=$trailer_search_fwd{$chr}{$threePrime}=~/^([^\;]+)\;(\d+)$/;
                        $trailer_counts_three_LSU{$gene}{$meta_pos}++;
                    }   
    
                    while ($cigar !~ /^$/){
                        if ($cigar =~ /^([0-9]+[MIDN])/){
                            my $cigar_part = $1;
                            if ($cigar_part =~ /(\d+)M/){   #alignment matching 
                                for my $pos ($leftMost .. ($leftMost+$1-1)){ #search though this position
    
                                    if (exists ($leader_search_fwd{$chr}{$pos})){
                                        my ($gene,$meta_pos)=$leader_search_fwd{$chr}{$pos}=~/^([^\;]+)\;(\d+)$/;
                                        $leader_counts_coverage_LSU{$gene}{$meta_pos}++;
                                    }
    
                                    if (exists ($CDS_search_fwd{$chr}{$pos})){
                                        my ($gene,$meta_pos)=$CDS_search_fwd{$chr}{$pos}=~/^([^\;]+)\;(\d+)$/;
                                        $CDS_counts_coverage_LSU{$gene}{$meta_pos}++;
                                    }
    
                                    if (exists ($trailer_search_fwd{$chr}{$pos})){
                                        my ($gene,$meta_pos)=$trailer_search_fwd{$chr}{$pos}=~/^([^\;]+)\;(\d+)$/;
                                        $trailer_counts_coverage_LSU{$gene}{$meta_pos}++;
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
}
close(LSU);

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
my $gene_matrix_5=$outDir."/".$prefix."_gene_matrix_5prime.csv";
my $gene_matrix_3=$outDir."/".$prefix."_gene_matrix_3prime.csv";
my $gene_matrix_c=$outDir."/".$prefix."_gene_matrix_coverage.csv";

open (OUT1,">$gene_matrix_5")  || die "can't open $gene_matrix_5\n";
open (OUT2,">$gene_matrix_3")  || die "can't open $gene_matrix_3\n";
open (OUT3,">$gene_matrix_c")  || die "can't open $gene_matrix_c\n";

print OUT1 "Gene,Fraction,Feature,Positon,Count\n";
print OUT2 "Gene,Fraction,Feature,Positon,Count\n";
print OUT3 "Gene,Fraction,Feature,Positon,Count\n";

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
#output the 5' plots

for my $gene (keys %length_CDS){
    my $chr=$gene_2_chr{$gene};
    my $lead_length=$length_leader{$gene};
    my $code_length=$length_CDS{$gene};
    my $tail_length=$length_trailer{$gene};

    if ( $overlaps_inframe_gene{$gene} == 0 && $leader_overlaps_upstream{$gene} == 0 && $gene_overlaps_downstream_leader{$gene} ==0  && $gene_overlaps_nc{$gene}==0 ){

        if ($lead_length >= 100 && $code_length >= 100 && $tail_length >= 100){
                                                                                                        
            my $SSU_aggregated_leader_ref=&aggregate_to_100($lead_length, \%{$leader_counts_five_SSU{$gene}});
            my %SSU_aggregated_leader=%{$SSU_aggregated_leader_ref};

            my $SSU_aggregated_CDS_ref=&aggregate_to_100($code_length, \%{$CDS_counts_five_SSU{$gene}}, $gene);
            my %SSU_aggregated_CDS=%{$SSU_aggregated_CDS_ref};

            my $SSU_aggregated_trailer_ref=&aggregate_to_100($tail_length, \%{$trailer_counts_five_SSU{$gene}});
            my %SSU_aggregated_trailer=%{$SSU_aggregated_trailer_ref};

            for my $pos (sort {$a <=> $b} keys %SSU_aggregated_leader){   print OUT1 "$gene,SSU,leader,$pos,$SSU_aggregated_leader{$pos}\n"; }
            for my $pos (sort {$a <=> $b} keys %SSU_aggregated_CDS){      print OUT1 "$gene,SSU,cds,$pos,$SSU_aggregated_CDS{$pos}\n"; }
            for my $pos (sort {$a <=> $b} keys %SSU_aggregated_trailer){  print OUT1 "$gene,SSU,trailer,$pos,$SSU_aggregated_trailer{$pos}\n"; } 

            my $LSU_aggregated_leader_ref=&aggregate_to_100($lead_length, \%{$leader_counts_five_LSU{$gene}});
            my %LSU_aggregated_leader=%{$LSU_aggregated_leader_ref};
 
            my $LSU_aggregated_CDS_ref=&aggregate_to_100($code_length, \%{$CDS_counts_five_LSU{$gene}});
            my %LSU_aggregated_CDS=%{$LSU_aggregated_CDS_ref};

            my $LSU_aggregated_trailer_ref=&aggregate_to_100($tail_length, \%{$trailer_counts_five_LSU{$gene}});
            my %LSU_aggregated_trailer=%{$LSU_aggregated_trailer_ref};

            for my $pos (sort {$a <=> $b} keys %LSU_aggregated_leader){   print OUT1 "$gene,LSU,leader,$pos,$LSU_aggregated_leader{$pos}\n"; }
            for my $pos (sort {$a <=> $b} keys %LSU_aggregated_CDS){      print OUT1 "$gene,LSU,cds,$pos,$LSU_aggregated_CDS{$pos}\n"; }
            for my $pos (sort {$a <=> $b} keys %LSU_aggregated_trailer){  print OUT1 "$gene,LSU,trailer,$pos,$LSU_aggregated_trailer{$pos}\n"; }
        }
    }
}

close (OUT1);

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#output the 3' plots

for my $gene (keys %length_CDS){
    my $chr=$gene_2_chr{$gene};
    my $lead_length=$length_leader{$gene};
    my $code_length=$length_CDS{$gene};
    my $tail_length=$length_trailer{$gene};

    if ( $overlaps_inframe_gene{$gene} == 0 && $leader_overlaps_upstream{$gene} == 0 && $gene_overlaps_downstream_leader{$gene} ==0 && $gene_overlaps_nc{$gene}==0 ){

        if ($lead_length >= 100 && $code_length >= 100 && $tail_length >= 100){

            my $SSU_aggregated_leader_ref=&aggregate_to_100($lead_length, \%{$leader_counts_three_SSU{$gene}});
            my %SSU_aggregated_leader=%{$SSU_aggregated_leader_ref};

            my $SSU_aggregated_CDS_ref=&aggregate_to_100($code_length, \%{$CDS_counts_three_SSU{$gene}});
            my %SSU_aggregated_CDS=%{$SSU_aggregated_CDS_ref};

            my $SSU_aggregated_trailer_ref=&aggregate_to_100($tail_length, \%{$trailer_counts_three_SSU{$gene}});
            my %SSU_aggregated_trailer=%{$SSU_aggregated_trailer_ref};

            for my $pos (sort {$a <=> $b} keys %SSU_aggregated_leader){   print OUT2 "$gene,SSU,leader,$pos,$SSU_aggregated_leader{$pos}\n"; }
            for my $pos (sort {$a <=> $b} keys %SSU_aggregated_CDS){      print OUT2 "$gene,SSU,cds,$pos,$SSU_aggregated_CDS{$pos}\n"; }
            for my $pos (sort {$a <=> $b} keys %SSU_aggregated_trailer){  print OUT2 "$gene,SSU,trailer,$pos,$SSU_aggregated_trailer{$pos}\n"; }

            my $LSU_aggregated_leader_ref=&aggregate_to_100($lead_length, \%{$leader_counts_three_LSU{$gene}});
            my %LSU_aggregated_leader=%{$LSU_aggregated_leader_ref};

            my $LSU_aggregated_CDS_ref=&aggregate_to_100($code_length, \%{$CDS_counts_three_LSU{$gene}});
            my %LSU_aggregated_CDS=%{$LSU_aggregated_CDS_ref};

            my $LSU_aggregated_trailer_ref=&aggregate_to_100($tail_length, \%{$trailer_counts_three_LSU{$gene}});
            my %LSU_aggregated_trailer=%{$LSU_aggregated_trailer_ref};

            for my $pos (sort {$a <=> $b} keys %LSU_aggregated_leader){   print OUT2 "$gene,LSU,leader,$pos,$LSU_aggregated_leader{$pos}\n"; }
            for my $pos (sort {$a <=> $b} keys %LSU_aggregated_CDS){      print OUT2 "$gene,LSU,cds,$pos,$LSU_aggregated_CDS{$pos}\n"; }
            for my $pos (sort {$a <=> $b} keys %LSU_aggregated_trailer){  print OUT2 "$gene,LSU,trailer,$pos,$LSU_aggregated_trailer{$pos}\n"; }
        }
    }
}

close (OUT2);

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#output the coverage plots

for my $gene (keys %length_CDS){
    my $chr=$gene_2_chr{$gene};
    my $lead_length=$length_leader{$gene};
    my $code_length=$length_CDS{$gene};
    my $tail_length=$length_trailer{$gene};

    if ( $overlaps_inframe_gene{$gene} == 0 && $leader_overlaps_upstream{$gene} == 0 && $gene_overlaps_downstream_leader{$gene} ==0  && $gene_overlaps_nc{$gene}==0 ){

        if ($lead_length >= 100 && $code_length >= 100 && $tail_length >= 100){

            my $SSU_aggregated_leader_ref=&aggregate_to_100($lead_length, \%{$leader_counts_coverage_SSU{$gene}});
            my %SSU_aggregated_leader=%{$SSU_aggregated_leader_ref};

            my $SSU_aggregated_CDS_ref=&aggregate_to_100($code_length, \%{$CDS_counts_coverage_SSU{$gene}});
            my %SSU_aggregated_CDS=%{$SSU_aggregated_CDS_ref};

            my $SSU_aggregated_trailer_ref=&aggregate_to_100($tail_length, \%{$trailer_counts_coverage_SSU{$gene}});
            my %SSU_aggregated_trailer=%{$SSU_aggregated_trailer_ref};
    
            for my $pos (sort {$a <=> $b} keys %SSU_aggregated_leader){   print OUT3 "$gene,SSU,leader,$pos,$SSU_aggregated_leader{$pos}\n"; }
            for my $pos (sort {$a <=> $b} keys %SSU_aggregated_CDS){      print OUT3 "$gene,SSU,cds,$pos,$SSU_aggregated_CDS{$pos}\n"; }
            for my $pos (sort {$a <=> $b} keys %SSU_aggregated_trailer){  print OUT3 "$gene,SSU,trailer,$pos,$SSU_aggregated_trailer{$pos}\n"; } 

            my $LSU_aggregated_leader_ref=&aggregate_to_100($lead_length, \%{$leader_counts_coverage_LSU{$gene}});
            my %LSU_aggregated_leader=%{$LSU_aggregated_leader_ref};
    
            my $LSU_aggregated_CDS_ref=&aggregate_to_100($code_length, \%{$CDS_counts_coverage_LSU{$gene}});
            my %LSU_aggregated_CDS=%{$LSU_aggregated_CDS_ref};
    
            my $LSU_aggregated_trailer_ref=&aggregate_to_100($tail_length, \%{$trailer_counts_coverage_LSU{$gene}});
            my %LSU_aggregated_trailer=%{$LSU_aggregated_trailer_ref};

            for my $pos (sort {$a <=> $b} keys %LSU_aggregated_leader){   print OUT3 "$gene,LSU,leader,$pos,$LSU_aggregated_leader{$pos}\n"; }
            for my $pos (sort {$a <=> $b} keys %LSU_aggregated_CDS){      print OUT3 "$gene,LSU,cds,$pos,$LSU_aggregated_CDS{$pos}\n"; }
            for my $pos (sort {$a <=> $b} keys %LSU_aggregated_trailer){  print OUT3 "$gene,LSU,trailer,$pos,$LSU_aggregated_trailer{$pos}\n"; } 
        }
    }
}

close (OUT3);

exit;

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#scale (aggrgate) the values in a positional hash to 100nt in length
sub aggregate_to_100{

    my $lead_length = shift;
    my $positional_hash_reference = shift;
    my %positional_hash = %$positional_hash_reference;
    my %tmp_hash;
    my %out_hash;
    my $lead_scaling_factor=$lead_length/100;

    for my $pos (sort {$a <=> $b} keys %positional_hash){ 

        #if ($gene eq "YDR395W"){ print "$pos,$positional_hash{$pos}\n"; } #checked
        my $updated_pos=int(($pos/$lead_scaling_factor)+0.99);  #round the fraction up
        #my $updated_pos=int($pos/$lead_scaling_factor);  #round the fraction down

        unless(exists $tmp_hash{$updated_pos} ){
            $tmp_hash{$updated_pos}=$positional_hash{$pos};
        }else{
            $tmp_hash{$updated_pos}.="_".$positional_hash{$pos};
        }
    }

    #aggregate
    for my $scaled_position (sort {$a <=> $b} keys %tmp_hash){
        my @poss=split("_",$tmp_hash{$scaled_position});
        my $count=@poss;
        my $sum=0;
        for (@poss){
            $sum+=$_;
        }
        my $aggregated_value=$sum/$count;
        $out_hash{$scaled_position}=$aggregated_value;
        #if ($gene eq "YML009W-B"){ print "$gene,$scaled_position,$aggregated_value\n"; } #checked
    } 
    return (\%out_hash);
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
