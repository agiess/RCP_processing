#!/usr/bin/perl -w
use strict;

#28/09/18
#script to assign calculate the count of fragment lengths mappong to '5 leader, TIS, CDS, stop codon and 3' trailer features.
#updated to include start of transcript (TSS) and end of transcript as features

my $inGtf=$ARGV[0]; 
my $fasta=$ARGV[1];
my $bam_TCPseq=$ARGV[2];
my $leaders=$ARGV[3]; #from cageR (cage_peaks_assignment_gtf_transcript_coords_for_TCPseq_v8.pl)
my $most_highly_expressed=$ARGV[4];
my $use_cage=$ARGV[5];  # unless($use_cage eq "exclude")

#excluding:
#genes where the TSS is downstream of the start codon
#genes without a detectable cage peak 
#genes that are annotated as protien_coding

#reads are assigned to leaders/CDS on a stand specific basis (including the RNA-seq) as follows:
#TCP-seq:   Assign fragments that overlap the TIS or stop codon to those feartures only.
#           Assign fragments that overlap the TSS or end of transcript to those features only (after start/stop codons)

#Flags:
#overlapping_gene: The longest transcript of this gene overlaps with the longest transcript of another gene
#leader_potentially_overlaps_upstream_gene: There is an upstream gene within 500nt of the start codon of this gene 
#gene_potentially_overlaps_downstream_leader: There is an downstream gene start codon withing 500nt of the 3' most position of this gene (in yeast this is the stop codon).

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
#second pass through the genome, find annotated start codons and setup transcript models for longest transcript of each gene
my %gene_start_codon_fwd;
my %gene_stop_codon_fwd;
my %gene_exons_fwd;
my %gene_start_codon_rev;
my %gene_stop_codon_rev;
my %gene_exons_rev;
my %gene_2_chr; #key = gene_id; value = chr

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

        if ($chr =~ /^chr(.*)/){
           $chr=$1; #if the chr name have a chr* prefix, remove it 
        }

        if ($gene_id && $transcript_id){

            if (exists ( $most_expressed_transcript{$gene_id} )){    #if the transcript is in the list of longest transcripts

                if ($transcript_id eq $most_expressed_transcript{$gene_id}){

                    $gene_2_chr{$gene_id}=$chr;

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
#open fasta for codon sequeces

my %fasta_sequences; #key = sequence_name, value=sequence
my $name;
open (FA, $fasta) || die "can't open $fasta";
while (<FA>){
    chomp;
    if (/^>([^\s]+)/){ #take header up to the first space
        $name=$1;
        if ($name =~ /^chr(.*)/){
           $name=$1; #if the chr name have a chr* prefix, remove it 
        }
    }else{
        $fasta_sequences{$name}.=$_;
    }
}
close(FA);

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#setup trascript models

my %gene_model_fwd;
my %start_coord_fwd;
my %stop_coord_fwd;

#5' is the #gene_model{$gene}{0}
#3' is the last coord
my %five_prime_most_coord_fwd;
my %three_prime_most_coord_fwd;

#gene_model{$gene}{0}=12345   #1st nt of start codon
#gene_model{$gene}{1}=12346
#gene_model{$gene}{2}=12347
#gene_model{$gene}{3}=12348
#gene_model{$gene}{4}=12349
#...
#to end of exons             #last nt of stop codon

for my $gene (keys %gene_exons_fwd){

    unless ($gene eq "ENSDARG00000077330" || $gene eq "ENSDARG00000102873" || $gene eq "ENSDARG00000089382"){ #probematic transcripts

    if ( (exists ($gene_start_codon_fwd{$gene})) && (exists ($gene_stop_codon_fwd{$gene})) ) { #restrict to genes with annotated start + stop codon

        my $model_pos=0;

        for my $exon_start (sort {$a <=> $b} keys %{ $gene_exons_fwd{$gene} } ){
            my $exon_end=$gene_exons_fwd{$gene}{$exon_start};

            #fwd exons are in ascending order
            # start(-1)-> 100958 100975
            #             101077 101715 <-end(+1)

            for ($exon_start .. $exon_end){
                $gene_model_fwd{$gene}{$model_pos}=$_;

                if ($model_pos == 0){                    #the default start to zero, can be extended by CAGE later
                    $five_prime_most_coord_fwd{$gene}=$model_pos;
                }

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
}

my %gene_model_rev;
my %start_coord_rev;
my %stop_coord_rev;

#5' is the #gene_model{$gene}{0}
#3' is the last coord
my %five_prime_most_coord_rev;
my %three_prime_most_coord_rev;

for my $gene (keys %gene_exons_rev){

    unless ($gene eq "ENSDARG00000077330" || $gene eq "ENSDARG00000102873" || $gene eq "ENSDARG00000089382"){ #probematic transcripts

    if ( (exists ($gene_start_codon_rev{$gene})) && (exists ($gene_stop_codon_rev{$gene})) ) { #restrict to genes with annotated start + stop codon

        my $model_pos=0;

        for my $exon_end (reverse (sort {$a <=> $b} keys %{ $gene_exons_rev{$gene} } )){
            my $exon_start=$gene_exons_rev{$gene}{$exon_end};

            #rev exons are sorted in decending order  
            #           447087 447794 <-start(+1)
            # end(-1)-> 446060 446254

            while ($exon_start >= $exon_end){
                $gene_model_rev{$gene}{$model_pos}=$exon_start;

                if ($model_pos == 0){                    #the default start to zero, can be extended by CAGE later
                    $five_prime_most_coord_rev{$gene}=$model_pos;
                }

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
}


#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#parse leaders
my %leader_positions_fwd;
my %leader_positions_rev;

#filters:
my %leader_length;
my %cage_peak_value;
my %overlaps_inframe_gene;
my %leader_overlaps_upstream;
my %gene_overlaps_downstream_leader;
my %gene_overlaps_nc;

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

        if ($chr =~ /^chr(.*)/){
           $chr=$1; #if the chr name have a chr* prefix, remove it 
        }

        if ($overlaps_inframe_gene eq "TRUE"){ $overlaps_inframe_gene{$gene}=1; }
        if ($leader_overlaps_upstream eq "TRUE"){ $leader_overlaps_upstream{$gene}=1; }
        if ($gene_overlaps_downstream_leader eq "TRUE"){ $gene_overlaps_downstream_leader{$gene}=1; }
        if ($gene_overlaps_nc eq "TRUE"){ $gene_overlaps_nc{$gene}=1; }

        unless ($leader_length eq "NaN"){  #only take genes that have a detectable cage peak

            unless ($leader_length < 0){  #exlude genes with negative leader sizes, as they cause problems with FPKM

                if ($dir eq "fwd"){ 
                    if (exists ($start_coord_fwd{$gene})){
                        unless ($highest_cage_peak >=  $gene_model_fwd{$gene}{$start_coord_fwd{$gene}}){  #exclude genes where the TSS is downstream of the start codon
                            $leader_positions_fwd{$gene}=$highest_cage_peak;
                            $cage_peak_value{$gene}=$count_at_highest_cage_peak;
                            $leader_length{$gene}=$leader_length;
                        } 
                    }
                }else{
                    if (exists ($start_coord_rev{$gene})){
                        unless ($highest_cage_peak <=  $gene_model_rev{$gene}{$start_coord_rev{$gene}}){  #exclude genes where the TSS is downstream of the start codon
                            $leader_positions_rev{$gene}=$highest_cage_peak;
                            $cage_peak_value{$gene}=$count_at_highest_cage_peak;
                            $leader_length{$gene}=$leader_length;
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
            my $five_prime_coord=$five_prime_most_coord_fwd{$gene}; 
            my $five_prime_pos=$gene_model_fwd{$gene}{$five_prime_coord};

            if ($leader_start >= $five_prime_pos){  #great, just find and save the coordinate
 
                for my $coord ($five_prime_coord .. $three_prime_coord){
                    my $pos=$gene_model_fwd{$gene}{$coord};       
                    if ($pos == $leader_start){
                       $five_prime_most_coord_fwd{$gene}=$coord;
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
            }
        }
    } 

    for my $gene (keys %leader_positions_rev){

        if (exists ($gene_2_chr{$gene})){  #to restict to protien_coding genes
            my $leader_start=$leader_positions_rev{$gene};
            my $three_prime_coord=$three_prime_most_coord_rev{$gene};
            my $three_prime_pos=$gene_model_rev{$gene}{$three_prime_coord};
            my $five_prime_coord=$five_prime_most_coord_rev{$gene};
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
#Transcript 5' coord   0                                  #1st nt of annotated transcript
#Transcript 3' coord   $three_prime_most_coord_???{$gene} #last nt of transcript
#Start codon coord     $start_coord_???{$gene}            #1st nt in start codon
#Stop codon coord      $stop_coord_???{$gene}             #1st nt in stop codon

#$gene_model_fwd{$gene}{$coord}==genomic position

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#loop though genes and assign leader, CDS and trailer regions to hashes for quick searching. (added start + stop codons).
my %leader_search_fwd;  #key1=chr, key2=pos, value=gene
my %leader_search_rev;  #key1=chr, key2=pos, value=gene
my %CDS_search_fwd;     #key1=chr, key2=pos, value=gene
my %CDS_search_rev;     #key1=chr, key2=pos, value=gene
my %trailer_search_fwd; #key1=chr, key2=pos, value=gene
my %trailer_search_rev; #key1=chr, key2=pos, value=gene

my %start_codons_search_fwd; #key1=chr, key2=pos, value=gene
my %start_codons_search_rev; #key1=chr, key2=pos, value=gene
my %stop_codons_search_fwd;  #key1=chr, key2=pos, value=gene
my %stop_codons_search_rev;  #key1=chr, key2=pos, value=gene

my %start_of_transcript_search_fwd;  #key1=chr, key2=pos, value=gene
my %start_of_transcript_search_rev;  #key1=chr, key2=pos, value=gene
my %end_of_transcript_search_fwd;    #key1=chr, key2=pos, value=gene
my %end_of_transcript_search_rev;    #key1=chr, key2=pos, value=gene
 
my %CDS_counts_TCPseq; 
my %leader_counts_TCPseq;
my %trailer_counts_TCPseq;
my %start_codon_counts_TCPseq;
my %stop_codon_counts_TCPseq;
my %start_of_transcript_counts_TCPseq;
my %end_of_transcript_counts_TCPseq;

for my $gene (keys %gene_model_fwd){
    my $chr=$gene_2_chr{$gene};
    my $start_coord=$start_coord_fwd{$gene};
    my $stop_coord=$stop_coord_fwd{$gene};
    my $transcript_start=$five_prime_most_coord_fwd{$gene};
    my $transcript_end=$three_prime_most_coord_fwd{$gene};

#    $CDS_counts_TCPseq{$gene}=0;
#    $leader_counts_TCPseq{$gene}=0;
#    $trailer_counts_TCPseq{$gene}=0;
#    $start_codon_counts_TCPseq{$gene}=0;
#    $stop_codon_counts_TCPseq{$gene}=0;

     unless (exists ($overlaps_inframe_gene{$gene}) || exists ( $leader_overlaps_upstream{$gene} ) || exists ( $gene_overlaps_downstream_leader{$gene} ) || exists ( $gene_overlaps_nc{$gene} ) ) {

        for my $coord (sort {$a <=> $b} keys %{ $gene_model_fwd{$gene} } ){
            my $pos=$gene_model_fwd{$gene}{$coord};

            #need to check thah the coord is passed the 5' most position (possibly updated by cage)
            if ($coord >= $five_prime_most_coord_fwd{$gene}){

                if ($coord == $start_coord || $coord == ($start_coord+1) || $coord == ($start_coord+2)){ 
                    $start_codons_search_fwd{$chr}{$pos}=$gene; 
                }elsif( $coord == $transcript_start){
                    $start_of_transcript_search_fwd{$chr}{$pos}=$gene;
                }

                elsif ($coord == $stop_coord || $coord == ($stop_coord+1) || $coord == ($stop_coord+2)){
                    $stop_codons_search_fwd{$chr}{$pos}=$gene;
                }elsif( $coord == $transcript_end){
                   $end_of_transcript_search_fwd{$chr}{$pos}=$gene;
                }
 
                elsif ($coord < $start_coord){
                    $leader_search_fwd{$chr}{$pos}=$gene;
                }elsif($coord <= ($stop_coord+2)){  #limit to the stop codon
                    $CDS_search_fwd{$chr}{$pos}=$gene;
                }elsif( $coord > ($stop_coord+2)){
                    $trailer_search_fwd{$chr}{$pos}=$gene;
                }
            }
        }
    }
}

for my $gene (keys %gene_model_rev){
    my $chr=$gene_2_chr{$gene};
    my $start_coord=$start_coord_rev{$gene};
    my $stop_coord=$stop_coord_rev{$gene};
    my $transcript_start=$five_prime_most_coord_rev{$gene};
    my $transcript_end=$three_prime_most_coord_rev{$gene};

#    $CDS_counts_TCPseq{$gene}=0;
#    $leader_counts_TCPseq{$gene}=0;
#    $trailer_counts_TCPseq{$gene}=0;
#    $start_codon_counts_TCPseq{$gene}=0;
#    $stop_codon_counts_TCPseq{$gene}=0;

     unless (exists ($overlaps_inframe_gene{$gene}) || exists ( $leader_overlaps_upstream{$gene} ) || exists ( $gene_overlaps_downstream_leader{$gene} ) || exists ( $gene_overlaps_nc{$gene} ) ) {
 
        for my $coord (sort {$a <=> $b} keys %{ $gene_model_rev{$gene} } ){
            my $pos=$gene_model_rev{$gene}{$coord};

            #need to check thah the coord is passed the 5' most position (possibly updated by cage)
            if ($coord >= $five_prime_most_coord_rev{$gene}){

                if ($coord == $start_coord || $coord == ($start_coord+1) || $coord == ($start_coord+2)){
                    $start_codons_search_rev{$chr}{$pos}=$gene; 
                }elsif( $coord == $transcript_start){
                    $start_of_transcript_search_rev{$chr}{$pos}=$gene;
                }

                #if ($coord == $start_coord){   $start_codons_search_rev{$chr}{$pos}=$gene; }
                #if ($coord == $start_coord+1){ $start_codons_search_rev{$chr}{$pos}=$gene; }
                #if ($coord == $start_coord+2){ $start_codons_search_rev{$chr}{$pos}=$gene; }

                elsif ($coord == $stop_coord || $coord == ($stop_coord+1) || $coord == ($stop_coord+2)){
                    $stop_codons_search_rev{$chr}{$pos}=$gene; 
                }elsif( $coord == $transcript_end){
                    $end_of_transcript_search_rev{$chr}{$pos}=$gene;
                }

                #if ($coord == $stop_coord){   $stop_codons_search_rev{$chr}{$pos}=$gene; }
                #if ($coord == $stop_coord+1){ $stop_codons_search_rev{$chr}{$pos}=$gene; }
                #if ($coord == $stop_coord+2){ $stop_codons_search_rev{$chr}{$pos}=$gene; }

                elsif ($coord < $start_coord){
                    $leader_search_rev{$chr}{$pos}=$gene;
                }elsif($coord <= ($stop_coord+2)){
                    $CDS_search_rev{$chr}{$pos}=$gene;
                }elsif( $coord > ($stop_coord+2)){
                    $trailer_search_rev{$chr}{$pos}=$gene;
                }
            }
        }
    }
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#open and store TCPseq counts
my $min_length=1000000;
my $max_length=0;

my $read_count=0;
my $matched_count=0;

open BAM,"samtools view $bam_TCPseq |";

while(<BAM>){

    next if(/^(\@)/);  ## skipping the header lines (if you used -h in the samools command)
    s/\n//;  s/\r//;  ## removing new line
    my @sam = split(/\t+/);  ## splitting SAM line into array

    my $leftMost=$sam[3]; #leftmost position of match 5' for fwd, 3' for rev
    my $flag=$sam[1];
    my $chr=$sam[2];
    my $mapq=$sam[4];
    my $cigar=$sam[5];
    my $seq=$sam[9];
    my $CDS_hit=0;
    my $leader_hit=0;
    my $trailer_hit=0;
    my $start_hit=0;
    my $stop_hit=0;
    my $start_of_transcript_hit=0;
    my $end_of_transcript_hit=0;

    if ($chr=~/chr(.*)/){ #ensembl chromosome names do not contain the "chr" prefix
        $chr=$1;
    }

    unless ($flag & 0x4){   #if aligned

        if ($mapq >= 10){     #mapping uniqnes filter

            $read_count++;
    
            if ($flag & 0x10){  #Reverse reads. Starting from the leftmost position parse the cigar and check if matching positions overlap leaders or cds's
    
                while ($cigar !~ /^$/){
                    if ($cigar =~ /^([0-9]+[MIDN])/){
                        my $cigar_part = $1;
                        if ($cigar_part =~ /(\d+)M/){   #alignment matching 
                            for my $pos ($leftMost .. ($leftMost+$1-1)){ #search though this position
    
                                #flag TCP-seq reads that overalap with TIS's or stop codons
                                if (exists ($start_codons_search_rev{$chr}{$pos})){
                                    $start_hit=$start_codons_search_rev{$chr}{$pos};
                                    last;
                                }elsif (exists ($stop_codons_search_rev{$chr}{$pos})){
                                    $stop_hit=$stop_codons_search_rev{$chr}{$pos};
                                    last;
                                }elsif (exists ($start_of_transcript_search_rev{$chr}{$pos})){
                                    $start_of_transcript_hit=$start_of_transcript_search_rev{$chr}{$pos};
                                    last;
                                }elsif (exists ($end_of_transcript_search_rev{$chr}{$pos})){
                                    $end_of_transcript_hit=$end_of_transcript_search_rev{$chr}{$pos};
                                    last;
                                }elsif (exists ( $leader_search_rev{$chr}{$pos})){                                 
                                    $leader_hit=$leader_search_rev{$chr}{$pos};
                                }elsif( exists($CDS_search_rev{$chr}{$pos})){
                                    $CDS_hit=$CDS_search_rev{$chr}{$pos};
                                }elsif( exists($trailer_search_rev{$chr}{$pos})){
                                    $trailer_hit=$trailer_search_rev{$chr}{$pos};
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
    
                                #flag TCP-seq reads that overalap with TIS's or stop codons
                                if (exists ($start_codons_search_fwd{$chr}{$pos})){
                                    $start_hit=$start_codons_search_fwd{$chr}{$pos};
                                    last;
                                }elsif (exists ($stop_codons_search_fwd{$chr}{$pos})){
                                    $stop_hit=$stop_codons_search_fwd{$chr}{$pos};
                                    last;
                                }elsif (exists ($start_of_transcript_search_fwd{$chr}{$pos})){
                                    $start_of_transcript_hit=$start_of_transcript_search_fwd{$chr}{$pos};
                                    last;
                                }elsif (exists ($end_of_transcript_search_fwd{$chr}{$pos})){
                                    $end_of_transcript_hit=$end_of_transcript_search_fwd{$chr}{$pos};
                                    last;
                                }elsif (exists ( $leader_search_fwd{$chr}{$pos})){                                 
                                    $leader_hit=$leader_search_fwd{$chr}{$pos};
                                }elsif( exists($CDS_search_fwd{$chr}{$pos})){
                                    $CDS_hit=$CDS_search_fwd{$chr}{$pos};
                                }elsif( exists($trailer_search_fwd{$chr}{$pos})){
                                    $trailer_hit=$trailer_search_fwd{$chr}{$pos};
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
    
            #reads should not be able to be assigned to multiple features, becuase I am seperating the reads that overlap start and stop codons (i.e. the boundaries between features)
    
            #change to length 
            my $length=length($seq);
            #store longest and smallest lengths, for printing
            if ($length < $min_length){ $min_length=$length; }
            if ($length > $max_length){ $max_length=$length; }
    
            if ($start_hit){
                $start_codon_counts_TCPseq{$length}++;
            }elsif($start_of_transcript_hit){
                $start_of_transcript_counts_TCPseq{$length}++;
            }elsif ($stop_hit){
                $stop_codon_counts_TCPseq{$length}++;
            }elsif($end_of_transcript_hit){
                $end_of_transcript_counts_TCPseq{$length}++;
            }else{
                if ($leader_hit){ #now contains gene name   
                    $leader_counts_TCPseq{$length}++;
                }
                if ($CDS_hit){ #now contains gene name
                     $CDS_counts_TCPseq{$length}++;
                }
                if ($trailer_hit){ #now contains gene name
                    $trailer_counts_TCPseq{$length}++;
                }
            }
        }
    }
}
close (BAM);

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#output

#header
print "#read_length,count,feature_type\n";

for my $read_length ($min_length .. $max_length){
    if (exists ($start_of_transcript_counts_TCPseq{$read_length})){
        print "$read_length,$start_of_transcript_counts_TCPseq{$read_length},start_of_transcript\n";
    }else{
        print "$read_length,0,start_of_transcript\n";
    }
}

for my $read_length ($min_length .. $max_length){
    if (exists ($leader_counts_TCPseq{$read_length})){
        print "$read_length,$leader_counts_TCPseq{$read_length},leader5\n";
    }else{
        print "$read_length,0,leader5\n";
    }
}

for my $read_length ($min_length .. $max_length){
    if (exists ($start_codon_counts_TCPseq{$read_length})){
        print "$read_length,$start_codon_counts_TCPseq{$read_length},start_codon\n";
    }else{
        print "$read_length,0,start_codon\n";
    } 
}

for my $read_length ($min_length .. $max_length){
    if (exists ($CDS_counts_TCPseq{$read_length})){
        print "$read_length,$CDS_counts_TCPseq{$read_length},CDS\n";
    }else{
        print "$read_length,0,CDS\n";
    }
}

for my $read_length ($min_length .. $max_length){
    if (exists ($stop_codon_counts_TCPseq{$read_length})){
        print "$read_length,$stop_codon_counts_TCPseq{$read_length},stop_codon\n";
    }else{
        print "$read_length,0,stop_codon\n";
    }
}

for my $read_length ($min_length .. $max_length){
    if (exists ($trailer_counts_TCPseq{$read_length})){
        print "$read_length,$trailer_counts_TCPseq{$read_length},trailer3\n";
    }else{
        print "$read_length,0,trailer3\n";
    }
}

for my $read_length ($min_length .. $max_length){
    if (exists ($end_of_transcript_counts_TCPseq{$read_length})){
        print "$read_length,$end_of_transcript_counts_TCPseq{$read_length},end_of_transcript\n";
    }else{
        print "$read_length,0,end_of_transcript\n";
    }
}

exit;

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
