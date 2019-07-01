#!/usr/bin/perl -w
use strict;

#14/01/19
#Script to count, SSU and RFP up and downstream of uORFs
#Upstream = TSS to uORF start-1
#Downstream = uORF start to TIS-1
#Output kozak scores and Up/downstream distances
#And stop codon, and distance from stop to TIS?

my $inGtf=$ARGV[0]; 
my $fasta=$ARGV[1];
my $most_highly_expressed=$ARGV[2];
my $leaders=$ARGV[3];
my $bam_SSU=$ARGV[4];  
my $bam_LSU=$ARGV[5];  
my $out_file=$ARGV[6];

#updated to calculate coverage over stop codons 

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#Zebrafish Zozak position frequency matrix (PFM) from Grzegorski et. al, 2015                   
#Position   -4  -3  -2  -1  4   5
#A  35  62  39  28  24  27
#C  32  5   23  36  12  42
#G  19  28  17  27  46  16
#T  14  5   21  10  17  15

my %raw_kozak =
   (
       "A" => [ 35, 62, 39, 28, 24, 27 ],
       "C" => [ 32, 5, 23, 36, 12, 42 ],
       "G" => [ 19, 28, 17, 27, 46, 16 ],
       "T" => [ 14, 5, 21, 10, 17, 15 ],
   );

my %PWM; #key1: position, key2: base, value: weight
for my $pos (0 .. 5){

    #0 == -4 #nucleotide position in relation to start codon
    #1 == -3
    #2 == -2
    #3 == -1
    #4 == +4
    #5 == +5

    my $pos_sum=0;
    my $pwm_sum=0;
    for my $base (keys %raw_kozak){ #sum the nucleotide frequencies per position
        $pos_sum+=$raw_kozak{$base}[$pos];
    }

    for my $base(keys %raw_kozak){ #score the PWM
        my $psudo_count= sqrt($pos_sum);
        my $background_probability=0.25; #no base preference
        my $pwm=&log2( ($raw_kozak{$base}[$pos] + $psudo_count * $background_probability) / ($pos_sum + $psudo_count * $background_probability));
        $PWM{$pos}{$base}=$pwm;
        $pwm_sum+=$pwm;
    }

    $PWM{$pos}{"N"}=($pwm_sum/4); #set "N" to be equal to the column mean. For genes with short leaders, missing upstream positions 
} 

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

my %overlaps_inframe_gene;
my %leader_overlaps_upstream;
my %gene_overlaps_downstream_leader;
my %gene_overlaps_ncRNA;
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
                    $gene_overlaps_ncRNA{$gene_id}=0;
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
        my $gene_overlaps_ncRNA=$b[7];
        my $highest_cage_peak=$b[8];
        my $count_at_highest_cage_peak=$b[9];
        my $leader_length=$b[10];

        if ($overlaps_inframe_gene eq "TRUE"){ $overlaps_inframe_gene{$gene}=1; }
        if ($leader_overlaps_upstream eq "TRUE"){ $leader_overlaps_upstream{$gene}=1; }
        if ($gene_overlaps_downstream_leader eq "TRUE"){ $gene_overlaps_downstream_leader{$gene}=1; }
        if ($gene_overlaps_ncRNA eq "TRUE"){ $gene_overlaps_ncRNA{$gene}=1; }

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

my %leader_search_fwd;
my %leader_search_rev;
my %start_codons_search_fwd; #key1=chr, key2=pos, value=gene
my %start_codons_search_rev; #key1=chr, key2=pos, value=gene
my %stop_codons_search_fwd; #key1=chr, key2=pos, value=gene
my %stop_codons_search_rev; #key1=chr, key2=pos, value=gene

my %CDS_search_fwd; #key1=chr, key2=pos, value=gene
my %CDS_search_rev; #key1=chr, key2=pos, value=gene
my %CDS_counts_SSU;
my %CDS_counts_LSU;
 
my %CDS_length;

for my $gene (keys %gene_model_fwd){
    my $chr=$gene_2_chr{$gene};
    my $start_coord=$start_coord_fwd{$gene};
    my $stop_coord=$stop_coord_fwd{$gene};

    $CDS_counts_SSU{$gene}=0;
    $CDS_counts_LSU{$gene}=0;

    $CDS_length{$gene}=(($stop_coord+2)-$start_coord)+1;

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

            if ($coord == $start_coord){   $start_codons_search_fwd{$chr}{$pos}=$gene; } 
            if ($coord == $start_coord+1){ $start_codons_search_fwd{$chr}{$pos}=$gene; }
            if ($coord == $start_coord+2){ $start_codons_search_fwd{$chr}{$pos}=$gene; }

            if ($coord == $stop_coord){   $stop_codons_search_fwd{$chr}{$pos}=$gene; }
            if ($coord == $stop_coord+1){ $stop_codons_search_fwd{$chr}{$pos}=$gene; }
            if ($coord == $stop_coord+2){ $stop_codons_search_fwd{$chr}{$pos}=$gene; }
           
            if ($coord < $start_coord){
#                $leader_search_fwd{$chr}{$pos}=$gene;
            }elsif($coord <= ($stop_coord+2)){  #limit to the stop codon
                if (exists($CDS_search_fwd{$chr}{$pos})){
                    $CDS_search_fwd{$chr}{$pos}.="_".$gene;
                }else{
                    $CDS_search_fwd{$chr}{$pos}=$gene;
                }
            }
        }
    }
}

for my $gene (keys %gene_model_rev){
    my $chr=$gene_2_chr{$gene};
    my $start_coord=$start_coord_rev{$gene};
    my $stop_coord=$stop_coord_rev{$gene};

    $CDS_counts_SSU{$gene}=0;
    $CDS_counts_LSU{$gene}=0;

    $CDS_length{$gene}=(($stop_coord+2)-$start_coord)+1;

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

            if ($coord == $start_coord){   $start_codons_search_rev{$chr}{$pos}=$gene; }
            if ($coord == $start_coord+1){ $start_codons_search_rev{$chr}{$pos}=$gene; }
            if ($coord == $start_coord+2){ $start_codons_search_rev{$chr}{$pos}=$gene; }

            if ($coord == $stop_coord){   $stop_codons_search_rev{$chr}{$pos}=$gene; }
            if ($coord == $stop_coord+1){ $stop_codons_search_rev{$chr}{$pos}=$gene; }
            if ($coord == $stop_coord+2){ $stop_codons_search_rev{$chr}{$pos}=$gene; }

            if ($coord < $start_coord){
#                $leader_search_rev{$chr}{$pos}=$gene;
            }elsif($coord <= ($stop_coord+2)){  #limit to the stop codon
                if (exists($CDS_search_rev{$chr}{$pos})){
                    $CDS_search_rev{$chr}{$pos}.="_".$gene;
                }else{
                    $CDS_search_rev{$chr}{$pos}=$gene;
                }
            }
        }
    }
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#assign uORFs

#output
#start codon
#stop codon
#kozak
#distance to TSS
#distance to TIS
#gene CDS FPKM #I can take this from the genes matricies

my %uORF_start_coord;
my %uORF_stop_coord;
my %uORF_start_codon;
my %uORF_stop_codon;
my %uORF_kozak;
my %uORF_to_TSS;
my %uORF_to_TIS;
my %uORF_stop_to_TIS;
my %uORF_length;
my %gene_to_uORF;

#my @uORF_break_distribution; #location or the uORF start as a proportion of leader length
#my %uORF_to_distrubution;

for my $gene (keys %gene_model_fwd){
    my $chr=$gene_2_chr{$gene};
    my $start_coord=$start_coord_fwd{$gene};
    my $stop_coord=$stop_coord_fwd{$gene};
    my $five_prime_coord=$five_prime_most_coord_fwd{$gene};
    my $three_prime_coord=$three_prime_most_coord_fwd{$gene};

    my $uORF_id_count=0;

    #This methods locates all uORFs in the leader region
    for my $coord (sort {$a <=> $b} keys %{ $gene_model_fwd{$gene} } ){

        if (($coord >= ($five_prime_coord+2)) && ($coord <= ($start_coord-1))){  
            my @TIS;
            for(-2 .. 0){
                if (exists ($gene_model_fwd{$gene}{ ($coord+$_) } )){
                    push (@TIS, substr($fasta_sequences{$chr}, ($gene_model_fwd{$gene}{ ($coord+$_) } -1), 1) );
                }else{
                    push (@TIS, "N");
                }
            }
            my $seq_TIS=join("", @TIS);

            if ( ($seq_TIS =~ /AT\w/) || ($seq_TIS =~ /A\wG/) || ($seq_TIS =~ /\wTG/) ){ #near/canoical
            #if ( $seq_TIS eq "ATG" ){

                my $uORF_start=$gene_model_fwd{$gene}{($coord-2)};
                my $uORF_start_coord=$coord-2;  #the 1st nt of the start coord
                my $search_coord=$coord+3;
                while ($search_coord < $three_prime_coord-2){ #-2 to allow enough room for the last codon

                    my @STOP;
                    for(-2 .. 0){
                        if (exists ($gene_model_fwd{$gene}{ ($search_coord+$_) } )){
                            push (@STOP, substr($fasta_sequences{$chr}, ($gene_model_fwd{$gene}{ ($search_coord+$_) } -1), 1) );
                        }else{
                            push (@STOP, "N");
                        }
                    }
                    my $seq_STOP=join("", @STOP);
                    if (($seq_STOP eq "TAA") || ($seq_STOP eq "TAG") || ($seq_STOP eq "TGA") ){

                        #exclude CDS extensions
                        unless ($search_coord eq ($stop_coord+2)){
                   
                            my $kozak_score=0; #calculate kozak score here
                            my @up;
                            my @down;
                            for (-6 .. -3){
                                if (exists ($gene_model_fwd{$gene}{ ($uORF_start_coord+$_) } )){
                                    if (exists ($five_prime_most_coord_fwd{$gene})){
                                        if (($coord+$_) >= ($five_prime_most_coord_fwd{$gene}-1)){ #check we're not extending behyond leader
                                            push (@up, substr($fasta_sequences{$chr}, ($gene_model_fwd{$gene}{ ($uORF_start_coord+$_) } -1), 1) );
                                        }else{
                                           push (@up, "N");
                                        }
                                    }
                                }else{
                                    push (@up, "N");
                                }
                            }
                            for (1 .. 2){
                                if (exists ($gene_model_fwd{$gene}{ ($uORF_start_coord+$_) } )){
                                    push (@down, substr($fasta_sequences{$chr}, ($gene_model_fwd{$gene}{ ($uORF_start_coord+$_) } -1), 1) );
                                }else{
                                    push (@down, "N");
                                }
                            }

                            my $seq_up=join("", @up);
                            my $seq_down=join("", @down);
                            $kozak_score=&score_kozak($seq_up,$seq_down);
                            
                            my $uORF_id=$gene."_".$uORF_id_count;
                            my $uORF_stop=$gene_model_fwd{$gene}{$search_coord};  #the last nt of the stop codon
                            my $leader_length=$start_coord-$five_prime_coord;
                            my $distance_from_TSS=($uORF_start_coord-1)-$five_prime_coord;
                            my $distance_to_TIS=$start_coord-$uORF_start_coord;
                            my $distance_stop_to_TIS=$start_coord-$search_coord;
                            my $proportional_location=$distance_from_TSS/$leader_length;

                            $uORF_start_coord{$gene}{$uORF_id}=$uORF_start_coord;
                            $uORF_stop_coord{$gene}{$uORF_id}=$search_coord;
                            $uORF_start_codon{$uORF_id}=$seq_TIS;
                            $uORF_stop_codon{$uORF_id}=$seq_STOP;
                            $uORF_kozak{$uORF_id}=$kozak_score;
                            $uORF_to_TSS{$uORF_id}=$distance_from_TSS;
                            $uORF_to_TIS{$uORF_id}=$distance_to_TIS;
                            $uORF_stop_to_TIS{$uORF_id}=$distance_stop_to_TIS;
                            $uORF_length{$uORF_id}=$search_coord-$uORF_start_coord;
                            $gene_to_uORF{$gene}{$uORF_id}=1;
                            #$uORF_to_distrubution{$uORF_id}=$proportional_location;
                            $uORF_id_count++;
                            last; #exit the stop codon search
                        }  
                    }
                    $search_coord+=3;
                }
            }
        }
    }
}


for my $gene (keys %gene_model_rev){
    my $chr=$gene_2_chr{$gene};
    my $start_coord=$start_coord_rev{$gene};
    my $stop_coord=$stop_coord_rev{$gene};
    my $five_prime_coord=$five_prime_most_coord_rev{$gene};
    my $three_prime_coord=$three_prime_most_coord_rev{$gene};

    my $uORF_id_count=0;

    #This methods locates all uORFs in the leader region
    for my $coord (sort {$a <=> $b} keys %{ $gene_model_rev{$gene} } ){
 
        if (($coord >= ($five_prime_coord+2)) && ($coord <= ($start_coord-1))){  #(+2) 
            my @TIS;
            for(-2 .. 0){
                if (exists ($gene_model_rev{$gene}{ ($coord+$_) } )){
                    push (@TIS, substr($fasta_sequences{$chr}, ($gene_model_rev{$gene}{ ($coord+$_) } -1), 1) );
                }else{
                     push (@TIS, "N");
                }
            }
            my $seq_TIS=join("", @TIS);
            $seq_TIS=~tr/ACGTacgt/TGCAtgca/;

            if ( ($seq_TIS =~ /AT\w/) || ($seq_TIS =~ /A\wG/) || ($seq_TIS =~ /\wTG/) ){ #near/canoical
            #if ( $seq_TIS eq "ATG" ){

                my $uORF_start=$gene_model_rev{$gene}{($coord-2)};
                my $uORF_start_coord=$coord-2;
                my $search_coord=$coord+3;
                while ($search_coord < $three_prime_coord-2){ #-2 to allow enough room for the last codon

                    my @STOP;
                    for(-2 .. 0){
                        if (exists ($gene_model_rev{$gene}{ ($search_coord+$_) } )){
                            push (@STOP, substr($fasta_sequences{$chr}, ($gene_model_rev{$gene}{ ($search_coord+$_) } -1), 1) );
                        }else{
                            push (@STOP, "N");
                        }
                    }
                    my $seq_STOP=join("", @STOP);
                    $seq_STOP=~tr/ACGTacgt/TGCAtgca/;
                    if (($seq_STOP eq "TAA") || ($seq_STOP eq "TAG") || ($seq_STOP eq "TGA") ){
 
                        #exclude extensions to CDSs
                        unless ($search_coord eq ($stop_coord+2)){

                            my $kozak_score=0; #calculate kozak score here
                            my @up;
                            my @down;
                            for (-6 .. -3){
                                if (exists ($gene_model_rev{$gene}{ ($uORF_start_coord+$_) } )){
                                    if (exists ($five_prime_most_coord_rev{$gene})){
                                        if (($coord+$_) >= ($five_prime_most_coord_rev{$gene}-1)){ #check we're not extending behyond leader
                                            push (@up, substr($fasta_sequences{$chr}, ($gene_model_rev{$gene}{ ($uORF_start_coord+$_) } -1), 1) );
                                        }else{
                                            push (@up, "N");
                                        }
                                    }
                                }else{
                                    push (@up, "N");
                                }
                            }
                            for (1 .. 2){
                                if (exists ($gene_model_rev{$gene}{ ($uORF_start_coord+$_) } )){
                                    push (@down, substr($fasta_sequences{$chr}, ($gene_model_rev{$gene}{ ($uORF_start_coord+$_) } -1), 1) );
                                }else{
                                    push (@down, "N");
                                }
                            }

                            my $seq_up=join("", @up);
                            my $seq_down=join("", @down);
                            $seq_up=~tr/ACGTacgt/TGCAtgca/;
                            $seq_down=~tr/ACGTacgt/TGCAtgca/;
                            $kozak_score=&score_kozak($seq_up,$seq_down);

                            my $uORF_id=$gene."_".$uORF_id_count;
                            my $uORF_stop=$gene_model_rev{$gene}{$search_coord};  #the last nt of the stop codon
                            my $leader_length=$start_coord-$five_prime_coord;
                            my $distance_from_TSS=($uORF_start_coord-1)-$five_prime_coord;
                            my $distance_to_TIS=$start_coord-$uORF_start_coord;
                            my $distance_stop_to_TIS=$start_coord-$search_coord;
                            my $proportional_location=$distance_from_TSS/$leader_length;
                            $uORF_kozak{$uORF_id}=$kozak_score;
                            $uORF_to_TSS{$uORF_id}=$distance_from_TSS;
                            $uORF_to_TIS{$uORF_id}=$distance_to_TIS;
                            $uORF_stop_to_TIS{$uORF_id}=$distance_stop_to_TIS;
                            $uORF_length{$uORF_id}=$search_coord-$uORF_start_coord;
                            $uORF_start_coord{$gene}{$uORF_id}=$uORF_start_coord;
                            $uORF_stop_coord{$gene}{$uORF_id}=$search_coord;
                            $uORF_start_codon{$uORF_id}=$seq_TIS;
                            $uORF_stop_codon{$uORF_id}=$seq_STOP;
                            $gene_to_uORF{$gene}{$uORF_id}=1;
                            #$uORF_to_distrubution{$uORF_id}=$proportional_location; 
                            $uORF_id_count++;
                            last; #exit the stop codon search
                        } 
                    }
                    $search_coord+=3;
                }
            }
        }
    }
}


#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#reorganised so that downstream regions are between the uORF stop codon and TIS regions
#For the better measuremment of the effect of uORF stop codons

#assign genes to catgories with leaders containing 0,1,2 or 3 ATG uORF (more than 30nt appart)
#setup coverage search hashes and initialize positional counts

my %uORF_TSS_TIS_SSU_coverage_by_positon; 
my %uORF_TSS_TIS_LSU_coverage_by_positon;
my %uORF_TSS_TIS_window_search_fwd; #update to allow overlapping
my %uORF_TSS_TIS_window_search_rev;

my %uORF_stop_codon_SSU_coverage;
my %uORF_stop_codon_LSU_coverage;
my %uORF_stop_codon_search_fwd;
my %uORF_stop_codon_search_rev;

#my %gene4background;
#my %bg_TSS_TIS_SSU_coverage_by_positon;
#my %bg_TSS_TIS_LSU_coverage_by_positon;
#my %bg_TSS_TIS_window_search_fwd;
#my %bg_TSS_TIS_window_search_rev;

for my $gene (keys %gene_model_fwd){
    my $chr=$gene_2_chr{$gene};
    my $gene_start_coord=$start_coord_fwd{$gene};
    my $gene_five_prime_coord=$five_prime_most_coord_fwd{$gene};

    if (exists ($gene_to_uORF{$gene})){ 

        for my $uORF_id (sort keys %{ $gene_to_uORF{$gene} }){
            
            if ( (exists ($uORF_start_coord{$gene}{$uORF_id})) && (exists ($uORF_stop_coord{$gene}{$uORF_id})) ){

                if ($uORF_stop_coord{$gene}{$uORF_id} < $gene_start_coord){ #add an additional constraint to only include uORFs that end before the CDS start codon.

                    my $uORF_start_coord=$uORF_start_coord{$gene}{$uORF_id};
                    my $uORF_stop_coord=$uORF_stop_coord{$gene}{$uORF_id};

                    #my $distribution=$uORF_to_distrubution{$uORF_id}; #the proportional distance from TSS to uORF TIS. 
                    #push (@uORF_break_distribution, $distribution);  #add to the distibution for background selction
    
                    $uORF_TSS_TIS_SSU_coverage_by_positon{$uORF_id}{"tss_to_uORF_start"}=0; #initialise
                    $uORF_TSS_TIS_LSU_coverage_by_positon{$uORF_id}{"tss_to_uORF_start"}=0; #initialise
    
                    for my $coord ($gene_five_prime_coord .. ($uORF_start_coord-1)){
                        my $pos=$gene_model_fwd{$gene}{$coord};
                        if (exists ($uORF_TSS_TIS_window_search_fwd{$chr}{$pos})){
                            $uORF_TSS_TIS_window_search_fwd{$chr}{$pos}.="-".$uORF_id."~tss_to_uORF_start";  
                        }else{
                            $uORF_TSS_TIS_window_search_fwd{$chr}{$pos}=$uORF_id."~tss_to_uORF_start"; 
                        }
                    }
    
                    $uORF_TSS_TIS_SSU_coverage_by_positon{$uORF_id}{"uORF_stop_to_tis"}=0; #initialise
                    $uORF_TSS_TIS_LSU_coverage_by_positon{$uORF_id}{"uORF_stop_to_tis"}=0; #initialise
    
                    for my $coord (($uORF_stop_coord+1) .. ($gene_start_coord-1)){
                        my $pos=$gene_model_fwd{$gene}{$coord};
                        if (exists ($uORF_TSS_TIS_window_search_fwd{$chr}{$pos})){
                            $uORF_TSS_TIS_window_search_fwd{$chr}{$pos}.="-".$uORF_id."~uORF_stop_to_tis";
                        }else{
                            $uORF_TSS_TIS_window_search_fwd{$chr}{$pos}=$uORF_id."~uORF_stop_to_tis";                    
                        }
                    }

                    $uORF_stop_codon_SSU_coverage{$uORF_id}=0;
                    $uORF_stop_codon_LSU_coverage{$uORF_id}=0;
      
                    for my $coord (($uORF_stop_coord-2) .. $uORF_stop_coord){
                        my $pos=$gene_model_fwd{$gene}{$coord};
 
                        #print "uORF_stop,fwd,$gene,$uORF_id,$chr,$pos\n";
                        
                        if (exists ($uORF_stop_codon_search_fwd{$chr}{$pos})){
                            $uORF_stop_codon_search_fwd{$chr}{$pos}.="-".$uORF_id;
                        }else{
                            $uORF_stop_codon_search_fwd{$chr}{$pos}=$uORF_id;
                        }
                    }
                }
            }
        }
    }else{ #build a background distibution for the TSS to TIS region of transcripts without a uORF

#            $gene4background{$gene}=1;

#           my $relational_position=1;
#           for my $coord ($gene_five_prime_coord .. ($gene_start_coord-1)){
#               $bg_TSS_TIS_SSU_coverage_by_positon{$gene}{$relational_position}=0; #initialise
#               $bg_TSS_TIS_LSU_coverage_by_positon{$gene}{$relational_position}=0; #initialise
#               my $pos=$gene_model_fwd{$gene}{$coord};
#               $bg_TSS_TIS_window_search_fwd{$chr}{$pos}=$gene."~".$relational_position; 
#               $relational_position++;
#           }
    }
}

for my $gene (keys %gene_model_rev){
    my $chr=$gene_2_chr{$gene};
    my $gene_start_coord=$start_coord_rev{$gene};
    my $gene_five_prime_coord=$five_prime_most_coord_rev{$gene};

    if (exists ($gene_to_uORF{$gene})){ #at least one ATG uORF  

        for my $uORF_id (sort keys %{ $gene_to_uORF{$gene} }){

            if ( (exists ($uORF_start_coord{$gene}{$uORF_id})) && (exists ($uORF_stop_coord{$gene}{$uORF_id})) ){

                if ($uORF_stop_coord{$gene}{$uORF_id} < $gene_start_coord){ #add an additional constraint to only include uORFs that end before the CDS start codon.

                    my $uORF_start_coord=$uORF_start_coord{$gene}{$uORF_id};
                    my $uORF_stop_coord=$uORF_stop_coord{$gene}{$uORF_id};

                    #my $distribution=$uORF_to_distrubution{$uORF_id}; #the proportional distance from TSS to uORF TIS. 
                    #push (@uORF_break_distribution, $distribution);  #add to the distibution for background selction
    
                    $uORF_TSS_TIS_SSU_coverage_by_positon{$uORF_id}{"tss_to_uORF_start"}=0; #initialise
                    $uORF_TSS_TIS_LSU_coverage_by_positon{$uORF_id}{"tss_to_uORF_start"}=0; #initialise
    
                    for my $coord ($gene_five_prime_coord .. ($uORF_start_coord-1)){
                        my $pos=$gene_model_rev{$gene}{$coord};
    
                        if (exists ($uORF_TSS_TIS_window_search_rev{$chr}{$pos})){
                            $uORF_TSS_TIS_window_search_rev{$chr}{$pos}.="-".$uORF_id."~tss_to_uORF_start";
                        }else{
                            $uORF_TSS_TIS_window_search_rev{$chr}{$pos}=$uORF_id."~tss_to_uORF_start"; 
                        }
                    }
    
                    $uORF_TSS_TIS_SSU_coverage_by_positon{$uORF_id}{"uORF_stop_to_tis"}=0; #initialise
                    $uORF_TSS_TIS_LSU_coverage_by_positon{$uORF_id}{"uORF_stop_to_tis"}=0; #initialise
    
                    for my $coord (($uORF_stop_coord+1) .. ($gene_start_coord-1)){
                        my $pos=$gene_model_rev{$gene}{$coord};
                        if (exists ($uORF_TSS_TIS_window_search_rev{$chr}{$pos})){
                            $uORF_TSS_TIS_window_search_rev{$chr}{$pos}.="-".$uORF_id."~uORF_stop_to_tis";
                        }else{
                            $uORF_TSS_TIS_window_search_rev{$chr}{$pos}=$uORF_id."~uORF_stop_to_tis";
                        }
                    }

                    $uORF_stop_codon_SSU_coverage{$uORF_id}=0;
                    $uORF_stop_codon_LSU_coverage{$uORF_id}=0;

                    for my $coord (($uORF_stop_coord-2) .. $uORF_stop_coord){
                        my $pos=$gene_model_rev{$gene}{$coord};

                        #print "uORF_stop,rev,$gene,$uORF_id,$chr,$pos\n";

                        if (exists ($uORF_stop_codon_search_rev{$chr}{$pos})){
                            $uORF_stop_codon_search_rev{$chr}{$pos}.="-".$uORF_id;
                        }else{
                            $uORF_stop_codon_search_rev{$chr}{$pos}=$uORF_id;
                        }
                    }
                }
            }
        }
    }else{ #build a background distibution for the TSS to TIS region of transcripts without a uORF

#            $gene4background{$gene}=1;

#           my $relational_position=1;
#           for my $coord ($gene_five_prime_coord .. ($gene_start_coord-1)){
#               $bg_TSS_TIS_SSU_coverage_by_positon{$gene}{$relational_position}=0; #initialise
#               $bg_TSS_TIS_LSU_coverage_by_positon{$gene}{$relational_position}=0; #initialise
#               my $pos=$gene_model_rev{$gene}{$coord};
#
#              $bg_TSS_TIS_window_search_rev{$chr}{$pos}=$gene."~".$relational_position;                        
#               $relational_position++;
#           }
    }
}



#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#assign psudo uORF split to background distibutions (sampled from proportional uORF start positions relaitive to TSS)

#for my $gene (keys %gene_model_fwd){
#    my $chr=$gene_2_chr{$gene};
#    my $gene_start_coord=$start_coord_fwd{$gene};
#    my $gene_five_prime_coord=$five_prime_most_coord_fwd{$gene};
#    my $leader_length=$gene_start_coord-$gene_five_prime_coord;
# 
#    if (exists ($gene4background{$gene})){
#
#        if ($leader_length >= 100){
#
#            my $resampled_proportion=$uORF_break_distribution[rand @uORF_break_distribution];
#            my $resampled_distance_from_TSS=int($leader_length*$resampled_proportion);
#            my $psudo_uORF_start=$gene_five_prime_coord+$resampled_distance_from_TSS;
#
#
#            if ( ($psudo_uORF_start > ($gene_five_prime_coord+50)) && ($psudo_uORF_start < ($gene_start_coord-50) ) ){
#
#                my $relational_position=1;
#                for my $coord ($gene_five_prime_coord .. ($psudo_uORF_start-1)){
#                    $bg_TSS_TIS_SSU_coverage_by_positon{$gene}{"tss_to_psudo_start"}{$relational_position}=0; #initialise
#                    $bg_TSS_TIS_LSU_coverage_by_positon{$gene}{"tss_to_psudo_start"}{$relational_position}=0; #initialise
#                    my $pos=$gene_model_fwd{$gene}{$coord};
#                    $bg_TSS_TIS_window_search_fwd{$chr}{$pos}=$gene."~tss_to_psudo_start~".$relational_position;
#                    $relational_position++;
#                }
#
#                $relational_position=1;
#                for my $coord ($psudo_uORF_start .. ($gene_start_coord-1)){
#                    $bg_TSS_TIS_SSU_coverage_by_positon{$gene}{"psudo_start_to_tis"}{$relational_position}=0; #initialise
#                    $bg_TSS_TIS_LSU_coverage_by_positon{$gene}{"psudo_start_to_tis"}{$relational_position}=0; #initialise
#                    my $pos=$gene_model_fwd{$gene}{$coord};
#                    $bg_TSS_TIS_window_search_fwd{$chr}{$pos}=$gene."~psudo_start_to_tis~".$relational_position;
#                    $relational_position++;
#                }
#            }
#        }
#    }   
#}
#
#
#for my $gene (keys %gene_model_rev){
#    my $chr=$gene_2_chr{$gene};
#    my $gene_start_coord=$start_coord_rev{$gene};
#    my $gene_five_prime_coord=$five_prime_most_coord_rev{$gene};
#    my $leader_length=$gene_start_coord-$gene_five_prime_coord;
#
#    if (exists ($gene4background{$gene})){
#
#        if ($leader_length >= 100){
#
#            my $resampled_proportion=$uORF_break_distribution[rand @uORF_break_distribution];
#            my $resampled_distance_from_TSS=int($leader_length*$resampled_proportion);
#            my $psudo_uORF_start=$gene_five_prime_coord+$resampled_distance_from_TSS;
#
##           print "$gene,$resampled_proportion,$resampled_distance_from_TSS,$psudo_uORF_start\n";
#
#            if ( ($psudo_uORF_start > ($gene_five_prime_coord+50)) && ($psudo_uORF_start < ($gene_start_coord-50) ) ){
#
##                print "$gene,background fits,$resampled_proportion,$resampled_distance_from_TSS,$gene_five_prime_coord,$psudo_uORF_start,$gene_start_coord\n";
#
#                my $relational_position=1;
#                for my $coord ($gene_five_prime_coord .. ($psudo_uORF_start-1)){
#                    $bg_TSS_TIS_SSU_coverage_by_positon{$gene}{"tss_to_psudo_start"}{$relational_position}=0; #initialise
#                    $bg_TSS_TIS_LSU_coverage_by_positon{$gene}{"tss_to_psudo_start"}{$relational_position}=0; #initialise
#                    my $pos=$gene_model_rev{$gene}{$coord};
#                    $bg_TSS_TIS_window_search_rev{$chr}{$pos}=$gene."~tss_to_psudo_start~".$relational_position;
#                    $relational_position++;
#                }
#    
#                $relational_position=1;
#                for my $coord ($psudo_uORF_start .. ($gene_start_coord-1)){
#                    $bg_TSS_TIS_SSU_coverage_by_positon{$gene}{"psudo_start_to_tis"}{$relational_position}=0; #initialise
#                    $bg_TSS_TIS_LSU_coverage_by_positon{$gene}{"psudo_start_to_tis"}{$relational_position}=0; #initialise
#                    my $pos=$gene_model_rev{$gene}{$coord};
#                    $bg_TSS_TIS_window_search_rev{$chr}{$pos}=$gene."~psudo_start_to_tis~".$relational_position;
#                    $relational_position++;
#                }
#
##            }else{
##                print "$gene,background out of bounds,$resampled_proportion,$resampled_distance_from_TSS,$gene_five_prime_coord,$psudo_uORF_start,$gene_start_coord\n";
#            }
#        }
#    }
#}

#exit;

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#open riboseq and assign
my $total_SSU_bam_count=0;

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

    my @CDS_hit;

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

                while ($cigar !~ /^$/){
                    if ($cigar =~ /^([0-9]+[MIDN])/){
                        my $cigar_part = $1;
                        if ($cigar_part =~ /(\d+)M/){   #alignment matching 
                            for my $pos ($leftMost .. ($leftMost+$1-1)){ #search though this position

                                if (exists ($uORF_TSS_TIS_window_search_rev{$chr}{$pos})){
                                    my @uORFs=split("-", $uORF_TSS_TIS_window_search_rev{$chr}{$pos});
                                    for (@uORFs){
                                        my ($u,$feature)=split("~", $_);
                                        $uORF_TSS_TIS_SSU_coverage_by_positon{$u}{$feature}++;
                                    }
                                }

                                if (exists ($uORF_stop_codon_search_rev{$chr}{$pos})){
                                    my @uORFs=split("-", $uORF_stop_codon_search_rev{$chr}{$pos});
                                    for (@uORFs){
                                        $uORF_stop_codon_SSU_coverage{$_}++;                     
                                    }
                                }

#                                if (exists ($bg_TSS_TIS_window_search_rev{$chr}{$pos})){
#                                    my ($gene,$feature,$meta_pos)=split("~", $bg_TSS_TIS_window_search_rev{$chr}{$pos});
#                                    $bg_TSS_TIS_SSU_coverage_by_positon{$gene}{$feature}{$meta_pos}++;
#                                }

                                if( exists($CDS_search_rev{$chr}{$pos})){
                                    my @genes=split("_", $CDS_search_rev{$chr}{$pos});
                                    for (@genes){
                                        push(@CDS_hit,$_);
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

                while ($cigar !~ /^$/){
                    if ($cigar =~ /^([0-9]+[MIDN])/){
                        my $cigar_part = $1;
                        if ($cigar_part =~ /(\d+)M/){   #alignment matching 
                            for my $pos ($leftMost .. ($leftMost+$1-1)){ #search though this position

                                if (exists ($uORF_TSS_TIS_window_search_fwd{$chr}{$pos})){
                                    my @uORFs=split("-", $uORF_TSS_TIS_window_search_fwd{$chr}{$pos});
                                    for (@uORFs){
                                        my ($u,$feature)=split("~", $_);
                                        $uORF_TSS_TIS_SSU_coverage_by_positon{$u}{$feature}++;
                                    }
                                }

                                if (exists ($uORF_stop_codon_search_fwd{$chr}{$pos})){
                                    my @uORFs=split("-", $uORF_stop_codon_search_fwd{$chr}{$pos});
                                    for (@uORFs){
                                        $uORF_stop_codon_SSU_coverage{$_}++;
                                    }
                                }

#                                if (exists ($bg_TSS_TIS_window_search_fwd{$chr}{$pos})){
#                                    my ($gene,$feature,$meta_pos)=split("~", $bg_TSS_TIS_window_search_fwd{$chr}{$pos});
#                                    $bg_TSS_TIS_SSU_coverage_by_positon{$gene}{$feature}{$meta_pos}++;
#                                }

                                if( exists($CDS_search_fwd{$chr}{$pos})){
                                    my @genes=split("_", $CDS_search_fwd{$chr}{$pos});
                                    for (@genes){
                                        push(@CDS_hit,$_);
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

            if (@CDS_hit){
                for my $gene (@CDS_hit){
                    $CDS_counts_SSU{$gene}+=1;
                }
            }
            $total_SSU_bam_count++;
        }
    }
}
close(SSU);

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#again for LSU
my $total_LSU_bam_count=0;

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

    my @CDS_hit;

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

                while ($cigar !~ /^$/){
                    if ($cigar =~ /^([0-9]+[MIDN])/){
                        my $cigar_part = $1;
                        if ($cigar_part =~ /(\d+)M/){   #alignment matching 
                            for my $pos ($leftMost .. ($leftMost+$1-1)){ #search though this position

                                if (exists ($uORF_TSS_TIS_window_search_rev{$chr}{$pos})){
                                    my @uORFs=split("-", $uORF_TSS_TIS_window_search_rev{$chr}{$pos});
                                    for (@uORFs){
                                        my ($u,$feature)=split("~", $_);
                                        $uORF_TSS_TIS_LSU_coverage_by_positon{$u}{$feature}++;
                                    }
                                }

                                if (exists ($uORF_stop_codon_search_rev{$chr}{$pos})){
                                    my @uORFs=split("-", $uORF_stop_codon_search_rev{$chr}{$pos});
                                    for (@uORFs){
                                        $uORF_stop_codon_LSU_coverage{$_}++;
                                    }
                                }

#                                if (exists ($bg_TSS_TIS_window_search_rev{$chr}{$pos})){
#                                    my ($gene,$feature,$meta_pos)=split("~", $bg_TSS_TIS_window_search_rev{$chr}{$pos});
#                                    $bg_TSS_TIS_LSU_coverage_by_positon{$gene}{$feature}{$meta_pos}++;
#                                }

                                if( exists($CDS_search_rev{$chr}{$pos})){
                                    my @genes=split("_", $CDS_search_rev{$chr}{$pos});
                                    for (@genes){
                                        push(@CDS_hit,$_);
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

                while ($cigar !~ /^$/){
                    if ($cigar =~ /^([0-9]+[MIDN])/){
                        my $cigar_part = $1;
                        if ($cigar_part =~ /(\d+)M/){   #alignment matching 
                            for my $pos ($leftMost .. ($leftMost+$1-1)){ #search though this position


                                if (exists ($uORF_TSS_TIS_window_search_fwd{$chr}{$pos})){
                                    my @uORFs=split("-", $uORF_TSS_TIS_window_search_fwd{$chr}{$pos});
                                    for (@uORFs){
                                        my ($u,$feature)=split("~", $_);
                                        $uORF_TSS_TIS_LSU_coverage_by_positon{$u}{$feature}++;
                                    }
                                }

                                if (exists ($uORF_stop_codon_search_fwd{$chr}{$pos})){
                                    my @uORFs=split("-", $uORF_stop_codon_search_fwd{$chr}{$pos});
                                    for (@uORFs){
                                        $uORF_stop_codon_LSU_coverage{$_}++;
                                    }
                                }

#                                if (exists ($bg_TSS_TIS_window_search_fwd{$chr}{$pos})){
#                                    my ($gene,$feature,$meta_pos)=split("~", $bg_TSS_TIS_window_search_fwd{$chr}{$pos});
#                                    $bg_TSS_TIS_LSU_coverage_by_positon{$gene}{$feature}{$meta_pos}++;
#                                }

                                if( exists($CDS_search_fwd{$chr}{$pos})){
                                    my @genes=split("_", $CDS_search_fwd{$chr}{$pos});
                                    for (@genes){
                                        push(@CDS_hit,$_);
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

            if (@CDS_hit){
                for my $gene (@CDS_hit){
                    $CDS_counts_LSU{$gene}+=1;
                }
            }
            $total_LSU_bam_count++;
        }
    }
}
close(LSU);

#for my $u (keys %LSU_coverage_by_positon){
#    for my $f (keys %{ $LSU_coverage_by_positon{$u} }){
#        print "$u,$f\n";
#        for my $p (sort keys %{ $LSU_coverage_by_positon{$u}{$f} }){
#            print "$u,$f,$p,$LSU_coverage_by_positon{$u}{$f}{$p}\n";
#        }
#    }
#}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Positional windows #
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
open (OUT1,">$out_file")  || die "can't open $out_file\n";
print OUT1 "gene_id,gene_cds_SSU_FPKM,gene_cds_LSU_FPKM,gene_overlaps_another_gene,leader_overlap_an_upstream_gene,trailer_overlaps_a_downstream_gene,gene_overlaps_ncRNA,uorf_kozak,uorf_length,uorf_distance_to_TSS,uorf_distance_to_TIS,uorf_stop_distance_to_TIS,uorf_start_codon,uorf_stop_codon,leader_SSU_count_upstream,leader_LSU_count_upstream,leader_SSU_count_downstream,leader_LSU_count_downstream,leader_SSU_stop_codon,leader_LSU_stop_codon\n";

for my $gene (keys %gene_model_fwd){
 
    my $chr=$gene_2_chr{$gene};
    my $gene_start_coord=$start_coord_fwd{$gene};
    my $gene_five_prime_coord=$five_prime_most_coord_fwd{$gene};

    my $gene_overlaps_another_gene=$overlaps_inframe_gene{$gene};
    my $leader_overlap_an_upstream_gene=$leader_overlaps_upstream{$gene};
    my $trailer_overlaps_a_downstream_gene=$gene_overlaps_downstream_leader{$gene};
    my $gene_overlaps_ncRNA=$gene_overlaps_ncRNA{$gene};

    my $CDS_TCPseq_SSU_FPKM=eval{ $CDS_counts_SSU{$gene}/($CDS_length{$gene}*$total_SSU_bam_count)*1000000000 } || 0;
    my $CDS_TCPseq_LSU_FPKM=eval{ $CDS_counts_LSU{$gene}/($CDS_length{$gene}*$total_LSU_bam_count)*1000000000 } || 0;

    my $gene_info="$gene,$CDS_TCPseq_SSU_FPKM,$CDS_TCPseq_LSU_FPKM,$gene_overlaps_another_gene,$leader_overlap_an_upstream_gene,$trailer_overlaps_a_downstream_gene,$gene_overlaps_ncRNA";

    if (exists ($gene_to_uORF{$gene})){

        for my $uORF (keys %{ $gene_to_uORF{$gene} } ) {

            my $u_kozak=$uORF_kozak{$uORF};
            my $u_up_distance=$uORF_to_TSS{$uORF};
            my $u_down_distance=$uORF_to_TIS{$uORF};
            my $u_stop_to_TIS=$uORF_stop_to_TIS{$uORF};
            my $u_length=$uORF_length{$uORF};
            my $u_start=$uORF_start_codon{$uORF};
            my $u_stop=$uORF_stop_codon{$uORF};
          
            my $u_SSU_up=0;
            my $u_LSU_up=0;
            my $u_SSU_down=0;
            my $u_LSU_down=0;
            my $u_SSU_stop=0;
            my $u_LSU_stop=0;

            if (exists ($uORF_TSS_TIS_SSU_coverage_by_positon{$uORF}{"tss_to_uORF_start"})){
                $u_SSU_up=$uORF_TSS_TIS_SSU_coverage_by_positon{$uORF}{"tss_to_uORF_start"};
            }
            #}else{ print "bad_uORF,fwd,SSU,from_TSS,$uORF,$gene\n"; }

            if (exists ($uORF_TSS_TIS_LSU_coverage_by_positon{$uORF}{"tss_to_uORF_start"})){
                $u_LSU_up=$uORF_TSS_TIS_LSU_coverage_by_positon{$uORF}{"tss_to_uORF_start"};
            }
            #}else{ print "bad_uORF,fwd,LSU,from_TSS,$uORF,$gene\n"; }

            if (exists ($uORF_TSS_TIS_SSU_coverage_by_positon{$uORF}{"uORF_stop_to_tis"})){
                $u_SSU_down=$uORF_TSS_TIS_SSU_coverage_by_positon{$uORF}{"uORF_stop_to_tis"};
            }
            #}else{ print "bad_uORF,fwd,SSU,to_TIS,$uORF,$gene\n"; }

            if (exists ($uORF_TSS_TIS_LSU_coverage_by_positon{$uORF}{"uORF_stop_to_tis"})){
                $u_LSU_down=$uORF_TSS_TIS_LSU_coverage_by_positon{$uORF}{"uORF_stop_to_tis"};
            }
            #}else{ print "bad_uORF,fwd,LSU,to_TIS,$uORF,$gene\n"; }

            if (exists ($uORF_stop_codon_SSU_coverage{$uORF})){
                $u_SSU_stop=$uORF_stop_codon_SSU_coverage{$uORF};
            }

            if (exists ($uORF_stop_codon_LSU_coverage{$uORF})){
                $u_LSU_stop=$uORF_stop_codon_LSU_coverage{$uORF};
            }

            print OUT1 "$gene_info,$u_kozak,$u_length,$u_up_distance,$u_down_distance,$u_stop_to_TIS,$u_start,$u_stop,$u_SSU_up,$u_LSU_up,$u_SSU_down,$u_LSU_down,$u_SSU_stop,$u_LSU_stop\n"; 
        }
    }
}

for my $gene (keys %gene_model_rev){
 
    my $chr=$gene_2_chr{$gene};
    my $gene_start_coord=$start_coord_rev{$gene};
    my $gene_five_prime_coord=$five_prime_most_coord_rev{$gene};

    my $gene_overlaps_another_gene=$overlaps_inframe_gene{$gene};
    my $leader_overlap_an_upstream_gene=$leader_overlaps_upstream{$gene};
    my $trailer_overlaps_a_downstream_gene=$gene_overlaps_downstream_leader{$gene};
    my $gene_overlaps_ncRNA=$gene_overlaps_ncRNA{$gene};

    my $CDS_TCPseq_SSU_FPKM=eval{ $CDS_counts_SSU{$gene}/($CDS_length{$gene}*$total_SSU_bam_count)*1000000000 } || 0;
    my $CDS_TCPseq_LSU_FPKM=eval{ $CDS_counts_LSU{$gene}/($CDS_length{$gene}*$total_LSU_bam_count)*1000000000 } || 0;

    my $gene_info="$gene,$CDS_TCPseq_SSU_FPKM,$CDS_TCPseq_LSU_FPKM,$gene_overlaps_another_gene,$leader_overlap_an_upstream_gene,$trailer_overlaps_a_downstream_gene,$gene_overlaps_ncRNA";

    if (exists ($gene_to_uORF{$gene})){

        for my $uORF (keys %{ $gene_to_uORF{$gene} } ) {

            my $u_kozak=$uORF_kozak{$uORF};
            my $u_up_distance=$uORF_to_TSS{$uORF};
            my $u_down_distance=$uORF_to_TIS{$uORF};
            my $u_stop_to_TIS=$uORF_stop_to_TIS{$uORF};
            my $u_length=$uORF_length{$uORF};
            my $u_start=$uORF_start_codon{$uORF};
            my $u_stop=$uORF_stop_codon{$uORF};

            my $u_SSU_up=0;
            my $u_LSU_up=0;
            my $u_SSU_down=0;
            my $u_LSU_down=0;
            my $u_SSU_stop=0;
            my $u_LSU_stop=0;

            if (exists ($uORF_TSS_TIS_SSU_coverage_by_positon{$uORF}{"tss_to_uORF_start"})){
                $u_SSU_up=$uORF_TSS_TIS_SSU_coverage_by_positon{$uORF}{"tss_to_uORF_start"};   #doesn't work
            }
            #}else{ print "bad_uORF,rev,SSU,from_TSS,$uORF,$gene\n"; }

            if (exists ($uORF_TSS_TIS_LSU_coverage_by_positon{$uORF}{"tss_to_uORF_start"})){
                $u_LSU_up=$uORF_TSS_TIS_LSU_coverage_by_positon{$uORF}{"tss_to_uORF_start"};   #doesn't work ~/700 start at TSS?
            }
            #}else{ print "bad_uORF,rev,LSU,from_TSS,$uORF,$gene\n"; }

            if (exists ($uORF_TSS_TIS_SSU_coverage_by_positon{$uORF}{"uORF_stop_to_tis"})){
                $u_SSU_down=$uORF_TSS_TIS_SSU_coverage_by_positon{$uORF}{"uORF_stop_to_tis"}; 
            }
            #}else{ print "bad_uORF,rev,SSU,to_TIS,$uORF,$gene\n"; }

            if (exists ($uORF_TSS_TIS_LSU_coverage_by_positon{$uORF}{"uORF_stop_to_tis"})){
                $u_LSU_down=$uORF_TSS_TIS_LSU_coverage_by_positon{$uORF}{"uORF_stop_to_tis"};
            }
            #}else{ print "bad_uORF,rev,LSU,to_TIS,$uORF,$gene\n"; }
          
            #my $u_SSU_up=$uORF_TSS_TIS_SSU_coverage_by_positon{$uORF}{"tss_to_uORF_start"};
            #my $u_LSU_up=$uORF_TSS_TIS_LSU_coverage_by_positon{$uORF}{"tss_to_uORF_start"};
            #my $u_SSU_down=$uORF_TSS_TIS_SSU_coverage_by_positon{$uORF}{"uORF_stop_to_tis"};
            #my $u_LSU_down=$uORF_TSS_TIS_LSU_coverage_by_positon{$uORF}{"uORF_stop_to_tis"};

            if (exists ($uORF_stop_codon_SSU_coverage{$uORF})){
                $u_SSU_stop=$uORF_stop_codon_SSU_coverage{$uORF};
            }

            if (exists ($uORF_stop_codon_LSU_coverage{$uORF})){
                $u_LSU_stop=$uORF_stop_codon_LSU_coverage{$uORF};
            }

            print OUT1 "$gene_info,$u_kozak,$u_length,$u_up_distance,$u_down_distance,$u_stop_to_TIS,$u_start,$u_stop,$u_SSU_up,$u_LSU_up,$u_SSU_down,$u_LSU_down,$u_SSU_stop,$u_LSU_stop\n"; 
        }
    }
}

close(OUT1);

exit;

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Background windows #
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#open (OUT2,">$out_background") || die "can't open $out_background\n";
#print OUT2 "gene_id,cds_SSU_FPKM,cds_LSU_FPKM,gene_overlaps_another_gene,leader_overlap_an_upstream_gene,trailer_overlaps_a_downstream_gene,gene_overlaps_ncRNA,fraction,region,position_in_region,count\n";
#
#for my $gene (keys %bg_TSS_TIS_SSU_coverage_by_positon){
#
#    my $CDS_TCPseq_SSU_FPKM=eval{ $CDS_counts_SSU{$gene}/($CDS_length{$gene}*$total_SSU_bam_count)*1000000000 } || 0;
#    my $CDS_TCPseq_LSU_FPKM=eval{ $CDS_counts_LSU{$gene}/($CDS_length{$gene}*$total_LSU_bam_count)*1000000000 } || 0;
#
#    my $gene_overlaps_another_gene=$overlaps_inframe_gene{$gene};
#    my $leader_overlap_an_upstream_gene=$leader_overlaps_upstream{$gene};
#    my $trailer_overlaps_a_downstream_gene=$gene_overlaps_downstream_leader{$gene};
#    my $gene_overlaps_ncRNA=$gene_overlaps_ncRNA{$gene};
#
#    my $gene_info="$gene,$CDS_TCPseq_SSU_FPKM,$CDS_TCPseq_LSU_FPKM,$gene_overlaps_another_gene,$leader_overlap_an_upstream_gene,$trailer_overlaps_a_downstream_gene,$gene_overlaps_ncRNA";
#
#    my $region="tss_to_psudo_start";
#    if (exists ($bg_TSS_TIS_SSU_coverage_by_positon{$gene})){ #SSU
#        my $feature_length=keys $bg_TSS_TIS_SSU_coverage_by_positon{$gene}{$region};
#        my $SSU_aggregated_leader_ref=&aggregate_to_50($feature_length, \%{ $bg_TSS_TIS_SSU_coverage_by_positon{$gene}{$region} });
#        my %SSU_aggregated_leader=%{$SSU_aggregated_leader_ref};
#        for my $pos (sort {$a <=> $b} keys %SSU_aggregated_leader){ print OUT2 "$gene_info,SSU,$region,$pos,$SSU_aggregated_leader{$pos}\n"; }
#    }
#
#    if (exists ($bg_TSS_TIS_LSU_coverage_by_positon{$gene})){ #LSU
#        my $feature_length=keys $bg_TSS_TIS_LSU_coverage_by_positon{$gene}{$region};
#        my $LSU_aggregated_leader_ref=&aggregate_to_50($feature_length, \%{ $bg_TSS_TIS_LSU_coverage_by_positon{$gene}{$region} });
#        my %LSU_aggregated_leader=%{$LSU_aggregated_leader_ref};
#        for my $pos (sort {$a <=> $b} keys %LSU_aggregated_leader){ print OUT2 "$gene_info,LSU,$region,$pos,$LSU_aggregated_leader{$pos}\n"; }
#    }
#
#    $region="psudo_start_to_tis";
#    if (exists ($bg_TSS_TIS_SSU_coverage_by_positon{$gene})){ #SSU
#        my $feature_length=keys $bg_TSS_TIS_SSU_coverage_by_positon{$gene}{$region};
#        my $SSU_aggregated_leader_ref=&aggregate_to_50($feature_length, \%{ $bg_TSS_TIS_SSU_coverage_by_positon{$gene}{$region} });
#        my %SSU_aggregated_leader=%{$SSU_aggregated_leader_ref};
#        for my $pos (sort {$a <=> $b} keys %SSU_aggregated_leader){ print OUT2 "$gene_info,SSU,$region,$pos,$SSU_aggregated_leader{$pos}\n"; }
#    }
#
#    if (exists ($bg_TSS_TIS_LSU_coverage_by_positon{$gene})){ #LSU
#        my $feature_length=keys $bg_TSS_TIS_LSU_coverage_by_positon{$gene}{$region};
#        my $LSU_aggregated_leader_ref=&aggregate_to_50($feature_length, \%{ $bg_TSS_TIS_LSU_coverage_by_positon{$gene}{$region} });
#        my %LSU_aggregated_leader=%{$LSU_aggregated_leader_ref};
#        for my $pos (sort {$a <=> $b} keys %LSU_aggregated_leader){ print OUT2 "$gene_info,LSU,$region,$pos,$LSU_aggregated_leader{$pos}\n"; }
#    }
#}
#
#close(OUT2);

exit;

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#caclulate log2
sub log2 {
    my $n = shift;
    my $l = log($n)/log(2);
    #$l = sprintf("%.2f",$l);
    return $l;
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#score the tis surrounding sequence against the Zebrafish kozak PWM
sub score_kozak{
    my $upstream = shift;
    my $downstream = shift;
    my $score=0;
    my $seq_to_score=uc($upstream.$downstream); #concaternate and set to uppercase   
    my @seq_to_score=split("",$seq_to_score);
    my $count=0;

    for my $base (@seq_to_score){
        if (exists ($PWM{$count}{$base} )){
            $score+=$PWM{$count}{$base};
        }
        $count++;
    }
    return $score;
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#scale (aggrgate) the values in a positional hash to 30nt in length
sub aggregate_to_30{

    my $lead_length = shift;
    my $positional_hash_reference = shift;
    #my $uorf = shift;
    my %positional_hash = %$positional_hash_reference;
    my %tmp_hash;
    my %out_hash;
    my $lead_scaling_factor=$lead_length/30;

    for my $pos (sort {$a <=> $b} keys %positional_hash){

        #if ($uorf eq "YGR234W_1"){ print "$pos,$positional_hash{$pos}\n"; } #checked
        my $updated_pos=int(($pos/$lead_scaling_factor)+0.999999);  #round the fraction up   #this should be robust to regions smaller than 30,000,000
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
        #if ($uorf eq "YGR234W_1"){ print "$uorf,$scaled_position,$aggregated_value\n"; } #checked
    }
    return (\%out_hash);
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#scale (aggrgate) the values in a positional hash to 30nt in length
sub aggregate_to_50{

    my $lead_length = shift;
    my $positional_hash_reference = shift;
    #my $uorf = shift;
    my %positional_hash = %$positional_hash_reference;
    my %tmp_hash;
    my %out_hash;
    my $lead_scaling_factor=$lead_length/50;

    for my $pos (sort {$a <=> $b} keys %positional_hash){

        #if ($uorf eq "YGR234W_1"){ print "$pos,$positional_hash{$pos}\n"; } #checked
        my $updated_pos=int(($pos/$lead_scaling_factor)+0.999999);  #round the fraction up   #this should be robust to regions smaller than 30,000,000
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
        #if ($uorf eq "YGR234W_1"){ print "$uorf,$scaled_position,$aggregated_value\n"; } #checked
    }
    return (\%out_hash);
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#scale (aggrgate) the values in a positional hash to 30nt in length
sub aggregate_to_200{

    my $lead_length = shift;
    my $positional_hash_reference = shift;
    my %positional_hash = %$positional_hash_reference;
    my %tmp_hash;
    my %out_hash;
    my $lead_scaling_factor=$lead_length/200;

    for my $pos (sort {$a <=> $b} keys %positional_hash){

        #if ($uorf eq "YGR234W_1"){ print "$pos,$positional_hash{$pos}\n"; } #checked
        my $updated_pos=int(($pos/$lead_scaling_factor)+0.999999);  #round the fraction up   #this should be robust to regions smaller than 30,000,000
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
        #if ($uorf eq "YGR234W_1"){ print "$uorf,$scaled_position,$aggregated_value\n"; } #checked
    }
    return (\%out_hash);
}
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

