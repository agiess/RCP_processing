#!/usr/bin/perl -w
use strict;

#15/01/19
#script to count count the SSU and RFP in fixed windows from the transcription start site

#plot 100 nt downstream of TSS for all genes with leaders greater than 100nt
#plot 200 nt downstream of TSS for all genes with leaders greater than 200nt
#plot 500 nt downstream of TSS for all genes with leaders greater than 500nt

my $most_highly_expressed=$ARGV[0];
my $inGtf=$ARGV[1]; 
my $fasta=$ARGV[2];
my $leaders=$ARGV[3];
my $bam_SSU=$ARGV[4];  
my $bam_LSU=$ARGV[5]; 
my $outdir=$ARGV[6]; 

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
my %stop_codons_search_fwd; #key1=chr, key2=pos, value=gene
my %stop_codons_search_rev; #key1=chr, key2=pos, value=gene

my %CDS_search_fwd; #key1=chr, key2=pos, value=gene
my %CDS_search_rev; #key1=chr, key2=pos, value=gene
my %CDS_counts_SSU;
my %CDS_counts_LSU;

my %transcript_search_fwd; #key1=chr, key2=pos, value=gene
my %transcript_search_rev; #key1=chr, key2=pos, value=gene
my %transcript_counts_SSU;
my %transcript_counts_LSU;

my %CDS_length;
my %transcript_length;

my %transcript_SSU_coverage_by_positon;
my %transcript_LSU_coverage_by_positon;
my %transcript_window_search_fwd;
my %transcript_window_search_rev;

my %transcript_SSU_3prime_counts_by_positon;
my %transcript_LSU_3prime_counts_by_positon;
my %transcript_SSU_5prime_counts_by_positon;
my %transcript_LSU_5prime_counts_by_positon;

my %TSS_100_window_search_fwd;
my %TSS_200_window_search_fwd;
my %TSS_500_window_search_fwd;
my %TSS_100_window_search_rev;
my %TSS_200_window_search_rev;
my %TSS_500_window_search_rev;

my %TSS_100_SSU_coverage_by_positon;
my %TSS_100_SSU_3prime_counts_by_positon;
my %TSS_100_SSU_5prime_counts_by_positon;

my %TSS_100_LSU_3prime_counts_by_positon;
my %TSS_100_LSU_5prime_counts_by_positon;
my %TSS_100_LSU_coverage_by_positon;

my %TSS_200_SSU_coverage_by_positon;
my %TSS_200_SSU_3prime_counts_by_positon;
my %TSS_200_SSU_5prime_counts_by_positon;

my %TSS_200_LSU_3prime_counts_by_positon;
my %TSS_200_LSU_5prime_counts_by_positon;
my %TSS_200_LSU_coverage_by_positon;

my %TSS_500_SSU_coverage_by_positon;
my %TSS_500_SSU_3prime_counts_by_positon;
my %TSS_500_SSU_5prime_counts_by_positon;

my %TSS_500_LSU_3prime_counts_by_positon;
my %TSS_500_LSU_5prime_counts_by_positon;
my %TSS_500_LSU_coverage_by_positon;


for my $gene (keys %gene_model_fwd){
    my $chr=$gene_2_chr{$gene};
    my $start_coord=$start_coord_fwd{$gene};
    my $stop_coord=$stop_coord_fwd{$gene};
    my $TSS_coord=$five_prime_most_coord_fwd{$gene};
    my $TTS_coord=$three_prime_most_coord_fwd{$gene};

    $CDS_counts_SSU{$gene}=0;
    $CDS_counts_LSU{$gene}=0;

    $transcript_counts_SSU{$gene}=0;
    $transcript_counts_LSU{$gene}=0;

    $CDS_length{$gene}=(($stop_coord+2)-$start_coord)+1;
    $transcript_length{$gene}=($TTS_coord-$TSS_coord)+1;

    my $relational_position=0;

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

            $transcript_search_fwd{$chr}{$pos}=$gene;

            if ($coord == $TSS_coord){   

                if ($start_coord-$TSS_coord > 100){ #check that there the leader is at least 100nt long
                    my $relational_position=0; #1 to 100 window upstream and downsteam
                    while ($relational_position < 100){
                        my $genomic_position=$gene_model_fwd{$gene}{($TSS_coord+$relational_position)};
                        $TSS_100_window_search_fwd{$chr}{$genomic_position}=$gene."~".$relational_position;
                        $TSS_100_SSU_coverage_by_positon{$gene}{$relational_position}=0; #initialise
                        $TSS_100_LSU_coverage_by_positon{$gene}{$relational_position}=0; #initialise
                        $TSS_100_SSU_3prime_counts_by_positon{$gene}{$relational_position}=0; #initialise
                        $TSS_100_LSU_3prime_counts_by_positon{$gene}{$relational_position}=0; #initialise
                        $TSS_100_SSU_5prime_counts_by_positon{$gene}{$relational_position}=0; #initialise
                        $TSS_100_LSU_5prime_counts_by_positon{$gene}{$relational_position}=0; #initialise
                        $relational_position++;
                    }
                }

                if ($start_coord-$TSS_coord > 200){ #check that there the leader is at least 100nt long
                    my $relational_position=0; #1 to 100 window upstream and downsteam
                    while ($relational_position < 200){
                        my $genomic_position=$gene_model_fwd{$gene}{($TSS_coord+$relational_position)};
                        $TSS_200_window_search_fwd{$chr}{$genomic_position}=$gene."~".$relational_position;
                        $TSS_200_SSU_coverage_by_positon{$gene}{$relational_position}=0; #initialise
                        $TSS_200_LSU_coverage_by_positon{$gene}{$relational_position}=0; #initialise
                        $TSS_200_SSU_3prime_counts_by_positon{$gene}{$relational_position}=0; #initialise
                        $TSS_200_LSU_3prime_counts_by_positon{$gene}{$relational_position}=0; #initialise
                        $TSS_200_SSU_5prime_counts_by_positon{$gene}{$relational_position}=0; #initialise
                        $TSS_200_LSU_5prime_counts_by_positon{$gene}{$relational_position}=0; #initialise
                        $relational_position++;
                    }
                }

                if ($start_coord-$TSS_coord > 500){ #check that there the leader is at least 100nt long
                    my $relational_position=0; #1 to 100 window upstream and downsteam
                    while ($relational_position < 500){
                        my $genomic_position=$gene_model_fwd{$gene}{($TSS_coord+$relational_position)};
                        $TSS_500_window_search_fwd{$chr}{$genomic_position}=$gene."~".$relational_position;
                        $TSS_500_SSU_coverage_by_positon{$gene}{$relational_position}=0; #initialise
                        $TSS_500_LSU_coverage_by_positon{$gene}{$relational_position}=0; #initialise
                        $TSS_500_SSU_3prime_counts_by_positon{$gene}{$relational_position}=0; #initialise
                        $TSS_500_LSU_3prime_counts_by_positon{$gene}{$relational_position}=0; #initialise
                        $TSS_500_SSU_5prime_counts_by_positon{$gene}{$relational_position}=0; #initialise
                        $TSS_500_LSU_5prime_counts_by_positon{$gene}{$relational_position}=0; #initialise
                        $relational_position++;
                    }
                }
            } 

            if ($coord == $stop_coord){   $stop_codons_search_fwd{$chr}{$pos}=$gene; }
            if ($coord == $stop_coord+1){ $stop_codons_search_fwd{$chr}{$pos}=$gene; }
            if ($coord == $stop_coord+2){ $stop_codons_search_fwd{$chr}{$pos}=$gene; }
           
            if ($coord < $start_coord){
                $leader_search_fwd{$chr}{$pos}=$gene;
            }elsif($coord <= ($stop_coord+2)){  #limit to the stop codon
                $CDS_search_fwd{$chr}{$pos}=$gene;
            }

            #overlapping genes are posible ENSDARG00000095174 ENSDARG00000062049
            if (exists ($transcript_window_search_fwd{$chr}{$pos})){
                $transcript_window_search_fwd{$chr}{$pos}.="_".$gene."~".$relational_position;          
                #print "overlapping $gene\n";
 
            }else{ 
                $transcript_window_search_fwd{$chr}{$pos}=$gene."~".$relational_position;
            }

            $transcript_window_search_fwd{$chr}{$pos}=$gene."~".$relational_position;
            $transcript_SSU_coverage_by_positon{$gene}{$relational_position}=0; #initialise                           
            $transcript_LSU_coverage_by_positon{$gene}{$relational_position}=0; #initialise
            $transcript_SSU_3prime_counts_by_positon{$gene}{$relational_position}=0; #initialise 
            $transcript_LSU_3prime_counts_by_positon{$gene}{$relational_position}=0; #initialise 
            $transcript_SSU_5prime_counts_by_positon{$gene}{$relational_position}=0; #initialise 
            $transcript_LSU_5prime_counts_by_positon{$gene}{$relational_position}=0; #initialise 
            $relational_position++;
        }
    }
}

for my $gene (keys %gene_model_rev){
    my $chr=$gene_2_chr{$gene};
    my $start_coord=$start_coord_rev{$gene};
    my $stop_coord=$stop_coord_rev{$gene};
    my $TSS_coord=$five_prime_most_coord_rev{$gene};    
    my $TTS_coord=$three_prime_most_coord_rev{$gene};

    $CDS_counts_SSU{$gene}=0;
    $CDS_counts_LSU{$gene}=0;

    $transcript_counts_SSU{$gene}=0;
    $transcript_counts_LSU{$gene}=0;

    $CDS_length{$gene}=(($stop_coord+2)-$start_coord)+1;
    $transcript_length{$gene}=($TTS_coord-$TSS_coord)+1;

    my $relational_position=0;

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


            $transcript_search_rev{$chr}{$pos}=$gene;

            if ($coord == $TSS_coord){

                if ($start_coord-$TSS_coord > 100){ #check that there the leader is at least 100nt long
                    my $relational_position=0; #1 to 100 window upstream and downsteam
                    while ($relational_position < 100){
                        my $genomic_position=$gene_model_rev{$gene}{($TSS_coord+$relational_position)};
                        $TSS_100_window_search_rev{$chr}{$genomic_position}=$gene."~".$relational_position;
                        $TSS_100_SSU_coverage_by_positon{$gene}{$relational_position}=0; #initialise
                        $TSS_100_LSU_coverage_by_positon{$gene}{$relational_position}=0; #initialise
                        $TSS_100_SSU_3prime_counts_by_positon{$gene}{$relational_position}=0; #initialise
                        $TSS_100_LSU_3prime_counts_by_positon{$gene}{$relational_position}=0; #initialise
                        $TSS_100_SSU_5prime_counts_by_positon{$gene}{$relational_position}=0; #initialise
                        $TSS_100_LSU_5prime_counts_by_positon{$gene}{$relational_position}=0; #initialise
                        $relational_position++;
                    }
                }

                if ($start_coord-$TSS_coord > 200){ #check that there the leader is at least 100nt long
                    my $relational_position=0; #1 to 100 window upstream and downsteam
                    while ($relational_position < 200){
                        my $genomic_position=$gene_model_rev{$gene}{($TSS_coord+$relational_position)};
                        $TSS_200_window_search_rev{$chr}{$genomic_position}=$gene."~".$relational_position;
                        $TSS_200_SSU_coverage_by_positon{$gene}{$relational_position}=0; #initialise
                        $TSS_200_LSU_coverage_by_positon{$gene}{$relational_position}=0; #initialise
                        $TSS_200_SSU_3prime_counts_by_positon{$gene}{$relational_position}=0; #initialise
                        $TSS_200_LSU_3prime_counts_by_positon{$gene}{$relational_position}=0; #initialise
                        $TSS_200_SSU_5prime_counts_by_positon{$gene}{$relational_position}=0; #initialise
                        $TSS_200_LSU_5prime_counts_by_positon{$gene}{$relational_position}=0; #initialise
                        $relational_position++;
                    }
                }

                if ($start_coord-$TSS_coord > 500){ #check that there the leader is at least 100nt long
                    my $relational_position=0; #1 to 100 window upstream and downsteam
                    while ($relational_position < 500){
                        my $genomic_position=$gene_model_rev{$gene}{($TSS_coord+$relational_position)};
                        $TSS_500_window_search_rev{$chr}{$genomic_position}=$gene."~".$relational_position;
                        $TSS_500_SSU_coverage_by_positon{$gene}{$relational_position}=0; #initialise
                        $TSS_500_LSU_coverage_by_positon{$gene}{$relational_position}=0; #initialise
                        $TSS_500_SSU_3prime_counts_by_positon{$gene}{$relational_position}=0; #initialise
                        $TSS_500_LSU_3prime_counts_by_positon{$gene}{$relational_position}=0; #initialise
                        $TSS_500_SSU_5prime_counts_by_positon{$gene}{$relational_position}=0; #initialise
                        $TSS_500_LSU_5prime_counts_by_positon{$gene}{$relational_position}=0; #initialise
                        $relational_position++;
                    }
                }
            }

            if ($coord == $stop_coord){   $stop_codons_search_rev{$chr}{$pos}=$gene; }
            if ($coord == $stop_coord+1){ $stop_codons_search_rev{$chr}{$pos}=$gene; }
            if ($coord == $stop_coord+2){ $stop_codons_search_rev{$chr}{$pos}=$gene; }

            if ($coord < $start_coord){
                $leader_search_rev{$chr}{$pos}=$gene;
            }elsif($coord <= ($stop_coord+2)){  #limit to the stop codon
                $CDS_search_rev{$chr}{$pos}=$gene;
            }

            #overlapping genes are posible ENSDARG00000095174 ENSDARG00000062049
            if (exists ($transcript_window_search_rev{$chr}{$pos})){
                $transcript_window_search_rev{$chr}{$pos}.="_".$gene."~".$relational_position;             
                #print "overlapping $gene\n";
            }else{ 
                $transcript_window_search_rev{$chr}{$pos}=$gene."~".$relational_position;
            }
            $transcript_SSU_coverage_by_positon{$gene}{$relational_position}=0; #initialise                           
            $transcript_LSU_coverage_by_positon{$gene}{$relational_position}=0; #initialise
            $transcript_SSU_3prime_counts_by_positon{$gene}{$relational_position}=0; #initialise 
            $transcript_LSU_3prime_counts_by_positon{$gene}{$relational_position}=0; #initialise 
            $transcript_SSU_5prime_counts_by_positon{$gene}{$relational_position}=0; #initialise 
            $transcript_LSU_5prime_counts_by_positon{$gene}{$relational_position}=0; #initialise 
            $relational_position++;
        }
    }
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#open riboseq and assign
my $SSU_total;
my $dORF_SSU_total;

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

    my $CDS_hit=0;
    my $transcript_hit=0;

    if ($chr=~/chr(.*)/){ #ensembl chromosome names do not contain the "chr" prefix
        $chr=$1;
    }

    unless ($flag & 0x4){   #if aligned

        if ($mapq >= 10){     #mapping uniqnes filter

        if ((length($seq) >= 26) && (length($seq) <= 30)){

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

                if (exists ($TSS_100_window_search_rev{$chr}{$threePrime})){
                    my ($gene,$meta_pos)=split("~", $TSS_100_window_search_rev{$chr}{$threePrime});
                    $TSS_100_SSU_3prime_counts_by_positon{$gene}{$meta_pos}++;
                }

                if (exists ($TSS_200_window_search_rev{$chr}{$threePrime})){
                    my ($gene,$meta_pos)=split("~", $TSS_200_window_search_rev{$chr}{$threePrime});
                    $TSS_200_SSU_3prime_counts_by_positon{$gene}{$meta_pos}++;
                }

                if (exists ($TSS_500_window_search_rev{$chr}{$threePrime})){
                    my ($gene,$meta_pos)=split("~", $TSS_500_window_search_rev{$chr}{$threePrime});
                    $TSS_500_SSU_3prime_counts_by_positon{$gene}{$meta_pos}++;
                }


                if (exists ($transcript_window_search_rev{$chr}{$threePrime})){
                    my @genes = split ("_", $transcript_window_search_rev{$chr}{$threePrime});
                    for (@genes){
                        my ($gene,$meta_pos)=split("~", $_);
                        $transcript_SSU_3prime_counts_by_positon{$gene}{$meta_pos}++;
                    }
                }
                                       
                if (exists ($TSS_100_window_search_rev{$chr}{$fivePrime})){
                    my ($gene,$meta_pos)=split("~", $TSS_100_window_search_rev{$chr}{$fivePrime});
                    $TSS_100_SSU_5prime_counts_by_positon{$gene}{$meta_pos}++;
                }

                if (exists ($TSS_200_window_search_rev{$chr}{$fivePrime})){
                    my ($gene,$meta_pos)=split("~", $TSS_200_window_search_rev{$chr}{$fivePrime});
                    $TSS_200_SSU_5prime_counts_by_positon{$gene}{$meta_pos}++;
                }

                if (exists ($TSS_500_window_search_rev{$chr}{$fivePrime})){
                    my ($gene,$meta_pos)=split("~", $TSS_500_window_search_rev{$chr}{$fivePrime});
                    $TSS_500_SSU_5prime_counts_by_positon{$gene}{$meta_pos}++;
                }


                if (exists ($transcript_window_search_rev{$chr}{$fivePrime})){
                    my @genes = split ("_", $transcript_window_search_rev{$chr}{$fivePrime});
                    for (@genes){
                        my ($gene,$meta_pos)=split("~", $_);
                        $transcript_SSU_5prime_counts_by_positon{$gene}{$meta_pos}++;
                    }
                }

                while ($cigar !~ /^$/){
                    if ($cigar =~ /^([0-9]+[MIDN])/){
                        my $cigar_part = $1;
                        if ($cigar_part =~ /(\d+)M/){   #alignment matching 
                            for my $pos ($leftMost .. ($leftMost+$1-1)){ #search though this position

                                if (exists ($TSS_100_window_search_rev{$chr}{$pos})){
                                    my ($gene,$meta_pos)=split("~", $TSS_100_window_search_rev{$chr}{$pos});
                                    $TSS_100_SSU_coverage_by_positon{$gene}{$meta_pos}++;
                                }

                                if (exists ($TSS_200_window_search_rev{$chr}{$pos})){
                                    my ($gene,$meta_pos)=split("~", $TSS_200_window_search_rev{$chr}{$pos});
                                    $TSS_200_SSU_coverage_by_positon{$gene}{$meta_pos}++;
                                }

                                if (exists ($TSS_500_window_search_rev{$chr}{$pos})){
                                    my ($gene,$meta_pos)=split("~", $TSS_500_window_search_rev{$chr}{$pos});
                                    $TSS_500_SSU_coverage_by_positon{$gene}{$meta_pos}++;
                                }

                                if (exists ($transcript_window_search_rev{$chr}{$pos})){
                                    my @genes = split ("_", $transcript_window_search_rev{$chr}{$pos});
                                    for (@genes){
                                        my ($gene,$meta_pos)=split("~", $_);
                                        $transcript_SSU_coverage_by_positon{$gene}{$meta_pos}++;
                                        #if ($gene eq "ENSDARG00000095174"){ print "ENSDARG00000095174 hit SSU rev\n"; }
                                    }
                                }

                                if( exists($CDS_search_rev{$chr}{$pos})){
                                    $CDS_hit=$CDS_search_rev{$chr}{$pos};
                                }

                                if (exists($transcript_search_rev{$chr}{$pos})){
                                    $transcript_hit=$transcript_search_rev{$chr}{$pos}++;
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

                if (exists ($TSS_100_window_search_fwd{$chr}{$threePrime})){
                     my ($gene,$meta_pos)=split("~", $TSS_100_window_search_fwd{$chr}{$threePrime});
                     $TSS_100_SSU_3prime_counts_by_positon{$gene}{$meta_pos}++;
                }

                if (exists ($TSS_200_window_search_fwd{$chr}{$threePrime})){
                     my ($gene,$meta_pos)=split("~", $TSS_200_window_search_fwd{$chr}{$threePrime});
                     $TSS_200_SSU_3prime_counts_by_positon{$gene}{$meta_pos}++;
                }

                if (exists ($TSS_500_window_search_fwd{$chr}{$threePrime})){
                     my ($gene,$meta_pos)=split("~", $TSS_500_window_search_fwd{$chr}{$threePrime});
                     $TSS_500_SSU_3prime_counts_by_positon{$gene}{$meta_pos}++;
                }


                if (exists ($transcript_window_search_fwd{$chr}{$threePrime})){
                    my @genes = split ("_", $transcript_window_search_fwd{$chr}{$threePrime});
                    for (@genes){
                        my ($gene,$meta_pos)=split("~", $_);
                        $transcript_SSU_3prime_counts_by_positon{$gene}{$meta_pos}++;
                    }
                }

                if (exists ($TSS_100_window_search_fwd{$chr}{$fivePrime})){
                     my ($gene,$meta_pos)=split("~", $TSS_100_window_search_fwd{$chr}{$fivePrime});
                     $TSS_100_SSU_5prime_counts_by_positon{$gene}{$meta_pos}++;
                }

                if (exists ($TSS_200_window_search_fwd{$chr}{$fivePrime})){
                     my ($gene,$meta_pos)=split("~", $TSS_200_window_search_fwd{$chr}{$fivePrime});
                     $TSS_200_SSU_5prime_counts_by_positon{$gene}{$meta_pos}++;
                }

                if (exists ($TSS_500_window_search_fwd{$chr}{$fivePrime})){
                     my ($gene,$meta_pos)=split("~", $TSS_500_window_search_fwd{$chr}{$fivePrime});
                     $TSS_500_SSU_5prime_counts_by_positon{$gene}{$meta_pos}++;
                }


                if (exists ($transcript_window_search_fwd{$chr}{$fivePrime})){
                    my @genes = split ("_", $transcript_window_search_fwd{$chr}{$fivePrime});
                    for (@genes){
                        my ($gene,$meta_pos)=split("~", $_);
                        $transcript_SSU_5prime_counts_by_positon{$gene}{$meta_pos}++;
                    }
                }

                while ($cigar !~ /^$/){
                    if ($cigar =~ /^([0-9]+[MIDN])/){
                        my $cigar_part = $1;
                        if ($cigar_part =~ /(\d+)M/){   #alignment matching 
                            for my $pos ($leftMost .. ($leftMost+$1-1)){ #search though this position

                                if (exists ($TSS_100_window_search_fwd{$chr}{$pos})){
                                    my ($gene,$meta_pos)=split("~", $TSS_100_window_search_fwd{$chr}{$pos});
                                    $TSS_100_SSU_coverage_by_positon{$gene}{$meta_pos}++;
                                }

                                if (exists ($TSS_200_window_search_fwd{$chr}{$pos})){
                                    my ($gene,$meta_pos)=split("~", $TSS_200_window_search_fwd{$chr}{$pos});
                                    $TSS_200_SSU_coverage_by_positon{$gene}{$meta_pos}++;
                                }                                

                                if (exists ($TSS_500_window_search_fwd{$chr}{$pos})){
                                    my ($gene,$meta_pos)=split("~", $TSS_500_window_search_fwd{$chr}{$pos});
                                    $TSS_500_SSU_coverage_by_positon{$gene}{$meta_pos}++;
                                }

                                if (exists ($transcript_window_search_fwd{$chr}{$pos})){
                                    my @genes = split ("_", $transcript_window_search_fwd{$chr}{$pos});
                                    for (@genes){
                                        my ($gene,$meta_pos)=split("~", $_);
                                        $transcript_SSU_coverage_by_positon{$gene}{$meta_pos}++;
                                        #if ($gene eq "ENSDARG00000095174"){ print "ENSDARG00000095174 hit SSU fwd\n"; }
                                    }
                                }   

                                if( exists($CDS_search_fwd{$chr}{$pos})){
                                    $CDS_hit=$CDS_search_fwd{$chr}{$pos};
                                }

                                if (exists($transcript_search_fwd{$chr}{$pos})){
                                    $transcript_hit=$transcript_search_fwd{$chr}{$pos}++;
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

            if ($transcript_hit){
                $transcript_counts_SSU{$transcript_hit}+=1;
            }

            if ($CDS_hit){
                $CDS_counts_SSU{$CDS_hit}+=1;
            }
            $total_SSU_bam_count++;
        }
        }
    }
}
close(SSU);

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#again for LSU
my $LSU_total;
my $dORF_LSU_total;
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

    my $CDS_hit=0;
    my $transcript_hit=0;

    if ($chr=~/chr(.*)/){ #ensembl chromosome names do not contain the "chr" prefix
        $chr=$1;
    }

    unless ($flag & 0x4){   #if aligned

        if ($mapq >= 10){     #mapping uniqnes filter

        if ((length($seq) >= 26) && (length($seq) <= 30)){

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

                if (exists ($TSS_100_window_search_rev{$chr}{$threePrime})){
                     my ($gene,$meta_pos)=split("~", $TSS_100_window_search_rev{$chr}{$threePrime});
                     $TSS_100_LSU_3prime_counts_by_positon{$gene}{$meta_pos}++;
                }

                if (exists ($TSS_200_window_search_rev{$chr}{$threePrime})){
                     my ($gene,$meta_pos)=split("~", $TSS_200_window_search_rev{$chr}{$threePrime});
                     $TSS_200_LSU_3prime_counts_by_positon{$gene}{$meta_pos}++;
                }

                if (exists ($TSS_500_window_search_rev{$chr}{$threePrime})){
                     my ($gene,$meta_pos)=split("~", $TSS_500_window_search_rev{$chr}{$threePrime});
                     $TSS_500_LSU_3prime_counts_by_positon{$gene}{$meta_pos}++;
                }

                if (exists ($transcript_window_search_rev{$chr}{$threePrime})){
                    my @genes = split ("_", $transcript_window_search_rev{$chr}{$threePrime});
                    for (@genes){
                        my ($gene,$meta_pos)=split("~", $_);
                        $transcript_LSU_3prime_counts_by_positon{$gene}{$meta_pos}++;
                    }
                }

                if (exists ($TSS_100_window_search_rev{$chr}{$fivePrime})){
                     my ($gene,$meta_pos)=split("~", $TSS_100_window_search_rev{$chr}{$fivePrime});
                     $TSS_100_LSU_5prime_counts_by_positon{$gene}{$meta_pos}++;
                }

                if (exists ($TSS_200_window_search_rev{$chr}{$fivePrime})){
                     my ($gene,$meta_pos)=split("~", $TSS_200_window_search_rev{$chr}{$fivePrime});
                     $TSS_200_LSU_5prime_counts_by_positon{$gene}{$meta_pos}++;
                }

                if (exists ($TSS_500_window_search_rev{$chr}{$fivePrime})){
                     my ($gene,$meta_pos)=split("~", $TSS_500_window_search_rev{$chr}{$fivePrime});
                     $TSS_500_LSU_5prime_counts_by_positon{$gene}{$meta_pos}++;
                }

                if (exists ($transcript_window_search_rev{$chr}{$fivePrime})){
                    my @genes = split ("_", $transcript_window_search_rev{$chr}{$fivePrime});
                    for (@genes){
                        my ($gene,$meta_pos)=split("~", $_);
                        $transcript_LSU_5prime_counts_by_positon{$gene}{$meta_pos}++;
                    }
                }

                while ($cigar !~ /^$/){
                    if ($cigar =~ /^([0-9]+[MIDN])/){
                        my $cigar_part = $1;
                        if ($cigar_part =~ /(\d+)M/){   #alignment matching 
                            for my $pos ($leftMost .. ($leftMost+$1-1)){ #search though this position

                                if (exists ($TSS_100_window_search_rev{$chr}{$pos})){
                                    my ($gene,$meta_pos)=split("~", $TSS_100_window_search_rev{$chr}{$pos});
                                    $TSS_100_LSU_coverage_by_positon{$gene}{$meta_pos}++;
                                }

                                if (exists ($TSS_200_window_search_rev{$chr}{$pos})){
                                    my ($gene,$meta_pos)=split("~", $TSS_200_window_search_rev{$chr}{$pos});
                                    $TSS_200_LSU_coverage_by_positon{$gene}{$meta_pos}++;
                                }

                                if (exists ($TSS_500_window_search_rev{$chr}{$pos})){
                                    my ($gene,$meta_pos)=split("~", $TSS_500_window_search_rev{$chr}{$pos});
                                    $TSS_500_LSU_coverage_by_positon{$gene}{$meta_pos}++;
                                }

                                if (exists ($transcript_window_search_rev{$chr}{$pos})){
                                    my @genes = split ("_", $transcript_window_search_rev{$chr}{$pos});
                                    for (@genes){
                                        my ($gene,$meta_pos)=split("~", $_);
                                        $transcript_LSU_coverage_by_positon{$gene}{$meta_pos}++;
                                        #if ($gene eq "ENSDARG00000095174"){ print "ENSDARG00000095174 hit LSU rev\n"; }
                                    }
                                }   

                                if( exists($CDS_search_rev{$chr}{$pos})){
                                    $CDS_hit=$CDS_search_rev{$chr}{$pos};
                                }

                                if (exists($transcript_search_rev{$chr}{$pos})){
                                    $transcript_hit=$transcript_search_rev{$chr}{$pos}++;
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

                if (exists ($TSS_100_window_search_fwd{$chr}{$threePrime})){
                     my ($gene,$meta_pos)=split("~", $TSS_100_window_search_fwd{$chr}{$threePrime});
                     $TSS_100_LSU_3prime_counts_by_positon{$gene}{$meta_pos}++;
                }

                if (exists ($TSS_200_window_search_fwd{$chr}{$threePrime})){
                     my ($gene,$meta_pos)=split("~", $TSS_200_window_search_fwd{$chr}{$threePrime});
                     $TSS_200_LSU_3prime_counts_by_positon{$gene}{$meta_pos}++;
                }

                if (exists ($TSS_500_window_search_fwd{$chr}{$threePrime})){
                     my ($gene,$meta_pos)=split("~", $TSS_500_window_search_fwd{$chr}{$threePrime});
                     $TSS_500_LSU_3prime_counts_by_positon{$gene}{$meta_pos}++;
                }

                if (exists ($transcript_window_search_fwd{$chr}{$threePrime})){
                    my @genes = split ("_", $transcript_window_search_fwd{$chr}{$threePrime});
                    for (@genes){
                        my ($gene,$meta_pos)=split("~", $_);
                        $transcript_LSU_3prime_counts_by_positon{$gene}{$meta_pos}++;
                    }
                }

                if (exists ($TSS_100_window_search_fwd{$chr}{$fivePrime})){
                     my ($gene,$meta_pos)=split("~", $TSS_100_window_search_fwd{$chr}{$fivePrime});
                     $TSS_100_LSU_5prime_counts_by_positon{$gene}{$meta_pos}++;
                }

                if (exists ($TSS_200_window_search_fwd{$chr}{$fivePrime})){
                     my ($gene,$meta_pos)=split("~", $TSS_200_window_search_fwd{$chr}{$fivePrime});
                     $TSS_200_LSU_5prime_counts_by_positon{$gene}{$meta_pos}++;
                }

                if (exists ($TSS_500_window_search_fwd{$chr}{$fivePrime})){
                     my ($gene,$meta_pos)=split("~", $TSS_500_window_search_fwd{$chr}{$fivePrime});
                     $TSS_500_LSU_5prime_counts_by_positon{$gene}{$meta_pos}++;
                }

                if (exists ($transcript_window_search_fwd{$chr}{$fivePrime})){
                    my @genes = split ("_", $transcript_window_search_fwd{$chr}{$fivePrime});
                    for (@genes){
                        my ($gene,$meta_pos)=split("~", $_);
                        $transcript_LSU_5prime_counts_by_positon{$gene}{$meta_pos}++;
                    }
                }

                while ($cigar !~ /^$/){
                    if ($cigar =~ /^([0-9]+[MIDN])/){
                        my $cigar_part = $1;
                        if ($cigar_part =~ /(\d+)M/){   #alignment matching 
                            for my $pos ($leftMost .. ($leftMost+$1-1)){ #search though this position

                                if (exists ($TSS_100_window_search_fwd{$chr}{$pos})){
                                    my ($gene,$meta_pos)=split("~", $TSS_100_window_search_fwd{$chr}{$pos});
                                    $TSS_100_LSU_coverage_by_positon{$gene}{$meta_pos}++;
                                }

                                if (exists ($TSS_200_window_search_fwd{$chr}{$pos})){
                                    my ($gene,$meta_pos)=split("~", $TSS_200_window_search_fwd{$chr}{$pos});
                                    $TSS_200_LSU_coverage_by_positon{$gene}{$meta_pos}++;
                                }

                                if (exists ($TSS_500_window_search_fwd{$chr}{$pos})){
                                    my ($gene,$meta_pos)=split("~", $TSS_500_window_search_fwd{$chr}{$pos});
                                    $TSS_500_LSU_coverage_by_positon{$gene}{$meta_pos}++;
                                }

                                if (exists ($transcript_window_search_fwd{$chr}{$pos})){
                                    my @genes = split ("_", $transcript_window_search_fwd{$chr}{$pos});
                                    for (@genes){
                                        my ($gene,$meta_pos)=split("~", $_);
                                        $transcript_LSU_coverage_by_positon{$gene}{$meta_pos}++;
                                        #if ($gene eq "ENSDARG00000095174"){ print "ENSDARG00000095174 hit LSU fwd\n"; }
                                    }
                                }   

                                if( exists($CDS_search_fwd{$chr}{$pos})){
                                    $CDS_hit=$CDS_search_fwd{$chr}{$pos};
                                }

                                if (exists($transcript_search_fwd{$chr}{$pos})){
                                    $transcript_hit=$transcript_search_fwd{$chr}{$pos}++;
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

            if ($transcript_hit){
                $transcript_counts_LSU{$transcript_hit}+=1;
            }

            if ($CDS_hit){
                $CDS_counts_LSU{$CDS_hit}+=1;
            }
            $total_LSU_bam_count++;
        }
        }
    }
}
close(LSU);

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#open files in outdir
my $out_file1 = $outdir."/TSS_coverage_100.csv";
my $out_file2 = $outdir."/TSS_coverage_200.csv";
my $out_file3 = $outdir."/TSS_coverage_500.csv";
my $out_file4 = $outdir."/transcript_coverage.csv";

my $out_file5 = $outdir."/TSS_3prime_100.csv";
my $out_file6 = $outdir."/TSS_3prime_200.csv";
my $out_file7 = $outdir."/TSS_3prime_500.csv";
my $out_file8 = $outdir."/transcript_3prime.csv";

my $out_file9  = $outdir."/TSS_5prime_100.csv";
my $out_file10 = $outdir."/TSS_5prime_200.csv";
my $out_file11 = $outdir."/TSS_5prime_500.csv";
my $out_file12 = $outdir."/transcript_5prime.csv";

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Positional windows for uORFS #
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
open (OUT1,">$out_file1")  || die "can't open $out_file1\n";
print OUT1 "gene_id,transcript_SSU_count,transcript_LSU_count,cds_SSU_count,cds_LSU_count,cds_SSU_FPKM,cds_LSU_FPKM,transcript_length,cds_length,gene_overlaps_another_gene,leader_overlap_an_upstream_gene,trailer_overlaps_a_downstream_gene,gene_overlaps_ncRNA,fraction,position_in_region,count\n";

for my $gene (keys %TSS_100_SSU_coverage_by_positon){

    my $CDS_TCPseq_SSU_FPKM=eval{ $CDS_counts_SSU{$gene}/($CDS_length{$gene}*$total_SSU_bam_count)*1000000000 } || 0;
    my $CDS_TCPseq_LSU_FPKM=eval{ $CDS_counts_LSU{$gene}/($CDS_length{$gene}*$total_LSU_bam_count)*1000000000 } || 0;
 
    my $transcript_length=$transcript_length{$gene};
    my $cds_length=$CDS_length{$gene};

    my $gene_overlaps_another_gene=$overlaps_inframe_gene{$gene};
    my $leader_overlap_an_upstream_gene=$leader_overlaps_upstream{$gene};
    my $trailer_overlaps_a_downstream_gene=$gene_overlaps_downstream_leader{$gene};
    my $gene_overlaps_ncRNA=$gene_overlaps_ncRNA{$gene};

    my $gene_info="$gene,$transcript_counts_SSU{$gene},$transcript_counts_LSU{$gene},$CDS_counts_SSU{$gene},$CDS_counts_LSU{$gene},$CDS_TCPseq_SSU_FPKM,$CDS_TCPseq_LSU_FPKM,$transcript_length,$cds_length,$gene_overlaps_another_gene,$leader_overlap_an_upstream_gene,$trailer_overlaps_a_downstream_gene,$gene_overlaps_ncRNA";

    for my $pos (sort {$a <=> $b} keys %{$TSS_100_SSU_coverage_by_positon{$gene}}){ 
        print OUT1 "$gene_info,SSU,$pos,$TSS_100_SSU_coverage_by_positon{$gene}{$pos}\n"; 
    }

    for my $pos (sort {$a <=> $b} keys %{$TSS_100_LSU_coverage_by_positon{$gene}}){
        print OUT1 "$gene_info,LSU,$pos,$TSS_100_LSU_coverage_by_positon{$gene}{$pos}\n"; 
    }
}
close(OUT1);


open (OUT2,">$out_file2")  || die "can't open $out_file2\n";
print OUT2 "gene_id,transcript_SSU_count,transcript_LSU_count,cds_SSU_count,cds_LSU_count,cds_SSU_FPKM,cds_LSU_FPKM,transcript_length,cds_length,gene_overlaps_another_gene,leader_overlap_an_upstream_gene,trailer_overlaps_a_downstream_gene,gene_overlaps_ncRNA,fraction,position_in_region,count\n";

for my $gene (keys %TSS_200_SSU_coverage_by_positon){

    my $CDS_TCPseq_SSU_FPKM=eval{ $CDS_counts_SSU{$gene}/($CDS_length{$gene}*$total_SSU_bam_count)*1000000000 } || 0;
    my $CDS_TCPseq_LSU_FPKM=eval{ $CDS_counts_LSU{$gene}/($CDS_length{$gene}*$total_LSU_bam_count)*1000000000 } || 0;

    my $transcript_length=$transcript_length{$gene};
    my $cds_length=$CDS_length{$gene};

    my $gene_overlaps_another_gene=$overlaps_inframe_gene{$gene};
    my $leader_overlap_an_upstream_gene=$leader_overlaps_upstream{$gene};
    my $trailer_overlaps_a_downstream_gene=$gene_overlaps_downstream_leader{$gene};
    my $gene_overlaps_ncRNA=$gene_overlaps_ncRNA{$gene};

    my $gene_info="$gene,$transcript_counts_SSU{$gene},$transcript_counts_LSU{$gene},$CDS_counts_SSU{$gene},$CDS_counts_LSU{$gene},$CDS_TCPseq_SSU_FPKM,$CDS_TCPseq_LSU_FPKM,$transcript_length,$cds_length,$gene_overlaps_another_gene,$leader_overlap_an_upstream_gene,$trailer_overlaps_a_downstream_gene,$gene_overlaps_ncRNA";

    for my $pos (sort {$a <=> $b} keys %{$TSS_200_SSU_coverage_by_positon{$gene}}){
        print OUT2 "$gene_info,SSU,$pos,$TSS_200_SSU_coverage_by_positon{$gene}{$pos}\n";
    }

    for my $pos (sort {$a <=> $b} keys %{$TSS_200_LSU_coverage_by_positon{$gene}}){
        print OUT2 "$gene_info,LSU,$pos,$TSS_200_LSU_coverage_by_positon{$gene}{$pos}\n";
    }
}
close(OUT2);


open (OUT3,">$out_file3")  || die "can't open $out_file3\n";
print OUT3 "gene_id,transcript_SSU_count,transcript_LSU_count,cds_SSU_count,cds_LSU_count,cds_SSU_FPKM,cds_LSU_FPKM,transcript_length,cds_length,gene_overlaps_another_gene,leader_overlap_an_upstream_gene,trailer_overlaps_a_downstream_gene,gene_overlaps_ncRNA,fraction,position_in_region,count\n";

for my $gene (keys %TSS_500_SSU_coverage_by_positon){

    my $CDS_TCPseq_SSU_FPKM=eval{ $CDS_counts_SSU{$gene}/($CDS_length{$gene}*$total_SSU_bam_count)*1000000000 } || 0;
    my $CDS_TCPseq_LSU_FPKM=eval{ $CDS_counts_LSU{$gene}/($CDS_length{$gene}*$total_LSU_bam_count)*1000000000 } || 0;

    my $transcript_length=$transcript_length{$gene};
    my $cds_length=$CDS_length{$gene};

    my $gene_overlaps_another_gene=$overlaps_inframe_gene{$gene};
    my $leader_overlap_an_upstream_gene=$leader_overlaps_upstream{$gene};
    my $trailer_overlaps_a_downstream_gene=$gene_overlaps_downstream_leader{$gene};
    my $gene_overlaps_ncRNA=$gene_overlaps_ncRNA{$gene};

    my $gene_info="$gene,$transcript_counts_SSU{$gene},$transcript_counts_LSU{$gene},$CDS_counts_SSU{$gene},$CDS_counts_LSU{$gene},$CDS_TCPseq_SSU_FPKM,$CDS_TCPseq_LSU_FPKM,$transcript_length,$cds_length,$gene_overlaps_another_gene,$leader_overlap_an_upstream_gene,$trailer_overlaps_a_downstream_gene,$gene_overlaps_ncRNA";

    for my $pos (sort {$a <=> $b} keys %{$TSS_500_SSU_coverage_by_positon{$gene}}){
        print OUT3 "$gene_info,SSU,$pos,$TSS_500_SSU_coverage_by_positon{$gene}{$pos}\n";
    }

    for my $pos (sort {$a <=> $b} keys %{$TSS_500_LSU_coverage_by_positon{$gene}}){
        print OUT3 "$gene_info,LSU,$pos,$TSS_500_LSU_coverage_by_positon{$gene}{$pos}\n";
    }
}
close(OUT3);

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Positional windows for whole regions transcripts #
##¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
open (OUT4,">$out_file4")  || die "can't open $out_file4\n";
print OUT4 "gene_id,fraction,position_in_region,count\n";

for my $gene (keys %transcript_SSU_coverage_by_positon){

    #SSU
    for my $pos (sort {$a <=> $b} keys %{$transcript_SSU_coverage_by_positon{$gene}}){
         print OUT4 "$gene,SSU,$pos,$transcript_SSU_coverage_by_positon{$gene}{$pos}\n";
    }
    
    #LSU
    for my $pos (sort {$a <=> $b} keys %{$transcript_LSU_coverage_by_positon{$gene}}){
         print OUT4 "$gene,LSU,$pos,$transcript_LSU_coverage_by_positon{$gene}{$pos}\n";
    }
}
close(OUT4);

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Positional windows for uORFS # 3'
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
open (OUT5,">$out_file5")  || die "can't open $out_file5\n";
print OUT5 "gene_id,transcript_SSU_count,transcript_LSU_count,cds_SSU_count,cds_LSU_count,cds_SSU_FPKM,cds_LSU_FPKM,transcript_length,cds_length,gene_overlaps_another_gene,leader_overlap_an_upstream_gene,trailer_overlaps_a_downstream_gene,gene_overlaps_ncRNA,fraction,position_in_region,count\n";

for my $gene (keys %TSS_100_SSU_3prime_counts_by_positon){

    my $CDS_TCPseq_SSU_FPKM=eval{ $CDS_counts_SSU{$gene}/($CDS_length{$gene}*$total_SSU_bam_count)*1000000000 } || 0;
    my $CDS_TCPseq_LSU_FPKM=eval{ $CDS_counts_LSU{$gene}/($CDS_length{$gene}*$total_LSU_bam_count)*1000000000 } || 0;

    my $transcript_length=$transcript_length{$gene};
    my $cds_length=$CDS_length{$gene};

    my $gene_overlaps_another_gene=$overlaps_inframe_gene{$gene};
    my $leader_overlap_an_upstream_gene=$leader_overlaps_upstream{$gene};
    my $trailer_overlaps_a_downstream_gene=$gene_overlaps_downstream_leader{$gene};
    my $gene_overlaps_ncRNA=$gene_overlaps_ncRNA{$gene};

    my $gene_info="$gene,$transcript_counts_SSU{$gene},$transcript_counts_LSU{$gene},$CDS_counts_SSU{$gene},$CDS_counts_LSU{$gene},$CDS_TCPseq_SSU_FPKM,$CDS_TCPseq_LSU_FPKM,$transcript_length,$cds_length,$gene_overlaps_another_gene,$leader_overlap_an_upstream_gene,$trailer_overlaps_a_downstream_gene,$gene_overlaps_ncRNA";

    for my $pos (sort {$a <=> $b} keys %{$TSS_100_SSU_3prime_counts_by_positon{$gene}}){
        print OUT5 "$gene_info,SSU,$pos,$TSS_100_SSU_3prime_counts_by_positon{$gene}{$pos}\n";
    }

    for my $pos (sort {$a <=> $b} keys %{$TSS_100_LSU_3prime_counts_by_positon{$gene}}){
        print OUT5 "$gene_info,LSU,$pos,$TSS_100_LSU_3prime_counts_by_positon{$gene}{$pos}\n";
    }
}
close(OUT5);


open (OUT6,">$out_file6")  || die "can't open $out_file6\n";
print OUT6 "gene_id,transcript_SSU_count,transcript_LSU_count,cds_SSU_count,cds_LSU_count,cds_SSU_FPKM,cds_LSU_FPKM,transcript_length,cds_length,gene_overlaps_another_gene,leader_overlap_an_upstream_gene,trailer_overlaps_a_downstream_gene,gene_overlaps_ncRNA,fraction,position_in_region,count\n";

for my $gene (keys %TSS_200_SSU_3prime_counts_by_positon){

    my $CDS_TCPseq_SSU_FPKM=eval{ $CDS_counts_SSU{$gene}/($CDS_length{$gene}*$total_SSU_bam_count)*1000000000 } || 0;
    my $CDS_TCPseq_LSU_FPKM=eval{ $CDS_counts_LSU{$gene}/($CDS_length{$gene}*$total_LSU_bam_count)*1000000000 } || 0;

    my $transcript_length=$transcript_length{$gene};
    my $cds_length=$CDS_length{$gene};

    my $gene_overlaps_another_gene=$overlaps_inframe_gene{$gene};
    my $leader_overlap_an_upstream_gene=$leader_overlaps_upstream{$gene};
    my $trailer_overlaps_a_downstream_gene=$gene_overlaps_downstream_leader{$gene};
    my $gene_overlaps_ncRNA=$gene_overlaps_ncRNA{$gene};

    my $gene_info="$gene,$transcript_counts_SSU{$gene},$transcript_counts_LSU{$gene},$CDS_counts_SSU{$gene},$CDS_counts_LSU{$gene},$CDS_TCPseq_SSU_FPKM,$CDS_TCPseq_LSU_FPKM,$transcript_length,$cds_length,$gene_overlaps_another_gene,$leader_overlap_an_upstream_gene,$trailer_overlaps_a_downstream_gene,$gene_overlaps_ncRNA";

    for my $pos (sort {$a <=> $b} keys %{$TSS_200_SSU_3prime_counts_by_positon{$gene}}){
        print OUT6 "$gene_info,SSU,$pos,$TSS_200_SSU_3prime_counts_by_positon{$gene}{$pos}\n";
    }

    for my $pos (sort {$a <=> $b} keys %{$TSS_200_LSU_3prime_counts_by_positon{$gene}}){
        print OUT6 "$gene_info,LSU,$pos,$TSS_200_LSU_3prime_counts_by_positon{$gene}{$pos}\n";
    }
}
close(OUT6);


open (OUT7,">$out_file7")  || die "can't open $out_file7\n";
print OUT7 "gene_id,transcript_SSU_count,transcript_LSU_count,cds_SSU_count,cds_LSU_count,cds_SSU_FPKM,cds_LSU_FPKM,transcript_length,cds_length,gene_overlaps_another_gene,leader_overlap_an_upstream_gene,trailer_overlaps_a_downstream_gene,gene_overlaps_ncRNA,fraction,position_in_region,count\n";

for my $gene (keys %TSS_500_SSU_3prime_counts_by_positon){

    my $CDS_TCPseq_SSU_FPKM=eval{ $CDS_counts_SSU{$gene}/($CDS_length{$gene}*$total_SSU_bam_count)*1000000000 } || 0;
    my $CDS_TCPseq_LSU_FPKM=eval{ $CDS_counts_LSU{$gene}/($CDS_length{$gene}*$total_LSU_bam_count)*1000000000 } || 0;

    my $transcript_length=$transcript_length{$gene};
    my $cds_length=$CDS_length{$gene};

    my $gene_overlaps_another_gene=$overlaps_inframe_gene{$gene};
    my $leader_overlap_an_upstream_gene=$leader_overlaps_upstream{$gene};
    my $trailer_overlaps_a_downstream_gene=$gene_overlaps_downstream_leader{$gene};
    my $gene_overlaps_ncRNA=$gene_overlaps_ncRNA{$gene};

    my $gene_info="$gene,$transcript_counts_SSU{$gene},$transcript_counts_LSU{$gene},$CDS_counts_SSU{$gene},$CDS_counts_LSU{$gene},$CDS_TCPseq_SSU_FPKM,$CDS_TCPseq_LSU_FPKM,$transcript_length,$cds_length,$gene_overlaps_another_gene,$leader_overlap_an_upstream_gene,$trailer_overlaps_a_downstream_gene,$gene_overlaps_ncRNA";

    for my $pos (sort {$a <=> $b} keys %{$TSS_500_SSU_3prime_counts_by_positon{$gene}}){
        print OUT7 "$gene_info,SSU,$pos,$TSS_500_SSU_3prime_counts_by_positon{$gene}{$pos}\n";
    }

    for my $pos (sort {$a <=> $b} keys %{$TSS_500_LSU_3prime_counts_by_positon{$gene}}){
        print OUT7 "$gene_info,LSU,$pos,$TSS_500_LSU_3prime_counts_by_positon{$gene}{$pos}\n";
    }
}
close(OUT7);


#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Positional windows for whole regions transcripts # 3'
##¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
open (OUT8,">$out_file8")  || die "can't open $out_file8\n";
print OUT8 "gene_id,fraction,position_in_region,count\n";

for my $gene (keys %transcript_SSU_3prime_counts_by_positon){

    #SSU
    for my $pos (sort {$a <=> $b} keys %{$transcript_SSU_3prime_counts_by_positon{$gene}}){
         print OUT8 "$gene,SSU,$pos,$transcript_SSU_3prime_counts_by_positon{$gene}{$pos}\n";
    }

    #LSU
    for my $pos (sort {$a <=> $b} keys %{$transcript_LSU_3prime_counts_by_positon{$gene}}){
         print OUT8 "$gene,LSU,$pos,$transcript_LSU_3prime_counts_by_positon{$gene}{$pos}\n";
    }
}
close(OUT8);

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Positional windows for uORFS # 5'
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
open (OUT9,">$out_file9")  || die "can't open $out_file9\n";
print OUT9 "gene_id,transcript_SSU_count,transcript_LSU_count,cds_SSU_count,cds_LSU_count,cds_SSU_FPKM,cds_LSU_FPKM,transcript_length,cds_length,gene_overlaps_another_gene,leader_overlap_an_upstream_gene,trailer_overlaps_a_downstream_gene,gene_overlaps_ncRNA,fraction,position_in_region,count\n";

for my $gene (keys %TSS_100_SSU_5prime_counts_by_positon){

    my $CDS_TCPseq_SSU_FPKM=eval{ $CDS_counts_SSU{$gene}/($CDS_length{$gene}*$total_SSU_bam_count)*1000000000 } || 0;
    my $CDS_TCPseq_LSU_FPKM=eval{ $CDS_counts_LSU{$gene}/($CDS_length{$gene}*$total_LSU_bam_count)*1000000000 } || 0;

    my $transcript_length=$transcript_length{$gene};
    my $cds_length=$CDS_length{$gene};

    my $gene_overlaps_another_gene=$overlaps_inframe_gene{$gene};
    my $leader_overlap_an_upstream_gene=$leader_overlaps_upstream{$gene};
    my $trailer_overlaps_a_downstream_gene=$gene_overlaps_downstream_leader{$gene};
    my $gene_overlaps_ncRNA=$gene_overlaps_ncRNA{$gene};

    my $gene_info="$gene,$transcript_counts_SSU{$gene},$transcript_counts_LSU{$gene},$CDS_counts_SSU{$gene},$CDS_counts_LSU{$gene},$CDS_TCPseq_SSU_FPKM,$CDS_TCPseq_LSU_FPKM,$transcript_length,$cds_length,$gene_overlaps_another_gene,$leader_overlap_an_upstream_gene,$trailer_overlaps_a_downstream_gene,$gene_overlaps_ncRNA";

    for my $pos (sort {$a <=> $b} keys %{$TSS_100_SSU_5prime_counts_by_positon{$gene}}){
        print OUT9 "$gene_info,SSU,$pos,$TSS_100_SSU_5prime_counts_by_positon{$gene}{$pos}\n";
    }

    for my $pos (sort {$a <=> $b} keys %{$TSS_100_LSU_5prime_counts_by_positon{$gene}}){
        print OUT9 "$gene_info,LSU,$pos,$TSS_100_LSU_5prime_counts_by_positon{$gene}{$pos}\n";
    }
}
close(OUT9);


open (OUT10,">$out_file10")  || die "can't open $out_file10\n";
print OUT10 "gene_id,transcript_SSU_count,transcript_LSU_count,cds_SSU_count,cds_LSU_count,cds_SSU_FPKM,cds_LSU_FPKM,transcript_length,cds_length,gene_overlaps_another_gene,leader_overlap_an_upstream_gene,trailer_overlaps_a_downstream_gene,gene_overlaps_ncRNA,fraction,position_in_region,count\n";

for my $gene (keys %TSS_200_SSU_5prime_counts_by_positon){

    my $CDS_TCPseq_SSU_FPKM=eval{ $CDS_counts_SSU{$gene}/($CDS_length{$gene}*$total_SSU_bam_count)*1000000000 } || 0;
    my $CDS_TCPseq_LSU_FPKM=eval{ $CDS_counts_LSU{$gene}/($CDS_length{$gene}*$total_LSU_bam_count)*1000000000 } || 0;

    my $transcript_length=$transcript_length{$gene};
    my $cds_length=$CDS_length{$gene};

    my $gene_overlaps_another_gene=$overlaps_inframe_gene{$gene};
    my $leader_overlap_an_upstream_gene=$leader_overlaps_upstream{$gene};
    my $trailer_overlaps_a_downstream_gene=$gene_overlaps_downstream_leader{$gene};
    my $gene_overlaps_ncRNA=$gene_overlaps_ncRNA{$gene};

    my $gene_info="$gene,$transcript_counts_SSU{$gene},$transcript_counts_LSU{$gene},$CDS_counts_SSU{$gene},$CDS_counts_LSU{$gene},$CDS_TCPseq_SSU_FPKM,$CDS_TCPseq_LSU_FPKM,$transcript_length,$cds_length,$gene_overlaps_another_gene,$leader_overlap_an_upstream_gene,$trailer_overlaps_a_downstream_gene,$gene_overlaps_ncRNA";

    for my $pos (sort {$a <=> $b} keys %{$TSS_200_SSU_5prime_counts_by_positon{$gene}}){
        print OUT10 "$gene_info,SSU,$pos,$TSS_200_SSU_5prime_counts_by_positon{$gene}{$pos}\n";
    }

    for my $pos (sort {$a <=> $b} keys %{$TSS_200_LSU_5prime_counts_by_positon{$gene}}){
        print OUT10 "$gene_info,LSU,$pos,$TSS_200_LSU_5prime_counts_by_positon{$gene}{$pos}\n";
    }
}
close(OUT10);


open (OUT11,">$out_file11")  || die "can't open $out_file11\n";
print OUT11 "gene_id,transcript_SSU_count,transcript_LSU_count,cds_SSU_count,cds_LSU_count,cds_SSU_FPKM,cds_LSU_FPKM,transcript_length,cds_length,gene_overlaps_another_gene,leader_overlap_an_upstream_gene,trailer_overlaps_a_downstream_gene,gene_overlaps_ncRNA,fraction,position_in_region,count\n";

for my $gene (keys %TSS_500_SSU_5prime_counts_by_positon){

    my $CDS_TCPseq_SSU_FPKM=eval{ $CDS_counts_SSU{$gene}/($CDS_length{$gene}*$total_SSU_bam_count)*1000000000 } || 0;
    my $CDS_TCPseq_LSU_FPKM=eval{ $CDS_counts_LSU{$gene}/($CDS_length{$gene}*$total_LSU_bam_count)*1000000000 } || 0;

    my $transcript_length=$transcript_length{$gene};
    my $cds_length=$CDS_length{$gene};

    my $gene_overlaps_another_gene=$overlaps_inframe_gene{$gene};
    my $leader_overlap_an_upstream_gene=$leader_overlaps_upstream{$gene};
    my $trailer_overlaps_a_downstream_gene=$gene_overlaps_downstream_leader{$gene};
    my $gene_overlaps_ncRNA=$gene_overlaps_ncRNA{$gene};

    my $gene_info="$gene,$transcript_counts_SSU{$gene},$transcript_counts_LSU{$gene},$CDS_counts_SSU{$gene},$CDS_counts_LSU{$gene},$CDS_TCPseq_SSU_FPKM,$CDS_TCPseq_LSU_FPKM,$transcript_length,$cds_length,$gene_overlaps_another_gene,$leader_overlap_an_upstream_gene,$trailer_overlaps_a_downstream_gene,$gene_overlaps_ncRNA";

    for my $pos (sort {$a <=> $b} keys %{$TSS_500_SSU_5prime_counts_by_positon{$gene}}){
        print OUT11 "$gene_info,SSU,$pos,$TSS_500_SSU_5prime_counts_by_positon{$gene}{$pos}\n";
    }

    for my $pos (sort {$a <=> $b} keys %{$TSS_500_LSU_5prime_counts_by_positon{$gene}}){
        print OUT11 "$gene_info,LSU,$pos,$TSS_500_LSU_5prime_counts_by_positon{$gene}{$pos}\n";
    }
}
close(OUT11);


#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Positional windows for whole regions transcripts # 5'
##¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

open (OUT12,">$out_file12")  || die "can't open $out_file12\n";
print OUT12 "gene_id,fraction,position_in_region,count\n";

for my $gene (keys %transcript_SSU_5prime_counts_by_positon){

    #SSU
    for my $pos (sort {$a <=> $b} keys %{$transcript_SSU_5prime_counts_by_positon{$gene}}){
         print OUT12 "$gene,SSU,$pos,$transcript_SSU_5prime_counts_by_positon{$gene}{$pos}\n";
    }

    #LSU
    for my $pos (sort {$a <=> $b} keys %{$transcript_LSU_5prime_counts_by_positon{$gene}}){
         print OUT12 "$gene,LSU,$pos,$transcript_LSU_5prime_counts_by_positon{$gene}{$pos}\n";
    }
}
close(OUT12);

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

exit;

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
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
sub aggregate_to_100{

    my $lead_length = shift;
    my $positional_hash_reference = shift;
    #my $uorf = shift;
    my %positional_hash = %$positional_hash_reference;
    my %tmp_hash;
    my %out_hash;
    my $lead_scaling_factor=$lead_length/100;

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
sub aggregate_to_300{

    my $lead_length = shift;
    my $positional_hash_reference = shift;
    my %positional_hash = %$positional_hash_reference;
    my %tmp_hash;
    my %out_hash;
    my $lead_scaling_factor=$lead_length/300;

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

