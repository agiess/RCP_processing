#!/usr/bin/perl -w
use strict;

#to do 06/06/2019
#script to prrocess TCP-seq reads for 5' and 3' fragment length plots aorund start and stop codons
##do not scale the genes for read lengths, we will do this in R
#scale the bar charts though
#where a read overlaps multiple in strand genes counts it towards all those genes

#input files:
#gtf exons, start, stop codons
#sam = aligned 
my $gtf=$ARGV[0];
my $leaders=$ARGV[1];
my $rna_bam=$ARGV[2];
my $bam=$ARGV[3];
my $outDir=$ARGV[4];
my $most_highly_expressed=$ARGV[5];
my $use_cage=$ARGV[6];  # unless($use_cage eq "exclude")

my ($prefix)=$bam=~/([^\/]+).bam$/;

my $START_UPSTREAM=-100; 
my $START_DOWNSTREAM=500;  
my $STOP_UPSTREAM=-500; 
my $STOP_DOWNSTREAM=100;  

#exclude genes with leaders shorter than START_UPSTREAM
#exclude genes with trailers shorter than STOP_DOWNSTREAM
#exclude genes with CDSs shorter than START_DOWNSTREAM or STOP_UPSTREAM

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

open(GENES2,$gtf) || die "can't open $gtf";      #gft is 1 based
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

        #strip the chr prefixc
        if ($chr=~/chr(.*)/){ #ensembl chromosome names do not contain the "chr" prefix
            $chr=$1;
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
my %five_prime_most_coord_rev;
my %three_prime_most_coord_rev;

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

        #strip the chr prefix
        if ($chr=~/chr(.*)/){ #ensembl chromosome names do not contain the "chr" prefix
            $chr=$1;
        }

	if ($overlaps_inframe_gene eq "TRUE"){ $overlaps_inframe_gene{$gene}=1; }
	if ($leader_overlaps_upstream eq "TRUE"){ $leader_overlaps_upstream{$gene}=1; }
	if ($gene_overlaps_downstream_leader eq "TRUE"){ $gene_overlaps_downstream_leader{$gene}=1; }
        if ($gene_overlaps_nc eq "TRUE"){ $gene_overlaps_nc{$gene}=1; }

        unless ($leader_length eq "NaN"){  #only take genes that have a detectable cage peak

            unless ($leader_length < 0){  #exlude genes with negative leader sizes, as they cause problems with FPKM

                if ($dir eq "fwd"){
                    if (exists ($start_coord_fwd{$gene})){
                        unless ($highest_cage_peak >= $gene_model_fwd{$gene}{$start_coord_fwd{$gene}}){  #exclude genes where the TSS is downstream of the start codon
                            $leader_positions_fwd{$gene}=$highest_cage_peak;
                            $cage_peak_value{$gene}=$count_at_highest_cage_peak;
                            $leader_length{$gene}=$leader_length;
                        }
                    }
                }else{
                    if (exists ($start_coord_rev{$gene})){
                        unless ($highest_cage_peak <= $gene_model_rev{$gene}{$start_coord_rev{$gene}}){  #exclude genes where the TSS is downstream of the start codon
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
my %start_codons_search_fwd; #key1 = chr, key2 = position, value = gene_id(s) ; window position
my %start_codons_search_rev; 
my %stop_codons_search_fwd; 
my %stop_codons_search_rev; 

#for initalising
my %start_codon_meta_positions_5; #key1 = gene_id, key2=metapos, key3=value
my %start_codon_meta_positions_3;
my %stop_codon_meta_positions_5;
my %stop_codon_meta_positions_3;

my %start_region_signal_5; #gene_id = sum signal_at_start_region
my %start_region_signal_3; 
my %stop_region_signal_5; 
my %stop_region_signal_3; 

my %cds_search_fwd;
my %cds_search_rev;
my %cds_signal_3;
my %cds_signal_5;

my %cds_total_rna;
my %cds_length;

for my $gene (keys %gene_model_fwd){
    my $chr=$gene_2_chr{$gene};
    my $start_coord=$start_coord_fwd{$gene};
    my $stop_coord=$stop_coord_fwd{$gene};

    my $five_prime_coord=$five_prime_most_coord_fwd{$gene};
    my $three_prime_coord=$three_prime_most_coord_fwd{$gene};

    my $include_leader=0;
    my $include_trailer=0;
 
    #check if there is 100nt of leader.
    if (exists ($gene_model_fwd{$gene}{($start_coord+$START_UPSTREAM)})){  #-100
        if (($start_coord+$START_UPSTREAM) >= $five_prime_most_coord_fwd{$gene}) { #restrict to gene with annoated or CAGE updated leader of a suitable length
             $include_leader=1;
        }
    }

    #check if there is 100nt of trailer.
    if (exists ($gene_model_fwd{$gene}{($stop_coord+$STOP_DOWNSTREAM)})){  #+100
         $include_trailer=1;
    }

    #should I also exclude CDS's shorter than 500 nt?
    if ($stop_coord-$start_coord < $START_DOWNSTREAM){
        $include_leader=0;
        $include_trailer=0;
    }

    if (exists ($overlaps_inframe_gene{$gene})){
        $include_leader=0;
        $include_trailer=0;
    }

    if (exists ($gene_overlaps_nc{$gene})){
        $include_leader=0;
        $include_trailer=0;
    }

    if (exists ($leader_overlaps_upstream{$gene})){
        $include_leader=0;
    }

    if (exists ($gene_overlaps_downstream_leader{$gene})){
        $include_trailer=0;
    }

    #if ($gene eq "ENSDARG00000036180"){  #gene with large -15 peak
    if (($gene eq "ENSDARG00000036180") || ($gene eq "ENSDARG00000014790")){  #gene with large -15 + -31 peak
    #if (($gene eq "ENSDARG00000036180") || ($gene eq "ENSDARG00000014790") || ($gene eq "ENSDARG00000075113")){   #-15, -31, -3 peaks
        $include_leader=0;
        $include_trailer=0;
    }

    if ($include_leader && $include_trailer){

        $cds_signal_5{$gene}=0;
        $cds_signal_3{$gene}=0;

        $cds_total_rna{$gene}=0;
        $cds_length{$gene}=$three_prime_coord-$five_prime_coord;

        #for my $coord ($start_coord .. $stop_coord){          #store CDS position for ranking on expression
        for my $coord ($five_prime_coord .. $three_prime_coord){ #store transcript position for ranking on expression
            if (exists ($gene_model_fwd{$gene}{$coord})){
                my $genomic_coord=$gene_model_fwd{$gene}{($coord)};
                if (exists ($cds_search_fwd{$chr}{$genomic_coord})){
                    $cds_search_fwd{$chr}{$genomic_coord}.=";".$gene;  #linker == ";"
                }else{
                    $cds_search_fwd{$chr}{$genomic_coord}=$gene;
                }
            }
        }

        for my $pos ($START_UPSTREAM .. $START_DOWNSTREAM){   #-100 to 500
          
            $start_codon_meta_positions_5{$gene}{$pos}=0;
            $start_codon_meta_positions_3{$gene}{$pos}=0;
            $start_region_signal_5{$gene}=0;
            $start_region_signal_3{$gene}=0;

            my $genomic_coord=$gene_model_fwd{$gene}{($start_coord+$pos)};
            if (exists ($start_codons_search_fwd{$chr}{$genomic_coord})){
                $start_codons_search_fwd{$chr}{$genomic_coord}.=",".$gene.";".$pos;  #linker == ";"
            }else{
                $start_codons_search_fwd{$chr}{$genomic_coord}=$gene.";".$pos;
            }
        }

        for my $pos ($STOP_UPSTREAM .. $STOP_DOWNSTREAM){   #-500 to 100

            $stop_codon_meta_positions_5{$gene}{$pos}=0;
            $stop_codon_meta_positions_3{$gene}{$pos}=0;
            $stop_region_signal_5{$gene}=0;
            $stop_region_signal_3{$gene}=0;

            my $genomic_coord=$gene_model_fwd{$gene}{($stop_coord+$pos)};
            if (exists ($stop_codons_search_fwd{$chr}{$genomic_coord})){
                 $stop_codons_search_fwd{$chr}{$genomic_coord}.=",".$gene.";".$pos;  #linker == ";"
            }else{
                 $stop_codons_search_fwd{$chr}{$genomic_coord}=$gene.";".$pos;
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

    my $include_leader=0;
    my $include_trailer=0;

    #check if there is 100nt of leader.
    if (exists ($gene_model_rev{$gene}{($start_coord+$START_UPSTREAM)})){  #-100
        if (($start_coord + $START_UPSTREAM) >= $five_prime_most_coord_rev{$gene}){ #restrict to gene with annoated or CAGE updated leader of a suitable length
            $include_leader=1;
        }
    }

    #check if there is 100nt of trailer.
    if (exists ($gene_model_rev{$gene}{($stop_coord+$STOP_DOWNSTREAM)})){  #+100
         $include_trailer=1;
    }

    #should I also exclude CDS's shorter than 500 nt?
    if ($stop_coord-$start_coord < $START_DOWNSTREAM){
        $include_leader=0;
        $include_trailer=0;
    }

    if (exists ($gene_overlaps_nc{$gene})){
        $include_leader=0;
        $include_trailer=0;
    }

    if (exists ($overlaps_inframe_gene{$gene})){
        $include_leader=0;
        $include_trailer=0;
    }

    if (exists ($leader_overlaps_upstream{$gene})){
        $include_leader=0;
    }

    if (exists ($gene_overlaps_downstream_leader{$gene})){
        $include_trailer=0;
    }

    #if ($gene eq "ENSDARG00000036180"){  #gene with large -15 peak
    if (($gene eq "ENSDARG00000036180") || ($gene eq "ENSDARG00000014790")){  #gene with large -15 + -31 peak
    #if (($gene eq "ENSDARG00000036180") || ($gene eq "ENSDARG00000014790") || ($gene eq "ENSDARG00000075113")){   #-15, -31, -3 peaks
        $include_leader=0;
        $include_trailer=0;        
    }

    if ($include_leader && $include_trailer){

        $cds_signal_5{$gene}=0;
        $cds_signal_3{$gene}=0;

        $cds_total_rna{$gene}=0;
        $cds_length{$gene}=$three_prime_coord-$five_prime_coord;

        #for my $coord ($start_coord .. $stop_coord){          #store CDS position for ranking on expression
        for my $coord ($five_prime_coord .. $three_prime_coord){ #store transcript position for ranking on expression
            if (exists ($gene_model_rev{$gene}{$coord})){
                my $genomic_coord=$gene_model_rev{$gene}{($coord)};
                if (exists ($cds_search_rev{$chr}{$genomic_coord})){
                    $cds_search_rev{$chr}{$genomic_coord}.=";".$gene;  #linker == ";"
                }else{
                    $cds_search_rev{$chr}{$genomic_coord}=$gene;
                }
            }
        }

        for my $pos ($START_UPSTREAM .. $START_DOWNSTREAM){   #-100 to 500

            $start_codon_meta_positions_5{$gene}{$pos}=0;
            $start_codon_meta_positions_3{$gene}{$pos}=0;
            $start_region_signal_5{$gene}=0;
            $start_region_signal_3{$gene}=0;

            my $genomic_coord=$gene_model_rev{$gene}{($start_coord+$pos)};

            if (exists ($start_codons_search_rev{$chr}{$genomic_coord})){
                $start_codons_search_rev{$chr}{$genomic_coord}.=",".$gene.";".$pos;  #linker == ";"
            }else{
                $start_codons_search_rev{$chr}{$genomic_coord}=$gene.";".$pos;
            }
        }

        for my $pos ($STOP_UPSTREAM .. $STOP_DOWNSTREAM){   #-500 to 100

            $stop_codon_meta_positions_5{$gene}{$pos}=0;
            $stop_codon_meta_positions_3{$gene}{$pos}=0;
            $stop_region_signal_5{$gene}=0;
            $stop_region_signal_3{$gene}=0;

            my $genomic_coord=$gene_model_rev{$gene}{($stop_coord+$pos)};
            if (exists ($stop_codons_search_rev{$chr}{$genomic_coord})){
                 $stop_codons_search_rev{$chr}{$genomic_coord}.=",".$gene.";".$pos;  #linker == ";"
            }else{
                 $stop_codons_search_rev{$chr}{$genomic_coord}=$gene.";".$pos;
            }
        }
    }
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#open and store totalRNA counts
my $total_totalRNA_count=0;

open BAM4,"samtools view $rna_bam |";

while(<BAM4>){
  
  next if(/^(\@)/);  ## skipping the header lines (if you used -h in the samools command)
  s/\n//;  s/\r//;  ## removing new line
  my @sam = split(/\t+/);  ## splitting SAM line into array
  
  my $leftMost=$sam[3]; #leftmost position of match 5' for fwd, 3' for rev
  my $flag=$sam[1];
  my $chr=$sam[2];
  my $mapq=$sam[4];
  my $cigar=$sam[5];
  
#  my $region_CDS_hit_in=0;
  my $region_CDS_hit_off=0;
  
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
                 
#                if (exists ($CDS_region_search_rev{$chr}{$pos})){
#                  $region_CDS_hit_in=1;
#                }
               
                if (exists ($cds_search_fwd{$chr}{$pos})){
                    my @over1=split(";",$cds_search_fwd{$chr}{$pos});
                    for my $gene (@over1){
                        $region_CDS_hit_off=$gene;
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
               
                #if (exists ($cds_search_rev{$chr}{$pos})){
                #  $region_CDS_hit_off=$cds_search_rev{$chr}{$pos};
                #}

                if (exists ($cds_search_rev{$chr}{$pos})){
                    my @over1=split(";",$cds_search_rev{$chr}{$pos});
                    for my $gene (@over1){
                        $region_CDS_hit_off=$gene;
                    }
                }

#                if (exists ($CDS_region_search_fwd{$chr}{$pos})){
#                  $region_CDS_hit_in=1;
#                }
                
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
      
      if($region_CDS_hit_off){ 
        $total_totalRNA_count++;
        $cds_total_rna{$region_CDS_hit_off}++;
        #assign hit to gene
      }
    }
  }
}
close (BAM4);


#calculate FPKMS
#loop through transcripts
#if ($TSS_beginning_count_RNAseq > 0){ $TSS_beginning_RNAseq_FPKM=eval{ ($TSS_beginning_count_RNAseq/($total_RNAseq_aligned_sum))*1000000              } || 0; }
#if ($leader_count_RNAseq        > 0){ $leader_RNAseq_FPKM=eval{        ($leader_count_RNAseq/($leader_length*$total_RNAseq_aligned_sum))*1000000000   } || 0; }
#if ($start_codon_count_RNAseq   > 0){ $start_codon_RNAseq_FPKM=eval{   ($start_codon_count_RNAseq/($total_RNAseq_aligned_sum))*1000000                } || 0; }
#if ($CDS_count_RNAseq           > 0){ $CDS_RNAseq_FPKM=eval{           ($CDS_count_RNAseq/($CDS_length*$total_RNAseq_aligned_sum))*1000000000         } || 0; }
#if ($stop_codon_count_RNAseq    > 0){ $stop_codon_RNAseq_FPKM=eval{    ($stop_codon_count_RNAseq/($total_RNAseq_aligned_sum))*1000000                 } || 0; }
#if ($trailer_count_RNAseq       > 0){ $trailer_RNAseq_FPKM=eval{       ($trailer_count_RNAseq/($trailer_length*$total_RNAseq_aligned_sum))*1000000000 } || 0; }


#$CDS_RNAseq_FPKM=eval{ ($cds_total_rna{$gene}/($cds_length{$gene}*$total_totalRNA_count))*1000000000 } || 0; } #actually full transcript region
#if ($CDS_RNAseq_FPKM >=10){

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#open and assign TCPseq counts
my %start_codon_meta_positions_by_length_5; #gene, meta position (-50 to 500), length = count
my %start_codon_meta_positions_by_length_3;
my %stop_codon_meta_positions_by_length_5;
my %stop_codon_meta_positions_by_length_3;

my $longest_length=0;
my $shortest_length=1000;

open BAM,"samtools view $bam |";
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
    my $threePrime;
    my $fivePrime;

    if ($chr=~/chr(.*)/){ #ensembl chromosome names do not contain the "chr" prefix
        $chr=$1;
    }

    unless ($flag & 0x4){   #if aligned

        if ($mapq >= 10){     #mapping uniqnes filter

            if ($flag & 0x10){  #Reverse reads. Starting from the leftmost position parse the cigar and check if matching positions overlap leaders or cds's

                #assign 5' amd 3' to positions
                $threePrime=$leftMost;

                #parse cigar for indels and adjust the length of the alignment             
                my $length=length($seq);
                while ($cigar =~/(\d+)I/g){   #add to length for insertions
                    $length+=$1;
                }
                while ($cigar =~/(\d+)D/g){   #substact from length for deletions
                    $length-=$1;
                }

                if ($length > $longest_length){ $longest_length=$length; }
                if ($length < $shortest_length){ $shortest_length=$length; }

                $fivePrime=$leftMost+($length-1);              #SAM is 1 based

                #assign to metaplots
                if (exists ($start_codons_search_rev{$chr}{$threePrime})){
                    my @over1=split(",",$start_codons_search_rev{$chr}{$threePrime});
                    for (@over1){
                        my ($gene,$meta_pos)=$_=~/(^[^;]+);(-?\w+)$/;  #/(\w+)_(-?\w+)/;
                        $start_region_signal_3{$gene}++;
                        $start_codon_meta_positions_3{$gene}{$meta_pos}+=1;
                        $start_codon_meta_positions_by_length_3{$gene}{$meta_pos}{length($seq)}+=1;
                    }
                }

                if (exists ($stop_codons_search_rev{$chr}{$threePrime})){
                    my @over1=split(",",$stop_codons_search_rev{$chr}{$threePrime});
                    for (@over1){
                        my ($gene,$meta_pos)=$_=~/(^[^;]+);(-?\w+)$/;  #/(\w+)_(-?\w+)/;
                        $stop_region_signal_3{$gene}++;
                        $stop_codon_meta_positions_3{$gene}{$meta_pos}+=1;
                        $stop_codon_meta_positions_by_length_3{$gene}{$meta_pos}{length($seq)}+=1;
                    }
                }

                if (exists ($cds_search_rev{$chr}{$threePrime})){
                    my @over1=split(";",$cds_search_rev{$chr}{$threePrime});
                    for my $gene (@over1){
                        $cds_signal_3{$gene}++;
                    }
                }

                #same again for the 5 prime
                if (exists ($start_codons_search_rev{$chr}{$fivePrime})){
                    my @over1=split(",",$start_codons_search_rev{$chr}{$fivePrime});
                    for (@over1){
                        my ($gene,$meta_pos)=$_=~/(^[^;]+);(-?\w+)$/;  #/(\w+)_(-?\w+)/;
                        $start_region_signal_5{$gene}++;
                        $start_codon_meta_positions_5{$gene}{$meta_pos}+=1;
                        $start_codon_meta_positions_by_length_5{$gene}{$meta_pos}{length($seq)}+=1;
                    }
                }

                if (exists ($stop_codons_search_rev{$chr}{$fivePrime})){
                    my @over1=split(",",$stop_codons_search_rev{$chr}{$fivePrime});
                    for (@over1){
                        my ($gene,$meta_pos)=$_=~/(^[^;]+);(-?\w+)$/;  #/(\w+)_(-?\w+)/;
                        $stop_region_signal_5{$gene}++;
                        $stop_codon_meta_positions_5{$gene}{$meta_pos}+=1;
                        $stop_codon_meta_positions_by_length_5{$gene}{$meta_pos}{length($seq)}+=1;
                    }
                }

               if (exists ($cds_search_rev{$chr}{$fivePrime})){
                    my @over1=split(";",$cds_search_rev{$chr}{$fivePrime});
                    for my $gene (@over1){
                        $cds_signal_5{$gene}++;
                    }
                }

            }else{ #if fwd 3' == sam coordinate (leftmost) + read length 

                #parse cigar for indels and adjust the length of the alignment             
                my $length=length($seq);
                while ($cigar =~/(\d+)I/g){   #add to length for insertions
                    $length+=$1;
                }
                while ($cigar =~/(\d+)D/g){   #substact from length for deletions
                    $length-=$1;
                }

                if ($length > $longest_length){ $longest_length=$length; }
                if ($length < $shortest_length){ $shortest_length=$length; }

                $threePrime=$leftMost+($length-1);              #SAM is 1 based
                $fivePrime=$leftMost;

                #assign 3' to metaplots
                if (exists ($start_codons_search_fwd{$chr}{$threePrime})){
                    my @over4=split(",",$start_codons_search_fwd{$chr}{$threePrime});
                    for (@over4){
                        my ($gene,$meta_pos)=$_=~/(^[^;]+);(-?\w+$)/;
                        $start_region_signal_3{$gene}++;
                        $start_codon_meta_positions_3{$gene}{$meta_pos}+=1;    
                        $start_codon_meta_positions_by_length_3{$gene}{$meta_pos}{length($seq)}+=1;
                    }
                }

                if (exists ($stop_codons_search_fwd{$chr}{$threePrime})){
                    my @over4=split(",",$stop_codons_search_fwd{$chr}{$threePrime});
                    for (@over4){
                        my ($gene,$meta_pos)=$_=~/(^[^;]+);(-?\w+$)/;
                        $stop_region_signal_3{$gene}++;
                        $stop_codon_meta_positions_3{$gene}{$meta_pos}+=1;
                        $stop_codon_meta_positions_by_length_3{$gene}{$meta_pos}{length($seq)}+=1;
                    }
                }

                if (exists ($cds_search_fwd{$chr}{$threePrime})){
                    my @over1=split(";",$cds_search_fwd{$chr}{$threePrime});
                    for my $gene (@over1){
                        $cds_signal_3{$gene}++;
                    }
                }
  
                #do that again for five prime
                if (exists ($start_codons_search_fwd{$chr}{$fivePrime})){
                    my @over4=split(",",$start_codons_search_fwd{$chr}{$fivePrime});
                    for (@over4){
                        my ($gene,$meta_pos)=$_=~/(^[^;]+);(-?\w+$)/;
                        $start_region_signal_5{$gene}++;
                        $start_codon_meta_positions_5{$gene}{$meta_pos}+=1;
                        $start_codon_meta_positions_by_length_5{$gene}{$meta_pos}{length($seq)}+=1;
                    }
                }

                if (exists ($stop_codons_search_fwd{$chr}{$fivePrime})){
                    my @over4=split(",",$stop_codons_search_fwd{$chr}{$fivePrime});
                    for (@over4){
                        my ($gene,$meta_pos)=$_=~/(^[^;]+);(-?\w+$)/;
                        $stop_region_signal_5{$gene}++;
                        $stop_codon_meta_positions_5{$gene}{$meta_pos}+=1;
                        $stop_codon_meta_positions_by_length_5{$gene}{$meta_pos}{length($seq)}+=1;
                    }
                } 

                if (exists ($cds_search_fwd{$chr}{$fivePrime})){
                    my @over1=split(";",$cds_search_fwd{$chr}{$fivePrime});
                    for my $gene (@over1){
                        $cds_signal_5{$gene}++;
                    }
                }
            }
        }
    }                    
}
close (BAM);

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#scale by the counts within the start or stop window by the total counts in the window, per gene

my %barchart_upstream_5_q1;
my %barchart_upstream_3_q1;
my %barchart_downstream_5_q1;
my %barchart_downstream_3_q1;

#find the gene with the highest SSU 3' -15
#   exclude gene

my $highest_peak="";
my $highest_value=0;
my $highest_so_far=0;


my $highest_peak_31="";
my $highest_value_31=0;
my $highest_so_far_31=0;


my $highest_peak_3="";
my $highest_value_3=0;
my $highest_so_far_3=0;



### ENSDARG00000036180

#transcripts with strong 3' UTR peaks:
#	ENSDARG00000077330, 
#	ENSDARG00000102873,
#	ENSDARG00000089382

#Add counts of passed transcriots
my $included_gene_5=0;
my $included_gene_3=0;

for my $gene_id (keys %start_codon_meta_positions_5){
    my $CDS_RNAseq_FPKM=eval{ ($cds_total_rna{$gene_id}/($cds_length{$gene_id}*$total_totalRNA_count))*1000000000 } || 0;  #actually full transcript region
    if ($CDS_RNAseq_FPKM >=10){
        $included_gene_5++;
        for my $pos (sort { $a <=> $b } keys %{ $start_codon_meta_positions_5{$gene_id}} ){
            $barchart_upstream_5_q1{$pos}+=$start_codon_meta_positions_5{$gene_id}{$pos};  #barchart positions are summed over gene
        }
    }
}

for my $gene_id (keys %start_codon_meta_positions_3){
    my $CDS_RNAseq_FPKM=eval{ ($cds_total_rna{$gene_id}/($cds_length{$gene_id}*$total_totalRNA_count))*1000000000 } || 0;  #actually full transcript region
    if ($CDS_RNAseq_FPKM >=10){
        $included_gene_3++;
        for my $pos (sort { $a <=> $b } keys %{ $start_codon_meta_positions_3{$gene_id}} ){

            if ($pos == -15){ #find genes resonosible for peak
                if ($start_codon_meta_positions_3{$gene_id}{$pos} > $highest_so_far){ 
                    $highest_peak=$gene_id;
                    $highest_value=$start_codon_meta_positions_3{$gene_id}{$pos};
                    $highest_so_far=$start_codon_meta_positions_3{$gene_id}{$pos};
                }
            }

            if ($pos == -3){ #-3 find genes resonosible for peak
                if ($start_codon_meta_positions_3{$gene_id}{$pos} > $highest_so_far_3){
                    $highest_peak_3=$gene_id;
                    $highest_value_3=$start_codon_meta_positions_3{$gene_id}{$pos};
                    $highest_so_far_3=$start_codon_meta_positions_3{$gene_id}{$pos};
                }
            }

            if ($pos == -31){ #-31 find genes resonosible for peak
                if ($start_codon_meta_positions_3{$gene_id}{$pos} > $highest_so_far_31){
                    $highest_peak_31=$gene_id;
                    $highest_value_31=$start_codon_meta_positions_3{$gene_id}{$pos};
                    $highest_so_far_31=$start_codon_meta_positions_3{$gene_id}{$pos};
                }
            }

            $barchart_upstream_3_q1{$pos}+=$start_codon_meta_positions_3{$gene_id}{$pos};  #barchart positions
        }
    }
}

for my $gene_id (keys %stop_codon_meta_positions_5){
    my $CDS_RNAseq_FPKM=eval{ ($cds_total_rna{$gene_id}/($cds_length{$gene_id}*$total_totalRNA_count))*1000000000 } || 0;  #actually full transcript region
    if ($CDS_RNAseq_FPKM >=10){
        for my $pos (sort { $a <=> $b } keys %{ $stop_codon_meta_positions_5{$gene_id}} ){
            $barchart_downstream_5_q1{$pos}+=$stop_codon_meta_positions_5{$gene_id}{$pos};  #barchart positions are summed over genes
        }
    }
}

for my $gene_id (keys %stop_codon_meta_positions_3){
    my $CDS_RNAseq_FPKM=eval{ ($cds_total_rna{$gene_id}/($cds_length{$gene_id}*$total_totalRNA_count))*1000000000 } || 0;  #actually full transcript region
    if ($CDS_RNAseq_FPKM >=10){
        for my $pos (sort { $a <=> $b } keys %{ $stop_codon_meta_positions_3{$gene_id}} ){
            $barchart_downstream_3_q1{$pos}+=$stop_codon_meta_positions_3{$gene_id}{$pos};  #barchart positions
        }
    }
}

#scale lengths
my %scaled_upstream_3_q1;
my %scaled_upstream_5_q1;
my %scaled_downstream_3_q1;
my %scaled_downstream_5_q1;

#for my $pos ($START_UPSTREAM .. $START_DOWNSTREAM){   #-100 to 500
#    for my $length ($shortest_length .. $longest_length){
#        $scaled_upstream_3_q1{$pos}{$length}=0;
#        $scaled_upstream_5_q1{$pos}{$length}=0;
#        $scaled_upstream_3_q2{$pos}{$length}=0;
#        $scaled_upstream_5_q2{$pos}{$length}=0;
#        $scaled_upstream_3_q3{$pos}{$length}=0;
#        $scaled_upstream_5_q3{$pos}{$length}=0;
#        $scaled_upstream_3_q4{$pos}{$length}=0;
#        $scaled_upstream_5_q4{$pos}{$length}=0;
#    }
#}

#for my $pos ($STOP_UPSTREAM .. $STOP_DOWNSTREAM){
#    for my $length ($shortest_length .. $longest_length){
#        $scaled_downstream_3_q1{$pos}{$length}=0;
#        $scaled_downstream_5_q1{$pos}{$length}=0;
#        $scaled_downstream_3_q2{$pos}{$length}=0;
#        $scaled_downstream_5_q2{$pos}{$length}=0;
#        $scaled_downstream_3_q3{$pos}{$length}=0;
#        $scaled_downstream_5_q3{$pos}{$length}=0;
#        $scaled_downstream_3_q4{$pos}{$length}=0;
#        $scaled_downstream_5_q4{$pos}{$length}=0;
#    }
#}

print "/n/n/nThere are $included_gene_5 genes in the final 5' plots/n";
print "There are $included_gene_3 genes in the final 3' plots/n/n/n";

print "\n\n\ngene with highest -15 3' peak: $highest_peak\n";
print "peak value: $highest_value\n\n\n";

print "\n\n\ngene with highest -3 3' peak: $highest_peak_3\n";
print "peak value: $highest_value_3\n\n\n";

print "\n\n\ngene with highest -31 3' peak: $highest_peak_31\n";
print "peak value: $highest_value_31\n\n\n";

for my $gene (sort keys %start_codon_meta_positions_by_length_5){
    my $CDS_RNAseq_FPKM=eval{ ($cds_total_rna{$gene}/($cds_length{$gene}*$total_totalRNA_count))*1000000000 } || 0;  #actually full transcript region
    if ($CDS_RNAseq_FPKM >=10){
        for my $pos (sort {$a <=> $b} keys %{$start_codon_meta_positions_by_length_5{$gene}}){
            for my $length (sort {$a <=> $b} keys %{$start_codon_meta_positions_by_length_5{$gene}{$pos}}){
                $scaled_upstream_5_q1{$pos}{$length}+=$start_codon_meta_positions_by_length_5{$gene}{$pos}{$length};
            }
        }
    }
}

for my $gene (sort keys %start_codon_meta_positions_by_length_3){   
    my $CDS_RNAseq_FPKM=eval{ ($cds_total_rna{$gene}/($cds_length{$gene}*$total_totalRNA_count))*1000000000 } || 0;  #actually full transcript region
    if ($CDS_RNAseq_FPKM >=10){
        for my $pos (sort {$a <=> $b} keys %{$start_codon_meta_positions_by_length_3{$gene}}){
            for my $length (sort {$a <=> $b} keys %{$start_codon_meta_positions_by_length_3{$gene}{$pos}}){
                $scaled_upstream_3_q1{$pos}{$length}+=$start_codon_meta_positions_by_length_3{$gene}{$pos}{$length}; 
            }
        }
    }
}

for my $gene (sort keys %stop_codon_meta_positions_by_length_5){
    my $CDS_RNAseq_FPKM=eval{ ($cds_total_rna{$gene}/($cds_length{$gene}*$total_totalRNA_count))*1000000000 } || 0;  #actually full transcript region
    if ($CDS_RNAseq_FPKM >=10){
        for my $pos (sort {$a <=> $b} keys %{$stop_codon_meta_positions_by_length_5{$gene}}){
            for my $length (sort {$a <=> $b} keys %{$stop_codon_meta_positions_by_length_5{$gene}{$pos}}){
                 $scaled_downstream_5_q1{$pos}{$length}+=$stop_codon_meta_positions_by_length_5{$gene}{$pos}{$length};
            }
        }
    }
}

for my $gene (sort keys %stop_codon_meta_positions_by_length_3){
    my $CDS_RNAseq_FPKM=eval{ ($cds_total_rna{$gene}/($cds_length{$gene}*$total_totalRNA_count))*1000000000 } || 0;  #actually full transcript region
    if ($CDS_RNAseq_FPKM >=10){
        for my $pos (sort {$a <=> $b} keys %{$stop_codon_meta_positions_by_length_3{$gene}}){
            for my $length (sort {$a <=> $b} keys %{$stop_codon_meta_positions_by_length_3{$gene}{$pos}}){
                $scaled_downstream_3_q1{$pos}{$length}+=$stop_codon_meta_positions_by_length_3{$gene}{$pos}{$length};
            }
        }
    }
}

#for my $gene (sort keys %stop_codon_meta_positions_by_length_3){
#    if (exists ($quantile1_3{$gene})){
#        for my $pos (sort {$a <=> $b} keys %{$stop_codon_meta_positions_by_length_3{$gene}}){
#            for my $length (sort {$a <=> $b} keys %{$stop_codon_meta_positions_by_length_3{$gene}{$pos}}){
#                my $scaled_count=eval { $stop_codon_meta_positions_by_length_3{$gene}{$pos}{$length}/$stop_region_signal_3{$gene}} || 0 ;
#                $scaled_downstream_3_q1{$pos}{$length}+=$scaled_count;
#            }
#        }
#    }
#}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
my $outStart_q1=$outDir."/".$prefix."_start_3prime_scale_q1.csv";
my $outStop_q1=$outDir."/".$prefix."_stop_3prime_scale_q1.csv";
my $outStartLengthsScale_q1=$outDir."/".$prefix."_start_lengths_scale_3prime_q1.csv";
my $outStopLengthsScale_q1=$outDir."/".$prefix."_stop_lengths_scale_3prime_q1.csv";

open (OUT1,">$outStart_q1") || die "can't open $outStart_q1\n";
open (OUT2,">$outStop_q1")  || die "can't open $outStop_q1\n";
open (OUT3,">$outStartLengthsScale_q1") || die "can't open $outStartLengthsScale_q1\n";
open (OUT4,">$outStopLengthsScale_q1")  || die "can't open $outStopLengthsScale_q1\n";

# meta plots start codons
for my $pos (sort {$a <=> $b} keys %barchart_upstream_3_q1) { print OUT1  "$pos,$barchart_upstream_3_q1{$pos}\n"; }

# meta_plots stop codons
for my $pos (sort {$a <=> $b} keys %barchart_downstream_3_q1){ print OUT2  "$pos,$barchart_downstream_3_q1{$pos}\n"; }

# output scaled length values
for my $pos (sort {$a <=> $b} keys %scaled_upstream_3_q1){ for my $length (sort {$a <=> $b} keys %{$scaled_upstream_3_q1{$pos}}){ print OUT3  "$pos\t$length\t$scaled_upstream_3_q1{$pos}{$length}\n"; } }

for my $pos (sort {$a <=> $b} keys %scaled_downstream_3_q1){ for my $length (sort {$a <=> $b} keys %{$scaled_downstream_3_q1{$pos}}){  print OUT4  "$pos\t$length\t$scaled_downstream_3_q1{$pos}{$length}\n"; } }

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
$outStart_q1=$outDir."/".$prefix."_start_5prime_scale_q1.csv";
$outStop_q1=$outDir."/".$prefix."_stop_5prime_scale_q1.csv";
$outStartLengthsScale_q1=$outDir."/".$prefix."_start_lengths_scale_5prime_q1.csv";
$outStopLengthsScale_q1=$outDir."/".$prefix."_stop_lengths_scale_5prime_q1.csv";

open (OUT17,">$outStart_q1") || die "can't open $outStart_q1\n";
open (OUT18,">$outStop_q1")  || die "can't open $outStop_q1\n";
open (OUT19,">$outStartLengthsScale_q1") || die "can't open $outStartLengthsScale_q1\n";
open (OUT20,">$outStopLengthsScale_q1")  || die "can't open $outStopLengthsScale_q1\n";

# meta plots start codons
for my $pos (sort {$a <=> $b} keys %barchart_upstream_5_q1) { print OUT17 "$pos,$barchart_upstream_5_q1{$pos}\n"; }

# meta_plots stop codons
for my $pos (sort {$a <=> $b} keys %barchart_downstream_5_q1){ print OUT18 "$pos,$barchart_downstream_5_q1{$pos}\n"; }

# output scaled length values
for my $pos (sort {$a <=> $b} keys %scaled_upstream_5_q1){ for my $length (sort {$a <=> $b} keys %{$scaled_upstream_5_q1{$pos}}){ print OUT19 "$pos\t$length\t$scaled_upstream_5_q1{$pos}{$length}\n"; } }

for my $pos (sort {$a <=> $b} keys %scaled_downstream_5_q1){  for my $length (sort {$a <=> $b} keys %{$scaled_downstream_5_q1{$pos}}){  print OUT20 "$pos\t$length\t$scaled_downstream_5_q1{$pos}{$length}\n"; } }

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

exit;
