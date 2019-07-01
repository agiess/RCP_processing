#!/usr/bin/perl -w
use strict;
#28/09/2018

#script to take normalised cage counts from cageR, and assign leader regions to the most highly expressed longest transcript of each gene.

my $inGtf=$ARGV[0];
my $cage_fwd=$ARGV[1];
my $cage_rev=$ARGV[2];
my $most_highly_expressed=$ARGV[3];
my $out_leaders=$ARGV[4];

my $UPSTEAM_DISTANCE=1000; #the distance to extend from the start codon, to search for a cage peak (zebrafish)
#my $UPSTEAM_DISTANCE=500; #the distance to extend from the start codon, to search for a cage peak (yeast)

#restrictions
#gene must have gene_id and transcript_id
#gene must have protien_coding biotype  	#possibly better to wait until after overlapps to filter this
#Longest transcript of each gene taken forward
#the highest peak is selected acros the whole transcript upto 500nt upstream (or until an upstream gene is encountered)

#Marking:
#genes whose longest transcript overlaps with another genes longest transcript, on the same strand
#genes whose leaders potentially run into upstream genes on the same strand (within 5000bp)
#genes that are potentially overlapped by the leaders of downstream genes

#Possible additions:
#distance to upstream genes

#Suggested filtering criteria
#permissive. 1 normalised cage tags (~5 raw cage tags)
#restrictive: 10 normalised cage tags (~50 raw cage tags)
#highest_cage_peak="NaN";
#leader_length="NaN";

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#open the most highly expressed trasncript list and store 
my %most_expressed_transcript; #key = gene_id, transcript_id, #value = sum_exon_lengths;

open(EXP, $most_highly_expressed) || die "can't open $most_highly_expressed";     
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
#open gtf and get transcript lengths
my %nc_exons_fwd; #key = gene_id, transcript_id, #value = sum_exon_lengths;
my %nc_exons_rev; #key = gene_id, transcript_id, #value = sum_exon_lengths;
my %ncgene_2_chr;

open(GENES1,$inGtf) || die "can't open $inGtf";      #gft is 1 based
while (<GENES1>){
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
            if ($gene_biotype ne "protein_coding"){ #store the locations of non protien coding exons for overlaps

                $ncgene_2_chr{$gene_id}=$chr;

                if ($class eq "exon"){
                    if ($dir eq "+"){
                        $nc_exons_fwd{$gene_id}{$start}=$end;
                    }else{
                        $nc_exons_rev{$gene_id}{$start}=$end;
                    }
                }
            }
        }
    }
}
close (GENES1);

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
#open and store cage counts
my %fwd_tags;

open (CTS1, $cage_fwd) || die "can't open $cage_fwd\n";
while(<CTS1>){
    unless (/^track/){
        my @c=split();
        my $chr=$c[0];
        my $start=$c[1]+1; #bedgraph start is zero based
        my $end=$c[2];
        my $value=$c[3];

        #This is gtf dependant (required for RG64)
        if ($chr=~/chr(.*)/){ #ensembl chromosome names do not contain the "chr" prefix
           $chr=$1;
        }

        for ($start .. $end){
            $fwd_tags{$chr}{$_}=$value;
        }
    }
}
close(CTS1);

my %rev_tags;

open (CTS2, $cage_rev) || die "can't open $cage_rev\n";
while(<CTS2>){
    unless (/^track/){
        my @c=split();
        my $chr=$c[0];
        my $start=$c[1]+1; #bedgraph start is zero based
        my $end=$c[2];
        my $value=$c[3];

        #This is gtf dependant (required for RG64)
        if ($chr=~/chr(.*)/){ #ensembl chromosome names do not contain the "chr" prefix
            $chr=$1;
        }

        if ($value =~ /^\-(.*)$/){  #remove the negative sign
            $value=$1;
        }

        for ($start .. $end){
            $rev_tags{$chr}{$_}=$value;
        }
    }
}
close(CTS2);


#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#calculate 5 and 3 prime positions of non coding genes
my %nc_three_prime_most_position_fwd;
my %nc_three_prime_most_position_rev;
my %nc_five_prime_most_position_fwd;
my %nc_five_prime_most_position_rev;

for my $gene (keys %nc_exons_fwd){

    my $chr=$ncgene_2_chr{$gene};
    my $three_prime_pos=0;
    my $five_prime_pos=0;
 
    for my $exon_start (sort {$a <=> $b} keys %{ $nc_exons_fwd{$gene} } ){

        unless ($five_prime_pos){ 
            $five_prime_pos=$exon_start; #take the first exon
        }
 
        my $exon_end=$nc_exons_fwd{$gene}{$exon_start};
        $three_prime_pos=$exon_end; #take the last exon
    }

    $nc_three_prime_most_position_fwd{$chr}{$three_prime_pos}=$gene;
    $nc_five_prime_most_position_fwd{$chr}{$five_prime_pos}=$gene;
}

for my $gene (keys %nc_exons_rev){

    my $chr=$ncgene_2_chr{$gene};
    my $three_prime_pos=0;
    my $five_prime_pos=0;
 
    for my $exon_end (reverse (sort {$a <=> $b} keys %{ $nc_exons_rev{$gene} } )){

        unless ($three_prime_pos){
            $three_prime_pos=$exon_end; #take the first exon
        }

        my $exon_start=$nc_exons_rev{$gene}{$exon_end};
        $five_prime_pos=$exon_start; #take the last exon
    }

    $nc_three_prime_most_position_rev{$chr}{$three_prime_pos}=$gene;
    $nc_five_prime_most_position_rev{$chr}{$five_prime_pos}=$gene;
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#setup trascript models
my %gene_model_fwd;
my %start_coord_fwd;
my %stop_coord_fwd;

my %gene_overlaps_nc;   

#5' is the #gene_model{$gene}{0}
#3' is the last coord
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

        my $chr=$gene_2_chr{$gene};
        my $model_pos=0;

        for my $exon_start (sort {$a <=> $b} keys %{ $gene_exons_fwd{$gene} } ){
            my $exon_end=$gene_exons_fwd{$gene}{$exon_start};

            #fwd exons are in ascending order
            # start(-1)-> 100958 100975
            #             101077 101715 <-end(+1)

            for ($exon_start .. $exon_end){

                ##check for direct overlaps with ncRNA here. Best performed after leader extension
                #if (exists ($nc_three_prime_most_position_fwd{$chr}{$_})){ $gene_overlaps_nc{$gene=1}; }
                #if (exists ($nc_five_prime_most_position_fwd{$chr}{$_})){ $gene_overlaps_nc{$gene=1}; }

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

for my $gene (keys %gene_exons_rev){
    if ( (exists ($gene_start_codon_rev{$gene})) && (exists ($gene_stop_codon_rev{$gene})) ) { #restrict to genes with annotated start + stop codon

        my $chr=$gene_2_chr{$gene};
        my $model_pos=0;

        for my $exon_end (reverse (sort {$a <=> $b} keys %{ $gene_exons_rev{$gene} } )){
            my $exon_start=$gene_exons_rev{$gene}{$exon_end};

            #rev exons are sorted in decending order  
            #           447087 447794 <-start(+1)
            # end(-1)-> 446060 446254

            while ($exon_start >= $exon_end){

                #check for direct overlaps with ncRNA here #best done after leader extension
                #if (exists ($nc_three_prime_most_position_rev{$chr}{$_})){ $gene_overlaps_nc{$gene=1}; }
                #if (exists ($nc_five_prime_most_position_rev{$chr}{$_})){ $gene_overlaps_nc{$gene=1}; }

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

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#find overlapping genes and mark them for removal
#TCP-seq is directional, so I only need to exlude genes that overlap with another gene on the same strand

#Classes:
#  This gene overlap with another longest trancript of a gene on the same strand
#  This gennes leader could overlap with an upstream gene (within 500nt search)
#  This Gene could overlap with downstream leader (within 500nt search window)

my %five_primes_fwd;
my %three_primes_fwd;  #store the three and five prime genomic coords
my %five_primes_rev;
my %three_primes_rev;
my %gene_overlaps;   #shares a 5 or 3 prime with another gene?

#I need to add other types of ncRNA here

for my $gene (keys %three_prime_most_coord_fwd){
    my $chr=$gene_2_chr{$gene};
    my $coord5=0;
    my $coord3=$three_prime_most_coord_fwd{$gene};
    my $pos5=$gene_model_fwd{$gene}{$coord5};
    my $pos3=$gene_model_fwd{$gene}{$coord3};
  
    unless (exists ($five_primes_fwd{$chr}{$pos5})){
        $five_primes_fwd{$chr}{$pos5}=$gene;
    }else{
         for my $pre_gene (split(";",$five_primes_fwd{$chr}{$pos5})){
             $gene_overlaps{$pre_gene}=1; #mark the genes that it matches with
         }
         $gene_overlaps{$gene}=1; #mark this gene
         $five_primes_fwd{$chr}{$pos5}.=";".$gene;  #setup search-able genomic coordinates
    }

    unless (exists ($three_primes_fwd{$chr}{$pos3})){
        $three_primes_fwd{$chr}{$pos3}=$gene;
    }else{
         for my $pre_gene (split(";",$three_primes_fwd{$chr}{$pos3})){
             $gene_overlaps{$pre_gene}=1; #mark the genes that it matches with
         }
         $gene_overlaps{$gene}=1; #mark this gene
         $three_primes_fwd{$chr}{$pos3}.=";".$gene;  #setup search-able genomic coordinates
    }
}

for my $gene (keys %three_prime_most_coord_rev){
    my $chr=$gene_2_chr{$gene};
    my $coord5=0;
    my $coord3=$three_prime_most_coord_rev{$gene};
    my $pos5=$gene_model_rev{$gene}{$coord5};
    my $pos3=$gene_model_rev{$gene}{$coord3};

    unless (exists ($five_primes_rev{$chr}{$pos5})){
        $five_primes_rev{$chr}{$pos5}=$gene;
    }else{
         for my $pre_gene (split(";",$five_primes_rev{$chr}{$pos5})){
             $gene_overlaps{$pre_gene}=1; #mark the genes that it matches with
         }
         $gene_overlaps{$gene}=1; #mark this gene
         $five_primes_rev{$chr}{$pos5}.=";".$gene;  #setup search-able genomic coordinates
    }

    unless (exists ($three_primes_rev{$chr}{$pos3})){
        $three_primes_rev{$chr}{$pos3}=$gene;
    }else{
         for my $pre_gene (split(";",$three_primes_rev{$chr}{$pos3})){
             $gene_overlaps{$pre_gene}=1; #mark the genes that it matches with
         }
         $gene_overlaps{$gene}=1; #mark this gene
         $three_primes_rev{$chr}{$pos3}.=";".$gene;  #setup search-able genomic coordinates
    }
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#assign leaders here  
#search regions is the whole transcript + 500bp upstream (stopping if an upstream gene is encountered)

my %cage_highest; #key = gene, value = highest cage peak
my %cage_peak; #key = gene, value = position of highest cage peak

my %cage_sum_leader;
my %rna_sum_leader;
my %cage_sum_CDS;
my %rna_sum_CDS;
my %leader_length;
my %CDS_length;
my %length_count;

my %leader_overlaps_upstream_gene;
my %gene_overlaps_downstream_leader;

for my $gene (keys %gene_model_fwd){
    my $coordFivePrime=$start_coord_fwd{$gene};   
    my $coordThreePrime=$stop_coord_fwd{$gene};
    my $chr=$gene_2_chr{$gene};
    $cage_highest{$gene}=0; #initialise
    my $coord=0;
    my $pos=0;
    my $length_count=0;
    my $extended_coord;

    #search through the transcript up until the start codon
    for $coord (reverse (sort {$a <=> $b} keys %{ $gene_model_fwd{$gene} } )){

        $pos=$gene_model_fwd{$gene}{$coord};
        $length_count++;

        #check for upstream gene (overlapping with CDS)
        if (exists ($five_primes_fwd{$chr}{$pos})){
            for my $alt_gene (split(";",$five_primes_fwd{$chr}{$pos})){
                unless ($alt_gene eq $gene){  #not the stop codon of this gene
                    $gene_overlaps{$gene}=1;
                    $gene_overlaps{$alt_gene}=1;
                }
            }
        }

        #check for upstream Gene (overlapping with CDS)
        if (exists ($three_primes_fwd{$chr}{$pos})){
            for my $alt_gene (split(";",$three_primes_fwd{$chr}{$pos})){
                unless ($alt_gene eq $gene){  #not the start codon of the current gene
                    $gene_overlaps{$gene}=1;
                    $gene_overlaps{$alt_gene}=1;
                }
            }
        }
 
        #check for nc overlap
        if (exists ($nc_three_prime_most_position_fwd{$chr}{$pos})){ $gene_overlaps_nc{$gene}=1; }
        if (exists ($nc_five_prime_most_position_fwd{$chr}{$pos})){ $gene_overlaps_nc{$gene}=1; }
                                                                                               
        if (exists ($fwd_tags{$chr}{$pos})){
            if ($fwd_tags{$chr}{$pos} >= 1){ #minimum cage count 
                if ($fwd_tags{$chr}{$pos} >= $cage_highest{$gene}){
                    $cage_highest{$gene}=$fwd_tags{$chr}{$pos};
                    $cage_peak{$gene}=$pos;
                    $length_count{$gene}=$length_count;
                }
            }
        }

        if ($coord <= $coordFivePrime){ #stop at the annotated start codon
            $extended_coord=$coord-1; #last will exit the loop and also set the value of coord to the last position in the loop (zero)
            last;                     #therefore I need to store the value of the coord in a new variable to continue the search
        }
    }

    #then proceed upstream (from the coord if there was an annotated leader)
    for (1..$UPSTEAM_DISTANCE){

        $length_count++;

        #check for upstream gene;
        if (exists ($three_primes_fwd{$chr}{$pos})){
            for my $alt_gene (split(";",$three_primes_fwd{$chr}{$pos})){
                unless ($alt_gene eq $gene){  #not the stop codon of this gene
                    $leader_overlaps_upstream_gene{$gene}=1;
                    $gene_overlaps_downstream_leader{$alt_gene}=1;
                }
            }
        }
 
        #check for nc overlap
        if (exists ($nc_three_prime_most_position_fwd{$chr}{$pos})){ $gene_overlaps_nc{$gene}=1; }
        if (exists ($nc_five_prime_most_position_fwd{$chr}{$pos})){ $gene_overlaps_nc{$gene}=1; }

        #the coord gets set to zero when the previous loop exits?
        if (exists ($gene_model_fwd{$gene}{$extended_coord})){
            $pos=$gene_model_fwd{$gene}{$extended_coord};
       
            if (exists ($fwd_tags{$chr}{$pos})){
                if ($fwd_tags{$chr}{$pos} >= 1){ #minimum cage count 
                    if ($fwd_tags{$chr}{$pos} >= $cage_highest{$gene}){
                        $cage_highest{$gene}=$fwd_tags{$chr}{$pos};
                        $cage_peak{$gene}=$pos;
                        $length_count{$gene}=$length_count;
                    }
                }
            }
            $extended_coord--;
        }else{  #otherwise move into genommic coordinate space
            $pos--;      

            if (exists ($fwd_tags{$chr}{$pos})){
                if ($fwd_tags{$chr}{$pos} >= 1){ #minimum cage count 
                    if ($fwd_tags{$chr}{$pos} >= $cage_highest{$gene}){
                        $cage_highest{$gene}=$fwd_tags{$chr}{$pos};
                        $cage_peak{$gene}=$pos;
                        $length_count{$gene}=$length_count;
                    }
                }
            }
        }
    }
}

for my $gene (keys %gene_model_rev){
    my $coordFivePrime=$start_coord_rev{$gene};
    my $coordThreePrime=$stop_coord_rev{$gene};
    my $chr=$gene_2_chr{$gene};
    $cage_highest{$gene}=0; #initialise
    my $coord=0;
    my $pos=0;
    my $length_count=0;
    my $extended_coord;

    #search through the transcript up until the start codon
    for $coord (reverse (sort {$a <=> $b} keys %{ $gene_model_rev{$gene} } )){

        $pos=$gene_model_rev{$gene}{$coord};
        $length_count++;

        #check for upstream gene (overlapping with CDS)
        if (exists ($five_primes_rev{$chr}{$pos})){
            for my $alt_gene (split(";",$five_primes_rev{$chr}{$pos})){
                unless ($alt_gene eq $gene){  #not the stop codon of this gene
                    $gene_overlaps{$gene}=1;
                    $gene_overlaps{$alt_gene}=1;
                }
            }
        }

        #check for upstream gene (overlapping with CDS)
        if (exists ($three_primes_rev{$chr}{$pos})){
            for my $alt_gene (split(";",$three_primes_rev{$chr}{$pos})){
                unless ($alt_gene eq $gene){  #not the start codon of the current gene
                    $gene_overlaps{$gene}=1;
                    $gene_overlaps{$alt_gene}=1;
                }
            }
        }

        #check for nc overlap
        if (exists ($nc_three_prime_most_position_rev{$chr}{$pos})){ $gene_overlaps_nc{$gene}=1; }
        if (exists ($nc_five_prime_most_position_rev{$chr}{$pos})){ $gene_overlaps_nc{$gene}=1; }
        
        if (exists ($rev_tags{$chr}{$pos})){
            if ($rev_tags{$chr}{$pos} >= 1){ #minimum cage count 
                if ($rev_tags{$chr}{$pos} >= $cage_highest{$gene}){
                    $cage_highest{$gene}=$rev_tags{$chr}{$pos};
                    $cage_peak{$gene}=$pos;
                    $length_count{$gene}=$length_count;
                }
            }
        }

        if ($coord <= $coordFivePrime){ #stop at the annotaated start codon
            $extended_coord=$coord-1;
            last;
        }
    }

    #then proceed upstream (from the coord if there was an annotated leader
    for (1..$UPSTEAM_DISTANCE){

        $length_count++;

        #check for upstream gene;
        if (exists ($three_primes_rev{$chr}{$pos})){
            for my $alt_gene (split(";",$three_primes_rev{$chr}{$pos})){
                unless ($alt_gene eq $gene){  #not the start codon of the current gene
                    $leader_overlaps_upstream_gene{$gene}=1;
                    $gene_overlaps_downstream_leader{$alt_gene}=1;
                }
            }
        }                                    
 
        #check for nc overlap
        if (exists ($nc_three_prime_most_position_rev{$chr}{$pos})){ $gene_overlaps_nc{$gene}=1; }
        if (exists ($nc_five_prime_most_position_rev{$chr}{$pos})){ $gene_overlaps_nc{$gene}=1; }

        if (exists ($gene_model_rev{$gene}{$extended_coord})) {   
            $pos=$gene_model_rev{$gene}{$extended_coord};
            if (exists ($rev_tags{$chr}{$pos})){
                if ($rev_tags{$chr}{$pos} >= 1){ #minimum cage count 
                    if ($rev_tags{$chr}{$pos} >= $cage_highest{$gene}){
                        $cage_highest{$gene}=$rev_tags{$chr}{$pos};
                        $cage_peak{$gene}=$pos;
                        $length_count{$gene}=$length_count;
                    }
                }
            }
            $extended_coord--;
        }else{
            $pos++;  #reverse genes genomic coordinate will increase as we move upstream
            if (exists ($rev_tags{$chr}{$pos})){
                if ($rev_tags{$chr}{$pos} >= 1){ #minimum cage count 
                    if ($rev_tags{$chr}{$pos} >= $cage_highest{$gene}){
                        $cage_highest{$gene}=$rev_tags{$chr}{$pos};
                        $cage_peak{$gene}=$pos;
                        $length_count{$gene}=$length_count;
                    }
                }
            }
        }
    }
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#output

open (OUT1, '>', $out_leaders) || die "can't open $out_leaders\n";

print OUT1 "#gene_id,transcript_id,chr,dir,overlaps_inframe_gene,leader_overlaps_upstream,gene_overlaps_downstream_leader,gene_overlaps_nc_transcript,highest_cage_peak,count_at_highest_cage_peak,leader_length\n";

#add:
#%gene_overlaps_nc

for my $gene (keys %gene_model_fwd){
    my $coordFivePrime=$start_coord_fwd{$gene};
    my $coordThreePrime=$stop_coord_fwd{$gene};
    my $pos5=$gene_model_fwd{$gene}{$coordFivePrime};
    my $pos3=$gene_model_fwd{$gene}{$coordThreePrime};

    my $transcript=$most_expressed_transcript{$gene};
    my $chr=$gene_2_chr{$gene};
    my $cage_peak_value=$cage_highest{$gene};

    my $leader_length="NaN";
    if (exists($length_count{$gene})){
        $leader_length=$length_count{$gene}-($three_prime_most_coord_fwd{$gene}+2)+$start_coord_fwd{$gene};
        #length_count{$gene}        #count in coordinate space from last position in transctipt to highest peak
        #start_coord_fwd{$gene}     #coordinate to start codon
        #three_prime_most_coord_fwd{$gene} #coordinate of the last position in the transcript
    }   

    my $cage_peak_position="NaN";
    if (exists($cage_peak{$gene})){ $cage_peak_position=$cage_peak{$gene}; }

    my $overlaps_inframe_gene="FALSE";
    my $leader_overlaps_upstream="FALSE";
    my $gene_overlaps_downstream_leader="FALSE";
    my $gene_overlaps_nc="FALSE";

    if (exists ($gene_overlaps{$gene})){ $overlaps_inframe_gene="TRUE"; }
    if (exists ($leader_overlaps_upstream_gene{$gene})){ $leader_overlaps_upstream="TRUE"; }
    if (exists ($gene_overlaps_downstream_leader{$gene})){ $gene_overlaps_downstream_leader="TRUE"; }
    if (exists ($gene_overlaps_nc{$gene})){ $gene_overlaps_nc="TRUE"; }

    print OUT1 "$gene,$transcript,$chr,fwd,$overlaps_inframe_gene,$leader_overlaps_upstream,$gene_overlaps_downstream_leader,$gene_overlaps_nc,$cage_peak_position,$cage_peak_value,$leader_length\n";
}

for my $gene (keys %gene_model_rev){
    my $coordFivePrime=$start_coord_rev{$gene};
    my $coordThreePrime=$stop_coord_rev{$gene};
    my $pos5=$gene_model_rev{$gene}{$coordFivePrime};
    my $pos3=$gene_model_rev{$gene}{$coordThreePrime};
   
    my $transcript=$most_expressed_transcript{$gene};
    my $chr=$gene_2_chr{$gene};
    my $cage_peak_value=$cage_highest{$gene};

    my $leader_length="NaN";
    if (exists($length_count{$gene})){
        $leader_length=$length_count{$gene}-($three_prime_most_coord_rev{$gene}+2)+$start_coord_rev{$gene};
    }

    my $cage_peak_position="NaN";
    if (exists($cage_peak{$gene})){ $cage_peak_position=$cage_peak{$gene}; }

    my $overlaps_inframe_gene="FALSE";
    my $leader_overlaps_upstream="FALSE";
    my $gene_overlaps_downstream_leader="FALSE";
    my $gene_overlaps_nc="FALSE";

    if (exists ($gene_overlaps{$gene})){ $overlaps_inframe_gene="TRUE"; }
    if (exists ($leader_overlaps_upstream_gene{$gene})){ $leader_overlaps_upstream="TRUE"; }
    if (exists ($gene_overlaps_downstream_leader{$gene})){ $gene_overlaps_downstream_leader="TRUE"; }
    if (exists ($gene_overlaps_nc{$gene})){ $gene_overlaps_nc="TRUE"; }

    print OUT1 "$gene,$transcript,$chr,rev,$overlaps_inframe_gene,$leader_overlaps_upstream,$gene_overlaps_downstream_leader,$gene_overlaps_nc,$cage_peak_position,$cage_peak_value,$leader_length\n";
}
close(OUT1);

exit;
