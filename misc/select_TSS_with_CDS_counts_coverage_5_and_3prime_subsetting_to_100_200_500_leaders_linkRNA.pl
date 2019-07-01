#!/usr/bin/perl -w
use strict;

#30/01/19
#script to count count the SSU and RFP in fixed windows from the transcription start site

#plot 100 nt downstream of TSS for all genes with leaders greater than 100nt
#plot 200 nt downstream of TSS for all genes with leaders greater than 200nt
#plot 500 nt downstream of TSS for all genes with leaders greater than 500nt

#for lincRNAs

my $inGtf=$ARGV[0]; 
my $fasta=$ARGV[1];
my $bam_SSU=$ARGV[2];  
my $bam_LSU=$ARGV[3]; 
my $outdir=$ARGV[4]; 

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#second pass through the genome, find annotated start codons and setup transcript models for longest transcript of each gene
my %gene_exons_fwd;
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
        my ($biotype) = $b[8] =~ /gene_biotype\s"([^\"]+)";/;

        if ($gene_id && $transcript_id){

            if ($biotype){    #if the transcript is in the list of longest transcripts

                if ($biotype eq "lincRNA"){

                    $gene_2_chr{$gene_id}=$chr;

                    if ($dir eq "+"){ #fwd cases. Use start positions as 5'

                        if ($class eq "exon"){
                            $gene_exons_fwd{$gene_id}{$start}=$end;
                        }

                    }else{ #revese cases use end as 5'

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
my %three_prime_most_coord_fwd;
my %five_prime_most_coord_fwd;

for my $gene (keys %gene_exons_fwd){

     my $model_pos=0;
     $five_prime_most_coord_fwd{$gene}=$model_pos;  #initalise the 5' to the first coord coord

     for my $exon_start (sort {$a <=> $b} keys %{ $gene_exons_fwd{$gene} } ){
         my $exon_end=$gene_exons_fwd{$gene}{$exon_start};

         #fwd exons are in ascending order
         # start(-1)-> 100958 100975
         #             101077 101715 <-end(+1)

         for ($exon_start .. $exon_end){
             $gene_model_fwd{$gene}{$model_pos}=$_;
             $model_pos++;
         }
    }           
    $three_prime_most_coord_fwd{$gene}=$model_pos-1; #store the 3 prime most position of each gene
}

my %gene_model_rev;
my %three_prime_most_coord_rev;
my %five_prime_most_coord_rev;

for my $gene (keys %gene_exons_rev){

    my $model_pos=0;
    $five_prime_most_coord_rev{$gene}=$model_pos;  #initalise the 5' to the first coord coord
 
    for my $exon_end (reverse (sort {$a <=> $b} keys %{ $gene_exons_rev{$gene} } )){
        my $exon_start=$gene_exons_rev{$gene}{$exon_end};

        #rev exons are sorted in decending order  
        #           447087 447794 <-start(+1)
        # end(-1)-> 446060 446254

        while ($exon_start >= $exon_end){
            $gene_model_rev{$gene}{$model_pos}=$exon_start;
            $model_pos++;
            $exon_start--;
        }
    }
    $three_prime_most_coord_rev{$gene}=$model_pos-1; #store the 3 prime most position of each gene
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#Transcript coord key
#Transcript 5' coord   $five_prime_most_coord_fwd{$gene}  #1st nt of annotated transcript
#Transcript 3' coord   $three_prime_most_coord_???{$gene} #last nt of transcript
#$gene_model_fwd{$gene}{$coord}==genomic position

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#loop though genes and assign leader, CDS and trailer regions to hashes for quick searching.
my %whole_transcripts_fwd;
my %whole_transcripts_rev;

my %transcript_search_fwd; #key1=chr, key2=pos, value=gene
my %transcript_search_rev; #key1=chr, key2=pos, value=gene
my %transcript_counts_SSU;
my %transcript_counts_LSU;

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

my %overlaps_inframe_gene;

for my $gene (keys %gene_model_fwd){
    my $chr=$gene_2_chr{$gene};
    my $TSS_coord=$five_prime_most_coord_fwd{$gene};
    my $TTS_coord=$three_prime_most_coord_fwd{$gene};

    $transcript_counts_SSU{$gene}=0;
    $transcript_counts_LSU{$gene}=0;

    $transcript_length{$gene}=($TTS_coord-$TSS_coord)+1;
    $overlaps_inframe_gene{$gene}=0;

    my $relational_position=0;

    for my $coord (sort {$a <=> $b} keys %{ $gene_model_fwd{$gene} } ){
        if ($coord >= $five_prime_most_coord_fwd{$gene}){

            my $pos=$gene_model_fwd{$gene}{$coord};

            if (exists ($whole_transcripts_fwd{$chr}{$pos})){ #mark overlapping genes
                #print "overlapping lincRNA,$gene\n";
                my @gene = split("_", $whole_transcripts_fwd{$chr}{$pos});
                for (@gene){ $overlaps_inframe_gene{$_}=1; }
                $overlaps_inframe_gene{$gene}=1;
                $whole_transcripts_fwd{$chr}{$pos}.="_".$gene;
            }else{
                $whole_transcripts_fwd{$chr}{$pos}=$gene;
            }

            $transcript_search_fwd{$chr}{$pos}=$gene;

            if ($coord == $TSS_coord){   

                if ($TTS_coord-$TSS_coord > 100){ #check that there the transcript is at least 100nt long
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

                if ($TTS_coord-$TSS_coord > 200){ 
                    my $relational_position=0; 
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

                if ($TTS_coord-$TSS_coord > 500){ 
                    my $relational_position=0;                     while ($relational_position < 500){
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
    my $TSS_coord=$five_prime_most_coord_rev{$gene};    
    my $TTS_coord=$three_prime_most_coord_rev{$gene};

    $transcript_counts_SSU{$gene}=0;
    $transcript_counts_LSU{$gene}=0;
    $transcript_length{$gene}=($TTS_coord-$TSS_coord)+1;
    $overlaps_inframe_gene{$gene}=0;

    my $relational_position=0;

    for my $coord (sort {$a <=> $b} keys %{ $gene_model_rev{$gene} } ){
        if ($coord >= $five_prime_most_coord_rev{$gene}){

            my $pos=$gene_model_rev{$gene}{$coord};

            if (exists ($whole_transcripts_rev{$chr}{$pos})){ #mark overlapping genes
                #print "overlapping lincRNA,$gene\n";
                my @gene = split("_", $whole_transcripts_rev{$chr}{$pos});
                for (@gene){ $overlaps_inframe_gene{$_}=1; }
                $overlaps_inframe_gene{$gene}=1;
                $whole_transcripts_rev{$chr}{$pos}.="_".$gene;
            }else{
                $whole_transcripts_rev{$chr}{$pos}=$gene;
            }

            $transcript_search_rev{$chr}{$pos}=$gene;

            if ($coord == $TSS_coord){

                if ($TTS_coord-$TSS_coord > 100){ #check that the transcript is at least 100nt long
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

                if ($TTS_coord-$TSS_coord > 200){ 
                    my $relational_position=0; 
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

                if ($TTS_coord-$TSS_coord > 500){ 
                    my $relational_position=0; 
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
#counts

my $link=0;
for (keys %transcript_counts_SSU){
    $link++;
}

print "there are $link linkRNA's\n";


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
    my $transcript_hit=0;

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
            $total_SSU_bam_count++;
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
    my $transcript_hit=0;

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
            $total_LSU_bam_count++;
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
print OUT1 "gene_id,transcript_SSU_count,transcript_LSU_count,transcript_length,gene_overlaps_another_gene,fraction,position_in_region,count\n";

for my $gene (keys %TSS_100_SSU_coverage_by_positon){
 
    my $transcript_length=$transcript_length{$gene};
    my $gene_info="$gene,$transcript_counts_SSU{$gene},$transcript_counts_LSU{$gene},$transcript_length,$overlaps_inframe_gene{$gene}";

    for my $pos (sort {$a <=> $b} keys %{$TSS_100_SSU_coverage_by_positon{$gene}}){ 
        print OUT1 "$gene_info,SSU,$pos,$TSS_100_SSU_coverage_by_positon{$gene}{$pos}\n"; 
    }

    for my $pos (sort {$a <=> $b} keys %{$TSS_100_LSU_coverage_by_positon{$gene}}){
        print OUT1 "$gene_info,LSU,$pos,$TSS_100_LSU_coverage_by_positon{$gene}{$pos}\n"; 
    }
}
close(OUT1);


open (OUT2,">$out_file2")  || die "can't open $out_file2\n";
print OUT2 "gene_id,transcript_SSU_count,transcript_LSU_count,transcript_length,gene_overlaps_another_gene,fraction,position_in_region,count\n";

for my $gene (keys %TSS_200_SSU_coverage_by_positon){

    my $transcript_length=$transcript_length{$gene};
    my $gene_info="$gene,$transcript_counts_SSU{$gene},$transcript_counts_LSU{$gene},$transcript_length,$overlaps_inframe_gene{$gene}";

    for my $pos (sort {$a <=> $b} keys %{$TSS_200_SSU_coverage_by_positon{$gene}}){
        print OUT2 "$gene_info,SSU,$pos,$TSS_200_SSU_coverage_by_positon{$gene}{$pos}\n";
    }

    for my $pos (sort {$a <=> $b} keys %{$TSS_200_LSU_coverage_by_positon{$gene}}){
        print OUT2 "$gene_info,LSU,$pos,$TSS_200_LSU_coverage_by_positon{$gene}{$pos}\n";
    }
}
close(OUT2);


open (OUT3,">$out_file3")  || die "can't open $out_file3\n";
print OUT3 "gene_id,transcript_SSU_count,transcript_LSU_count,transcript_length,gene_overlaps_another_gene,fraction,position_in_region,count\n";

for my $gene (keys %TSS_500_SSU_coverage_by_positon){

    my $transcript_length=$transcript_length{$gene};
    my $gene_info="$gene,$transcript_counts_SSU{$gene},$transcript_counts_LSU{$gene},$transcript_length,$overlaps_inframe_gene{$gene}";

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
print OUT5 "gene_id,transcript_SSU_count,transcript_LSU_count,transcript_length,gene_overlaps_another_gene,fraction,position_in_region,count\n";

for my $gene (keys %TSS_100_SSU_3prime_counts_by_positon){

    my $transcript_length=$transcript_length{$gene};
    my $gene_info="$gene,$transcript_counts_SSU{$gene},$transcript_counts_LSU{$gene},$transcript_length,$overlaps_inframe_gene{$gene}";

    for my $pos (sort {$a <=> $b} keys %{$TSS_100_SSU_3prime_counts_by_positon{$gene}}){
        print OUT5 "$gene_info,SSU,$pos,$TSS_100_SSU_3prime_counts_by_positon{$gene}{$pos}\n";
    }

    for my $pos (sort {$a <=> $b} keys %{$TSS_100_LSU_3prime_counts_by_positon{$gene}}){
        print OUT5 "$gene_info,LSU,$pos,$TSS_100_LSU_3prime_counts_by_positon{$gene}{$pos}\n";
    }
}
close(OUT5);


open (OUT6,">$out_file6")  || die "can't open $out_file6\n";
print OUT6 "gene_id,transcript_SSU_count,transcript_LSU_count,transcript_length,gene_overlaps_another_gene,fraction,position_in_region,count\n";

for my $gene (keys %TSS_200_SSU_3prime_counts_by_positon){

    my $transcript_length=$transcript_length{$gene};
    my $gene_info="$gene,$transcript_counts_SSU{$gene},$transcript_counts_LSU{$gene},$transcript_length,$overlaps_inframe_gene{$gene}";

    for my $pos (sort {$a <=> $b} keys %{$TSS_200_SSU_3prime_counts_by_positon{$gene}}){
        print OUT6 "$gene_info,SSU,$pos,$TSS_200_SSU_3prime_counts_by_positon{$gene}{$pos}\n";
    }

    for my $pos (sort {$a <=> $b} keys %{$TSS_200_LSU_3prime_counts_by_positon{$gene}}){
        print OUT6 "$gene_info,LSU,$pos,$TSS_200_LSU_3prime_counts_by_positon{$gene}{$pos}\n";
    }
}
close(OUT6);


open (OUT7,">$out_file7")  || die "can't open $out_file7\n";
print OUT7 "gene_id,transcript_SSU_count,transcript_LSU_count,transcript_length,gene_overlaps_another_gene,fraction,position_in_region,count\n";

for my $gene (keys %TSS_500_SSU_3prime_counts_by_positon){

    my $transcript_length=$transcript_length{$gene};
    my $gene_info="$gene,$transcript_counts_SSU{$gene},$transcript_counts_LSU{$gene},$transcript_length,$overlaps_inframe_gene{$gene}";

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
print OUT9 "gene_id,transcript_SSU_count,transcript_LSU_count,transcript_length,gene_overlaps_another_gene,fraction,position_in_region,count\n";

for my $gene (keys %TSS_100_SSU_5prime_counts_by_positon){

    my $transcript_length=$transcript_length{$gene};
    my $gene_info="$gene,$transcript_counts_SSU{$gene},$transcript_counts_LSU{$gene},$transcript_length,$overlaps_inframe_gene{$gene}";

    for my $pos (sort {$a <=> $b} keys %{$TSS_100_SSU_5prime_counts_by_positon{$gene}}){
        print OUT9 "$gene_info,SSU,$pos,$TSS_100_SSU_5prime_counts_by_positon{$gene}{$pos}\n";
    }

    for my $pos (sort {$a <=> $b} keys %{$TSS_100_LSU_5prime_counts_by_positon{$gene}}){
        print OUT9 "$gene_info,LSU,$pos,$TSS_100_LSU_5prime_counts_by_positon{$gene}{$pos}\n";
    }
}
close(OUT9);


open (OUT10,">$out_file10")  || die "can't open $out_file10\n";
print OUT10 "gene_id,transcript_SSU_count,transcript_LSU_count,transcript_length,gene_overlaps_another_gene,fraction,position_in_region,count\n";

for my $gene (keys %TSS_200_SSU_5prime_counts_by_positon){

    my $transcript_length=$transcript_length{$gene};
    my $gene_info="$gene,$transcript_counts_SSU{$gene},$transcript_counts_LSU{$gene},$transcript_length,$overlaps_inframe_gene{$gene}";

    for my $pos (sort {$a <=> $b} keys %{$TSS_200_SSU_5prime_counts_by_positon{$gene}}){
        print OUT10 "$gene_info,SSU,$pos,$TSS_200_SSU_5prime_counts_by_positon{$gene}{$pos}\n";
    }

    for my $pos (sort {$a <=> $b} keys %{$TSS_200_LSU_5prime_counts_by_positon{$gene}}){
        print OUT10 "$gene_info,LSU,$pos,$TSS_200_LSU_5prime_counts_by_positon{$gene}{$pos}\n";
    }
}
close(OUT10);


open (OUT11,">$out_file11")  || die "can't open $out_file11\n";
print OUT11 "gene_id,transcript_SSU_count,transcript_LSU_count,transcript_length,gene_overlaps_another_gene,fraction,position_in_region,count\n";

for my $gene (keys %TSS_500_SSU_5prime_counts_by_positon){

    my $transcript_length=$transcript_length{$gene};
    my $gene_info="$gene,$transcript_counts_SSU{$gene},$transcript_counts_LSU{$gene},$transcript_length,$overlaps_inframe_gene{$gene}";

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
