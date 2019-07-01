#!/usr/bin/perl -w
use strict;

#29/04/18
##script to parse am aligned bam file that has been poly A trimmed and to add in sequences that were cropped by which match the genome reference (incorrectly trimmed), also remove 3' A's that do not match the reference

#I want to only look read overlapping transcript?
#    so that I can account for reads mapping introns

my $inGtf=$ARGV[0]; 
my $fasta=$ARGV[1];
my $leaders=$ARGV[2]; #from cageR (cage_peaks_assignment_gtf_transcript_coords_for_TCPseq_v8.pl)
my $fastq=$ARGV[3];
my $bam_TCPseq=$ARGV[4];
my $most_highly_expressed=$ARGV[5];

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
        if ($chr =~ /^chr(.*)/){ $chr=$1; } #if the chr name have a chr* prefix, remove it 
      
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
my %five_prime_most_coord_fwd;
my %three_prime_most_coord_fwd;

for my $gene (keys %gene_exons_fwd){
    if ( (exists ($gene_start_codon_fwd{$gene})) && (exists ($gene_stop_codon_fwd{$gene})) ) { #restrict to genes with annotated start + stop codon
        my $model_pos=0;
        for my $exon_start (sort {$a <=> $b} keys %{ $gene_exons_fwd{$gene} } ){
            my $exon_end=$gene_exons_fwd{$gene}{$exon_start};
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

my %gene_model_rev;
my %start_coord_rev;
my %stop_coord_rev;
my %five_prime_most_coord_rev;
my %three_prime_most_coord_rev;

for my $gene (keys %gene_exons_rev){
    if ( (exists ($gene_start_codon_rev{$gene})) && (exists ($gene_stop_codon_rev{$gene})) ) { #restrict to genes with annotated start + stop codon
        my $model_pos=0;
        for my $exon_end (reverse (sort {$a <=> $b} keys %{ $gene_exons_rev{$gene} } )){
            my $exon_start=$gene_exons_rev{$gene}{$exon_end};
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

        if ($chr =~ /^chr(.*)/){ $chr=$1; } #if the chr name have a chr* prefix, remove it 
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
my %leader_start_coord; #key=gene, value=coord

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
                    $leader_start_coord{$gene}=$coord;
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
            $leader_start_coord{$gene}=$extended_coord;
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
                    $leader_start_coord{$gene}=$coord;
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
            $leader_start_coord{$gene}=$extended_coord;
            $five_prime_most_coord_rev{$gene}=$extended_coord;
        }
    }
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#Transcript coord key
#Transcript 5' coord   0                                  #1st nt of annotated transcript
#Transcript 3' coord   $three_prime_most_coord_???{$gene} #last nt of transcript
#Start codon coord     $start_coord_???{$gene}            #1st nt in start codon
#Stop codon coord      $stop_coord_???{$gene}             #1st nt in stop codon
#Leader start coord    $leader_start_coord{$gene}         #1st nt of cage defined leader

#$gene_model_fwd{$gene}{$coord}==genomic position

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#loop though genes and assign leader, CDS and trailer regions to hashes for quick searching. (added start + stop codons).
my %transcript_search_fwd;
my %transcript_search_rev;

for my $gene (keys %gene_model_fwd){
    my $chr=$gene_2_chr{$gene};
    my $start_coord=$start_coord_fwd{$gene};
    my $stop_coord=$stop_coord_fwd{$gene};
    my $transcript_start=$five_prime_most_coord_fwd{$gene};
    my $transcript_end=$three_prime_most_coord_fwd{$gene};
    unless (exists ($overlaps_inframe_gene{$gene}) || exists ( $leader_overlaps_upstream{$gene} ) || exists ( $gene_overlaps_downstream_leader{$gene} ) || exists ( $gene_overlaps_nc{$gene}) ) {
        for my $coord (sort {$a <=> $b} keys %{ $gene_model_fwd{$gene} } ){
            my $pos=$gene_model_fwd{$gene}{$coord};
            if ($coord >= $five_prime_most_coord_fwd{$gene}){ #check if we are passed the cage assigned 5' most position
#               $transcript_search_fwd{$chr}{$pos}=$gene;
                $transcript_search_fwd{$chr}{$pos}=$gene.";"."$coord";
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
    unless (exists ($overlaps_inframe_gene{$gene}) || exists ( $leader_overlaps_upstream{$gene} ) || exists ( $gene_overlaps_downstream_leader{$gene} ) || exists ( $gene_overlaps_nc{$gene}) ) {
        for my $coord (sort {$a <=> $b} keys %{ $gene_model_rev{$gene} } ){
            my $pos=$gene_model_rev{$gene}{$coord};
            if ($coord >= $five_prime_most_coord_rev{$gene}){ #check if we are passed the cage assigned 5' most position
#               $transcript_search_rev{$chr}{$pos}=$gene;
                $transcript_search_rev{$chr}{$pos}=$gene.";"."$coord";
            }
        }
    }
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#open the bam file and store the read identifiers to reduce the number of fastq records that need to be saved
my %read_id;
my $aligned_count=0;

open BAM1, "samtools view $bam_TCPseq |";
while(<BAM1>){
    next if(/^(\@)/);        # skipping the header lines (if you used -h in the samools command)
    s/\n//;  s/\r//;         # removing new line
    my @sam = split(/\t+/);  # splitting SAM line into array
    unless ($sam[1] & 0x4){  #if the read is aligned
        $read_id{$sam[0]}=1; #store the read id
        $aligned_count++;
    }
}
close (BAM1);

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#of the fastq file and store the sequences and qualities
my %fastq_sequence;
my %fastq_qualities;
my $total_fastq=0;
my $matched_fastq=0;

my $header;
my $seq;

open(FASTQ, "gunzip -c $fastq |") || die "can’t open pipe to $fastq";
while (<FASTQ>) {
    chomp();


    #for zebrafish I need to add the option to trim 3nt?
    #or just write somehting fancy to get the match coord?
   
    if ($.%4==1){        #header1
 #       print "header=$_\n";
        if (/^\@([^\s]+)\s/){
            $header=$1;
            $total_fastq++;
        }
    }elsif($.%4==2){     #sequence
#        print "sequence=$_\n";
        $seq=$_;
    }elsif($.%4==3){     #header2 
    }else{               #quality scores
#        print "$header,$seq,$_\n";
        if (exists ($read_id{$header})){
             $fastq_sequence{$header}=$seq;
             $fastq_qualities{$header}=$_;
             $matched_fastq++;
        }
    }
}
close(FASTQ);

#print "Aligned_reads: $aligned_count\n";
#print "Fastq_reads: $total_fastq\n";
#print "matched_reads: $matched_fastq\n";

#exit;

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

#fastq example
#@SRR3458592.1 D9S4KXP1:312:HFK2FADXX:1:1101:1114:2142 length=151
#TATCTGGTTGATCCTGCCAGTAGTCATATGCTTGTCTCAAAGATTAAGCCATGCATGTCTAAGTATAAGCAATAAAAAAAAAAAAAAAAAAATGGAATTCTCGGGTGCCAAGGAACTCCAGTCACATTCCTATCTCGTATGCCGTCTTCTG
#+SRR3458592.1 D9S4KXP1:312:HFK2FADXX:1:1101:1114:2142 length=151
#CCCFFFFDFHHHHJJJIGIIIJJJJJJJIJJIJJHJJJIJJIIJJGIJJJJIJJIIIIIIJIIEHIJIIJIJIIJJJIHFDDDDDDDDDDDD?>@@CCC>>:?BD89BC>@??BCD>C@3>:>::3@>@>>:::3:<<<?CD@555<@:@(

#bam example
#SRR3458591.35877335	16	II	395979	50	26M	*	0	0	CGGCCTTGGCGGCGTCGACAATTTTT	IIJGGDFJJJJIJHHHHHFFDFFCC@	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:26	YT:Z:UU	XS:A:-	NH:i:1

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#assign read to genomic regions/transcripts
#count reads in transcripts
#calculate longest fragment length

open BAM2, "samtools view -h $bam_TCPseq |";
while(<BAM2>){

    #next if(/^(\@)/);  ## skipping the header lines (if you used -h in the samools command)
    if(/^(\@)/){
          print $_; #headers
    }else{

        #s/\n//;  s/\r//;  ## removing new line
        my @sam = split(/\t+/);  ## splitting SAM line into array

        my $leftMost=$sam[3]; #leftmost position of match 5' for fwd, 3' for rev
        my $readName=$sam[0];
        my $flag=$sam[1];
        my $chr=$sam[2];
        my $mapq=$sam[4];
        my $cigar=$sam[5];
        my $seq=$sam[9];
        my $qal=$sam[10];
        my $fivePrime;
        my $threePrime;

        if ($chr=~/chr(.*)/){ #ensembl chromosome names do not contain the "chr" prefix
            $chr=$1;
        }

        unless ($flag & 0x4){   #if aligned

            if ($mapq >= 10){     #mapping uniqnes filter

                my $length=length($seq);

                if ($flag & 0x10){  #Reverse reads. PolyA adaptors are expected in the leftmost positions encoded as "T"

                    my $length=length($seq);
                    while ($cigar =~/(\d+)I/g){   #add to length for insertions
                        $length+=$1;
                    }
                    while ($cigar =~/(\d+)D/g){   #substact from length for deletions
                       $length-=$1;
                    }
                    my $rightMost=$leftMost+($length-1);    
  
                    if ($cigar =~ /^(\d+)M/){ #if the cigar starts with a match
 
                        if (exists ($fasta_sequences{$chr})){

                            my $bases_removed=0;
                            my $bases_added=0;
                            my $genomic_seq=substr($fasta_sequences{$chr}, ($leftMost-1), $1);
 
                            #compare the sequences and remove 5' T's that do not match the sequence 
                            while ((substr ($seq, 0, 1)) eq "T"){
                                if ((substr($genomic_seq, 0, 1)) ne "T"){
                                    $seq=substr($seq, 1);
                                    $qal=substr($qal, 1);
                                    $bases_removed++;
                                    if (length $genomic_seq > 1){ 
                                        $genomic_seq=substr($genomic_seq, 1);
                                    }else{
                                        $genomic_seq="N";
                                    }
                                }else{
                                    last; #no need to trim if the 5' end matches the genomic position
                                }
                            }

                            if ($bases_removed){ #print out the updated alignments
      
                                my $left_over=""; #update the cigar string
                                while ($bases_removed >= 1){
                                   if ($cigar =~ /(^([0-9]+)([MIDN]))/){ #go though each part of the cigar, starting from the rightmost segment
                                        my $cigar_part = $1;
                                        unless ($3 eq "N" || $3 eq "D" ){ #if this cigar segment is skipped from the reference, just remove the whole segment (Insertions however need to be counted)
                                            if ($2 > $bases_removed){ #if this part of the cigar pattern is longer than the number of bases we removed, save the extras
                                                 $left_over=$2-$bases_removed."".$3;
                                                 $leftMost+=$bases_removed;  #update the lestMost position. Add the remaning length to the leftmost coord
                                            }else{
                                                $leftMost+=$2;  #otherwise add the length of the cigar segment that we just process to the leftmost coord
                                            }
                                            $bases_removed-=$2;        #update the count
                                        }
                                        $cigar =~ s/^$cigar_part//; #remove the part we just processed (anchor substitution to the start of the cigar string)
                                    }
                                }
                                $cigar=$left_over.$cigar; #append the saved pattern (if any)  
                                print "$sam[0]\t$sam[1]\t$sam[2]\t$leftMost\t$sam[4]\t$cigar\t$sam[6]\t$sam[7]\t$sam[8]\t$seq\t$qal\n"; #I don't include the optional alignment fields

                            }else{ #see if we can add some bases to the 3'

                                if (exists ($fastq_sequence{$readName})){
                                     my $fastq_seq=reverse($fastq_sequence{$readName}); #get the fastq sequence
                                     my $qual_seq=reverse($fastq_qualities{$readName});
                                     $fastq_seq=~tr/ACGTacgt/TGCAtgca/;
                                     $qual_seq=~tr/ACGTacgt/TGCAtgca/;

                                    if ($fastq_seq =~ /$seq/ ){ #find the end coordindate of the sequence in the bam file
                                        my $seq_to_check=substr($fastq_seq, 0, $-[0]);  #the leftmost portion of the fastq (3' or the match)
                                        my $qual_to_check=substr($qual_seq, 0, $-[0]);
                                        my $match_pos=$leftMost-1;

                                        if (exists ($transcript_search_rev{$chr}{$match_pos})){
                                            my ($gene, $coord) = split(";", $transcript_search_rev{$chr}{$match_pos});
 
                                           while ($seq_to_check){ #check the transcript sequence upstream of the read one by one until a missmatch is found
                                               if (exists ($gene_model_rev{$gene}{$coord})){
                                                   my $base_fastq=substr($seq_to_check, -1);   #take last character
                                                   my $base_qual=substr($qual_to_check, -1);
                                                   my $base_genome=substr($fasta_sequences{$chr},($gene_model_rev{$gene}{$coord}-1),1);

                                                   if ($base_fastq eq $base_genome){
                                                       $leftMost=$gene_model_rev{$gene}{$coord}; #update the new leftmost coordinate (bases on the extension
                                                       $seq=$base_fastq.$seq; #add to the leftmost position 
                                                       $qal=$base_qual.$qal;  #add to the leftmost position
                                                       $bases_added++;
                                                       $coord++;
                                                       chop($seq_to_check); # chop the last charater 
                                                       chop($qual_to_check); # chop the last charater
                                                   }else{
                                                       last;
                                                   }
                                               }else{
                                                    last;
                                               }
                                            }
                                        }
                                    }
                                }
                            }
                    
                            if ($bases_added){  #print out the updated alignments. Add the extended bases at the leftmost end or the reverse read
                                 if ($cigar=~/(^\d+)M(\d+.*$)/){
                                      $cigar=($1+$bases_added)."M".($2);
                                 }elsif($cigar=~/^(\d+)M$/){
                                     $cigar=($1+$bases_added)."M";
                                 }
 
                                 print "$sam[0]\t$sam[1]\t$sam[2]\t$leftMost\t$sam[4]\t$cigar\t$sam[6]\t$sam[7]\t$sam[8]\t$seq\t$qal\n"; #I don't include the optional alignment fields 

                            }else{
                                unless ($bases_removed){ #if we have neither added nor removed bases, print the original sam line
                                   print $_;
                                }
                            }
                        }
                    }else{ #print the unmodified sam line
                        print $_;     
                    }

                }else{ #fwd cases

                    my $length=length($seq);
                    while ($cigar =~/(\d+)I/g){   #add to length for insertions
                        $length+=$1;
                    }
                    while ($cigar =~/(\d+)D/g){   #subtract from length for deletions
                       $length-=$1;
                    }
                    my $rightMost=$leftMost+($length-1);
        
                    if ($cigar =~ /(\d+)M$/){ #if the cigar finishes with a match
                        if (exists ($fasta_sequences{$chr})){
                            my $bases_removed=0;
                            my $bases_added=0;
                            my $genomic_seq=substr($fasta_sequences{$chr}, ($rightMost-$1), $1);
                  
                            #compare the sequences and remove 3' A's that do not match the sequence 
                            while ((substr ($seq, -1)) eq "A"){
                                if ((substr($genomic_seq, -1)) ne "A"){
                                    $seq=substr($seq, 0, -1);
                                    $qal=substr($qal, 0, -1);
                                    $genomic_seq=substr($genomic_seq, 0, -1);
                                    $bases_removed++;
                                }else{
                                    last; #no need to trim if the 3' end matches the genomic position
                                }
                            }
        
                            if($bases_removed){  #print out the updated sequence
        
                                my $left_over=""; #update the cigar string
                                while ($bases_removed >= 1){
                                   if ($cigar =~ /(([0-9]+)([MIDN])$)/){ #go though each part of the cigar, starting from the rightmost segment
                                        my $cigar_part = $1;
                                        unless ($3 eq "N" || $3 eq "D" ){ #if this cigar segment is skipped from the reference, just remove the whole segment (Insertions however need to be counted)
                                            if ($2 > $bases_removed){ #if this part of the cigar pattern is longer than the number of bases we removed, save the extras
                                                $left_over=$2-$bases_removed."".$3;
                                            }
                                            $bases_removed-=$2;        #update the count
                                        }
                                        $cigar =~ s/$cigar_part$//; #remove the part we just processed (anchor substitution to the end of the cigar string)
                                    }
                                }
                                $cigar.=$left_over; #append the saved pattern (if any)
                                print "$sam[0]\t$sam[1]\t$sam[2]\t$sam[3]\t$sam[4]\t$cigar\t$sam[6]\t$sam[7]\t$sam[8]\t$seq\t$qal\n"; #I don't include the optional alignment fields
        
                            }else{  #seee if we can add in some bases
         
                                if (exists ($fastq_sequence{$readName})){
                                    my $fastq_seq=$fastq_sequence{$readName}; #get the fastq sequence
                                    my $qual_seq=$fastq_qualities{$readName};
                                    if ($fastq_seq =~ /$seq/ ){ #find the end coordindate of the sequence in the bam file
                                        my $seq_to_check=substr($fastq_seq, $+[0], -1);
                                        my $qual_to_check=substr($qual_seq, $+[0], -1);
                                        my $match_pos=$leftMost+$+[0]; #could I not also use rightmnost + here?
                                                                        
                                        if (exists ($transcript_search_fwd{$chr}{$match_pos})){
                                            my ($gene, $coord) = split(";", $transcript_search_fwd{$chr}{$match_pos});
                                            while ($seq_to_check){ #check the transcript sequence upstream of the read one by one until a missmatch is found
                                                if (exists ($gene_model_fwd{$gene}{$coord})){
                                                    my $base_fastq=substr($seq_to_check, 0, 1);
                                                    my $base_qual=substr($qual_to_check, 0, 1);
                                                    my $base_genome=substr($fasta_sequences{$chr},($gene_model_fwd{$gene}{$coord}-1), 1);
                                                    if ($base_fastq eq $base_genome){
                                                        $seq.=$base_fastq;
                                                        $qal.=$base_qual;
                                                        $bases_added++;
                                                        $coord++;
                                                        $seq_to_check = substr($seq_to_check, 1); # chop the first charater
                                                        $qual_to_check = substr($qual_to_check, 1); # chop the first charater 
                                                    }else{
                                                        last;
                                                    }
                                                }else{
                                                    last;
                                                }
                                            }                                
                                        }
                                    }
                                } 
                            }
                            if($bases_added){  #print out the updated alignments, add the extended bases to the righmost position of the fwd read cigar
                                if ($cigar=~/(^.*[IDMN])(\d+)M$/){
                                    $cigar=$1."".($2+$bases_added)."M";
                                }elsif($cigar=~/^(\d+)M$/){
                                    $cigar=($1+$bases_added)."M";
                                }
                                print "$sam[0]\t$sam[1]\t$sam[2]\t$sam[3]\t$sam[4]\t$cigar\t$sam[6]\t$sam[7]\t$sam[8]\t$seq\t$qal\n"; #I don't include the optional alignment fields 
                            }else{
                                unless ($bases_removed){ #if we have neither added nor removed bases, print the original sam line
                                    print $_;
                                }
                            }
                        }
                    }else{ #print the unmodified sam line
                       print $_;
                    }
                }
            }
        }
    }
}
close(BAM2);

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

#print "There were $total_trimmed bases trimmed in total\n";
#print "There were $total_added bases trimmed in total\n"; 

#print "There were $total_reads_trimmed reads trimmed\n";
#print "There were $total_reads_extended reads extended\n"; 

exit;

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
