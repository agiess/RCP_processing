#!/usr/bin/perl -w
use strict;

#28/02/19
#script to assign FMPK counts to leader (CAGE updated), start_codon, CDS, stop_codon, 3'trailer
#calculates kozak scores for zebrafish and yeast
#mark histone genes based on biomart GO terms
#exclude the 3' most end of each transcript where polyA enrichment 

#identify overlapping transcripts

#updated to count coverage in:
# TSS proximal regions (+70 to +120)
# TIS proximal regions (-100 to -50)

my $inGtf=$ARGV[0]; 
my $fasta=$ARGV[1];
my $bam_RIBOseq=$ARGV[2];
my $bam_RNAseq=$ARGV[3]; #polyA selected RNA
my $bam_totalRNA=$ARGV[4];
my $bam_TCPseq_SSU=$ARGV[5];
my $bam_TCPseq_LSU=$ARGV[6];
my $leaders=$ARGV[7]; #from cageR (cage_peaks_assignment_gtf_transcript_coords_for_TCPseq_v8.pl)
my $go=$ARGV[8]; #biomart genes with GO term name and WikiGene description
my $most_highly_expressed=$ARGV[9];

my $organism=$ARGV[10];  # if ($organism eq "yeast")
my $use_cage=$ARGV[11];  # unless($use_cage eq "exclude")

#excluding:
#genes where the TSS is downstream of the start codon
#genes without a detectable cage peak 
#genes that are annotated as protien_coding

#reads are assigned to leaders/CDS on a stand specific basis (including the RNA-seq) as follows:
#Ribo-seq:  Using a shoelaces shifted bed, exclude reads that overlap TIS and stop codons 
#RNA-seq:   Assign fragments that overlap leader and CDS to both.
#TCP-seq_SSU:   Exclude fragments that overlap the TIS or stop codon.
#TCP-seq_LSU:   Exclude fragments that overlap the TIS or stop codon.

#Reads are normalised by the total number of reads that were found to map to leader, CDS and trailer regions of genes as defined above.

#Flags:
#overlapping_gene: The longest transcript of this gene overlaps with the longest transcript of another gene
#leader_potentially_overlaps_upstream_gene: There is an upstream gene within 500nt of the start codon of this gene 
#gene_potentially_overlaps_downstream_leader: There is an downstream gene start codon withing 500nt of the 3' most position of this gene (in yeast this is the stop codon).

#TCT motif - #The TCT motif, a key component of an RNA polymerase II transcription system for the translational machinery
#YC+1TTTY central 6 nt of the TCT motif (YCTTTY) TCT motif core   (Y=T/C)

my $MAX_READ_LENGTH=76; #Zv10

#my $CERT_motif=  [CG]\w[CG]\w[CG]\wC\wCCGC\w\w[CG]C;
#my $polyU = TTTTTTT;

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
#upstream kozak

my %raw_kozak_upstream =
   (
       "A" => [ 35, 62, 39, 28 ],
       "C" => [ 32, 5, 23, 36 ],
       "G" => [ 19, 28, 17, 27 ],
       "T" => [ 14, 5, 21, 10 ],
   );

my %PWM_upstream; #key1: position, key2: base, value: weight
for my $pos (0 .. 3){

    #0 == -4 #nucleotide position in relation to start codon
    #1 == -3
    #2 == -2
    #3 == -1

    my $pos_sum=0;
    my $pwm_sum=0;
    for my $base (keys %raw_kozak_upstream){ #sum the nucleotide frequencies per position
        $pos_sum+=$raw_kozak_upstream{$base}[$pos];
    }

    for my $base(keys %raw_kozak_upstream){ #score the PWM
        my $psudo_count= sqrt($pos_sum);
        my $background_probability=0.25; #no base preference
        my $pwm=&log2( ($raw_kozak_upstream{$base}[$pos] + $psudo_count * $background_probability) / ($pos_sum + $psudo_count * $background_probability));
        $PWM_upstream{$pos}{$base}=$pwm;
        $pwm_sum+=$pwm;
    }

    $PWM_upstream{$pos}{"N"}=($pwm_sum/4); #set "N" to be equal to the column mean. For genes with short leaders, missing upstream positions 
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

#Yeast 2011	The mRNA landscape at yeast translation initiation sites				
#my %yeast_matrix =
#   (            #-4   -3    -2    -1     4     5
#       "A" => [ 0.53, 1.02, 0.47, 0.58, -0.11, -0.6 ],
#       "C" => [ 0.16, -1.16, 0.16, -0.12, -0.47, 1.26 ],
#       "G" => [ -0.3, 0.2, -0.52, -0.15, 0.77, -0.11 ],
#       "T" => [ -0.74, -2.1, -0.49, -0.71, -0.18, -0.75 ],
#   );

my %PWM_yeast =
(
    0 => { 
        "A" => 0.53,
        "C" => 0.16,
        "G" => -0.3,
        "T" => -0.74
    },
    1 => {
        "A" => 1.02,
        "C" => -1.16,
        "G" => 0.2,
        "T" => -2.1
    },
    2 => { 
        "A" => 0.47,
        "C" => 0.16,
        "G" => -0.52,
        "T" => -0.49
    },
    3 => { 
        "A" => 0.58,
        "C" => -0.12,
        "G" => -0.15,
        "T" => -0.71
    },
    4 => {
        "A" => -0.11,
        "C" => -0.47,
        "G" => 0.77,
        "T" => -0.18
    },
    5 => {
        "A" => -0.6,
        "C" => 1.26,
        "G" => -0.11,
        "T" => -0.75
    }
);

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#TISU
#
#GAAGATGGCGGC
#CAAGATGGCGGC
#GAACATGGCGGC
#CAACATGGCGGC

my %raw_TISU =
   (
       "A" => [ 0, 100, 100, 0, 100, 0, 0, 0, 0, 0, 0, 0 ],
       "C" => [ 50, 0, 0, 50, 0, 0, 0, 0, 100, 0, 0, 100 ],
       "G" => [ 50, 0, 0, 50, 0, 0, 100, 100, 0, 100, 100, 0 ],
       "T" => [ 0, 0, 0,  0, 0, 100, 0, 0, 0, 0, 0, 0 ],
   );

my %PWM_TISU; #key1: position, key2: base, value: weight
for my $pos (0 .. 11){

    my $pos_sum=0;
    my $pwm_sum=0;
    for my $base (keys %raw_TISU){ #sum the nucleotide frequencies per position
        $pos_sum+=$raw_TISU{$base}[$pos];
    }

    for my $base(keys %raw_TISU){ #score the PWM
        my $psudo_count= sqrt($pos_sum);
        my $background_probability=0.25; #no base preference
        my $pwm=&log2( ($raw_TISU{$base}[$pos] + $psudo_count * $background_probability) / ($pos_sum + $psudo_count * $background_probability));
        $PWM_TISU{$pos}{$base}=$pwm;
        $pwm_sum+=$pwm;
    }

    $PWM_TISU{$pos}{"N"}=($pwm_sum/4); #set "N" to be equal to the column mean. For genes with short leaders, missing upstream positions 
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#open gtf and get transcript lengths function
#my %transcripts; #key = gene_id, transcript_id, #value = sum_exon_lengths;
#get_transcript_lengths($inGtf, \%transcripts);

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
#Store histone genes
my %histone; #key=gene_id; value = 1

open (GO, $go) || die "can't open $go";
while (<GO>){
   chomp();
   my @lines=split(",");

    my $gene=$lines[0];
    my $go_terms=$lines[3];
    my $wiki_desc=$lines[13];

    if ($gene){
        if ($go_terms){
            if ($go_terms =~ /^nucleosome\sassembly$/i){
                $histone{$gene}=1;
            }elsif ($go_terms =~ /^nucleosome$/i){
                $histone{$gene}=1;    
            }
        }
        if ($wiki_desc){
            if ($wiki_desc =~ /^histone\sfamily$/i){
                $histone{$gene}=1;
            }
        }
    }
}
close(GO);

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#setup trascript models

my %gene_model_fwd;
my %start_coord_fwd;
my %stop_coord_fwd;

my %gene_exon_exon_in_leader;

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

        $five_prime_most_coord_fwd{$gene}=0;
        my $model_pos=0;
        my $exon_count=0; 

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
                    $gene_exon_exon_in_leader{$gene}=$exon_count;
                }
                $model_pos++;
            }
            $exon_count++;
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

        $five_prime_most_coord_rev{$gene}=0;
        my $model_pos=0;
        my $exon_count=0;

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
                    $gene_exon_exon_in_leader{$gene}=$exon_count;
                }
                $model_pos++;
                $exon_start--;
            }
            $exon_count++;

        }
        $three_prime_most_coord_rev{$gene}=$model_pos-1; #store the 3 prime most position of each gene
    }
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#parse leaders
my %leader_positions_fwd;
my %leader_positions_rev;

#filters:
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

        if ($overlaps_inframe_gene eq "TRUE"){ $overlaps_inframe_gene{$gene}=1; }
        if ($leader_overlaps_upstream eq "TRUE"){ $leader_overlaps_upstream{$gene}=1; }
        if ($gene_overlaps_downstream_leader eq "TRUE"){ $gene_overlaps_downstream_leader{$gene}=1; }
        if ($gene_overlaps_nc eq "TRUE"){ $gene_overlaps_nc{$gene}=1; }

        unless ($leader_length eq "NaN"){  #only take genes that have a detectable cage peak

            unless ($leader_length < 0){  #exlude genes with negative leader sizes, as they cuase problems with FPKM

                if ($dir eq "fwd"){ 
                    if (exists ($start_coord_fwd{$gene})){
                        unless ($highest_cage_peak >=  $gene_model_fwd{$gene}{$start_coord_fwd{$gene}}){  #exclude genes where the TSS is downstream of the start codon
                            $leader_positions_fwd{$gene}=$highest_cage_peak;
                            $cage_peak_value{$gene}=$count_at_highest_cage_peak;
                        } 
                    }
                }else{
                    if (exists ($start_coord_rev{$gene})){
                        unless ($highest_cage_peak <=  $gene_model_rev{$gene}{$start_coord_rev{$gene}}){  #exclude genes where the TSS is downstream of the start codon
                            $leader_positions_rev{$gene}=$highest_cage_peak;
                            $cage_peak_value{$gene}=$count_at_highest_cage_peak;
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
my %whole_transcripts_fwd;
my %whole_transcripts_rev;

my %leader_search_fwd; #key1=chr, key2=pos, value=gene
my %leader_search_rev; #key1=chr, key2=pos, value=gene
my %CDS_search_fwd; #key1=chr, key2=pos, value=gene
my %CDS_search_rev; #key1=chr, key2=pos, value=gene
my %trailer_search_fwd; #key1=chr, key2=pos, value=gene
my %trailer_search_rev; #key1=chr, key2=pos, value=gene

my %start_codons_search_fwd; #key1=chr, key2=pos, value=gene
my %start_codons_search_rev; #key1=chr, key2=pos, value=gene
my %stop_codons_search_fwd; #key1=chr, key2=pos, value=gene
my %stop_codons_search_rev; #key1=chr, key2=pos, value=gene

my %TSS_beginning_fwd; #key1=chr, key2=pos, value=gene
my %TSS_beginning_rev; #key1=chr, key2=pos, value=gene

my %TTS_ending_fwd; #key1=chr, key2=pos, value=gene
my %TTS_ending_rev; #key1=chr, key2=pos, value=gene
 
my %leader_region_search_fwd;
my %leader_region_search_rev;
my %CDS_region_search_fwd;
my %CDS_region_search_rev;
my %trailer_region_search_fwd;
my %trailer_region_search_rev;

my %CDS_counts_RIBOseq; #key = gene, value = count sum
my %leader_counts_RIBOseq; #key, value = count sum
my %trailer_counts_RIBOseq; #key, value = count sum
my %start_codon_counts_RIBOseq;
my %stop_codon_counts_RIBOseq;

my %CDS_counts_RNAseq; 
my %leader_counts_RNAseq; 
my %trailer_counts_RNAseq;
my %start_codon_counts_RNAseq;
my %stop_codon_counts_RNAseq;

my %CDS_counts_totalRNA;
my %leader_counts_totalRNA;
my %trailer_counts_totalRNA;
my %start_codon_counts_totalRNA;
my %stop_codon_counts_totalRNA;

my %CDS_counts_TCPseq_SSU; 
my %leader_counts_TCPseq_SSU;
my %trailer_counts_TCPseq_SSU;
my %start_codon_counts_TCPseq_SSU;
my %stop_codon_counts_TCPseq_SSU;

my %CDS_counts_TCPseq_LSU;
my %leader_counts_TCPseq_LSU;
my %trailer_counts_TCPseq_LSU;
my %start_codon_counts_TCPseq_LSU;
my %stop_codon_counts_TCPseq_LSU;

my %leader_length;

my %TSS_beginning_count_TCPseq_SSU; #use to count reads that begin at the TSS (5' most part of transcript)
my %TSS_beginning_count_TCPseq_LSU;
my %TSS_beginning_count_RNAseq;
my %TSS_beginning_count_totalRNA;
my %TSS_beginning_count_RIBOseq;

#my %TTS_ending_count_TCPseq_SSU; #use to count reads that begin at the TSS (5' most part of transcript)
#my %TTS_ending_count_TCPseq_LSU;
#my %TTS_ending_count_RNAseq;
#my %TTS_ending_count_totalRNA;
#my %TTS_ending_count_RIBOseq;

my %leader_region_count_TCPseq_SSU; 
my %leader_region_count_TCPseq_LSU;
my %leader_region_count_RNAseq;
my %leader_region_count_totalRNA;
my %leader_region_count_RIBOseq;

my %CDS_region_count_TCPseq_SSU;
my %CDS_region_count_TCPseq_LSU;
my %CDS_region_count_RNAseq;
my %CDS_region_count_totalRNA;
my %CDS_region_count_RIBOseq;

my %trailer_region_count_TCPseq_SSU;
my %trailer_region_count_TCPseq_LSU;
my %trailer_region_count_RNAseq;
my %trailer_region_count_totalRNA;
my %trailer_region_count_RIBOseq;

my %TSS_proximal_region_search_fwd;
my %TSS_proximal_region_search_rev;
my %TIS_proximal_region_search_fwd;
my %TIS_proximal_region_search_rev;

my %TSS_proximal_region_count_TCPseq_SSU;
my %TSS_proximal_region_count_TCPseq_LSU;
my %TSS_proximal_region_count_RNAseq;
my %TSS_proximal_region_count_totalRNA;
my %TSS_proximal_region_count_RIBOseq;

my %TIS_proximal_region_count_TCPseq_SSU;
my %TIS_proximal_region_count_TCPseq_LSU;
my %TIS_proximal_region_count_RNAseq;
my %TIS_proximal_region_count_totalRNA;
my %TIS_proximal_region_count_RIBOseq;

for my $gene (keys %gene_model_fwd){
    my $chr=$gene_2_chr{$gene};
    my $start_coord=$start_coord_fwd{$gene};
    my $stop_coord=$stop_coord_fwd{$gene};
    my $five_prime_coord=$five_prime_most_coord_fwd{$gene};
    my $three_prime_coord=$three_prime_most_coord_fwd{$gene};

    $leader_length{$gene}=$start_coord-$five_prime_coord;

    $CDS_counts_RIBOseq{$gene}=0;
    $leader_counts_RIBOseq{$gene}=0;
    $trailer_counts_RIBOseq{$gene}=0; 
    $start_codon_counts_RIBOseq{$gene}=0;
    $stop_codon_counts_RIBOseq{$gene}=0;

    $CDS_counts_RNAseq{$gene}=0; 
    $leader_counts_RNAseq{$gene}=0; 
    $trailer_counts_RNAseq{$gene}=0;
    $start_codon_counts_RNAseq{$gene}=0;
    $stop_codon_counts_RNAseq{$gene}=0;

    $CDS_counts_totalRNA{$gene}=0;
    $leader_counts_totalRNA{$gene}=0;
    $trailer_counts_totalRNA{$gene}=0;
    $start_codon_counts_totalRNA{$gene}=0;
    $stop_codon_counts_totalRNA{$gene}=0;

    $CDS_counts_TCPseq_SSU{$gene}=0;
    $leader_counts_TCPseq_SSU{$gene}=0;
    $trailer_counts_TCPseq_SSU{$gene}=0;
    $start_codon_counts_TCPseq_SSU{$gene}=0;
    $stop_codon_counts_TCPseq_SSU{$gene}=0;

    $CDS_counts_TCPseq_LSU{$gene}=0;
    $leader_counts_TCPseq_LSU{$gene}=0;
    $trailer_counts_TCPseq_LSU{$gene}=0;
    $start_codon_counts_TCPseq_LSU{$gene}=0;
    $stop_codon_counts_TCPseq_LSU{$gene}=0;

    $TSS_beginning_count_TCPseq_SSU{$gene}=0;
    $TSS_beginning_count_TCPseq_LSU{$gene}=0;
    $TSS_beginning_count_RNAseq{$gene}=0;
    $TSS_beginning_count_totalRNA{$gene}=0;
    $TSS_beginning_count_RIBOseq{$gene}=0;

#    $TTS_ending_count_TCPseq_SSU{$gene}=0;
#    $TTS_ending_count_TCPseq_LSU{$gene}=0;
#    $TTS_ending_count_RNAseq{$gene}=0;
#    $TTS_ending_count_totalRNA{$gene}=0;
#    $TTS_ending_count_RIBOseq{$gene}=0;
 
    $leader_region_count_TCPseq_SSU{$gene}=0;
    $leader_region_count_TCPseq_LSU{$gene}=0;
    $leader_region_count_RNAseq{$gene}=0;
    $leader_region_count_totalRNA{$gene}=0;
    $leader_region_count_RIBOseq{$gene}=0;

    $CDS_region_count_TCPseq_SSU{$gene}=0;
    $CDS_region_count_TCPseq_LSU{$gene}=0;
    $CDS_region_count_RNAseq{$gene}=0;
    $CDS_region_count_totalRNA{$gene}=0;
    $CDS_region_count_RIBOseq{$gene}=0;

    $trailer_region_count_TCPseq_SSU{$gene}=0;
    $trailer_region_count_TCPseq_LSU{$gene}=0;
    $trailer_region_count_RNAseq{$gene}=0;
    $trailer_region_count_totalRNA{$gene}=0;
    $trailer_region_count_RIBOseq{$gene}=0;

    $TSS_proximal_region_count_TCPseq_SSU{$gene}=0;
    $TSS_proximal_region_count_TCPseq_LSU{$gene}=0;
    $TSS_proximal_region_count_RNAseq{$gene}=0;
    $TSS_proximal_region_count_totalRNA{$gene}=0;
    $TSS_proximal_region_count_RIBOseq{$gene}=0;

    $TIS_proximal_region_count_TCPseq_SSU{$gene}=0;
    $TIS_proximal_region_count_TCPseq_LSU{$gene}=0;
    $TIS_proximal_region_count_RNAseq{$gene}=0;
    $TIS_proximal_region_count_totalRNA{$gene}=0;
    $TIS_proximal_region_count_RIBOseq{$gene}=0;

    for my $coord (sort {$a <=> $b} keys %{ $gene_model_fwd{$gene} } ){
        my $pos=$gene_model_fwd{$gene}{$coord};

        #I need to check that the coords are within the reannoatted leader!
        if ($coord >= $five_prime_coord){

            if (exists ($whole_transcripts_fwd{$chr}{$pos})){ #mark overlapping genes
                my @gene = split("_", $whole_transcripts_fwd{$chr}{$pos});
                for (@gene){ $overlaps_inframe_gene{$_}=1; }
                $overlaps_inframe_gene{$gene}=1;
                $whole_transcripts_fwd{$chr}{$pos}.="_".$gene;
            }else{
                $whole_transcripts_fwd{$chr}{$pos}=$gene;
            }

            if ( ($coord >= ($five_prime_coord+69)) && ($coord <= ($five_prime_coord+119) ) ){
                $TSS_proximal_region_search_fwd{$chr}{$pos}=$gene; 
            }  

            if ( ($coord >= ($start_coord-99)) && (($coord <= $start_coord-49)) ){
                $TIS_proximal_region_search_fwd{$chr}{$pos}=$gene; 
            }

            if ($coord <= $five_prime_coord+9){  $TSS_beginning_fwd{$chr}{$pos}=$gene; }
            if ($coord >= $three_prime_coord-9){ $TTS_ending_fwd{$chr}{$pos}=$gene;    }

            if ($coord == $start_coord){   $start_codons_search_fwd{$chr}{$pos}=$gene; } 
            if ($coord == $start_coord+1){ $start_codons_search_fwd{$chr}{$pos}=$gene; }
            if ($coord == $start_coord+2){ $start_codons_search_fwd{$chr}{$pos}=$gene; }

            if ($coord == $stop_coord){   $stop_codons_search_fwd{$chr}{$pos}=$gene; }
            if ($coord == $stop_coord+1){ $stop_codons_search_fwd{$chr}{$pos}=$gene; }
            if ($coord == $stop_coord+2){ $stop_codons_search_fwd{$chr}{$pos}=$gene; }
           
            if ($coord < $start_coord){
                $leader_search_fwd{$chr}{$pos}=$gene;
            }elsif($coord <= ($stop_coord+2)){  #limit to the stop codon
                $CDS_search_fwd{$chr}{$pos}=$gene;
            }elsif( $coord > ($stop_coord+2)){
                $trailer_search_fwd{$chr}{$pos}=$gene;
            }

            if ($coord < $start_coord){ #to count thw values across the whole of these regions
                $leader_region_search_fwd{$chr}{$pos}=$gene;
            }elsif($coord < ($stop_coord+2)){
                $CDS_region_search_fwd{$chr}{$pos}=$gene;
            }else{
                $trailer_region_search_fwd{$chr}{$pos}=$gene;
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

    $leader_length{$gene}=$start_coord-$five_prime_coord;

    $CDS_counts_RIBOseq{$gene}=0;    
    $leader_counts_RIBOseq{$gene}=0;
    $trailer_counts_RIBOseq{$gene}=0;
    $start_codon_counts_RIBOseq{$gene}=0;
    $stop_codon_counts_RIBOseq{$gene}=0;

    $CDS_counts_RNAseq{$gene}=0;
    $leader_counts_RNAseq{$gene}=0;
    $trailer_counts_RNAseq{$gene}=0;
    $start_codon_counts_RNAseq{$gene}=0;
    $stop_codon_counts_RNAseq{$gene}=0;

    $CDS_counts_totalRNA{$gene}=0;
    $leader_counts_totalRNA{$gene}=0;
    $trailer_counts_totalRNA{$gene}=0;
    $start_codon_counts_totalRNA{$gene}=0;
    $stop_codon_counts_totalRNA{$gene}=0;

    $CDS_counts_TCPseq_SSU{$gene}=0;
    $leader_counts_TCPseq_SSU{$gene}=0;
    $trailer_counts_TCPseq_SSU{$gene}=0;
    $start_codon_counts_TCPseq_SSU{$gene}=0;
    $stop_codon_counts_TCPseq_SSU{$gene}=0;

    $CDS_counts_TCPseq_LSU{$gene}=0;
    $leader_counts_TCPseq_LSU{$gene}=0;
    $trailer_counts_TCPseq_LSU{$gene}=0;
    $start_codon_counts_TCPseq_LSU{$gene}=0;
    $stop_codon_counts_TCPseq_LSU{$gene}=0;

    $TSS_beginning_count_TCPseq_SSU{$gene}=0;
    $TSS_beginning_count_TCPseq_LSU{$gene}=0;
    $TSS_beginning_count_RNAseq{$gene}=0;
    $TSS_beginning_count_totalRNA{$gene}=0;
    $TSS_beginning_count_RIBOseq{$gene}=0;

#    $TTS_ending_count_TCPseq_SSU{$gene}=0;
#    $TTS_ending_count_TCPseq_LSU{$gene}=0;
#    $TTS_ending_count_RNAseq{$gene}=0;
#    $TTS_ending_count_totalRNA{$gene}=0;
#    $TTS_ending_count_RIBOseq{$gene}=0;

    $leader_region_count_TCPseq_SSU{$gene}=0;
    $leader_region_count_TCPseq_LSU{$gene}=0;
    $leader_region_count_RNAseq{$gene}=0;
    $leader_region_count_totalRNA{$gene}=0;
    $leader_region_count_RIBOseq{$gene}=0;

    $CDS_region_count_TCPseq_SSU{$gene}=0;
    $CDS_region_count_TCPseq_LSU{$gene}=0;
    $CDS_region_count_RNAseq{$gene}=0;
    $CDS_region_count_totalRNA{$gene}=0;
    $CDS_region_count_RIBOseq{$gene}=0;

    $trailer_region_count_TCPseq_SSU{$gene}=0;
    $trailer_region_count_TCPseq_LSU{$gene}=0;
    $trailer_region_count_RNAseq{$gene}=0;
    $trailer_region_count_totalRNA{$gene}=0;
    $trailer_region_count_RIBOseq{$gene}=0;

    $TSS_proximal_region_count_TCPseq_SSU{$gene}=0;
    $TSS_proximal_region_count_TCPseq_LSU{$gene}=0;
    $TSS_proximal_region_count_RNAseq{$gene}=0;
    $TSS_proximal_region_count_totalRNA{$gene}=0;
    $TSS_proximal_region_count_RIBOseq{$gene}=0;

    $TIS_proximal_region_count_TCPseq_SSU{$gene}=0;
    $TIS_proximal_region_count_TCPseq_LSU{$gene}=0;
    $TIS_proximal_region_count_RNAseq{$gene}=0;
    $TIS_proximal_region_count_totalRNA{$gene}=0;
    $TIS_proximal_region_count_RIBOseq{$gene}=0;

    for my $coord (sort {$a <=> $b} keys %{ $gene_model_rev{$gene} } ){
        my $pos=$gene_model_rev{$gene}{$coord};

        #I need to check that the coords are within the reannoatted leader!
        if ($coord >= $five_prime_coord){
        
            if (exists ($whole_transcripts_rev{$chr}{$pos})){ #mark overlapping genes
                my @gene = split("_", $whole_transcripts_rev{$chr}{$pos});
                for (@gene){ $overlaps_inframe_gene{$_}=1; }
                $overlaps_inframe_gene{$gene}=1;
                $whole_transcripts_rev{$chr}{$pos}.="_".$gene;
            }else{
                $whole_transcripts_rev{$chr}{$pos}=$gene;
            }

            if ( ($coord >= ($five_prime_coord+69)) && ($coord <= ($five_prime_coord+119) ) ){
                $TSS_proximal_region_search_rev{$chr}{$pos}=$gene;
            }

            if ( ($coord >= ($start_coord-99)) && (($coord <= $start_coord-49)) ){
                $TIS_proximal_region_search_rev{$chr}{$pos}=$gene;
            }

            if ($coord <= $five_prime_coord+9){  $TSS_beginning_rev{$chr}{$pos}=$gene; }
            if ($coord >= $three_prime_coord-9){ $TTS_ending_rev{$chr}{$pos}=$gene;    }

            if ($coord == $start_coord){   $start_codons_search_rev{$chr}{$pos}=$gene; }
            if ($coord == $start_coord+1){ $start_codons_search_rev{$chr}{$pos}=$gene; }
            if ($coord == $start_coord+2){ $start_codons_search_rev{$chr}{$pos}=$gene; }

            if ($coord == $stop_coord){   $stop_codons_search_rev{$chr}{$pos}=$gene; }
            if ($coord == $stop_coord+1){ $stop_codons_search_rev{$chr}{$pos}=$gene; }
            if ($coord == $stop_coord+2){ $stop_codons_search_rev{$chr}{$pos}=$gene; }

            if ($coord < $start_coord){
                $leader_search_rev{$chr}{$pos}=$gene;
            }elsif($coord <= ($stop_coord+2)){
                $CDS_search_rev{$chr}{$pos}=$gene;
            }elsif( $coord > ($stop_coord+2)){
                $trailer_search_rev{$chr}{$pos}=$gene;
            }

            if ($coord < $start_coord){ #to count thw values across the whole of these regions
                $leader_region_search_rev{$chr}{$pos}=$gene;
            }elsif($coord < ($stop_coord+2)){
                $CDS_region_search_rev{$chr}{$pos}=$gene;
            }else{
                $trailer_region_search_rev{$chr}{$pos}=$gene;
            }
        }
    }
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#leader region motif scan 
my %cert_count;
my %cert_short_count;
my %polyU_count;

my $longC_fwd=0;
my $shortC_fwd=0;
my $u_fwd=0;

for my $gene (keys %gene_model_fwd){
    my $chr=$gene_2_chr{$gene};
    my $start_coord=$start_coord_fwd{$gene};
    my $five_prime_coord=$five_prime_most_coord_fwd{$gene};

    $cert_count{$gene}=0;
    $cert_short_count{$gene}=0;
    $polyU_count{$gene}=0;

    my $leader_length=$start_coord-$five_prime_coord;
    my $leader_seq=substr($fasta_sequences{$chr}, ($gene_model_fwd{$gene}{$five_prime_coord} -1), $leader_length);

    while ($leader_seq =~ /C\wCCGC\w\w[CG]C/g){
        $cert_short_count{$gene}++;
    #    print "$gene,$leader_seq\nCERT start:$-[0]\nCERT stop:$+[0]\n";
        $shortC_fwd++;
    }

    while ($leader_seq =~ /[CG]\w[CG][CG]\wC\wCCGC\w\w[CG]C/g){
        $cert_count{$gene}++;
        #print "$gene,$leader_seq\nCERT start:$-[0]\nCERT stop:$+[0]\n";
        $longC_fwd++;
    }

    while ($leader_seq =~ /TTTTTTT/g){   #by default this looks for non overlapping matches
        $polyU_count{$gene}++;
        #print "$gene,$leader_seq\nPolyU start:$-[0]\npolyU stop:$+[0]\n";
        $u_fwd++;
    }
}

#print "fwd short cert $shortC_fwd\n";
#print "fwd long cert $longC_fwd\n";
#print "fwd poly U $u_fwd\n";

my $longC_rev=0;
my $shortC_rev=0;
my $u_rev=0;

for my $gene (keys %gene_model_rev){
    my $chr=$gene_2_chr{$gene};
    my $start_coord=$start_coord_rev{$gene};
    my $five_prime_coord=$five_prime_most_coord_rev{$gene};

    $cert_count{$gene}=0;
    $cert_short_count{$gene}=0;
    $polyU_count{$gene}=0;

    my $leader_length=$start_coord-$five_prime_coord;
    my $leader_seq=substr($fasta_sequences{$chr}, ($gene_model_rev{$gene}{$five_prime_coord} -1), $leader_length);
    $leader_seq=~tr/ACGTacgt/TGCAtgca/;

    while ($leader_seq =~ /C\wCCGC\w\w[CG]C/g){
        $cert_short_count{$gene}++;
    #    print "$gene,$leader_seq\nCERT start:$-[0]\nCERT stop:$+[0]\n";
        $shortC_rev++;
    }

    while ($leader_seq =~ /[CG]\w[CG][CG]\wC\wCCGC\w\w[CG]C/g){
        $cert_count{$gene}++;
        #print "$gene,$leader_seq\nCERT start:$-[0]\nCERT stop:$+[0]\n";
        $longC_rev++;
    }

    while ($leader_seq =~ /TTTTTTT/g){   #can I make this non overlapping? (split?)
        $polyU_count{$gene}++;
        #print "$gene,$leader_seq\nPolyU start:$-[0]\npolyU stop:$+[0]\n";
        $u_rev++;
    }
}

#print "rev short cert $shortC_rev\n";
#print "rev long cert $longC_rev\n";
#print "rev poly U $u_rev\n";

#exit;

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#open riboseq and assign
my $total_RIBOseq_bam_count=0;

open BAM0,"samtools view $bam_RIBOseq |";
while(<BAM0>){

    next if(/^(\@)/);  ## skipping the header lines (if you used -h in the samools command)
    s/\n//;  s/\r//;  ## removing new line
    my @sam = split(/\t+/);  ## splitting SAM line into array

    my $leftMost=$sam[3]; #leftmost position of match 5' for fwd, 3' for rev
    my $flag=$sam[1];
    my $chr=$sam[2];
    my $mapq=$sam[4];
    my $cigar=$sam[5];
    my $CDS_hit=0;
    my $leader_hit=0;
    my $trailer_hit=0;
    my $start_hit=0;
    my $stop_hit=0;
    my $TSS_beginning_hit=0;
    my $TTS_ending_hit=0;

    my $region_leader_hit=0;
    my $region_CDS_hit=0;
    my $region_trailer_hit=0;

    my $proximal_TSS_hit=0;
    my $proximal_TIS_hit=0;

    if ($chr=~/chr(.*)/){ #ensembl chromosome names do not contain the "chr" prefix
        $chr=$1;
    }

    unless ($flag & 0x4){   #if aligned

        if ($mapq >= 10){     #mapping uniqnes filter

            #both chew_2013 and subtelney_2014 riboseq are directional
            if ($flag & 0x10){  #Reverse reads. Starting from the leftmost position parse the cigar and check if matching positions overlap leaders or cds's
    
                while ($cigar !~ /^$/){
                    if ($cigar =~ /^([0-9]+[MIDN])/){
                        my $cigar_part = $1;
                        if ($cigar_part =~ /(\d+)M/){   #alignment matching 
                            for my $pos ($leftMost .. ($leftMost+$1-1)){ #search though this position
    
                                if (exists ($TSS_beginning_rev{$chr}{$pos})){
                                    $TSS_beginning_hit=$TSS_beginning_rev{$chr}{$pos};
                                }
                                if (exists ($leader_search_rev{$chr}{$pos})){
                                    $leader_hit=$leader_search_rev{$chr}{$pos};
                                }
                                if (exists ($start_codons_search_rev{$chr}{$pos})){
                                    $start_hit=$start_codons_search_rev{$chr}{$pos};
                                }
                                if( exists($CDS_search_rev{$chr}{$pos})){
                                    $CDS_hit=$CDS_search_rev{$chr}{$pos};
                                }
                                if( exists($trailer_search_rev{$chr}{$pos})){
                                    $trailer_hit=$trailer_search_rev{$chr}{$pos};
                                }
                                if (exists ($stop_codons_search_rev{$chr}{$pos})){
                                    $stop_hit=$stop_codons_search_rev{$chr}{$pos};
                                }
                                if (exists ($TTS_ending_rev{$chr}{$pos})){
                                    $TTS_ending_hit=$TTS_ending_rev{$chr}{$pos};
                                }

                                if (exists ($leader_region_search_rev{$chr}{$pos})){
                                    $region_leader_hit=$leader_region_search_rev{$chr}{$pos};
                                }
                                if (exists ($CDS_region_search_rev{$chr}{$pos})){
                                    $region_CDS_hit=$CDS_region_search_rev{$chr}{$pos};
                                }
                                if (exists ($trailer_region_search_rev{$chr}{$pos})){
                                    $region_trailer_hit=$trailer_region_search_rev{$chr}{$pos};
                                }
                      
                                if (exists ($TSS_proximal_region_search_rev{$chr}{$pos})){
                                    $proximal_TSS_hit=$TSS_proximal_region_search_rev{$chr}{$pos};
                                }

                                if (exists ($TIS_proximal_region_search_rev{$chr}{$pos})){
                                    $proximal_TIS_hit=$TIS_proximal_region_search_rev{$chr}{$pos};
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
   
                                if (exists ($TSS_beginning_fwd{$chr}{$pos})){
                                    $TSS_beginning_hit=$TSS_beginning_fwd{$chr}{$pos};
                                }
                                if (exists ($leader_search_fwd{$chr}{$pos})){
                                    $leader_hit=$leader_search_fwd{$chr}{$pos};
                                }
                                if (exists ($start_codons_search_fwd{$chr}{$pos})){
                                    $start_hit=$start_codons_search_fwd{$chr}{$pos};
                                }
                                if( exists($CDS_search_fwd{$chr}{$pos})){
                                    $CDS_hit=$CDS_search_fwd{$chr}{$pos};
                                }
                                if( exists($trailer_search_fwd{$chr}{$pos})){
                                    $trailer_hit=$trailer_search_fwd{$chr}{$pos};
                                }
                                if (exists ($stop_codons_search_fwd{$chr}{$pos})){
                                    $stop_hit=$stop_codons_search_fwd{$chr}{$pos};
                                }
                                if (exists ($TTS_ending_fwd{$chr}{$pos})){
                                    $TTS_ending_hit=$TTS_ending_fwd{$chr}{$pos};
                                }

                                if (exists ($leader_region_search_fwd{$chr}{$pos})){
                                    $region_leader_hit=$leader_region_search_fwd{$chr}{$pos};
                                }
                                if (exists ($CDS_region_search_fwd{$chr}{$pos})){
                                    $region_CDS_hit=$CDS_region_search_fwd{$chr}{$pos};
                                }
                                if (exists ($trailer_region_search_fwd{$chr}{$pos})){
                                    $region_trailer_hit=$trailer_region_search_fwd{$chr}{$pos};
                                }

                                if (exists ($TSS_proximal_region_search_fwd{$chr}{$pos})){
                                    $proximal_TSS_hit=$TSS_proximal_region_search_fwd{$chr}{$pos};
                                }

                                if (exists ($TIS_proximal_region_search_fwd{$chr}{$pos})){
                                    $proximal_TIS_hit=$TIS_proximal_region_search_fwd{$chr}{$pos};
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

            #assignment preference:
            #TSS > TIS > Start codon > Stop codon > CDS > Leader > Trailer

            if ($TTS_ending_hit){
                #exclude 3' peaks caused by polyA selection
                #$TTS_ending_count_RIBOseq{$TTS_ending_hit}+=1;
                #$total_RIBOseq_bam_count++;
            }elsif ($TSS_beginning_hit){
                $TSS_beginning_count_RIBOseq{$TSS_beginning_hit}+=1;
                $total_RIBOseq_bam_count++;
            }elsif ($start_hit){
                $start_codon_counts_RIBOseq{$start_hit}+=1;
                $total_RIBOseq_bam_count++;
            }elsif ($stop_hit){
                $stop_codon_counts_RIBOseq{$stop_hit}+=1;
                $total_RIBOseq_bam_count++;
            }elsif ($CDS_hit){
                $CDS_counts_RIBOseq{$CDS_hit}+=1;
                $total_RIBOseq_bam_count++;
            }elsif($leader_hit){
                $leader_counts_RIBOseq{$leader_hit}+=1;
                $total_RIBOseq_bam_count++;
            }elsif($trailer_hit){
                $trailer_counts_RIBOseq{$trailer_hit}+=1;
                $total_RIBOseq_bam_count++;
            }


            #assignment preference:
            #TTS > CDS > Leader > Trailer 

            if ($TTS_ending_hit){ #do nothing. Exlcude the 3' polyA selected peaks
            }elsif ($region_CDS_hit){ #assign to CDS preferencially
                $CDS_region_count_RIBOseq{$region_CDS_hit}+=1;
            }elsif ($region_leader_hit){
                $leader_region_count_RIBOseq{$region_leader_hit}+=1;
            }elsif ($region_trailer_hit){
                $trailer_region_count_RIBOseq{$region_trailer_hit}+=1;
            }
         
            #assignment of proximal regions
            if ($proximal_TSS_hit){
                $TSS_proximal_region_count_RIBOseq{$proximal_TSS_hit}+=1;
            }
 
            if ($proximal_TIS_hit){
                $TIS_proximal_region_count_RIBOseq{$proximal_TIS_hit}+=1;
            }
        }
    }
}
close(BAM0);

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤t¤
#open and store RNAseq counts
my $total_RNAseq_count=0;

open BAM1,"samtools view $bam_RNAseq |";

while(<BAM1>){

    next if(/^(\@)/);  ## skipping the header lines (if you used -h in the samools command)
    s/\n//;  s/\r//;  ## removing new line
    my @sam = split(/\t+/);  ## splitting SAM line into array

    my $leftMost=$sam[3]; #leftmost position of match 5' for fwd, 3' for rev
    my $flag=$sam[1];
    my $chr=$sam[2];
    my $mapq=$sam[4];
    my $cigar=$sam[5];
    my $CDS_hit=0;
    my $start_hit=0;
    my $stop_hit=0;
    my $leader_hit=0;
    my $trailer_hit=0;
    my $TSS_beginning_hit=0;
    my $TTS_ending_hit=0;

    my $region_leader_hit=0;
    my $region_CDS_hit=0;
    my $region_trailer_hit=0;

    my $proximal_TSS_hit=0;
    my $proximal_TIS_hit=0;

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

                                if (exists ($TSS_beginning_rev{$chr}{$pos})){
                                    $TSS_beginning_hit=$TSS_beginning_rev{$chr}{$pos};
                                }
                                if (exists ($leader_search_rev{$chr}{$pos})){
                                    $leader_hit=$leader_search_rev{$chr}{$pos};
                                }
                                if (exists ($start_codons_search_rev{$chr}{$pos})){
                                    $start_hit=$start_codons_search_rev{$chr}{$pos};
                                }
                                if( exists($CDS_search_rev{$chr}{$pos})){
                                    $CDS_hit=$CDS_search_rev{$chr}{$pos};
                                }
                                if( exists($trailer_search_rev{$chr}{$pos})){
                                    $trailer_hit=$trailer_search_rev{$chr}{$pos};
                                }
                                if (exists ($stop_codons_search_rev{$chr}{$pos})){
                                    $stop_hit=$stop_codons_search_rev{$chr}{$pos};
                                }
                                if (exists ($TTS_ending_rev{$chr}{$pos})){
                                    $TTS_ending_hit=$TTS_ending_rev{$chr}{$pos};
                                }

                                if (exists ($leader_region_search_rev{$chr}{$pos})){
                                    $region_leader_hit=$leader_region_search_rev{$chr}{$pos};
                                }
                                if (exists ($CDS_region_search_rev{$chr}{$pos})){
                                    $region_CDS_hit=$CDS_region_search_rev{$chr}{$pos};
                                }
                                if (exists ($trailer_region_search_rev{$chr}{$pos})){
                                    $region_trailer_hit=$trailer_region_search_rev{$chr}{$pos};
                                }

                                if (exists ($TSS_proximal_region_search_rev{$chr}{$pos})){
                                    $proximal_TSS_hit=$TSS_proximal_region_search_rev{$chr}{$pos};
                                }

                                if (exists ($TIS_proximal_region_search_rev{$chr}{$pos})){
                                    $proximal_TIS_hit=$TIS_proximal_region_search_rev{$chr}{$pos};
                                }

                                if (exists ($TSS_beginning_fwd{$chr}{$pos})){
                                    $TSS_beginning_hit=$TSS_beginning_fwd{$chr}{$pos};
                                }
                                if (exists ($leader_search_fwd{$chr}{$pos})){
                                    $leader_hit=$leader_search_fwd{$chr}{$pos};
                                }
                                if (exists ($start_codons_search_fwd{$chr}{$pos})){
                                    $start_hit=$start_codons_search_fwd{$chr}{$pos};
                                }
                                if( exists($CDS_search_fwd{$chr}{$pos})){
                                    $CDS_hit=$CDS_search_fwd{$chr}{$pos};
                                }
                                if( exists($trailer_search_fwd{$chr}{$pos})){
                                    $trailer_hit=$trailer_search_fwd{$chr}{$pos};
                                }
                                if (exists ($stop_codons_search_fwd{$chr}{$pos})){
                                    $stop_hit=$stop_codons_search_fwd{$chr}{$pos};
                                }
                                if (exists ($TTS_ending_fwd{$chr}{$pos})){
                                    $TTS_ending_hit=$TTS_ending_fwd{$chr}{$pos};
                                }

                                if (exists ($leader_region_search_fwd{$chr}{$pos})){
                                    $region_leader_hit=$leader_region_search_fwd{$chr}{$pos};
                                }
                                if (exists ($CDS_region_search_fwd{$chr}{$pos})){
                                    $region_CDS_hit=$CDS_region_search_fwd{$chr}{$pos};
                                }
                                if (exists ($trailer_region_search_fwd{$chr}{$pos})){
                                    $region_trailer_hit=$trailer_region_search_fwd{$chr}{$pos};
                                }   

                                if (exists ($TSS_proximal_region_search_fwd{$chr}{$pos})){
                                    $proximal_TSS_hit=$TSS_proximal_region_search_fwd{$chr}{$pos};
                                }

                                if (exists ($TIS_proximal_region_search_fwd{$chr}{$pos})){
                                    $proximal_TIS_hit=$TIS_proximal_region_search_fwd{$chr}{$pos};
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
    
                                if (exists ($TSS_beginning_rev{$chr}{$pos})){
                                    $TSS_beginning_hit=$TSS_beginning_rev{$chr}{$pos};
                                }
                                if (exists ($leader_search_rev{$chr}{$pos})){
                                    $leader_hit=$leader_search_rev{$chr}{$pos};
                                }
                                if (exists ($start_codons_search_rev{$chr}{$pos})){
                                    $start_hit=$start_codons_search_rev{$chr}{$pos};
                                }
                                if( exists($CDS_search_rev{$chr}{$pos})){
                                    $CDS_hit=$CDS_search_rev{$chr}{$pos};
                                }
                                if( exists($trailer_search_rev{$chr}{$pos})){
                                    $trailer_hit=$trailer_search_rev{$chr}{$pos};
                                }
                                if (exists ($stop_codons_search_rev{$chr}{$pos})){
                                    $stop_hit=$stop_codons_search_rev{$chr}{$pos};
                                }
                                if (exists ($TTS_ending_rev{$chr}{$pos})){
                                    $TTS_ending_hit=$TTS_ending_rev{$chr}{$pos};
                                }

                                if (exists ($leader_region_search_rev{$chr}{$pos})){
                                    $region_leader_hit=$leader_region_search_rev{$chr}{$pos};
                                }
                                if (exists ($CDS_region_search_rev{$chr}{$pos})){
                                    $region_CDS_hit=$CDS_region_search_rev{$chr}{$pos};
                                }
                                if (exists ($trailer_region_search_rev{$chr}{$pos})){
                                    $region_trailer_hit=$trailer_region_search_rev{$chr}{$pos};
                                }

                                if (exists ($TSS_proximal_region_search_rev{$chr}{$pos})){
                                    $proximal_TSS_hit=$TSS_proximal_region_search_rev{$chr}{$pos};
                                }

                                if (exists ($TIS_proximal_region_search_rev{$chr}{$pos})){
                                    $proximal_TIS_hit=$TIS_proximal_region_search_rev{$chr}{$pos};
                                }

                                if (exists ($TSS_beginning_fwd{$chr}{$pos})){
                                    $TSS_beginning_hit=$TSS_beginning_fwd{$chr}{$pos};
                                }
                                if (exists ($leader_search_fwd{$chr}{$pos})){
                                    $leader_hit=$leader_search_fwd{$chr}{$pos};
                                }
                                if (exists ($start_codons_search_fwd{$chr}{$pos})){
                                    $start_hit=$start_codons_search_fwd{$chr}{$pos};
                                }
                                if( exists($CDS_search_fwd{$chr}{$pos})){
                                    $CDS_hit=$CDS_search_fwd{$chr}{$pos};
                                }
                                if( exists($trailer_search_fwd{$chr}{$pos})){
                                    $trailer_hit=$trailer_search_fwd{$chr}{$pos};
                                }
                                if (exists ($stop_codons_search_fwd{$chr}{$pos})){
                                    $stop_hit=$stop_codons_search_fwd{$chr}{$pos};
                                }
                                if (exists ($TTS_ending_fwd{$chr}{$pos})){
                                    $TTS_ending_hit=$TTS_ending_fwd{$chr}{$pos};
                                }

                                if (exists ($leader_region_search_fwd{$chr}{$pos})){
                                    $region_leader_hit=$leader_region_search_fwd{$chr}{$pos};
                                }
                                if (exists ($CDS_region_search_fwd{$chr}{$pos})){
                                    $region_CDS_hit=$CDS_region_search_fwd{$chr}{$pos};
                                }
                                if (exists ($trailer_region_search_fwd{$chr}{$pos})){
                                    $region_trailer_hit=$trailer_region_search_fwd{$chr}{$pos};
                                }

                                if (exists ($TSS_proximal_region_search_fwd{$chr}{$pos})){
                                    $proximal_TSS_hit=$TSS_proximal_region_search_fwd{$chr}{$pos};
                                }

                                if (exists ($TIS_proximal_region_search_fwd{$chr}{$pos})){
                                    $proximal_TIS_hit=$TIS_proximal_region_search_fwd{$chr}{$pos};
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
    
            #assignment preference:
            #TSS > TIS > Start codon > Stop codon > CDS > Leader > Trailer

            if ($TTS_ending_hit){
                #exclude 3' peaks caused by polyA selection
                #$TTS_ending_count_RNAseq{$TTS_ending_hit}+=1;
                #$total_RNAseq_bam_count++;
            }elsif ($TSS_beginning_hit){
                $TSS_beginning_count_RNAseq{$TSS_beginning_hit}+=1;
                $total_RNAseq_count++;
            }elsif ($start_hit){
                $start_codon_counts_RNAseq{$start_hit}+=1;
                $total_RNAseq_count++;
            }elsif ($stop_hit){
                $stop_codon_counts_RNAseq{$stop_hit}+=1;
                $total_RNAseq_count++;
            }elsif ($CDS_hit){
                $CDS_counts_RNAseq{$CDS_hit}+=1;
                $total_RNAseq_count++;
            }elsif($leader_hit){
                $leader_counts_RNAseq{$leader_hit}+=1;
                $total_RNAseq_count++;
            }elsif($trailer_hit){
                $trailer_counts_RNAseq{$trailer_hit}+=1;
                $total_RNAseq_count++;
            }      


            #assignment preference:
            #TTS > CDS > Leader > Trailer 

            if ($TTS_ending_hit){ #do nothing. Exlcude the 3' polyA selected peaks
            }elsif ($region_CDS_hit){ #assign to CDS preferencially
                $CDS_region_count_RNAseq{$region_CDS_hit}+=1;
            }elsif ($region_leader_hit){
                $leader_region_count_RNAseq{$region_leader_hit}+=1;
            }elsif ($region_trailer_hit){
                $trailer_region_count_RNAseq{$region_trailer_hit}+=1;
            }

            if ($proximal_TSS_hit){
                $TSS_proximal_region_count_RNAseq{$proximal_TSS_hit}+=1;
            }

            if ($proximal_TIS_hit){
                $TIS_proximal_region_count_RNAseq{$proximal_TIS_hit}+=1;
            }
        }
    }
}
close (BAM1);

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤t¤
#open and store totalRNA counts
my $total_totalRNA_count=0;

open BAM4,"samtools view $bam_totalRNA |";

while(<BAM4>){

    next if(/^(\@)/);  ## skipping the header lines (if you used -h in the samools command)
    s/\n//;  s/\r//;  ## removing new line
    my @sam = split(/\t+/);  ## splitting SAM line into array

    my $leftMost=$sam[3]; #leftmost position of match 5' for fwd, 3' for rev
    my $flag=$sam[1];
    my $chr=$sam[2];
    my $mapq=$sam[4];
    my $cigar=$sam[5];
    my $CDS_hit=0;
    my $start_hit=0;
    my $stop_hit=0;
    my $leader_hit=0;
    my $trailer_hit=0;
    my $TSS_beginning_hit=0;
    my $TTS_ending_hit=0;

    my $region_leader_hit=0;
    my $region_CDS_hit=0;
    my $region_trailer_hit=0;

    my $proximal_TSS_hit=0;
    my $proximal_TIS_hit=0;

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

          #                      if (exists ($TSS_beginning_rev{$chr}{$pos})){
          #                          $TSS_beginning_hit=$TSS_beginning_rev{$chr}{$pos};
          #                      }
          #                      if (exists ($leader_search_rev{$chr}{$pos})){
          #                          $leader_hit=$leader_search_rev{$chr}{$pos};
          #                      }
          #                      if (exists ($start_codons_search_rev{$chr}{$pos})){
          #                          $start_hit=$start_codons_search_rev{$chr}{$pos};
          #                      }
          #                      if( exists($CDS_search_rev{$chr}{$pos})){
          #                          $CDS_hit=$CDS_search_rev{$chr}{$pos};
          #                      }
          #                      if( exists($trailer_search_rev{$chr}{$pos})){
          #                          $trailer_hit=$trailer_search_rev{$chr}{$pos};
          #                      }
          #                      if (exists ($stop_codons_search_rev{$chr}{$pos})){
          #                          $stop_hit=$stop_codons_search_rev{$chr}{$pos};
          #                      }
          #                      if (exists ($TTS_ending_rev{$chr}{$pos})){
          #                          $TTS_ending_hit=$TTS_ending_rev{$chr}{$pos};
          #                      }

          #                      if (exists ($leader_region_search_rev{$chr}{$pos})){
          #                          $region_leader_hit=$leader_region_search_rev{$chr}{$pos};
          #                      }
          #                      if (exists ($CDS_region_search_rev{$chr}{$pos})){
          #                          $region_CDS_hit=$CDS_region_search_rev{$chr}{$pos};
          #                      }
          #                      if (exists ($trailer_region_search_rev{$chr}{$pos})){
          #                          $region_trailer_hit=$trailer_region_search_rev{$chr}{$pos};
          #                      }

           #                     if (exists ($TSS_proximal_region_search_rev{$chr}{$pos})){
           #                         $proximal_TSS_hit=$TSS_proximal_region_search_rev{$chr}{$pos};
           #                     }

           #                     if (exists ($TIS_proximal_region_search_rev{$chr}{$pos})){
           #                         $proximal_TIS_hit=$TIS_proximal_region_search_rev{$chr}{$pos};
           #                     }


                                if (exists ($TSS_beginning_fwd{$chr}{$pos})){
                                    $TSS_beginning_hit=$TSS_beginning_fwd{$chr}{$pos};
                                }
                                if (exists ($leader_search_fwd{$chr}{$pos})){
                                    $leader_hit=$leader_search_fwd{$chr}{$pos};
                                }
                                if (exists ($start_codons_search_fwd{$chr}{$pos})){
                                    $start_hit=$start_codons_search_fwd{$chr}{$pos};
                                }
                                if( exists($CDS_search_fwd{$chr}{$pos})){
                                    $CDS_hit=$CDS_search_fwd{$chr}{$pos};
                                }
                                if( exists($trailer_search_fwd{$chr}{$pos})){
                                    $trailer_hit=$trailer_search_fwd{$chr}{$pos};
                                }
                                if (exists ($stop_codons_search_fwd{$chr}{$pos})){
                                    $stop_hit=$stop_codons_search_fwd{$chr}{$pos};
                                }
                                if (exists ($TTS_ending_fwd{$chr}{$pos})){
                                    $TTS_ending_hit=$TTS_ending_fwd{$chr}{$pos};
                                }

                                if (exists ($leader_region_search_fwd{$chr}{$pos})){
                                    $region_leader_hit=$leader_region_search_fwd{$chr}{$pos};
                                }
                                if (exists ($CDS_region_search_fwd{$chr}{$pos})){
                                    $region_CDS_hit=$CDS_region_search_fwd{$chr}{$pos};
                                }
                                if (exists ($trailer_region_search_fwd{$chr}{$pos})){
                                    $region_trailer_hit=$trailer_region_search_fwd{$chr}{$pos};
                                }

                                if (exists ($TSS_proximal_region_search_fwd{$chr}{$pos})){
                                    $proximal_TSS_hit=$TSS_proximal_region_search_fwd{$chr}{$pos};
                                }

                                if (exists ($TIS_proximal_region_search_fwd{$chr}{$pos})){
                                    $proximal_TIS_hit=$TIS_proximal_region_search_fwd{$chr}{$pos};
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

                                if (exists ($TSS_beginning_rev{$chr}{$pos})){
                                    $TSS_beginning_hit=$TSS_beginning_rev{$chr}{$pos};
                                }
                                if (exists ($leader_search_rev{$chr}{$pos})){
                                    $leader_hit=$leader_search_rev{$chr}{$pos};
                                }
                                if (exists ($start_codons_search_rev{$chr}{$pos})){
                                    $start_hit=$start_codons_search_rev{$chr}{$pos};
                                }
                                if( exists($CDS_search_rev{$chr}{$pos})){
                                    $CDS_hit=$CDS_search_rev{$chr}{$pos};
                                }
                                if( exists($trailer_search_rev{$chr}{$pos})){
                                    $trailer_hit=$trailer_search_rev{$chr}{$pos};
                                }
                                if (exists ($stop_codons_search_rev{$chr}{$pos})){
                                    $stop_hit=$stop_codons_search_rev{$chr}{$pos};
                                }
                                if (exists ($TTS_ending_rev{$chr}{$pos})){
                                    $TTS_ending_hit=$TTS_ending_rev{$chr}{$pos};
                                }

                                if (exists ($leader_region_search_rev{$chr}{$pos})){
                                    $region_leader_hit=$leader_region_search_rev{$chr}{$pos};
                                }
                                if (exists ($CDS_region_search_rev{$chr}{$pos})){
                                    $region_CDS_hit=$CDS_region_search_rev{$chr}{$pos};
                                }
                                if (exists ($trailer_region_search_rev{$chr}{$pos})){
                                    $region_trailer_hit=$trailer_region_search_rev{$chr}{$pos};
                                }

                                if (exists ($TSS_proximal_region_search_rev{$chr}{$pos})){
                                    $proximal_TSS_hit=$TSS_proximal_region_search_rev{$chr}{$pos};
                                }

                                if (exists ($TIS_proximal_region_search_rev{$chr}{$pos})){
                                    $proximal_TIS_hit=$TIS_proximal_region_search_rev{$chr}{$pos};
                                }

            #                    if (exists ($TSS_beginning_fwd{$chr}{$pos})){
            #                        $TSS_beginning_hit=$TSS_beginning_fwd{$chr}{$pos};
            #                    }
            #                    if (exists ($leader_search_fwd{$chr}{$pos})){
            #                        $leader_hit=$leader_search_fwd{$chr}{$pos};
            #                    }
            #                    if (exists ($start_codons_search_fwd{$chr}{$pos})){
            #                        $start_hit=$start_codons_search_fwd{$chr}{$pos};
            #                    }
            #                    if( exists($CDS_search_fwd{$chr}{$pos})){
            #                        $CDS_hit=$CDS_search_fwd{$chr}{$pos};
            #                    }
            #                    if( exists($trailer_search_fwd{$chr}{$pos})){
            #                        $trailer_hit=$trailer_search_fwd{$chr}{$pos};
            #                    }
            #                    if (exists ($stop_codons_search_fwd{$chr}{$pos})){
            #                        $stop_hit=$stop_codons_search_fwd{$chr}{$pos};
            #                    }
            #                    if (exists ($TTS_ending_fwd{$chr}{$pos})){
            #                        $TTS_ending_hit=$TTS_ending_fwd{$chr}{$pos};
            #                    }

            #                    if (exists ($leader_region_search_fwd{$chr}{$pos})){
            #                        $region_leader_hit=$leader_region_search_fwd{$chr}{$pos};
            #                    }
            #                    if (exists ($CDS_region_search_fwd{$chr}{$pos})){
            #                        $region_CDS_hit=$CDS_region_search_fwd{$chr}{$pos};
            #                    }
            #                    if (exists ($trailer_region_search_fwd{$chr}{$pos})){
            #                        $region_trailer_hit=$trailer_region_search_fwd{$chr}{$pos};
            #                    }

            #                    if (exists ($TSS_proximal_region_search_fwd{$chr}{$pos})){
            #                        $proximal_TSS_hit=$TSS_proximal_region_search_fwd{$chr}{$pos};
            #                    }

            #                    if (exists ($TIS_proximal_region_search_fwd{$chr}{$pos})){
            #                        $proximal_TIS_hit=$TIS_proximal_region_search_fwd{$chr}{$pos};
            #                    }

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

            #assignment preference:
            #TSS > TIS > Start codon > Stop codon > CDS > Leader > Trailer

            if ($TTS_ending_hit){
                #exclude 3' peaks caused by polyA selection
                #$TTS_ending_count_totalRNA{$TTS_ending_hit}+=1;
                #$total_totalRNA_bam_count++;
            }elsif ($TSS_beginning_hit){
                $TSS_beginning_count_totalRNA{$TSS_beginning_hit}+=1;
                $total_totalRNA_count++;
            }elsif ($start_hit){
                $start_codon_counts_totalRNA{$start_hit}+=1;
                $total_totalRNA_count++;
            }elsif ($stop_hit){
                $stop_codon_counts_totalRNA{$stop_hit}+=1;
                $total_totalRNA_count++;
            }elsif ($CDS_hit){
                $CDS_counts_totalRNA{$CDS_hit}+=1;
                $total_totalRNA_count++;
            }elsif($leader_hit){
                $leader_counts_totalRNA{$leader_hit}+=1;
                $total_totalRNA_count++;
            }elsif($trailer_hit){
                $trailer_counts_totalRNA{$trailer_hit}+=1;
                $total_totalRNA_count++;
            }

            #assignment preference:
            #TTS > CDS > Leader > Trailer 

            if ($TTS_ending_hit){ #do nothing. Exlcude the 3' polyA selected peaks
            }elsif ($region_CDS_hit){ #assign to CDS preferencially
                $CDS_region_count_totalRNA{$region_CDS_hit}+=1;
            }elsif ($region_leader_hit){
                $leader_region_count_totalRNA{$region_leader_hit}+=1;
            }elsif ($region_trailer_hit){
                $trailer_region_count_totalRNA{$region_trailer_hit}+=1;
            }

            if ($proximal_TSS_hit){
                $TSS_proximal_region_count_totalRNA{$proximal_TSS_hit}+=1;
            }

            if ($proximal_TIS_hit){
                $TIS_proximal_region_count_totalRNA{$proximal_TIS_hit}+=1;
            }
        }
    }
}
close (BAM4);

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#open and store TCPseq SSU counts
my $total_TCPseq_SSU_count=0;

open BAM2,"samtools view $bam_TCPseq_SSU |";

while(<BAM2>){

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
    my $TSS_beginning_hit=0;
    my $TTS_ending_hit=0;

    my $region_leader_hit=0;
    my $region_CDS_hit=0;
    my $region_trailer_hit=0;

    my $proximal_TSS_hit=0;
    my $proximal_TIS_hit=0;

    if ($chr=~/chr(.*)/){ #ensembl chromosome names do not contain the "chr" prefix
        $chr=$1;
    }

    unless (length($seq) < 17){

        unless ($flag & 0x4){   #if aligned

            if ($mapq >= 10){     #mapping uniqnes filter

                if ($flag & 0x10){  #Reverse reads. Starting from the leftmost position parse the cigar and check if matching positions overlap leaders or cds's
        
                    while ($cigar !~ /^$/){
                        if ($cigar =~ /^([0-9]+[MIDN])/){
                            my $cigar_part = $1;
                            if ($cigar_part =~ /(\d+)M/){   #alignment matching 
                                for my $pos ($leftMost .. ($leftMost+$1-1)){ #search though this position
        
                                    if (exists ($TSS_beginning_rev{$chr}{$pos})){
                                        $TSS_beginning_hit=$TSS_beginning_rev{$chr}{$pos};
                                    }
                                    if (exists ($leader_search_rev{$chr}{$pos})){
                                        $leader_hit=$leader_search_rev{$chr}{$pos};
                                    }
                                    if (exists ($start_codons_search_rev{$chr}{$pos})){
                                        $start_hit=$start_codons_search_rev{$chr}{$pos};
                                    }
                                    if( exists($CDS_search_rev{$chr}{$pos})){
                                        $CDS_hit=$CDS_search_rev{$chr}{$pos};
                                    }
                                    if( exists($trailer_search_rev{$chr}{$pos})){
                                        $trailer_hit=$trailer_search_rev{$chr}{$pos};
                                    }
                                    if (exists ($stop_codons_search_rev{$chr}{$pos})){
                                        $stop_hit=$stop_codons_search_rev{$chr}{$pos};
                                    }
                                    if (exists ($TTS_ending_rev{$chr}{$pos})){
                                        $TTS_ending_hit=$TTS_ending_rev{$chr}{$pos};
                                    }
    
                                    if (exists ($leader_region_search_rev{$chr}{$pos})){
                                        $region_leader_hit=$leader_region_search_rev{$chr}{$pos};
                                    }
                                    if (exists ($CDS_region_search_rev{$chr}{$pos})){
                                        $region_CDS_hit=$CDS_region_search_rev{$chr}{$pos};
                                    }
                                    if (exists ($trailer_region_search_rev{$chr}{$pos})){
                                        $region_trailer_hit=$trailer_region_search_rev{$chr}{$pos};
                                    }

                                    if (exists ($TSS_proximal_region_search_rev{$chr}{$pos})){
                                        $proximal_TSS_hit=$TSS_proximal_region_search_rev{$chr}{$pos};
                                    }

                                    if (exists ($TIS_proximal_region_search_rev{$chr}{$pos})){
                                        $proximal_TIS_hit=$TIS_proximal_region_search_rev{$chr}{$pos};
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
 
                                    if (exists ($TSS_beginning_fwd{$chr}{$pos})){
                                        $TSS_beginning_hit=$TSS_beginning_fwd{$chr}{$pos};
                                    }
                                    if (exists ($leader_search_fwd{$chr}{$pos})){
                                        $leader_hit=$leader_search_fwd{$chr}{$pos};
                                    }
                                    if (exists ($start_codons_search_fwd{$chr}{$pos})){
                                        $start_hit=$start_codons_search_fwd{$chr}{$pos};
                                    }
                                    if( exists($CDS_search_fwd{$chr}{$pos})){
                                        $CDS_hit=$CDS_search_fwd{$chr}{$pos};
                                    }
                                    if( exists($trailer_search_fwd{$chr}{$pos})){
                                        $trailer_hit=$trailer_search_fwd{$chr}{$pos};
                                    }
                                    if (exists ($stop_codons_search_fwd{$chr}{$pos})){
                                        $stop_hit=$stop_codons_search_fwd{$chr}{$pos};
                                    }
                                    if (exists ($TTS_ending_fwd{$chr}{$pos})){
                                        $TTS_ending_hit=$TTS_ending_fwd{$chr}{$pos};
                                    }
    
                                    if (exists ($leader_region_search_fwd{$chr}{$pos})){
                                        $region_leader_hit=$leader_region_search_fwd{$chr}{$pos};
                                    }
                                    if (exists ($CDS_region_search_fwd{$chr}{$pos})){
                                        $region_CDS_hit=$CDS_region_search_fwd{$chr}{$pos};
                                    }
                                    if (exists ($trailer_region_search_fwd{$chr}{$pos})){
                                        $region_trailer_hit=$trailer_region_search_fwd{$chr}{$pos};
                                    }       

                                    if (exists ($TSS_proximal_region_search_fwd{$chr}{$pos})){
                                        $proximal_TSS_hit=$TSS_proximal_region_search_fwd{$chr}{$pos};
                                    }

                                    if (exists ($TIS_proximal_region_search_fwd{$chr}{$pos})){
                                        $proximal_TIS_hit=$TIS_proximal_region_search_fwd{$chr}{$pos};
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

                #assignment preference:
                #TSS > TIS > Start codon > Stop codon > CDS > Leader > Trailer

                if ($TTS_ending_hit){
                    #exclude 3' peaks caused by polyA selection
                    #$TTS_ending_count_TCPseq_SSU{$TTS_ending_hit}+=1;
                    #$total_TCPseq_SSU_bam_count++;
                }elsif ($TSS_beginning_hit){
                    $TSS_beginning_count_TCPseq_SSU{$TSS_beginning_hit}+=1;
                    $total_TCPseq_SSU_count++;
                }elsif ($start_hit){
                    $start_codon_counts_TCPseq_SSU{$start_hit}+=1;
                    $total_TCPseq_SSU_count++;
                }elsif ($stop_hit){
                    $stop_codon_counts_TCPseq_SSU{$stop_hit}+=1;
                    $total_TCPseq_SSU_count++;
                }elsif ($CDS_hit){
                    $CDS_counts_TCPseq_SSU{$CDS_hit}+=1;
                    $total_TCPseq_SSU_count++;
                }elsif($leader_hit){
                    $leader_counts_TCPseq_SSU{$leader_hit}+=1;
                    $total_TCPseq_SSU_count++;
                }elsif($trailer_hit){
                    $trailer_counts_TCPseq_SSU{$trailer_hit}+=1;
                    $total_TCPseq_SSU_count++;
                }    


                #assignment preference:
                #TTS > CDS > Leader > Trailer 

                if ($TTS_ending_hit){ #do nothing. Exlcude the 3' polyA selected peaks
                }elsif ($region_CDS_hit){ #assign to CDS preferencially
                    $CDS_region_count_TCPseq_SSU{$region_CDS_hit}+=1;
                }elsif ($region_leader_hit){
                    $leader_region_count_TCPseq_SSU{$region_leader_hit}+=1;
                }elsif ($region_trailer_hit){
                    $trailer_region_count_TCPseq_SSU{$region_trailer_hit}+=1;
                }

                if ($proximal_TSS_hit){
                    $TSS_proximal_region_count_TCPseq_SSU{$proximal_TSS_hit}+=1;
                }

                if ($proximal_TIS_hit){
                    $TIS_proximal_region_count_TCPseq_SSU{$proximal_TIS_hit}+=1;
                }
            }
        }
    }
}
close (BAM2);

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
#
#open and store TPCseq LSU counts
my $total_TCPseq_LSU_count=0;

open BAM3,"samtools view $bam_TCPseq_LSU |";

while(<BAM3>){

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
    my $TSS_beginning_hit=0;
    my $TTS_ending_hit=0;

    my $region_leader_hit=0;
    my $region_CDS_hit=0;
    my $region_trailer_hit=0;

    my $proximal_TSS_hit=0;
    my $proximal_TIS_hit=0;

    if ($chr=~/chr(.*)/){ #ensembl chromosome names do not contain the "chr" prefix
        $chr=$1;
    }

    unless (length($seq) < 17){

        unless ($flag & 0x4){   #if aligned
         
            if ($mapq >= 10){     #mapping uniqnes filter

                if ($flag & 0x10){  #Reverse reads. Starting from the leftmost position parse the cigar and check if matching positions overlap leaders or cds's
                    
                    while ($cigar !~ /^$/){
                        if ($cigar =~ /^([0-9]+[MIDN])/){
                            my $cigar_part = $1;
                            if ($cigar_part =~ /(\d+)M/){   #alignment matching 
                                for my $pos ($leftMost .. ($leftMost+$1-1)){ #search though this position
        
                                    if (exists ($TSS_beginning_rev{$chr}{$pos})){
                                        $TSS_beginning_hit=$TSS_beginning_rev{$chr}{$pos};
                                    }
                                    if (exists ($leader_search_rev{$chr}{$pos})){
                                        $leader_hit=$leader_search_rev{$chr}{$pos};
                                    }
                                    if (exists ($start_codons_search_rev{$chr}{$pos})){
                                        $start_hit=$start_codons_search_rev{$chr}{$pos};
                                    }
                                    if( exists($CDS_search_rev{$chr}{$pos})){
                                        $CDS_hit=$CDS_search_rev{$chr}{$pos};
                                    }
                                    if( exists($trailer_search_rev{$chr}{$pos})){
                                        $trailer_hit=$trailer_search_rev{$chr}{$pos};
                                    }
                                    if (exists ($stop_codons_search_rev{$chr}{$pos})){
                                        $stop_hit=$stop_codons_search_rev{$chr}{$pos};
                                    }
                                    if (exists ($TTS_ending_rev{$chr}{$pos})){
                                        $TTS_ending_hit=$TTS_ending_rev{$chr}{$pos};
                                    }
    
                                    if (exists ($leader_region_search_rev{$chr}{$pos})){
                                        $region_leader_hit=$leader_region_search_rev{$chr}{$pos};
                                    }
                                    if (exists ($CDS_region_search_rev{$chr}{$pos})){
                                        $region_CDS_hit=$CDS_region_search_rev{$chr}{$pos};
                                    }
                                    if (exists ($trailer_region_search_rev{$chr}{$pos})){
                                        $region_trailer_hit=$trailer_region_search_rev{$chr}{$pos};
                                    }    

                                    if (exists ($TSS_proximal_region_search_rev{$chr}{$pos})){
                                        $proximal_TSS_hit=$TSS_proximal_region_search_rev{$chr}{$pos};
                                    }

                                    if (exists ($TIS_proximal_region_search_rev{$chr}{$pos})){
                                        $proximal_TIS_hit=$TIS_proximal_region_search_rev{$chr}{$pos};
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
            
                                    if (exists ($TSS_beginning_fwd{$chr}{$pos})){
                                        $TSS_beginning_hit=$TSS_beginning_fwd{$chr}{$pos};
                                    }
                                    if (exists ($leader_search_fwd{$chr}{$pos})){
                                        $leader_hit=$leader_search_fwd{$chr}{$pos};
                                    }
                                    if (exists ($start_codons_search_fwd{$chr}{$pos})){
                                        $start_hit=$start_codons_search_fwd{$chr}{$pos};
                                    }
                                    if( exists($CDS_search_fwd{$chr}{$pos})){
                                        $CDS_hit=$CDS_search_fwd{$chr}{$pos};
                                    }
                                    if( exists($trailer_search_fwd{$chr}{$pos})){
                                        $trailer_hit=$trailer_search_fwd{$chr}{$pos};
                                    }
                                    if (exists ($stop_codons_search_fwd{$chr}{$pos})){
                                        $stop_hit=$stop_codons_search_fwd{$chr}{$pos};
                                    }
                                    if (exists ($TTS_ending_fwd{$chr}{$pos})){
                                        $TTS_ending_hit=$TTS_ending_fwd{$chr}{$pos};
                                    }
    
                                    if (exists ($leader_region_search_fwd{$chr}{$pos})){
                                        $region_leader_hit=$leader_region_search_fwd{$chr}{$pos};
                                    }
                                    if (exists ($CDS_region_search_fwd{$chr}{$pos})){
                                        $region_CDS_hit=$CDS_region_search_fwd{$chr}{$pos};
                                    }
                                    if (exists ($trailer_region_search_fwd{$chr}{$pos})){
                                        $region_trailer_hit=$trailer_region_search_fwd{$chr}{$pos};
                                    }

                                    if (exists ($TSS_proximal_region_search_fwd{$chr}{$pos})){
                                        $proximal_TSS_hit=$TSS_proximal_region_search_fwd{$chr}{$pos};
                                    }

                                    if (exists ($TIS_proximal_region_search_fwd{$chr}{$pos})){
                                        $proximal_TIS_hit=$TIS_proximal_region_search_fwd{$chr}{$pos};
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
        
                #assignment preference:
                #TTS > TSS > Start codon > Stop codon > CDS > Leader > Trailer

                if ($TTS_ending_hit){
                    #exclude 3' peaks caused by polyA selection
                    #$TTS_ending_count_TCPseq_LSU{$TTS_ending_hit}+=1;
                    #$total_TCPseq_LSU_bam_count++;
                }elsif ($TSS_beginning_hit){
                    $TSS_beginning_count_TCPseq_LSU{$TSS_beginning_hit}+=1;
                    $total_TCPseq_LSU_count++;
                }elsif ($start_hit){
                    $start_codon_counts_TCPseq_LSU{$start_hit}+=1;
                    $total_TCPseq_LSU_count++;
                }elsif ($stop_hit){
                    $stop_codon_counts_TCPseq_LSU{$stop_hit}+=1;
                    $total_TCPseq_LSU_count++;
                }elsif ($CDS_hit){
                    $CDS_counts_TCPseq_LSU{$CDS_hit}+=1;
                    $total_TCPseq_LSU_count++;
                }elsif($leader_hit){
                    $leader_counts_TCPseq_LSU{$leader_hit}+=1;
                    $total_TCPseq_LSU_count++;
                }elsif($trailer_hit){
                    $trailer_counts_TCPseq_LSU{$trailer_hit}+=1;
                    $total_TCPseq_LSU_count++;
                }


                #assignment preference:
                #TTS > CDS > Leader > Trailer 

                if ($TTS_ending_hit){ #do nothing. Exlcude the 3' polyA selected peaks
                }elsif ($region_CDS_hit){ #assign to CDS preferencially
                    $CDS_region_count_TCPseq_LSU{$region_CDS_hit}+=1;
                }elsif ($region_leader_hit){
                    $leader_region_count_TCPseq_LSU{$region_leader_hit}+=1;
                }elsif ($region_trailer_hit){
                    $trailer_region_count_TCPseq_LSU{$region_trailer_hit}+=1;
                }

                if ($proximal_TSS_hit){
                    $TSS_proximal_region_count_TCPseq_LSU{$proximal_TSS_hit}+=1;
                }

                if ($proximal_TIS_hit){
                    $TIS_proximal_region_count_TCPseq_LSU{$proximal_TIS_hit}+=1;
                }
            }
        }
    }
}
close (BAM3);

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
#output
my $total_RIBOseq_aligned_sum=$total_RIBOseq_bam_count;
my $total_RNAseq_aligned_sum=$total_RNAseq_count;
my $total_totalRNA_aligned_sum=$total_totalRNA_count;
my $total_TCPseq_SSU_aligned_sum=$total_TCPseq_SSU_count;
my $total_TCPseq_LSU_aligned_sum=$total_TCPseq_LSU_count;

#header
print "#gene_id,transcript_id,chr,cage_tags_at_leader_peak,leader_length,CDS_length,trailer_length,TSS_beginning_RNA_FPKM,TSS_beginning_totalRNA_FPKM,TSS_beginning_RFP_FPKM,TSS_beginning_SSU_FPKM,TSS_beginning_LSU_FPKM,leader_RNA_FPKM,leader_totalRNA_FPKM,leader_RFP_FPKM,leader_SSU_FPKM,leader_LSU_FPKM,Start_codon_RNA_FPKM,Start_codon_totalRNA_FPKM,Start_codon_RFP_FPKM,Start_codon_SSU_FPKM,Start_codon_LSU_FPKM,CDS_RNA_FPKM,CDS_totalRNA_FPKM,CDS_RFP_FPKM,CDS_SSU_FPKM,CDS_LSU_FPKM,Stop_codon_RNA_FPKM,Stop_codon_totalRNA_FPKM,Stop_codon_RFP_FPKM,Stop_codon_SSU_FPKM,Stop_codon_LSU_FPKM,trailer_RNA_FPKM,trailer_totalRNA_FPKM,trailer_RFP_FPKM,trailer_SSU_FPKM,trailer_LSU_FPKM,complete_leader_RNA_FPKM,complete_leader_totalRNA_FPKM,complete_leader_RFP_FPKM,complete_leader_SSU_FPKM,complete_leader_LSU_FPKM,complete_CDS_RNA_FPKM,complete_CDS_totalRNA_FPKM,complete_CDS_RFP_FPKM,complete_CDS_SSU_FPKM,complete_CDS_LSU_FPKM,complete_trailer_RNA_FPKM,complete_trailer_totalRNA_FPKM,complete_trailer_RFP_FPKM,complete_trailer_SSU_FPKM,complete_trailer_LSU_FPKM,TSS_proximal_region_RNA_FPKM,TSS_proximal_region_totalRNA_FPKM,TSS_proximal_region_RFP_FPKM,TSS_proximal_region_SSU_FPKM,TSS_proximal_region_LSU_FPKM,TIS_proximal_region_RNA_FPKM,TIS_proximal_region_totalRNA_FPKM,TIS_proximal_region_RFP_FPKM,TIS_proximal_region_SSU_FPKM,TIS_proximal_region_LSU_FPKM,overlapping_gene,leader_potentially_overlaps_upstream_gene,gene_potentially_overlaps_downstream_leader,gene_overlaps_non_coding_transcript,Pre_TSS_start_sequence,TSS_start_sequence,Post_TSS_start_sequence,upstream_sequence,TIS_seqeuence,downstream_sequence,pre_stop_sequence,stop_codon_sequence,post_stop_sequence,exon_exon_junction_count_leader,kozak_strength,upstream_kozak_strength,TISU,histone,TOP,TCT,leader_cert_count,leader_short_cert_count,leader_polyU_count\n";

for my $gene (keys %gene_model_fwd){

    if (exists ($gene_2_chr{$gene})){  #restict to protien_coding genes
        my $chr=$gene_2_chr{$gene};
        my $transcript=$most_expressed_transcript{$gene};
        my $start_coord=$start_coord_fwd{$gene};
        my $stop_coord=$stop_coord_fwd{$gene};
        my $CDS_length=(($stop_coord+2)-$start_coord)+1;
        my $trailer_length=($three_prime_most_coord_fwd{$gene}-($stop_coord+3))+1;
        my $leader_start=$leader_positions_fwd{$gene};
        my $TSS_coord=$five_prime_most_coord_fwd{$gene};
        my $leader_length=$start_coord-$TSS_coord;
        my $highest_cage_peak_value=0;
        my $overlapping_gene="FALSE";
        my $leader_potentially_overlaps_upstream_gene="FALSE";
        my $gene_potentially_overlaps_downstream_leader="FALSE";
        my $gene_overlaps_nc="FALSE";
        my $is_histone="FALSE";
        my $exon_exon_junction_count_leader=$gene_exon_exon_in_leader{$gene};
        my $is_TOP="FALSE";
        my $is_TCT="FALSE";

        if (exists ($cage_peak_value{$gene})){ $highest_cage_peak_value=$cage_peak_value{$gene}; }
        if (exists ($overlaps_inframe_gene{$gene})){           $overlapping_gene="TRUE"; }
        if (exists ($leader_overlaps_upstream{$gene})){        $leader_potentially_overlaps_upstream_gene="TRUE"; }
        if (exists ($gene_overlaps_downstream_leader{$gene})){ $gene_potentially_overlaps_downstream_leader="TRUE"; }
        if (exists ($gene_overlaps_nc{$gene})){                $gene_overlaps_nc="TRUE"; }
        if (exists ($histone{$gene})){                         $is_histone="TRUE"; }

        my $TSS_beginning_count_RIBOseq=$TSS_beginning_count_RIBOseq{$gene};
        my $leader_count_RIBOseq=$leader_counts_RIBOseq{$gene};
        my $start_codon_count_RIBOseq=$start_codon_counts_RIBOseq{$gene};
        my $CDS_count_RIBOseq=$CDS_counts_RIBOseq{$gene};
        my $stop_codon_count_RIBOseq=$stop_codon_counts_RIBOseq{$gene};
        my $trailer_count_RIBOseq=$trailer_counts_RIBOseq{$gene};
#       my $TTS_ending_count_RIBOseq=$TTS_ending_count_RIBOseq{$gene};

        my $TSS_beginning_count_RNAseq=$TSS_beginning_count_RNAseq{$gene};
        my $leader_count_RNAseq=$leader_counts_RNAseq{$gene};
        my $start_codon_count_RNAseq=$start_codon_counts_RNAseq{$gene};
        my $CDS_count_RNAseq=$CDS_counts_RNAseq{$gene};
        my $stop_codon_count_RNAseq=$stop_codon_counts_RNAseq{$gene};
        my $trailer_count_RNAseq=$trailer_counts_RNAseq{$gene};
#       my $TTS_ending_count_RNAseq=$TTS_ending_count_RNAseq{$gene};

        my $TSS_beginning_count_totalRNA=$TSS_beginning_count_totalRNA{$gene};
        my $leader_count_totalRNA=$leader_counts_totalRNA{$gene};
        my $start_codon_count_totalRNA=$start_codon_counts_totalRNA{$gene};
        my $CDS_count_totalRNA=$CDS_counts_totalRNA{$gene};
        my $stop_codon_count_totalRNA=$stop_codon_counts_totalRNA{$gene};
        my $trailer_count_totalRNA=$trailer_counts_totalRNA{$gene};
#       my $TTS_ending_count_totalRNA=$TTS_ending_count_totalRNA{$gene};

        my $TSS_beginning_count_TCPseq_SSU=$TSS_beginning_count_TCPseq_SSU{$gene};
        my $leader_count_TCPseq_SSU=$leader_counts_TCPseq_SSU{$gene};
        my $start_codon_count_TCPseq_SSU=$start_codon_counts_TCPseq_SSU{$gene};
        my $CDS_count_TCPseq_SSU=$CDS_counts_TCPseq_SSU{$gene};
        my $stop_codon_count_TCPseq_SSU=$stop_codon_counts_TCPseq_SSU{$gene};
        my $trailer_count_TCPseq_SSU=$trailer_counts_TCPseq_SSU{$gene};
#       my $TTS_ending_count_TCPseq_SSU=$TTS_ending_count_TCPseq_SSU{$gene};

        my $TSS_beginning_count_TCPseq_LSU=$TSS_beginning_count_TCPseq_LSU{$gene};
        my $leader_count_TCPseq_LSU=$leader_counts_TCPseq_LSU{$gene};
        my $start_codon_count_TCPseq_LSU=$start_codon_counts_TCPseq_LSU{$gene};
        my $CDS_count_TCPseq_LSU=$CDS_counts_TCPseq_LSU{$gene};
        my $stop_codon_count_TCPseq_LSU=$stop_codon_counts_TCPseq_LSU{$gene};
        my $trailer_count_TCPseq_LSU=$trailer_counts_TCPseq_LSU{$gene};
#       my $TTS_ending_count_TCPseq_LSU=$TTS_ending_count_TCPseq_LSU{$gene};

        my $leader_region_count_TCPseq_SSU=$leader_region_count_TCPseq_SSU{$gene};
        my $leader_region_count_TCPseq_LSU=$leader_region_count_TCPseq_LSU{$gene};
        my $leader_region_count_RNAseq=$leader_region_count_RNAseq{$gene};
        my $leader_region_count_totalRNA=$leader_region_count_totalRNA{$gene};
        my $leader_region_count_RIBOseq=$leader_region_count_RIBOseq{$gene};

        my $CDS_region_count_TCPseq_SSU=$CDS_region_count_TCPseq_SSU{$gene};
        my $CDS_region_count_TCPseq_LSU=$CDS_region_count_TCPseq_LSU{$gene};
        my $CDS_region_count_RNAseq=$CDS_region_count_RNAseq{$gene};
        my $CDS_region_count_totalRNA=$CDS_region_count_totalRNA{$gene}; 
        my $CDS_region_count_RIBOseq=$CDS_region_count_RIBOseq{$gene};

        my $trailer_region_count_TCPseq_SSU=$trailer_region_count_TCPseq_SSU{$gene};
        my $trailer_region_count_TCPseq_LSU=$trailer_region_count_TCPseq_LSU{$gene};
        my $trailer_region_count_RNAseq=$trailer_region_count_RNAseq{$gene};
        my $trailer_region_count_totalRNA=$trailer_region_count_totalRNA{$gene};
        my $trailer_region_count_RIBOseq=$trailer_region_count_RIBOseq{$gene};

        my $TSS_proximal_region_count_TCPseq_SSU=$TSS_proximal_region_count_TCPseq_SSU{$gene};
        my $TSS_proximal_region_count_TCPseq_LSU=$TSS_proximal_region_count_TCPseq_LSU{$gene};
        my $TSS_proximal_region_count_RNAseq=$TSS_proximal_region_count_RNAseq{$gene};
        my $TSS_proximal_region_count_totalRNA=$TSS_proximal_region_count_totalRNA{$gene};
        my $TSS_proximal_region_count_RIBOseq=$TSS_proximal_region_count_RIBOseq{$gene};

        my $TIS_proximal_region_count_TCPseq_SSU=$TIS_proximal_region_count_TCPseq_SSU{$gene};
        my $TIS_proximal_region_count_TCPseq_LSU=$TIS_proximal_region_count_TCPseq_LSU{$gene};
        my $TIS_proximal_region_count_RNAseq=$TIS_proximal_region_count_RNAseq{$gene};
        my $TIS_proximal_region_count_totalRNA=$TIS_proximal_region_count_totalRNA{$gene};
        my $TIS_proximal_region_count_RIBOseq=$TIS_proximal_region_count_RIBOseq{$gene};

        #for TSS beginning, start codons and stop codons there is no need to adjust for length as they are all fixed at 1 or 3nt
        my $TSS_beginning_RIBOseq_FPKM=0;
	my $leader_RIBOseq_FPKM=0;
	my $start_codon_RIBOseq_FPKM=0;
        my $CDS_RIBOseq_FPKM=0;
        my $stop_codon_RIBOseq_FPKM=0;
        my $trailer_RIBOseq_FPKM=0;
#        my $TTS_ending_RIBOseq_FPKM=0;

        if ($TSS_beginning_count_RIBOseq > 0){ $TSS_beginning_RIBOseq_FPKM=eval{ ($TSS_beginning_count_RIBOseq/($total_RIBOseq_aligned_sum))*1000000              } || 0; }
        if ($leader_count_RIBOseq        > 0){ $leader_RIBOseq_FPKM=eval{        ($leader_count_RIBOseq/($leader_length*$total_RIBOseq_aligned_sum))*1000000000   } || 0; }
        if ($start_codon_count_RIBOseq   > 0){ $start_codon_RIBOseq_FPKM=eval{   ($start_codon_count_RIBOseq/($total_RIBOseq_aligned_sum))*1000000                } || 0; }
        if ($CDS_count_RIBOseq           > 0){ $CDS_RIBOseq_FPKM=eval{           ($CDS_count_RIBOseq/($CDS_length*$total_RIBOseq_aligned_sum))*1000000000         } || 0; }
        if ($stop_codon_count_RIBOseq    > 0){ $stop_codon_RIBOseq_FPKM=eval{    ($stop_codon_count_RIBOseq/($total_RIBOseq_aligned_sum))*1000000                 } || 0; }
        if ($trailer_count_RIBOseq       > 0){ $trailer_RIBOseq_FPKM=eval{       ($trailer_count_RIBOseq/($trailer_length*$total_RIBOseq_aligned_sum))*1000000000 } || 0; }
#        if ($TTS_ending_count_RIBOseq    > 0){ $TTS_ending_RIBOseq_FPKM=eval{    ($TTS_ending_count_RIBOseq/($total_RIBOseq_aligned_sum))*1000000                 } || 0; }

        my $TSS_beginning_RNAseq_FPKM=0;
	my $leader_RNAseq_FPKM=0;
	my $start_codon_RNAseq_FPKM=0;
        my $CDS_RNAseq_FPKM=0;
        my $stop_codon_RNAseq_FPKM=0;
        my $trailer_RNAseq_FPKM=0;
#        my $TTS_ending_RNAseq_FPKM=0;

        if ($TSS_beginning_count_RNAseq > 0){ $TSS_beginning_RNAseq_FPKM=eval{ ($TSS_beginning_count_RNAseq/($total_RNAseq_aligned_sum))*1000000              } || 0; }
        if ($leader_count_RNAseq        > 0){ $leader_RNAseq_FPKM=eval{        ($leader_count_RNAseq/($leader_length*$total_RNAseq_aligned_sum))*1000000000   } || 0; }
        if ($start_codon_count_RNAseq   > 0){ $start_codon_RNAseq_FPKM=eval{   ($start_codon_count_RNAseq/($total_RNAseq_aligned_sum))*1000000                } || 0; }
        if ($CDS_count_RNAseq           > 0){ $CDS_RNAseq_FPKM=eval{           ($CDS_count_RNAseq/($CDS_length*$total_RNAseq_aligned_sum))*1000000000         } || 0; }
        if ($stop_codon_count_RNAseq    > 0){ $stop_codon_RNAseq_FPKM=eval{    ($stop_codon_count_RNAseq/($total_RNAseq_aligned_sum))*1000000                 } || 0; }
        if ($trailer_count_RNAseq       > 0){ $trailer_RNAseq_FPKM=eval{       ($trailer_count_RNAseq/($trailer_length*$total_RNAseq_aligned_sum))*1000000000 } || 0; }
#        if ($TTS_ending_count_RNAseq    > 0){ $TTS_ending_RNAseq_FPKM=eval{    ($TTS_ending_count_RNAseq/($total_RNAseq_aligned_sum))*1000000                 } || 0; }

        my $TSS_beginning_totalRNA_FPKM=0;
        my $leader_totalRNA_FPKM=0;
        my $start_codon_totalRNA_FPKM=0;
        my $CDS_totalRNA_FPKM=0;
        my $stop_codon_totalRNA_FPKM=0;
        my $trailer_totalRNA_FPKM=0;
#        my $TTS_ending_totalRNA_FPKM=0;

        if ($TSS_beginning_count_totalRNA > 0){ $TSS_beginning_totalRNA_FPKM=eval{ ($TSS_beginning_count_totalRNA/($total_totalRNA_aligned_sum))*1000000              } || 0; }
        if ($leader_count_totalRNA        > 0){ $leader_totalRNA_FPKM=eval{        ($leader_count_totalRNA/($leader_length*$total_totalRNA_aligned_sum))*1000000000   } || 0; }
        if ($start_codon_count_totalRNA   > 0){ $start_codon_totalRNA_FPKM=eval{   ($start_codon_count_totalRNA/($total_totalRNA_aligned_sum))*1000000                } || 0; }
        if ($CDS_count_totalRNA           > 0){ $CDS_totalRNA_FPKM=eval{           ($CDS_count_totalRNA/($CDS_length*$total_totalRNA_aligned_sum))*1000000000         } || 0; }
        if ($stop_codon_count_totalRNA    > 0){ $stop_codon_totalRNA_FPKM=eval{    ($stop_codon_count_totalRNA/($total_totalRNA_aligned_sum))*1000000                 } || 0; }
        if ($trailer_count_totalRNA       > 0){ $trailer_totalRNA_FPKM=eval{       ($trailer_count_totalRNA/($trailer_length*$total_totalRNA_aligned_sum))*1000000000 } || 0; }
#        if ($TTS_ending_count_totalRNA    > 0){ $TTS_ending_totalRNA_FPKM=eval{    ($TTS_ending_count_totalRNA/($total_totalRNA_aligned_sum))*1000000                 } || 0; }

        my $TSS_beginning_TCPseq_SSU_FPKM=0;
	my $leader_TCPseq_SSU_FPKM=0;
	my $start_codon_TCPseq_SSU_FPKM=0;
        my $CDS_TCPseq_SSU_FPKM=0;
        my $stop_codon_TCPseq_SSU_FPKM=0;
        my $trailer_TCPseq_SSU_FPKM=0;
#        my $TTS_ending_TCPseq_SSU_FPKM=0;

        if ($TSS_beginning_count_TCPseq_SSU > 0){ $TSS_beginning_TCPseq_SSU_FPKM=eval{ ($TSS_beginning_count_TCPseq_SSU/($total_TCPseq_SSU_aligned_sum))*1000000              } || 0; }
        if ($leader_count_TCPseq_SSU        > 0){ $leader_TCPseq_SSU_FPKM=eval{        ($leader_count_TCPseq_SSU/($leader_length*$total_TCPseq_SSU_aligned_sum))*1000000000   } || 0; }
        if ($start_codon_count_TCPseq_SSU   > 0){ $start_codon_TCPseq_SSU_FPKM=eval{   ($start_codon_count_TCPseq_SSU/($total_TCPseq_SSU_aligned_sum))*1000000                } || 0; }
        if ($CDS_count_TCPseq_SSU           > 0){ $CDS_TCPseq_SSU_FPKM=eval{           ($CDS_count_TCPseq_SSU/($CDS_length*$total_TCPseq_SSU_aligned_sum))*1000000000         } || 0; }
        if ($stop_codon_count_TCPseq_SSU    > 0){ $stop_codon_TCPseq_SSU_FPKM=eval{    ($stop_codon_count_TCPseq_SSU/($total_TCPseq_SSU_aligned_sum))*1000000                 } || 0; }
        if ($trailer_count_TCPseq_SSU       > 0){ $trailer_TCPseq_SSU_FPKM=eval{       ($trailer_count_TCPseq_SSU/($trailer_length*$total_TCPseq_SSU_aligned_sum))*1000000000 } || 0; }
#        if ($TTS_ending_count_TCPseq_SSU    > 0){ $TTS_ending_TCPseq_SSU_FPKM=eval{    ($TTS_ending_count_TCPseq_SSU/($total_TCPseq_SSU_aligned_sum))*1000000                 } || 0; }

        my $TSS_beginning_TCPseq_LSU_FPKM=0;
	my $leader_TCPseq_LSU_FPKM=0;
	my $start_codon_TCPseq_LSU_FPKM=0;
        my $CDS_TCPseq_LSU_FPKM=0;
        my $stop_codon_TCPseq_LSU_FPKM=0;
        my $trailer_TCPseq_LSU_FPKM=0;
#        my $TTS_ending_TCPseq_LSU_FPKM=0;

        if ($TSS_beginning_count_TCPseq_LSU > 0){ $TSS_beginning_TCPseq_LSU_FPKM=eval{ ($TSS_beginning_count_TCPseq_LSU/($total_TCPseq_LSU_aligned_sum))*1000000              } || 0; }
        if ($leader_count_TCPseq_LSU        > 0){ $leader_TCPseq_LSU_FPKM=eval{        ($leader_count_TCPseq_LSU/($leader_length*$total_TCPseq_LSU_aligned_sum))*1000000000   } || 0; }
        if ($start_codon_count_TCPseq_LSU   > 0){ $start_codon_TCPseq_LSU_FPKM=eval{   ($start_codon_count_TCPseq_LSU/($total_TCPseq_LSU_aligned_sum))*1000000                } || 0; }
        if ($CDS_count_TCPseq_LSU           > 0){ $CDS_TCPseq_LSU_FPKM=eval{           ($CDS_count_TCPseq_LSU/($CDS_length*$total_TCPseq_LSU_aligned_sum))*1000000000         } || 0; }
        if ($stop_codon_count_TCPseq_LSU    > 0){ $stop_codon_TCPseq_LSU_FPKM=eval{    ($stop_codon_count_TCPseq_LSU/($total_TCPseq_LSU_aligned_sum))*1000000                 } || 0; }
        if ($trailer_count_TCPseq_LSU       > 0){ $trailer_TCPseq_LSU_FPKM=eval{       ($trailer_count_TCPseq_LSU/($trailer_length*$total_TCPseq_LSU_aligned_sum))*1000000000 } || 0; }
#        if ($TTS_ending_count_TCPseq_LSU    > 0){ $TTS_ending_TCPseq_LSU_FPKM=eval{    ($TTS_ending_count_TCPseq_LSU/($total_TCPseq_LSU_aligned_sum))*1000000                 } || 0; }

        my $complete_leader_RIBOseq_FPKM=0;
        my $complete_leader_RNAseq_FPKM=0;
        my $complete_leader_totalRNA_FPKM=0;
        my $complete_leader_TCPseq_SSU_FPKM=0;
        my $complete_leader_TCPseq_LSU_FPKM=0;

        if ($leader_region_count_RIBOseq     > 0){ $complete_leader_RIBOseq_FPKM=eval{     ($leader_region_count_RIBOseq    / ($leader_length*$total_RIBOseq_aligned_sum))    *1000000000   } || 0; }
        if ($leader_region_count_RNAseq      > 0){ $complete_leader_RNAseq_FPKM=eval{      ($leader_region_count_RNAseq     / ($leader_length*$total_RNAseq_aligned_sum))     *1000000000   } || 0; }
        if ($leader_region_count_totalRNA    > 0){ $complete_leader_totalRNA_FPKM=eval{    ($leader_region_count_totalRNA   / ($leader_length*$total_totalRNA_aligned_sum))   *1000000000   } || 0; }
        if ($leader_region_count_TCPseq_SSU  > 0){ $complete_leader_TCPseq_SSU_FPKM=eval{  ($leader_region_count_TCPseq_SSU / ($leader_length*$total_TCPseq_SSU_aligned_sum)) *1000000000   } || 0; }
        if ($leader_region_count_TCPseq_LSU  > 0){ $complete_leader_TCPseq_LSU_FPKM=eval{  ($leader_region_count_TCPseq_LSU / ($leader_length*$total_TCPseq_LSU_aligned_sum)) *1000000000   } || 0; }

        my $complete_CDS_RIBOseq_FPKM=0;
        my $complete_CDS_RNAseq_FPKM=0;
        my $complete_CDS_totalRNA_FPKM=0;
        my $complete_CDS_TCPseq_SSU_FPKM=0;
        my $complete_CDS_TCPseq_LSU_FPKM=0;

        if ($CDS_region_count_RIBOseq     > 0){ $complete_CDS_RIBOseq_FPKM=eval{     ($CDS_region_count_RIBOseq    / ($CDS_length*$total_RIBOseq_aligned_sum))    *1000000000   } || 0; }
        if ($CDS_region_count_RNAseq      > 0){ $complete_CDS_RNAseq_FPKM=eval{      ($CDS_region_count_RNAseq     / ($CDS_length*$total_RNAseq_aligned_sum))     *1000000000   } || 0; }
        if ($CDS_region_count_totalRNA    > 0){ $complete_CDS_totalRNA_FPKM=eval{    ($CDS_region_count_totalRNA   / ($CDS_length*$total_totalRNA_aligned_sum))   *1000000000   } || 0; }
        if ($CDS_region_count_TCPseq_SSU  > 0){ $complete_CDS_TCPseq_SSU_FPKM=eval{  ($CDS_region_count_TCPseq_SSU / ($CDS_length*$total_TCPseq_SSU_aligned_sum)) *1000000000   } || 0; }
        if ($CDS_region_count_TCPseq_LSU  > 0){ $complete_CDS_TCPseq_LSU_FPKM=eval{  ($CDS_region_count_TCPseq_LSU / ($CDS_length*$total_TCPseq_LSU_aligned_sum)) *1000000000   } || 0; }

        my $complete_trailer_RIBOseq_FPKM=0;
        my $complete_trailer_RNAseq_FPKM=0;
        my $complete_trailer_totalRNA_FPKM=0;
        my $complete_trailer_TCPseq_SSU_FPKM=0;
        my $complete_trailer_TCPseq_LSU_FPKM=0;

        if ($trailer_region_count_RIBOseq     > 0){ $complete_trailer_RIBOseq_FPKM=eval{     ($trailer_region_count_RIBOseq    / ($trailer_length*$total_RIBOseq_aligned_sum))    *1000000000   } || 0; }
        if ($trailer_region_count_RNAseq      > 0){ $complete_trailer_RNAseq_FPKM=eval{      ($trailer_region_count_RNAseq     / ($trailer_length*$total_RNAseq_aligned_sum))     *1000000000   } || 0; }
        if ($trailer_region_count_totalRNA    > 0){ $complete_trailer_totalRNA_FPKM=eval{    ($trailer_region_count_totalRNA   / ($trailer_length*$total_totalRNA_aligned_sum))   *1000000000   } || 0; }
        if ($trailer_region_count_TCPseq_SSU  > 0){ $complete_trailer_TCPseq_SSU_FPKM=eval{  ($trailer_region_count_TCPseq_SSU / ($trailer_length*$total_TCPseq_SSU_aligned_sum)) *1000000000   } || 0; }
        if ($trailer_region_count_TCPseq_LSU  > 0){ $complete_trailer_TCPseq_LSU_FPKM=eval{  ($trailer_region_count_TCPseq_LSU / ($trailer_length*$total_TCPseq_LSU_aligned_sum)) *1000000000   } || 0; }


        my $tss_proximal_region_RIBOseq_FPKM=0;
        my $tss_proximal_region_RNAseq_FPKM=0;
        my $tss_proximal_region_totalRNA_FPKM=0;
        my $tss_proximal_region_TCPseq_SSU_FPKM=0;
        my $tss_proximal_region_TCPseq_LSU_FPKM=0;

        if ($TSS_proximal_region_count_RIBOseq     > 0){ $tss_proximal_region_RIBOseq_FPKM=eval{     ($TSS_proximal_region_count_RIBOseq    / (50*$total_RIBOseq_aligned_sum))    *1000000000   } || 0; }
        if ($TSS_proximal_region_count_RNAseq      > 0){ $tss_proximal_region_RNAseq_FPKM=eval{      ($TSS_proximal_region_count_RNAseq     / (50*$total_RNAseq_aligned_sum))     *1000000000   } || 0; }
        if ($TSS_proximal_region_count_totalRNA    > 0){ $tss_proximal_region_totalRNA_FPKM=eval{    ($TSS_proximal_region_count_totalRNA   / (50*$total_totalRNA_aligned_sum))   *1000000000   } || 0; }
        if ($TSS_proximal_region_count_TCPseq_SSU  > 0){ $tss_proximal_region_TCPseq_SSU_FPKM=eval{  ($TSS_proximal_region_count_TCPseq_SSU / (50*$total_TCPseq_SSU_aligned_sum)) *1000000000   } || 0; }
        if ($TSS_proximal_region_count_TCPseq_LSU  > 0){ $tss_proximal_region_TCPseq_LSU_FPKM=eval{  ($TSS_proximal_region_count_TCPseq_LSU / (50*$total_TCPseq_LSU_aligned_sum)) *1000000000   } || 0; }


        my $tis_proximal_region_RIBOseq_FPKM=0;
        my $tis_proximal_region_RNAseq_FPKM=0;
        my $tis_proximal_region_totalRNA_FPKM=0;
        my $tis_proximal_region_TCPseq_SSU_FPKM=0;
        my $tis_proximal_region_TCPseq_LSU_FPKM=0;
        
        if ($TIS_proximal_region_count_RIBOseq     > 0){ $tis_proximal_region_RIBOseq_FPKM=eval{     ($TIS_proximal_region_count_RIBOseq    / (50*$total_RIBOseq_aligned_sum))    *1000000000   } || 0; }
        if ($TIS_proximal_region_count_RNAseq      > 0){ $tis_proximal_region_RNAseq_FPKM=eval{      ($TIS_proximal_region_count_RNAseq     / (50*$total_RNAseq_aligned_sum))     *1000000000   } || 0; }
        if ($TIS_proximal_region_count_totalRNA    > 0){ $tis_proximal_region_totalRNA_FPKM=eval{    ($TIS_proximal_region_count_totalRNA   / (50*$total_totalRNA_aligned_sum))   *1000000000   } || 0; }
        if ($TIS_proximal_region_count_TCPseq_SSU  > 0){ $tis_proximal_region_TCPseq_SSU_FPKM=eval{  ($TIS_proximal_region_count_TCPseq_SSU / (50*$total_TCPseq_SSU_aligned_sum)) *1000000000   } || 0; }
        if ($TIS_proximal_region_count_TCPseq_LSU  > 0){ $tis_proximal_region_TCPseq_LSU_FPKM=eval{  ($TIS_proximal_region_count_TCPseq_LSU / (50*$total_TCPseq_LSU_aligned_sum)) *1000000000   } || 0; }

        #find the downstream, TIS and upstream sequences, if there is no anotated leader, the sequences will be called "N"
        my @TIS;
        my @up;
        my @down;
        my @stop; 
        my @pre_stop;
        my @post_stop;
        my @TIS_window;
        my @top;
        my $pre_pre_TSS_start_nt="N";
        my $pre_TSS_start_nt="N";
        my $TSS_start_nt="N";
        my $post_TSS_start_nt="N";

        if (exists ($gene_model_fwd{$gene}{($TSS_coord-1)})){
            $pre_TSS_start_nt=substr($fasta_sequences{$chr}, ($gene_model_fwd{$gene}{($TSS_coord-1)} -1), 1);
        }else{
            #switch to genomic coords
            my $genomic_coord=$gene_model_fwd{$gene}{$TSS_coord};
            $pre_TSS_start_nt=substr($fasta_sequences{$chr}, ($genomic_coord-2), 1);
        }

        if (exists ($gene_model_fwd{$gene}{($TSS_coord-2)})){
            $pre_pre_TSS_start_nt=substr($fasta_sequences{$chr}, ($gene_model_fwd{$gene}{($TSS_coord-2)} -1), 1);
        }else{ #switch to genomic coords
            my $genomic_coord=$gene_model_fwd{$gene}{$TSS_coord};
            $pre_pre_TSS_start_nt=substr($fasta_sequences{$chr}, ($genomic_coord-3), 1);
        }
            
        if (exists ($gene_model_fwd{$gene}{$TSS_coord})){
            $TSS_start_nt=substr($fasta_sequences{$chr}, ($gene_model_fwd{$gene}{$TSS_coord} -1), 1);
        }

        if (exists ($gene_model_fwd{$gene}{($TSS_coord+1)})){
            $post_TSS_start_nt=substr($fasta_sequences{$chr}, ($gene_model_fwd{$gene}{($TSS_coord+1)} -1), 1);
        }

        for (0 .. 4){   #top gene. starts with a C + at least 4 T or C
            if (exists ($gene_model_fwd{$gene}{ ($TSS_coord+$_) } )){
                push (@top, substr($fasta_sequences{$chr}, ($gene_model_fwd{$gene}{ ($TSS_coord+$_) } -1), 1) );
            }else{
                push (@top, "N");
            }
        }

        for(-2 .. -1){
            if (exists ($gene_model_fwd{$gene}{ ($stop_coord+$_) } )){
                push (@pre_stop, substr($fasta_sequences{$chr}, ($gene_model_fwd{$gene}{ ($stop_coord+$_) } -1), 1) );
            }else{
                push (@pre_stop, "N");
            }
        }

        for(0 .. 2){
            if (exists ($gene_model_fwd{$gene}{ ($stop_coord+$_) } )){
                push (@stop, substr($fasta_sequences{$chr}, ($gene_model_fwd{$gene}{ ($stop_coord+$_) } -1), 1) );
            }else{
                push (@stop, "N");
            }
        }

        for(3 .. 4){
            if (exists ($gene_model_fwd{$gene}{ ($stop_coord+$_) } )){
                push (@post_stop, substr($fasta_sequences{$chr}, ($gene_model_fwd{$gene}{ ($stop_coord+$_) } -1), 1) );
            }else{
                push (@post_stop, "N");
            }
        }

        for(0 .. 2){
            if (exists ($gene_model_fwd{$gene}{ ($start_coord+$_) } )){
                push (@TIS, substr($fasta_sequences{$chr}, ($gene_model_fwd{$gene}{ ($start_coord+$_) } -1), 1) );
            }else{
                push (@TIS, "N");
            }
        }

        #further upstream sequence needs to be added in a seperate loop and then concaternated (becuase the PWM sequences neecds to start at -4) 
        for (-4 .. -1){
            if (exists ($gene_model_fwd{$gene}{ ($start_coord+$_) } )){

                if (exists ($five_prime_most_coord_fwd{$gene})){ #for genes with CAGE defined (updated) leaders
                    if (($start_coord+$_) >= ($five_prime_most_coord_fwd{$gene}-1)){ #check we're not going behyond the CAGE defined leader
                        push (@up, substr($fasta_sequences{$chr}, ($gene_model_fwd{$gene}{ ($start_coord+$_) } -1), 1) );
                    }else{
                        push (@up, "N");
                    } 
                }
            }else{
                push (@up, "N");
            }
            
        }

        #this can be freely extended in the 3' direction
        for (3 .. 4){
            if (exists ($gene_model_fwd{$gene}{ ($start_coord+$_) } )){
                push (@down, substr($fasta_sequences{$chr}, ($gene_model_fwd{$gene}{ ($start_coord+$_) } -1), 1) );
            }else{
                push (@down, "N");
            }
        }

        #score the kozak context of the gene
        #if a gene has no upstream sequence an "N" will be assigned with the mean PWM score of A,C,G,T  for that position

        for(-4 .. 7){  #4 downstream, 5 after TIS
            if (exists ($gene_model_fwd{$gene}{ ($start_coord+$_) } )){
                if (exists ($five_prime_most_coord_fwd{$gene})){ #for genes with CAGE defined (updated) leaders
                    if (($start_coord+$_) >= ($five_prime_most_coord_fwd{$gene}-1)){ #check we're not going behyond the CAGE defined leader
                        push (@TIS_window, substr($fasta_sequences{$chr}, ($gene_model_fwd{$gene}{ ($start_coord+$_) } -1), 1) );
                    }else{
                        push (@TIS_window, "N");
                    }
                }
            }else{
                push (@TIS_window, "N");
            }
        }

        my $seq_up=join("", @up);
        my $seq_TIS=join("", @TIS);
        my $seq_down=join("", @down);
        my $seq_stop=join("", @stop);
        my $seq_pre_stop=join("", @pre_stop);
        my $seq_post_stop=join("", @post_stop);
        my $seq_TIS_window=join("", @TIS_window);
        my $seq_top=join("", @top);
        my $seq_tct=$pre_pre_TSS_start_nt.$pre_TSS_start_nt.$seq_top; #this contains one extra

        my $kozak_strength=&score_kozak($seq_up,$seq_down); 
        my $tisu_strength=&score_TISU($seq_TIS_window);
        my $upstream_kozak_strength=&score_kozak_upstream($seq_up);

        if ($seq_top =~ /^C[TC]{4}$/){   #top gene. starts with a C + at least 4 T or C
            $is_TOP="TRUE";
        }

        if ($seq_tct =~ /^[TC]CTTT[TC]\w$/){  # TCT core motif (YCTTTY) Y=C/T
            $is_TCT="TRUE";
        } 

        my $cert=0; 
        if (exists ($cert_count{$gene})){ 
            $cert=$cert_count{$gene}; 
        }

        my $cert_short=0; 
        if (exists ($cert_short_count{$gene})){ 
            $cert_short=$cert_short_count{$gene}; 
        }

        my $polyU=0; 
        if (exists ($polyU_count{$gene})){ 
           $polyU=$polyU_count{$gene}; 
        }

        print "$gene,$transcript,$chr,$highest_cage_peak_value,$leader_length,$CDS_length,$trailer_length,$TSS_beginning_RNAseq_FPKM,$TSS_beginning_totalRNA_FPKM,$TSS_beginning_RIBOseq_FPKM,$TSS_beginning_TCPseq_SSU_FPKM,$TSS_beginning_TCPseq_LSU_FPKM,$leader_RNAseq_FPKM,$leader_totalRNA_FPKM,$leader_RIBOseq_FPKM,$leader_TCPseq_SSU_FPKM,$leader_TCPseq_LSU_FPKM,$start_codon_RNAseq_FPKM,$start_codon_totalRNA_FPKM,$start_codon_RIBOseq_FPKM,$start_codon_TCPseq_SSU_FPKM,$start_codon_TCPseq_LSU_FPKM,$CDS_RNAseq_FPKM,$CDS_totalRNA_FPKM,$CDS_RIBOseq_FPKM,$CDS_TCPseq_SSU_FPKM,$CDS_TCPseq_LSU_FPKM,$stop_codon_RNAseq_FPKM,$stop_codon_totalRNA_FPKM,$stop_codon_RIBOseq_FPKM,$stop_codon_TCPseq_SSU_FPKM,$stop_codon_TCPseq_LSU_FPKM,$trailer_RNAseq_FPKM,$trailer_totalRNA_FPKM,$trailer_RIBOseq_FPKM,$trailer_TCPseq_SSU_FPKM,$trailer_TCPseq_LSU_FPKM,$complete_leader_RNAseq_FPKM,$complete_leader_totalRNA_FPKM,$complete_leader_RIBOseq_FPKM,$complete_leader_TCPseq_SSU_FPKM,$complete_leader_TCPseq_LSU_FPKM,$complete_CDS_RNAseq_FPKM,$complete_CDS_totalRNA_FPKM,$complete_CDS_RIBOseq_FPKM,$complete_CDS_TCPseq_SSU_FPKM,$complete_CDS_TCPseq_LSU_FPKM,$complete_trailer_RNAseq_FPKM,$complete_trailer_totalRNA_FPKM,$complete_trailer_RIBOseq_FPKM,$complete_trailer_TCPseq_SSU_FPKM,$complete_trailer_TCPseq_LSU_FPKM,$tss_proximal_region_RNAseq_FPKM,$tss_proximal_region_totalRNA_FPKM,$tss_proximal_region_RIBOseq_FPKM,$tss_proximal_region_TCPseq_SSU_FPKM,$tss_proximal_region_TCPseq_LSU_FPKM,$tis_proximal_region_RNAseq_FPKM,$tis_proximal_region_totalRNA_FPKM,$tis_proximal_region_RIBOseq_FPKM,$tis_proximal_region_TCPseq_SSU_FPKM,$tis_proximal_region_TCPseq_LSU_FPKM,$overlapping_gene,$leader_potentially_overlaps_upstream_gene,$gene_potentially_overlaps_downstream_leader,$gene_overlaps_nc,$pre_TSS_start_nt,$TSS_start_nt,$post_TSS_start_nt,$seq_up,$seq_TIS,$seq_down,$seq_pre_stop,$seq_stop,$seq_post_stop,$exon_exon_junction_count_leader,$kozak_strength,$upstream_kozak_strength,$tisu_strength,$is_histone,$is_TOP,$is_TCT,$cert,$cert_short,$polyU\n";
    }
}

for my $gene (keys %gene_model_rev){ 

    if (exists ($gene_2_chr{$gene})){  #restict to protien_coding genes
        my $chr=$gene_2_chr{$gene};
        my $transcript=$most_expressed_transcript{$gene};
        my $start_coord=$start_coord_rev{$gene};
        my $stop_coord=$stop_coord_rev{$gene};
        my $CDS_length=(($stop_coord+2)-$start_coord)+1;
        my $trailer_length=($three_prime_most_coord_rev{$gene}-($stop_coord+3))+1;
        my $leader_start=$leader_positions_rev{$gene};
        my $TSS_coord=$five_prime_most_coord_rev{$gene};
        my $leader_length=$start_coord-$TSS_coord;
        my $highest_cage_peak_value=0;
        my $overlapping_gene="FALSE";
        my $leader_potentially_overlaps_upstream_gene="FALSE";
        my $gene_potentially_overlaps_downstream_leader="FALSE";
        my $gene_overlaps_nc="FALSE";
        my $is_histone="FALSE";
        my $exon_exon_junction_count_leader=$gene_exon_exon_in_leader{$gene};
        my $is_TOP="FALSE";
        my $is_TCT="FALSE";

        if (exists ($cage_peak_value{$gene})){ $highest_cage_peak_value=$cage_peak_value{$gene}; }
        if (exists ($overlaps_inframe_gene{$gene})){           $overlapping_gene="TRUE"; }
        if (exists ($leader_overlaps_upstream{$gene})){        $leader_potentially_overlaps_upstream_gene="TRUE"; }
        if (exists ($gene_overlaps_downstream_leader{$gene})){ $gene_potentially_overlaps_downstream_leader="TRUE"; }
        if (exists ($gene_overlaps_nc{$gene})){                $gene_overlaps_nc="TRUE"; }
        if (exists ($histone{$gene})){                         $is_histone="TRUE"; }

        my $TSS_beginning_count_RIBOseq=$TSS_beginning_count_RIBOseq{$gene};
        my $leader_count_RIBOseq=$leader_counts_RIBOseq{$gene};
        my $start_codon_count_RIBOseq=$start_codon_counts_RIBOseq{$gene};
        my $CDS_count_RIBOseq=$CDS_counts_RIBOseq{$gene};
        my $stop_codon_count_RIBOseq=$stop_codon_counts_RIBOseq{$gene};
        my $trailer_count_RIBOseq=$trailer_counts_RIBOseq{$gene};
#        my $TTS_ending_count_RIBOseq=$TTS_ending_count_RIBOseq{$gene};

        my $TSS_beginning_count_RNAseq=$TSS_beginning_count_RNAseq{$gene};
        my $leader_count_RNAseq=$leader_counts_RNAseq{$gene};
        my $start_codon_count_RNAseq=$start_codon_counts_RNAseq{$gene};
        my $CDS_count_RNAseq=$CDS_counts_RNAseq{$gene};
        my $stop_codon_count_RNAseq=$stop_codon_counts_RNAseq{$gene};
        my $trailer_count_RNAseq=$trailer_counts_RNAseq{$gene};
#        my $TTS_ending_count_RNAseq=$TTS_ending_count_RNAseq{$gene};

        my $TSS_beginning_count_totalRNA=$TSS_beginning_count_totalRNA{$gene};
        my $leader_count_totalRNA=$leader_counts_totalRNA{$gene};
        my $start_codon_count_totalRNA=$start_codon_counts_totalRNA{$gene};
        my $CDS_count_totalRNA=$CDS_counts_totalRNA{$gene};
        my $stop_codon_count_totalRNA=$stop_codon_counts_totalRNA{$gene};
        my $trailer_count_totalRNA=$trailer_counts_totalRNA{$gene};
#        my $TTS_ending_count_totalRNA=$TTS_ending_count_totalRNA{$gene};

        my $TSS_beginning_count_TCPseq_SSU=$TSS_beginning_count_TCPseq_SSU{$gene};
        my $leader_count_TCPseq_SSU=$leader_counts_TCPseq_SSU{$gene};
        my $start_codon_count_TCPseq_SSU=$start_codon_counts_TCPseq_SSU{$gene};
        my $CDS_count_TCPseq_SSU=$CDS_counts_TCPseq_SSU{$gene};
        my $stop_codon_count_TCPseq_SSU=$stop_codon_counts_TCPseq_SSU{$gene};
        my $trailer_count_TCPseq_SSU=$trailer_counts_TCPseq_SSU{$gene};
#        my $TTS_ending_count_TCPseq_SSU=$TTS_ending_count_TCPseq_SSU{$gene};

        my $TSS_beginning_count_TCPseq_LSU=$TSS_beginning_count_TCPseq_LSU{$gene};
        my $leader_count_TCPseq_LSU=$leader_counts_TCPseq_LSU{$gene};
        my $start_codon_count_TCPseq_LSU=$start_codon_counts_TCPseq_LSU{$gene};
        my $CDS_count_TCPseq_LSU=$CDS_counts_TCPseq_LSU{$gene};
        my $stop_codon_count_TCPseq_LSU=$stop_codon_counts_TCPseq_LSU{$gene};
        my $trailer_count_TCPseq_LSU=$trailer_counts_TCPseq_LSU{$gene};
#        my $TTS_ending_count_TCPseq_LSU=$TTS_ending_count_TCPseq_LSU{$gene};

        my $leader_region_count_TCPseq_SSU=$leader_region_count_TCPseq_SSU{$gene};
        my $leader_region_count_TCPseq_LSU=$leader_region_count_TCPseq_LSU{$gene};
        my $leader_region_count_RNAseq=$leader_region_count_RNAseq{$gene};
        my $leader_region_count_totalRNA=$leader_region_count_totalRNA{$gene};
        my $leader_region_count_RIBOseq=$leader_region_count_RIBOseq{$gene};

        my $CDS_region_count_TCPseq_SSU=$CDS_region_count_TCPseq_SSU{$gene};
        my $CDS_region_count_TCPseq_LSU=$CDS_region_count_TCPseq_LSU{$gene};
        my $CDS_region_count_RNAseq=$CDS_region_count_RNAseq{$gene};
        my $CDS_region_count_totalRNA=$CDS_region_count_totalRNA{$gene};      
        my $CDS_region_count_RIBOseq=$CDS_region_count_RIBOseq{$gene};

        my $trailer_region_count_TCPseq_SSU=$trailer_region_count_TCPseq_SSU{$gene};
        my $trailer_region_count_TCPseq_LSU=$trailer_region_count_TCPseq_LSU{$gene};
        my $trailer_region_count_RNAseq=$trailer_region_count_RNAseq{$gene};
        my $trailer_region_count_totalRNA=$trailer_region_count_totalRNA{$gene};
        my $trailer_region_count_RIBOseq=$trailer_region_count_RIBOseq{$gene};

        my $TSS_proximal_region_count_TCPseq_SSU=$TSS_proximal_region_count_TCPseq_SSU{$gene};
        my $TSS_proximal_region_count_TCPseq_LSU=$TSS_proximal_region_count_TCPseq_LSU{$gene};
        my $TSS_proximal_region_count_RNAseq=$TSS_proximal_region_count_RNAseq{$gene};
        my $TSS_proximal_region_count_totalRNA=$TSS_proximal_region_count_totalRNA{$gene};
        my $TSS_proximal_region_count_RIBOseq=$TSS_proximal_region_count_RIBOseq{$gene};

        my $TIS_proximal_region_count_TCPseq_SSU=$TIS_proximal_region_count_TCPseq_SSU{$gene};
        my $TIS_proximal_region_count_TCPseq_LSU=$TIS_proximal_region_count_TCPseq_LSU{$gene};
        my $TIS_proximal_region_count_RNAseq=$TIS_proximal_region_count_RNAseq{$gene};
        my $TIS_proximal_region_count_totalRNA=$TIS_proximal_region_count_totalRNA{$gene};
        my $TIS_proximal_region_count_RIBOseq=$TIS_proximal_region_count_RIBOseq{$gene};

        my $TSS_beginning_RIBOseq_FPKM=0;
	my $leader_RIBOseq_FPKM=0;
	my $start_codon_RIBOseq_FPKM=0;
        my $CDS_RIBOseq_FPKM=0;
        my $stop_codon_RIBOseq_FPKM=0;
        my $trailer_RIBOseq_FPKM=0;
#        my $TTS_ending_RIBOseq_FPKM=0;

        if ($TSS_beginning_count_RIBOseq > 0){ $TSS_beginning_RIBOseq_FPKM=eval{ ($TSS_beginning_count_RIBOseq/($total_RIBOseq_aligned_sum))*1000000              } || 0; }
        if ($leader_count_RIBOseq        > 0){ $leader_RIBOseq_FPKM=eval{        ($leader_count_RIBOseq/($leader_length*$total_RIBOseq_aligned_sum))*1000000000   } || 0; }
        if ($start_codon_count_RIBOseq   > 0){ $start_codon_RIBOseq_FPKM=eval{   ($start_codon_count_RIBOseq/($total_RIBOseq_aligned_sum))*1000000                } || 0; }
        if ($CDS_count_RIBOseq           > 0){ $CDS_RIBOseq_FPKM=eval{           ($CDS_count_RIBOseq/($CDS_length*$total_RIBOseq_aligned_sum))*1000000000         } || 0; }
        if ($stop_codon_count_RIBOseq    > 0){ $stop_codon_RIBOseq_FPKM=eval{    ($stop_codon_count_RIBOseq/($total_RIBOseq_aligned_sum))*1000000                 } || 0; }
        if ($trailer_count_RIBOseq       > 0){ $trailer_RIBOseq_FPKM=eval{       ($trailer_count_RIBOseq/($trailer_length*$total_RIBOseq_aligned_sum))*1000000000 } || 0; }
#        if ($TTS_ending_count_RIBOseq    > 0){ $TTS_ending_RIBOseq_FPKM=eval{    ($TTS_ending_count_RIBOseq/($total_RIBOseq_aligned_sum))*1000000                 } || 0; }

        my $TSS_beginning_RNAseq_FPKM=0;
	my $leader_RNAseq_FPKM=0;
	my $start_codon_RNAseq_FPKM=0;
        my $CDS_RNAseq_FPKM=0;
        my $stop_codon_RNAseq_FPKM=0;
        my $trailer_RNAseq_FPKM=0;
#        my $TTS_ending_RNAseq_FPKM=0;

        if ($TSS_beginning_count_RNAseq > 0){ $TSS_beginning_RNAseq_FPKM=eval{ ($TSS_beginning_count_RNAseq/($total_RNAseq_aligned_sum))*1000000              } || 0; }
        if ($leader_count_RNAseq        > 0){ $leader_RNAseq_FPKM=eval{        ($leader_count_RNAseq/($leader_length*$total_RNAseq_aligned_sum))*1000000000   } || 0; }
        if ($start_codon_count_RNAseq   > 0){ $start_codon_RNAseq_FPKM=eval{   ($start_codon_count_RNAseq/($total_RNAseq_aligned_sum))*1000000                } || 0; }
        if ($CDS_count_RNAseq           > 0){ $CDS_RNAseq_FPKM=eval{           ($CDS_count_RNAseq/($CDS_length*$total_RNAseq_aligned_sum))*1000000000         } || 0; }
        if ($stop_codon_count_RNAseq    > 0){ $stop_codon_RNAseq_FPKM=eval{    ($stop_codon_count_RNAseq/($total_RNAseq_aligned_sum))*1000000                 } || 0; }
        if ($trailer_count_RNAseq       > 0){ $trailer_RNAseq_FPKM=eval{       ($trailer_count_RNAseq/($trailer_length*$total_RNAseq_aligned_sum))*1000000000 } || 0; }
#        if ($TTS_ending_count_RNAseq    > 0){ $TTS_ending_RNAseq_FPKM=eval{    ($TTS_ending_count_RNAseq/($total_RNAseq_aligned_sum))*1000000                 } || 0; }

        my $TSS_beginning_totalRNA_FPKM=0;
        my $leader_totalRNA_FPKM=0;
        my $start_codon_totalRNA_FPKM=0;
        my $CDS_totalRNA_FPKM=0;
        my $stop_codon_totalRNA_FPKM=0;
        my $trailer_totalRNA_FPKM=0;
#        my $TTS_ending_totalRNA_FPKM=0;

        if ($TSS_beginning_count_totalRNA > 0){ $TSS_beginning_totalRNA_FPKM=eval{ ($TSS_beginning_count_totalRNA/($total_totalRNA_aligned_sum))*1000000              } || 0; }
        if ($leader_count_totalRNA        > 0){ $leader_totalRNA_FPKM=eval{        ($leader_count_totalRNA/($leader_length*$total_totalRNA_aligned_sum))*1000000000   } || 0; }
        if ($start_codon_count_totalRNA   > 0){ $start_codon_totalRNA_FPKM=eval{   ($start_codon_count_totalRNA/($total_totalRNA_aligned_sum))*1000000                } || 0; }
        if ($CDS_count_totalRNA           > 0){ $CDS_totalRNA_FPKM=eval{           ($CDS_count_totalRNA/($CDS_length*$total_totalRNA_aligned_sum))*1000000000         } || 0; }
        if ($stop_codon_count_totalRNA    > 0){ $stop_codon_totalRNA_FPKM=eval{    ($stop_codon_count_totalRNA/($total_totalRNA_aligned_sum))*1000000                 } || 0; }
        if ($trailer_count_totalRNA       > 0){ $trailer_totalRNA_FPKM=eval{       ($trailer_count_totalRNA/($trailer_length*$total_totalRNA_aligned_sum))*1000000000 } || 0; }
#        if ($TTS_ending_count_totalRNA    > 0){ $TTS_ending_totalRNA_FPKM=eval{    ($TTS_ending_count_totalRNA/($total_totalRNA_aligned_sum))*1000000                 } || 0; }

        my $TSS_beginning_TCPseq_SSU_FPKM=0;
	my $leader_TCPseq_SSU_FPKM=0;
	my $start_codon_TCPseq_SSU_FPKM=0;
        my $CDS_TCPseq_SSU_FPKM=0;
        my $stop_codon_TCPseq_SSU_FPKM=0;
        my $trailer_TCPseq_SSU_FPKM=0;
#        my $TTS_ending_TCPseq_SSU_FPKM=0;

        if ($TSS_beginning_count_TCPseq_SSU > 0){ $TSS_beginning_TCPseq_SSU_FPKM=eval{ ($TSS_beginning_count_TCPseq_SSU/($total_TCPseq_SSU_aligned_sum))*1000000              } || 0; }
        if ($leader_count_TCPseq_SSU        > 0){ $leader_TCPseq_SSU_FPKM=eval{        ($leader_count_TCPseq_SSU/($leader_length*$total_TCPseq_SSU_aligned_sum))*1000000000   } || 0; }
        if ($start_codon_count_TCPseq_SSU   > 0){ $start_codon_TCPseq_SSU_FPKM=eval{   ($start_codon_count_TCPseq_SSU/($total_TCPseq_SSU_aligned_sum))*1000000                } || 0; }
        if ($CDS_count_TCPseq_SSU           > 0){ $CDS_TCPseq_SSU_FPKM=eval{           ($CDS_count_TCPseq_SSU/($CDS_length*$total_TCPseq_SSU_aligned_sum))*1000000000         } || 0; }
        if ($stop_codon_count_TCPseq_SSU    > 0){ $stop_codon_TCPseq_SSU_FPKM=eval{    ($stop_codon_count_TCPseq_SSU/($total_TCPseq_SSU_aligned_sum))*1000000                 } || 0; }
        if ($trailer_count_TCPseq_SSU       > 0){ $trailer_TCPseq_SSU_FPKM=eval{       ($trailer_count_TCPseq_SSU/($trailer_length*$total_TCPseq_SSU_aligned_sum))*1000000000 } || 0; }
#        if ($TTS_ending_count_TCPseq_SSU    > 0){ $TTS_ending_TCPseq_SSU_FPKM=eval{    ($TTS_ending_count_TCPseq_SSU/($total_TCPseq_SSU_aligned_sum))*1000000                 } || 0; }

        my $TSS_beginning_TCPseq_LSU_FPKM=0;
	my $leader_TCPseq_LSU_FPKM=0;
	my $start_codon_TCPseq_LSU_FPKM=0;
        my $CDS_TCPseq_LSU_FPKM=0;
        my $stop_codon_TCPseq_LSU_FPKM=0;
        my $trailer_TCPseq_LSU_FPKM=0;
#        my $TTS_ending_TCPseq_LSU_FPKM=0;

        if ($TSS_beginning_count_TCPseq_LSU > 0){ $TSS_beginning_TCPseq_LSU_FPKM=eval{ ($TSS_beginning_count_TCPseq_LSU/($total_TCPseq_LSU_aligned_sum))*1000000              } || 0; }
        if ($leader_count_TCPseq_LSU        > 0){ $leader_TCPseq_LSU_FPKM=eval{        ($leader_count_TCPseq_LSU/($leader_length*$total_TCPseq_LSU_aligned_sum))*1000000000   } || 0; }
        if ($start_codon_count_TCPseq_LSU   > 0){ $start_codon_TCPseq_LSU_FPKM=eval{   ($start_codon_count_TCPseq_LSU/($total_TCPseq_LSU_aligned_sum))*1000000                } || 0; }
        if ($CDS_count_TCPseq_LSU           > 0){ $CDS_TCPseq_LSU_FPKM=eval{           ($CDS_count_TCPseq_LSU/($CDS_length*$total_TCPseq_LSU_aligned_sum))*1000000000         } || 0; }
        if ($stop_codon_count_TCPseq_LSU    > 0){ $stop_codon_TCPseq_LSU_FPKM=eval{    ($stop_codon_count_TCPseq_LSU/($total_TCPseq_LSU_aligned_sum))*1000000                 } || 0; }
        if ($trailer_count_TCPseq_LSU       > 0){ $trailer_TCPseq_LSU_FPKM=eval{       ($trailer_count_TCPseq_LSU/($trailer_length*$total_TCPseq_LSU_aligned_sum))*1000000000 } || 0; }
#        if ($TTS_ending_count_TCPseq_LSU    > 0){ $TTS_ending_TCPseq_LSU_FPKM=eval{    ($TTS_ending_count_TCPseq_LSU/($total_TCPseq_LSU_aligned_sum))*1000000                 } || 0; }

        my $complete_leader_RIBOseq_FPKM=0;
        my $complete_leader_RNAseq_FPKM=0;
        my $complete_leader_totalRNA_FPKM=0;
        my $complete_leader_TCPseq_SSU_FPKM=0;
        my $complete_leader_TCPseq_LSU_FPKM=0;

        if ($leader_region_count_RIBOseq     > 0){ $complete_leader_RIBOseq_FPKM=eval{     ($leader_region_count_RIBOseq    / ($leader_length*$total_RIBOseq_aligned_sum))    *1000000000   } || 0; }
        if ($leader_region_count_RNAseq      > 0){ $complete_leader_RNAseq_FPKM=eval{      ($leader_region_count_RNAseq     / ($leader_length*$total_RNAseq_aligned_sum))     *1000000000   } || 0; }
        if ($leader_region_count_totalRNA    > 0){ $complete_leader_totalRNA_FPKM=eval{    ($leader_region_count_totalRNA   / ($leader_length*$total_totalRNA_aligned_sum))   *1000000000   } || 0; }
        if ($leader_region_count_TCPseq_SSU  > 0){ $complete_leader_TCPseq_SSU_FPKM=eval{  ($leader_region_count_TCPseq_SSU / ($leader_length*$total_TCPseq_SSU_aligned_sum)) *1000000000   } || 0; }
        if ($leader_region_count_TCPseq_LSU  > 0){ $complete_leader_TCPseq_LSU_FPKM=eval{  ($leader_region_count_TCPseq_LSU / ($leader_length*$total_TCPseq_LSU_aligned_sum)) *1000000000   } || 0; }

        my $complete_CDS_RIBOseq_FPKM=0;
        my $complete_CDS_RNAseq_FPKM=0;
        my $complete_CDS_totalRNA_FPKM=0;
        my $complete_CDS_TCPseq_SSU_FPKM=0;
        my $complete_CDS_TCPseq_LSU_FPKM=0;

        if ($CDS_region_count_RIBOseq     > 0){ $complete_CDS_RIBOseq_FPKM=eval{     ($CDS_region_count_RIBOseq    / ($CDS_length*$total_RIBOseq_aligned_sum))    *1000000000   } || 0; }
        if ($CDS_region_count_RNAseq      > 0){ $complete_CDS_RNAseq_FPKM=eval{      ($CDS_region_count_RNAseq     / ($CDS_length*$total_RNAseq_aligned_sum))     *1000000000   } || 0; }
        if ($CDS_region_count_totalRNA    > 0){ $complete_CDS_totalRNA_FPKM=eval{    ($CDS_region_count_totalRNA   / ($CDS_length*$total_totalRNA_aligned_sum))   *1000000000   } || 0; }
        if ($CDS_region_count_TCPseq_SSU  > 0){ $complete_CDS_TCPseq_SSU_FPKM=eval{  ($CDS_region_count_TCPseq_SSU / ($CDS_length*$total_TCPseq_SSU_aligned_sum)) *1000000000   } || 0; }
        if ($CDS_region_count_TCPseq_LSU  > 0){ $complete_CDS_TCPseq_LSU_FPKM=eval{  ($CDS_region_count_TCPseq_LSU / ($CDS_length*$total_TCPseq_LSU_aligned_sum)) *1000000000   } || 0; }

        my $complete_trailer_RIBOseq_FPKM=0;
        my $complete_trailer_RNAseq_FPKM=0;
        my $complete_trailer_totalRNA_FPKM=0;
        my $complete_trailer_TCPseq_SSU_FPKM=0;
        my $complete_trailer_TCPseq_LSU_FPKM=0;

        if ($trailer_region_count_RIBOseq     > 0){ $complete_trailer_RIBOseq_FPKM=eval{     ($trailer_region_count_RIBOseq    / ($trailer_length*$total_RIBOseq_aligned_sum))    *1000000000   } || 0; }
        if ($trailer_region_count_RNAseq      > 0){ $complete_trailer_RNAseq_FPKM=eval{      ($trailer_region_count_RNAseq     / ($trailer_length*$total_RNAseq_aligned_sum))     *1000000000   } || 0; }
        if ($trailer_region_count_totalRNA    > 0){ $complete_trailer_totalRNA_FPKM=eval{    ($trailer_region_count_totalRNA   / ($trailer_length*$total_totalRNA_aligned_sum))   *1000000000   } || 0; }
        if ($trailer_region_count_TCPseq_SSU  > 0){ $complete_trailer_TCPseq_SSU_FPKM=eval{  ($trailer_region_count_TCPseq_SSU / ($trailer_length*$total_TCPseq_SSU_aligned_sum)) *1000000000   } || 0; }
        if ($trailer_region_count_TCPseq_LSU  > 0){ $complete_trailer_TCPseq_LSU_FPKM=eval{  ($trailer_region_count_TCPseq_LSU / ($trailer_length*$total_TCPseq_LSU_aligned_sum)) *1000000000   } || 0; }


        my $tss_proximal_region_RIBOseq_FPKM=0;
        my $tss_proximal_region_RNAseq_FPKM=0;
        my $tss_proximal_region_totalRNA_FPKM=0;
        my $tss_proximal_region_TCPseq_SSU_FPKM=0;
        my $tss_proximal_region_TCPseq_LSU_FPKM=0;

        if ($TSS_proximal_region_count_RIBOseq     > 0){ $tss_proximal_region_RIBOseq_FPKM=eval{     ($TSS_proximal_region_count_RIBOseq    / (50*$total_RIBOseq_aligned_sum))    *1000000000   } || 0; }
        if ($TSS_proximal_region_count_RNAseq      > 0){ $tss_proximal_region_RNAseq_FPKM=eval{      ($TSS_proximal_region_count_RNAseq     / (50*$total_RNAseq_aligned_sum))     *1000000000   } || 0; }
        if ($TSS_proximal_region_count_totalRNA    > 0){ $tss_proximal_region_totalRNA_FPKM=eval{    ($TSS_proximal_region_count_totalRNA   / (50*$total_totalRNA_aligned_sum))   *1000000000   } || 0; }
        if ($TSS_proximal_region_count_TCPseq_SSU  > 0){ $tss_proximal_region_TCPseq_SSU_FPKM=eval{  ($TSS_proximal_region_count_TCPseq_SSU / (50*$total_TCPseq_SSU_aligned_sum)) *1000000000   } || 0; }
        if ($TSS_proximal_region_count_TCPseq_LSU  > 0){ $tss_proximal_region_TCPseq_LSU_FPKM=eval{  ($TSS_proximal_region_count_TCPseq_LSU / (50*$total_TCPseq_LSU_aligned_sum)) *1000000000   } || 0; }


        my $tis_proximal_region_RIBOseq_FPKM=0;
        my $tis_proximal_region_RNAseq_FPKM=0;
        my $tis_proximal_region_totalRNA_FPKM=0;
        my $tis_proximal_region_TCPseq_SSU_FPKM=0;
        my $tis_proximal_region_TCPseq_LSU_FPKM=0;
        
        if ($TIS_proximal_region_count_RIBOseq      > 0){ $tis_proximal_region_RIBOseq_FPKM=eval{      ($TIS_proximal_region_count_RIBOseq    / (50*$total_RIBOseq_aligned_sum))    *1000000000   } || 0; }
        if ($TIS_proximal_region_count_RNAseq      > 0){ $tis_proximal_region_RNAseq_FPKM=eval{      ($TIS_proximal_region_count_RNAseq     / (50*$total_RNAseq_aligned_sum))     *1000000000   } || 0; }
        if ($TIS_proximal_region_count_totalRNA    > 0){ $tis_proximal_region_totalRNA_FPKM=eval{    ($TIS_proximal_region_count_totalRNA   / (50*$total_totalRNA_aligned_sum))   *1000000000   } || 0; }
        if ($TIS_proximal_region_count_TCPseq_SSU  > 0){ $tis_proximal_region_TCPseq_SSU_FPKM=eval{  ($TIS_proximal_region_count_TCPseq_SSU / (50*$total_TCPseq_SSU_aligned_sum)) *1000000000   } || 0; }
        if ($TIS_proximal_region_count_TCPseq_LSU  > 0){ $tis_proximal_region_TCPseq_LSU_FPKM=eval{  ($TIS_proximal_region_count_TCPseq_LSU / (50*$total_TCPseq_LSU_aligned_sum)) *1000000000   } || 0; }


        #Find the downstream, TIS and upstream squences, if there is no anotated leader, the sequences will be called "N"
        my @TIS;
        my @up;
        my @down;
        my @stop;
        my @pre_stop;
        my @post_stop;
        my @TIS_window;
        my @top;
        my $pre_pre_TSS_start_nt="N";
        my $pre_TSS_start_nt="N";
        my $TSS_start_nt="N";
        my $post_TSS_start_nt="N";

        if (exists ($gene_model_rev{$gene}{($TSS_coord-1)})){
            $pre_TSS_start_nt=substr($fasta_sequences{$chr}, ($gene_model_rev{$gene}{($TSS_coord-1)} -1), 1);
        }else{
            #switch to genomic coords
            my $genomic_coord=$gene_model_rev{$gene}{$TSS_coord};
            $pre_TSS_start_nt=substr($fasta_sequences{$chr}, ($genomic_coord), 1); #equivialnt to +1?
        }

        if (exists ($gene_model_fwd{$gene}{($TSS_coord-2)})){
            $pre_pre_TSS_start_nt=substr($fasta_sequences{$chr}, ($gene_model_rev{$gene}{($TSS_coord-2)} -1), 1);
        }else{ #switch to genomic coords
            my $genomic_coord=$gene_model_rev{$gene}{$TSS_coord};
            $pre_pre_TSS_start_nt=substr($fasta_sequences{$chr}, ($genomic_coord+1), 1);
        }

        if (exists ($gene_model_rev{$gene}{$TSS_coord})){
            $TSS_start_nt=substr($fasta_sequences{$chr}, ($gene_model_rev{$gene}{$TSS_coord} -1), 1);
        }

        if (exists ($gene_model_rev{$gene}{($TSS_coord+1)})){
            $post_TSS_start_nt=substr($fasta_sequences{$chr}, ($gene_model_rev{$gene}{($TSS_coord+1)} -1), 1);
        }

        for (0 .. 4){   #top gene. starts with a C + at least 4 T or C
            if (exists ($gene_model_rev{$gene}{ ($TSS_coord+$_) } )){
                push (@top, substr($fasta_sequences{$chr}, ($gene_model_rev{$gene}{ ($TSS_coord+$_) } -1), 1) );
            }else{
                push (@top, "N");
            }
        }

        for(-2 .. -1){
            if (exists ($gene_model_rev{$gene}{ ($stop_coord+$_) } )){
                push (@pre_stop, substr($fasta_sequences{$chr}, ($gene_model_rev{$gene}{ ($stop_coord+$_) } -1), 1) );
            }else{
                push (@pre_stop, "N");
            }
        }

        for(0 .. 2){
            if (exists ($gene_model_rev{$gene}{ ($stop_coord+$_) } )){
                push (@stop, substr($fasta_sequences{$chr}, ($gene_model_rev{$gene}{ ($stop_coord+$_) } -1), 1) );
            }else{
                push (@stop, "N");
            }
        }

        for(3 .. 4){
            if (exists ($gene_model_rev{$gene}{ ($stop_coord+$_) } )){
                push (@post_stop, substr($fasta_sequences{$chr}, ($gene_model_rev{$gene}{ ($stop_coord+$_) } -1), 1) );
            }else{
                push (@post_stop, "N");
            }
        }

        for(0 .. 2){
            if (exists ($gene_model_rev{$gene}{ ($start_coord+$_) } )){
                push (@TIS, substr($fasta_sequences{$chr}, ($gene_model_rev{$gene}{ ($start_coord+$_) } -1), 1) );
            }else{
                push (@TIS, "N");
            }
        }

        #further upstream sequence needs to be added in a seperate loop and then concaternated (becuase the PWM sequences neecds to start at -4) 
        my $upstream_coord=$start_coord;
        for (-4 .. -1){
            if (exists ($gene_model_rev{$gene}{ ($start_coord+$_) } )){
                if (exists ($five_prime_most_coord_rev{$gene})){ #for genes with CAGE defined (updated) leaders
                    if (($start_coord+$_) >= ($five_prime_most_coord_rev{$gene}-1)){ #check we're not going behyond the CAGE defined leader
                        push (@up, substr($fasta_sequences{$chr}, ($gene_model_rev{$gene}{ ($start_coord+$_) } -1), 1) );
                    }else{
                        push (@up, "N");
                    }
                }
            }else{
                push (@up, "N");
            }
        }

        #this can be freely extended in the 3' direction
        for (3 .. 4){
            if (exists ($gene_model_rev{$gene}{ ($start_coord+$_) } )){
                push (@down, substr($fasta_sequences{$chr}, ($gene_model_rev{$gene}{ ($start_coord+$_) } -1), 1) );
            }else{
                push (@down, "N");
            }
        }

        for(-4 .. 7){  #4 downstream, 5 after TIS
            if (exists ($gene_model_rev{$gene}{ ($start_coord+$_) } )){
                if (exists ($five_prime_most_coord_rev{$gene})){ #for genes with CAGE defined (updated) leaders
                    if (($start_coord+$_) >= ($five_prime_most_coord_rev{$gene}-1)){ #check we're not going behyond the CAGE defined leader
                        push (@TIS_window, substr($fasta_sequences{$chr}, ($gene_model_rev{$gene}{ ($start_coord+$_) } -1), 1) );
                    }else{
                        push (@TIS_window, "N");
                    }
                }
            }else{
                push (@TIS_window, "N");
            }
        }

        my $seq_up=join("", @up);
        my $seq_TIS=join("", @TIS);
        my $seq_down=join("", @down);
        my $seq_stop=join("", @stop);
        my $seq_pre_stop=join("", @pre_stop);
        my $seq_post_stop=join("", @post_stop);
        my $seq_TIS_window=join("", @TIS_window);
        my $seq_top=join("", @top);
        my $seq_tct=$pre_pre_TSS_start_nt.$pre_TSS_start_nt.$seq_top; #this contains one extra

        $seq_up=~tr/ACGTacgt/TGCAtgca/;
        $seq_TIS=~tr/ACGTacgt/TGCAtgca/;
        $seq_down=~tr/ACGTacgt/TGCAtgca/;
        $seq_stop=~tr/ACGTacgt/TGCAtgca/;
        $seq_pre_stop=~tr/ACGTacgt/TGCAtgca/;
        $seq_post_stop=~tr/ACGTacgt/TGCAtgca/;
        $pre_pre_TSS_start_nt=~tr/ACGTacgt/TGCAtgca/;
        $pre_TSS_start_nt=~tr/ACGTacgt/TGCAtgca/; 
        $TSS_start_nt=~tr/ACGTacgt/TGCAtgca/;  
        $post_TSS_start_nt=~tr/ACGTacgt/TGCAtgca/;  
        $seq_TIS_window=~tr/ACGTacgt/TGCAtgca/;  
        $seq_top=~tr/ACGTacgt/TGCAtgca/;
        $seq_tct=~tr/ACGTacgt/TGCAtgca/;

        if ($seq_top =~ /^C[TC]{4}$/){ 
            $is_TOP="TRUE";
        }

        if ($seq_tct =~ /^[TC]CTTT[TC]\w$/){  # TCT core motif (YCTTTY) Y=C/T
            $is_TCT="TRUE";
        }

        my $kozak_strength=&score_kozak($seq_up,$seq_down);
        my $tisu_strength=&score_TISU($seq_TIS_window);
        my $upstream_kozak_strength=&score_kozak_upstream($seq_up);

        my $cert=0; 
        if (exists ($cert_count{$gene})){ 
            $cert=$cert_count{$gene}; 
        }
        
        my $cert_short=0; 
        if (exists ($cert_short_count{$gene})){ 
            $cert_short=$cert_short_count{$gene}; 
        }

        my $polyU=0; 
        if (exists ($polyU_count{$gene})){ 
            $polyU=$polyU_count{$gene}; 
        }

        print "$gene,$transcript,$chr,$highest_cage_peak_value,$leader_length,$CDS_length,$trailer_length,$TSS_beginning_RNAseq_FPKM,$TSS_beginning_totalRNA_FPKM,$TSS_beginning_RIBOseq_FPKM,$TSS_beginning_TCPseq_SSU_FPKM,$TSS_beginning_TCPseq_LSU_FPKM,$leader_RNAseq_FPKM,$leader_totalRNA_FPKM,$leader_RIBOseq_FPKM,$leader_TCPseq_SSU_FPKM,$leader_TCPseq_LSU_FPKM,$start_codon_RNAseq_FPKM,$start_codon_totalRNA_FPKM,$start_codon_RIBOseq_FPKM,$start_codon_TCPseq_SSU_FPKM,$start_codon_TCPseq_LSU_FPKM,$CDS_RNAseq_FPKM,$CDS_totalRNA_FPKM,$CDS_RIBOseq_FPKM,$CDS_TCPseq_SSU_FPKM,$CDS_TCPseq_LSU_FPKM,$stop_codon_RNAseq_FPKM,$stop_codon_totalRNA_FPKM,$stop_codon_RIBOseq_FPKM,$stop_codon_TCPseq_SSU_FPKM,$stop_codon_TCPseq_LSU_FPKM,$trailer_RNAseq_FPKM,$trailer_totalRNA_FPKM,$trailer_RIBOseq_FPKM,$trailer_TCPseq_SSU_FPKM,$trailer_TCPseq_LSU_FPKM,$complete_leader_RNAseq_FPKM,$complete_leader_totalRNA_FPKM,$complete_leader_RIBOseq_FPKM,$complete_leader_TCPseq_SSU_FPKM,$complete_leader_TCPseq_LSU_FPKM,$complete_CDS_RNAseq_FPKM,$complete_CDS_totalRNA_FPKM,$complete_CDS_RIBOseq_FPKM,$complete_CDS_TCPseq_SSU_FPKM,$complete_CDS_TCPseq_LSU_FPKM,$complete_trailer_RNAseq_FPKM,$complete_trailer_totalRNA_FPKM,$complete_trailer_RIBOseq_FPKM,$complete_trailer_TCPseq_SSU_FPKM,$complete_trailer_TCPseq_LSU_FPKM,$tss_proximal_region_RNAseq_FPKM,$tss_proximal_region_totalRNA_FPKM,$tss_proximal_region_RIBOseq_FPKM,$tss_proximal_region_TCPseq_SSU_FPKM,$tss_proximal_region_TCPseq_LSU_FPKM,$tis_proximal_region_RNAseq_FPKM,$tis_proximal_region_totalRNA_FPKM,$tis_proximal_region_RIBOseq_FPKM,$tis_proximal_region_TCPseq_SSU_FPKM,$tis_proximal_region_TCPseq_LSU_FPKM,$overlapping_gene,$leader_potentially_overlaps_upstream_gene,$gene_potentially_overlaps_downstream_leader,$gene_overlaps_nc,$pre_TSS_start_nt,$TSS_start_nt,$post_TSS_start_nt,$seq_up,$seq_TIS,$seq_down,$seq_pre_stop,$seq_stop,$seq_post_stop,$exon_exon_junction_count_leader,$kozak_strength,$upstream_kozak_strength,$tisu_strength,$is_histone,$is_TOP,$is_TCT,$cert,$cert_short,$polyU\n";
    }
}

exit;

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#open a gtf file and get the transcript lengths of the protien coding genes
#sub transcript_lengths{
#
#    my $in_gtf = shift;
#    my %transcripts_ref = shift;

#    my $in_gtf = $_[0];
#    my $transcripts_ref = $_[1];            #key = gene_id, transcript_id, #value = sum_exon_lengths;
#
#    #to dereferece (makes a new copy of the hash)
#    my %transcripts=%{$transcripts_ref};
#
#    open(GENES1,$inGtf) || die "can't open $inGtf";      #gft is 1 based
#    while (<GENES1>){
#        unless(/^#/){
#            my @b=split("\t");
#            my $chr=$b[0];
#            my $class=$b[2];
#            my $start=$b[3];
#            my $end=$b[4];
#            my $dir=$b[6];
#            my ($gene_id) = $b[8] =~ /gene_id\s"([^\"]+)";/;
#            my ($transcript_id) = $b[8] =~ /transcript_id\s"([^\"]+)";/;
#            my $gene_biotype="NA";
#            if ($b[8] =~ /gene_biotype\s"([^\"]+)";/){
#               $gene_biotype=$1;
#            }
#
#            if ($gene_id && $transcript_id){
#                if ($gene_biotype eq "protein_coding"){ #restrict to protien coding genes (control for ncRNAs)
#                    if ($class eq "exon"){
#                        if ($dir eq "+"){
#                            for ($start .. $end){
#                                $transcripts{$gene_id}{$transcript_id}++;
#                           }
#                        }else{
#                            for ($start .. $end){
#                                $transcripts{$gene_id}{$transcript_id}++;
#                            }
#                        }
#                    }
#                }
#            }
#        }
#    }
#    close (GENES1);
#
#    return #I don't need to 
#    return (\%in_frame_fwd, \%in_frame_rev);
#}

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

        if ($organism eq "yeast"){
            if (exists ($PWM_yeast{$count}{$base} )){
                $score+=$PWM_yeast{$count}{$base};
            }
            $count++;
        }else{
            if (exists ($PWM{$count}{$base} )){
                $score+=$PWM{$count}{$base};
            }
            $count++;
       } 
    }
    return $score;
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#score the upstream tis sequence against the upstream Zebrafish kozak PWM
sub score_kozak_upstream{
    my $seq_to_score = shift;
    my $score=0;
    my @seq_to_score=split("",$seq_to_score);
    my $count=0;

   for my $base (@seq_to_score){
       if (exists ($PWM_upstream{$count}{$base} )){
           $score+=$PWM_upstream{$count}{$base};
       }
       $count++;
   }
   return $score;
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#score a 12nt sequence against the TISU motif
sub score_TISU{
    my $seq = shift;
    my $score=0;
    my $seq_to_score=uc($seq); #set to uppercase   
    my @seq_to_score=split("",$seq_to_score);
    my $count=0;

    for my $base (@seq_to_score){
        if (exists ($PWM_TISU{$count}{$base} )){
            $score+=$PWM_TISU{$count}{$base};
        }
        $count++;
    }
    return $score;
}
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#


