#!/bin/bash

#AG 11/03/19
#script to produce heatmaps of 5' and 3' length distributions from a bam file
#1 bam to sam
#2 get 5' and 3' reads
        #up/downstream of annotated start codon
        #up/downstream of annotated stop codons    
   
#dependanices (must be in path)
#samtools

usage(){
cat << EOF
usage: $0 options

script to trim and align riboseq fastq datasets 

OPTIONS:
    -s    path to the SSU BAM file
    -r    path to the LSU BAM file
    -n    path to total RNA BAM file
    -o    path to output dir
    -g    gtf file (matched to the bam gemone)
    -l    csv of leader definitions
    -t    csv or the most highly expressed transcripts per gene
    -c    flag to indicate if cage should be used to update leaders
    -h    this help message

example usage: run_transcript_coords_LSU_SSU.sh -b <in.bam> -g in.gtf -l in.leader -o <out_dir> -c use_cage

EOF
}

while getopts ":s:r:o:g:l:t:c:n:h" opt; do
    case $opt in 
        s)
            in_SSU=$OPTARG
            echo "-s SSU input bam $OPTARG"
            ;;
        r)
            in_LSU=$OPTARG
            echo "-r LSU input bam $OPTARG"
            ;;
        o)
            out_dir=$OPTARG
            echo "-o output folder $OPTARG"
            ;;
        g)
            in_gtf=$OPTARG
            echo "-g gtf file $OPTARG"
            ;;
        l)  
            in_leader=$OPTARG
            echo "-l leader file $OPTARG"
            ;;
        t)
            high_transcript=$OPTARG
            echo "-t highly expressed transcripts $OPTARG"
            ;;
        c)  
            use_cage=$OPTARG
            echo "-c use_leaders $OPTARG"
            ;;
        n)
            rna_bam=$OPTARG
            echo "-n rna bam $OPTARG"
            ;;
        h)
            usage
            exit
            ;;
        ?) 
            echo "Invalid option: -$OPTARG"
            usage
            exit 1
            ;;
    esac
done

if [ ! $in_SSU ] || [ ! $in_LSU ] || [ ! $rna_bam ] || [ ! $out_dir ] || [ ! $in_gtf ] || [ ! $in_leader ] || [ ! $high_transcript ]  || [ ! $use_cage ] ; then
        echo "something's missing - check your command lines args"
        usage
        exit 1
fi

if [ ! -d $out_dir ]; then
    mkdir $out_dir
fi

if [ ! -d ${out_dir}/tmp ]; then
    mkdir ${out_dir}/tmp
fi

sname=$(basename $in_SSU)
sprefix=${sname%.bam}

lname=$(basename $in_LSU)
lprefix=${lname%.bam}

PATH=$PATH:/net/apps/cbu/src/samtools-1.2/

#------------------------------------------------------------------------------------------
#1 Count 5' ends of reads
#------------------------------------------------------------------------------------------

#rank on the leader+CDS+trailer region. scale on the window region
nice -n 10 perl /export/valenfs/projects/adam/final_results/scripts/git/TCPeaking/positional_metaplots_TIS_Stop_codon/heatmap-matrix_by_transcript_raw_10rna.pl $in_gtf $in_leader $rna_bam $in_SSU $out_dir $high_transcript $use_cage

nice -n 10 perl /export/valenfs/projects/adam/final_results/scripts/git/TCPeaking/positional_metaplots_TIS_Stop_codon/heatmap-matrix_by_transcript_raw_10rna.pl $in_gtf $in_leader $rna_bam $in_LSU $out_dir $high_transcript $use_cage

#------------------------------------------------------------------------------------------ 
#2 Plot heatmaps
#------------------------------------------------------------------------------------------

Rscript /export/valenfs/projects/adam/final_results/scripts/git/TCPeaking/positional_metaplots_TIS_Stop_codon/TCP_heatmaps_TIS_raw_cal_bars.R ${out_dir}/${sprefix}_start_lengths_scale_5prime_q1.csv ${out_dir}/${sprefix}_start_lengths_scale_3prime_q1.csv ${out_dir}/${sprefix}_start_5prime_scale_q1.csv ${out_dir}/${sprefix}_start_3prime_scale_q1.csv ${out_dir}/${sprefix}_start_heatmaps_raw.pdf

Rscript /export/valenfs/projects/adam/final_results/scripts/git/TCPeaking/positional_metaplots_TIS_Stop_codon/TCP_heatmaps_TIS_raw_cal_bars.R ${out_dir}/${lprefix}_start_lengths_scale_5prime_q1.csv ${out_dir}/${lprefix}_start_lengths_scale_3prime_q1.csv ${out_dir}/${lprefix}_start_5prime_scale_q1.csv ${out_dir}/${lprefix}_start_3prime_scale_q1.csv ${out_dir}/${lprefix}_start_heatmaps_raw.pdf

#spcific colouring
Rscript /export/valenfs/projects/adam/final_results/scripts/git/TCPeaking/positional_metaplots_TIS_Stop_codon/TCP_heatmaps_TIS_raw_cal_bars_alt.R ${out_dir}/${lprefix}_start_lengths_scale_5prime_q1.csv ${out_dir}/${lprefix}_start_lengths_scale_3prime_q1.csv ${out_dir}/${lprefix}_start_5prime_scale_q1.csv ${out_dir}/${lprefix}_start_3prime_scale_q1.csv ${out_dir}/${lprefix}_start_heatmaps_raw_alt.pdf

#LSU 20 40
Rscript /export/valenfs/projects/adam/final_results/scripts/git/TCPeaking/positional_metaplots_TIS_Stop_codon/TCP_heatmaps_TIS_raw_cal_bars_20_40.R ${out_dir}/${lprefix}_start_lengths_scale_5prime_q1.csv ${out_dir}/${lprefix}_start_lengths_scale_3prime_q1.csv ${out_dir}/${lprefix}_start_5prime_scale_q1.csv ${out_dir}/${lprefix}_start_3prime_scale_q1.csv ${out_dir}/${lprefix}_start_heatmaps_raw_20_40.pdf
