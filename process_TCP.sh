#!/bin/bash

#remove peaks from TCP-seq bams
#reads along transcripts with the sam 5' and 3', that are more than 200X the average for the mezn expression of that transcript (sum_hits/transcript length)
#plot read length distributions
#plot scaled coverage over trasncript features
#plot reads length heatmaps over start and stop codons 

usage(){
cat << EOF
usage: $0 options

script to trim and align riboseq fastq datasets 

OPTIONS:
    -s  path to SSU bam
    -l  path to LSU bam
    -f  path to SSU fastq
    -r  path to LSU fastq
    -c  path to fwd cage TAGS
    -w  path to rev cage TAGs
    -n  path to RNA-seq bam
    -t  path to totalRNA-seq bam
    -i  path to ribo-seq bam
    -o  path to output dir
    -h  this help message
EOF
}

while getopts "s:l:f:r:c:w:o:n:t:i:h" opt; do
    case $opt in
        f)
            SSU_fastq=$OPTARG
            echo "-f SSU fastq $OPTARG"
            ;;

        r)
            LSU_fastq=$OPTARG
            echo "-r LSU fastq $OPTARG"
            ;;

        s)
            SSU_bam=$OPTARG
            echo "-s SSU bam $OPTARG":wq
            ;;

        l)
            LSU_bam=$OPTARG
            echo "-l LSU bam $OPTARG"
            ;;

        c)  
            CAGE_fwd=$OPTARG
            echo "-c CAGE fwd tags $OPTARG"
            ;;

        w)
            CAGE_rev=$OPTARG
            echo "-w CAGE rev tags $OPTARG"
            ;;

        n)  
            rna_bam=$OPTARG
            echo "-n rna bam $OPTARG"
            ;;

        t) 
            total_rna_bam=$OPTARG
            echo "-n total rna bam $OPTARG"
            ;;

        i)  
            ribo_bam=$OPTARG
            echo "-i ribo bam $OPTARG"
            ;;

        o)
            out_dir=$OPTARG
            echo "-o output folder $OPTARG"
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

if [ ! $SSU_bam ] || [ ! $LSU_bam ] || [ ! $SSU_fastq ] || [ ! $LSU_fastq ] || [ ! $CAGE_fwd ] || [ ! $CAGE_rev ] || [ ! $rna_bam ] || [ ! $total_rna_bam ] || [ ! $ribo_bam ] || [ ! $out_dir ] ; then
        echo "something's missing - check your command lines args"
        usage
        exit 1
fi

GTF=/export/valenfs/data/references/Zv10_zebrafish/Danio_rerio.GRCz10.81.gtf 
FA=/export/valenfs/data/references/Zv10_zebrafish/Danio_rerio.GRCz10.dna.toplevel.fa
GO=/export/valenfs/projects/adam/TCP_seq/data_files/GRCz10_GO_terms_ensembl_mart_export.csv

echo "GTF:$GTF"
echo "FA:$FA"

if [ ! -d $out_dir ]; then
    mkdir $out_dir
fi


#-----------------------------------------------------------------------------------------
# 1 Find the most highly expressed transcript of each gene
#-----------------------------------------------------------------------------------------

nice -n 10 perl /export/valenfs/projects/adam/final_results/scripts/git/TCPeaking/select_transcripts_based_on_totalRNA_coverage.pl $GTF $total_rna_bam ${out_dir}/most_highly_expressed_transcripts.csv


#-----------------------------------------------------------------------------------------
# 2 calculate leaders
#-----------------------------------------------------------------------------------------

nice -n 10 perl /export/valenfs/projects/adam/final_results/scripts/git/TCPeaking/cage_peak_assignment.pl $GTF $CAGE_fwd $CAGE_rev ${out_dir}/most_highly_expressed_transcripts.csv ${out_dir}/leaders.csv


#------------------------------------------------------------------------------------------
# 3 Tidy up polyA adaptor trimming
#------------------------------------------------------------------------------------------

if [ ! -d ${out_dir}/tidy_bams ]; then 
   mkdir ${out_dir}/tidy_bams
fi

nice -n 10 perl /export/valenfs/projects/adam/final_results/scripts/git/TCPeaking/extend_polyA_trimmed_alignments.pl $GTF $FA ${out_dir}/leaders.csv $SSU_fastq $SSU_bam ${out_dir}/most_highly_expressed_transcripts.csv | samtools view -Sb - | samtools sort -T sort_SSU_tmp -O bam -o ${out_dir}/tidy_bams/SSU_updated.bam -

nice -n 10 perl /export/valenfs/projects/adam/final_results/scripts/git/TCPeaking/extend_polyA_trimmed_alignments.pl $GTF $FA ${out_dir}/leaders.csv $LSU_fastq $LSU_bam ${out_dir}/most_highly_expressed_transcripts.csv | samtools view -Sb - | samtools sort -T sort_LSU_tmp -O bam -o ${out_dir}/tidy_bams/LSU_updated.bam -


#------------------------------------------------------------------------------------------
# 4 Remove peaks
#------------------------------------------------------------------------------------------

if [ ! -d ${out_dir}/tidy_bams ]; then 
   mkdir ${out_dir}/tidy_bams
fi

nice -n 10 perl /export/valenfs/projects/adam/final_results/scripts/git/TCPeaking/remove_peaks_from_TCP_seq_bam.pl $GTF $FA ${out_dir}/tidy_bams/SSU_updated.bam ${out_dir}/leaders.csv ${out_dir}/most_highly_expressed_transcripts.csv 200 | samtools view -o ${out_dir}/tidy_bams/SSU_peaks_removed_200.bam -Sb - 

nice -n 10 perl /export/valenfs/projects/adam/final_results/scripts/git/TCPeaking/remove_peaks_from_TCP_seq_bam.pl $GTF $FA ${out_dir}/tidy_bams/LSU_updated.bam ${out_dir}/leaders.csv ${out_dir}/most_highly_expressed_transcripts.csv 200 | samtools view -o ${out_dir}/tidy_bams/LSU_peaks_removed_200.bam -Sb -


#------------------------------------------------------------------------------------------
# 5 Filter LSU and SSU for translating ribosomes lengths
#------------------------------------------------------------------------------------------

#store the names of the translating filtered SSU + LSU
BAM_SSU_200=${out_dir}/tidy_bams/SSU_peaks_removed_200_removed_translating_lengths_25_35.bam
BAM_LSU_200=${out_dir}/tidy_bams/LSU_peaks_removed_200_selected_translating_lengths_25_35.bam

#select only the fragment length corroposoding to the translating ribosomes
nice -n 10 perl /export/valenfs/projects/adam/final_results/scripts/git/TCPeaking/subset_fragment_lengths_in_LSU_bam.pl ${out_dir}/tidy_bams/LSU_peaks_removed_200.bam 25 35 | samtools view -o $BAM_LSU_200 -Sb -

#remove the fragment lengths corroposnding to the translating ribosomes
nice -n 10 perl /export/valenfs/projects/adam/final_results/scripts/git/TCPeaking/subset_fragment_lengths_in_SSU_bam.pl ${out_dir}/tidy_bams/SSU_peaks_removed_200.bam 25 35 | samtools view -o $BAM_SSU_200 -Sb -

echo "BAM SSU 200:$BAM_SSU_200"
echo "BAM LSU 200:$BAM_LSU_200"


#------------------------------------------------------------------------------------------
# 6 get read length distributions per feature
#------------------------------------------------------------------------------------------

if [ ! -d ${out_dir}/lengths_per_feature_TSS ]; then
   mkdir ${out_dir}/lengths_per_feature_TSS
fi

nice -n 10 perl /export/valenfs/projects/adam/final_results/scripts/git/TCPeaking/length_distribution_plot/length_distributions_per_feature.pl $GTF $FA ${out_dir}/tidy_bams/SSU_peaks_removed_200.bam ${out_dir}/leaders.csv ${out_dir}/most_highly_expressed_transcripts.csv use_cage > ${out_dir}/lengths_per_feature_TSS/SSU_read_length_distributions_per_feature_200.csv 

nice -n 10 perl /export/valenfs/projects/adam/final_results/scripts/git/TCPeaking/length_distribution_plot/length_distributions_per_feature.pl $GTF $FA ${out_dir}/tidy_bams/LSU_peaks_removed_200.bam ${out_dir}/leaders.csv ${out_dir}/most_highly_expressed_transcripts.csv use_cage > ${out_dir}/lengths_per_feature_TSS/LSU_read_length_distributions_per_feature_200.csv

R CMD /export/valenfs/projects/adam/final_results/scripts/git/TCPeaking/length_distribution_plot/plot_fragment_length_distibutions_per_feature.R ${out_dir}/lengths_per_feature_TSS/ read_length_distributions_per_feature_200 ${out_dir}/lengths_per_feature_TSS/SSU_read_length_distributions_per_feature_200.csv ${out_dir}/lengths_per_feature_TSS/LSU_read_length_distributions_per_feature_200.csv


#------------------------------------------------------------------------------------------
# 7 scaled feature plots
#------------------------------------------------------------------------------------------

if [ ! -d ${out_dir}/scaled_feature_plots_peaks_removed ]; then
   mkdir ${out_dir}/scaled_feature_plots_peaks_removed
fi

name=$(basename {out_dir}/tidy_bams/SSU_peaks_removed_200.bam)
prefix=${name%.bam};

nice -n 10 perl /export/valenfs/projects/adam/final_results/scripts/git/TCPeaking/scaled_coverage_metaplot/scaled_window_plots_TCPseq_excluding_3prime_peak.pl $GTF ${out_dir}/leaders.csv ${out_dir}/tidy_bams/SSU_peaks_removed_200.bam ${out_dir}/tidy_bams/LSU_peaks_removed_200.bam ${out_dir}/scaled_feature_plots_peaks_removed/ ${out_dir}/most_highly_expressed_transcripts.csv use_cage

#zscore
R CMD /export/valenfs/projects/adam/final_results/scripts/git/TCPeaking/scaled_coverage_metaplot/plot_read_coverage_per_feature_zscore.R ${out_dir}/scaled_feature_plots_peaks_removed/ ${prefix} ${out_dir}/scaled_feature_plots_peaks_removed/${prefix}_gene_matrix_coverage.csv


if [ ! -d ${out_dir}/scaled_feature_plots_translating_filter ]; then
   mkdir ${out_dir}/scaled_feature_plots_translating_filter
fi

name=$(basename {out_dir}/tidy_bams/SSU_peaks_removed_200_removed_translating_lengths_25_35.bam)
prefix=${name%.bam};

nice -n 10 perl /export/valenfs/projects/adam/final_results/scripts/git/TCPeaking/scaled_coverage_metaplot/scaled_window_plots_TCPseq_excluding_3prime_peak.pl $GTF ${out_dir}/leaders.csv $BAM_SSU_200 $BAM_LSU_200 ${out_dir}/scaled_feature_plots_translating_filter/ ${out_dir}/most_highly_expressed_transcripts.csv use_cage

#zscore
R CMD /export/valenfs/projects/adam/final_results/scripts/git/TCPeaking/scaled_coverage_metaplot/plot_read_coverage_per_feature_zscore.R ${out_dir}/scaled_feature_plots_translating_filter/ ${prefix} ${out_dir}/scaled_feature_plots_translating_filter/${prefix}_gene_matrix_coverage.csv


#------------------------------------------------------------------------------------------
# 8 heatmaps of start and stop codons
#------------------------------------------------------------------------------------------

if [ ! -d ${out_dir}/start_stop_heatmaps_200_selected ]; then
   mkdir ${out_dir}/start_stop_heatmaps_200_selected
fi

nice -n 10 bash /export/valenfs/projects/adam/final_results/scripts/git/TCPeaking/positional_metaplots_TIS_Stop_codon/run_transcript_coord_LSU_SSU_raw_10RNA.sh -s ${out_dir}/tidy_bams/SSU_peaks_removed_200.bam -r ${out_dir}/tidy_bams/LSU_peaks_removed_200.bam -o ${out_dir}/start_stop_heatmaps_200_selected/ -n $total_rna_bam -g $GTF -l ${out_dir}/leaders.csv -t ${out_dir}/most_highly_expressed_transcripts.csv -c use_cage


#-----------------------------------------------------------------------------------------
# 9 per transcripot plots
#-----------------------------------------------------------------------------------------

if [ ! -d ${out_dir}/per_transcript_plots_peaks_removed ]; then
   mkdir ${out_dir}/per_transcript_plots_peaks_removed
fi
                                
name=$(basename {out_dir}/tidy_bams/SSU_peaks_removed_200.bam)
prefix=${name%.bam};

#store the names of the translating filtered SSU + LSU
BAM_SSU_peaks_200=${out_dir}/tidy_bams/SSU_peaks_removed_200.bam
BAM_LSU_peaks_200=${out_dir}/tidy_bams/LSU_peaks_removed_200.bam

nice -n 10 perl /export/valenfs/projects/adam/final_results/scripts/git/TCPeaking/TCPseq_matrix_leader_motif.pl $GTF $FA $ribo_bam $rna_bam $total_rna_bam $BAM_SSU_peaks_200 $BAM_LSU_peaks_200 ${out_dir}/leaders.csv $GO ${out_dir}/most_highly_expressed_transcripts.csv zebrafish use_cage > ${out_dir}/per_transcript_plots_peaks_removed/${prefix}_matrix.csv


if [ ! -d ${out_dir}/per_transcript_plots_200_translating_filter ]; then
   mkdir ${out_dir}/per_transcript_plots_200_translating_filter
fi

name=$(basename {out_dir}/tidy_bams/SSU_peaks_removed_200_removed_translating_lengths_25_35.bam)
prefix=${name%.bam};

nice -n 10 perl /export/valenfs/projects/adam/final_results/scripts/git/TCPeaking/TCPseq_matrix_leader_motif.pl $GTF $FA $ribo_bam $rna_bam $total_rna_bam $BAM_SSU_200 $BAM_LSU_200 ${out_dir}/leaders.csv $GO ${out_dir}/most_highly_expressed_transcripts.csv zebrafish use_cage > ${out_dir}/per_transcript_plots_200_translating_filter/${prefix}_matrix.csv
