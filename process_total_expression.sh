#!/bin/bash
#SBATCH --time=20:00:00 --mem=12G



input_bam_dir="$1"
exon_file="$2"
processed_total_expression_dir="$3"
time_series_quantile_normalized="$4"
covariate_dir="$5"

if false; then
Rscript preprocess_total_expression.R $processed_total_expression_dir $exon_file $input_bam_dir $time_series_quantile_normalized
fi


Rscript prepare_covariate_files.R $processed_total_expression_dir $covariate_dir