#!/bin/bash
#SBATCH --time=05:00:00 --mem=6GB

dosage_genotype_file="$1"
quantile_normalized_expression="$2"
gencode_gene_annotation_file="$3"
output_dir="$4"
eqtl_distance="$5"
maf_cutoff="$6"
normalization_method="$7"
data_prep_version="$8"
num_pcs="$9"
pca_loading_file="${10}"

output_file_prefix=$output_dir"eqtl_prepare_eqtl_distance_"$eqtl_distance"_maf_cutoff_"$maf_cutoff"_normalization_meth_"$normalization_method"_data_prep_"$data_prep_version"_num_pcs_"$num_pcs

# Perorm qtl analysis in each chromosome seperately (okay because doing cis eqtls)
for chrom_num in $(seq 1 22); do

    echo $chrom_num
    output_root=$output_file_prefix"_chrom_"$chrom_num"_"
    
    # Prepare input files to matrix eqtl
    python prepare_independent_time_step_matrix_eqtl_files.py $chrom_num $dosage_genotype_file $quantile_normalized_expression $gencode_gene_annotation_file $output_root $maf_cutoff $normalization_method $data_prep_version $num_pcs $pca_loading_file
    
    # Temporary Input file names
    matrix_eqtl_genotype_file=$output_root"genotype.txt"
    matrix_eqtl_variant_loc_file=$output_root"variant_location.txt"
    matrix_eqtl_expression_file=$output_root"expression.txt"
    matrix_eqtl_gene_location_file=$output_root"gene_location.txt"
    matrix_eqtl_covariate_file=$output_root"covariate.txt"
    matrix_eqtl_output_file=$output_root"matrix_eqtl_out.txt"


    # Run matrix eqtl
    Rscript matrix_eqtl_wrapper.R $matrix_eqtl_genotype_file $matrix_eqtl_variant_loc_file $matrix_eqtl_expression_file $matrix_eqtl_gene_location_file $eqtl_distance $matrix_eqtl_output_file $matrix_eqtl_covariate_file $data_prep_version
    # Convert matrix eqtl output to easier to use format
    python convert_matrix_eqtl_output_to_standard_association_output.py $matrix_eqtl_output_file $output_root"eqtl_results.txt" $matrix_eqtl_genotype_file $matrix_eqtl_variant_loc_file $matrix_eqtl_gene_location_file $chrom_num

    # Delete temporary files
    rm $matrix_eqtl_genotype_file
    rm $matrix_eqtl_variant_loc_file
    rm $matrix_eqtl_expression_file
    rm $matrix_eqtl_gene_location_file
    rm $matrix_eqtl_output_file
    rm $matrix_eqtl_covariate_file

done

# Multiple testing correction for genome wide significance
python independent_time_step_multiple_testing_correction.py $output_file_prefix"_" "eqtl_results.txt"


