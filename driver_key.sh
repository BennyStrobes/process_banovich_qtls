



##############################################################################
# Input Data
##############################################################################

# Directory containing bams created by Nick Banovich
# Files have naming system "RNAseq_"$sampleId".merged.sort.bam"
input_bam_dir="/project2/gilad/katie/nbanovich_iPSC/YRI/RNA-seq/Final_bams/"

# File created by "https://github.com/BennyStrobes/ipsc_preprocess_pipeline/"
# Contains exon annotations
exon_file="/project2/gilad/bstrober/ipsc_differentiation/preprocess/genome/exons.saf"


# File created by "https://github.com/BennyStrobes/ipsc_preprocess_pipeline/"
# Contains quantile normalized expression data from our (time series ipsc) analysis
time_series_quantile_normalized="/project2/gilad/bstrober/ipsc_differentiation/preprocess/processed_total_expression/quantile_normalized.txt"


# Dosage genotype for all cell lines in our analysis
dosage_genotype_file="/project2/gilad/bstrober/ipsc_differentiation/preprocess_input_data/genotypesYRI.gen.txt.gz"

# Gencode hg19 gene annotation file
gencode_gene_annotation_file="/project2/gilad/bstrober/ipsc_differentiation/preprocess_input_data/gencode.v19.annotation.gtf.gz"





##############################################################################
# Output Directories (script assumes these directories are created prior to running scripts)
##############################################################################

processed_total_expression_dir="/project2/gilad/bstrober/ipsc_differentiation/process_banovich_qtls/processsed_total_expression/"

covariate_dir="/project2/gilad/bstrober/ipsc_differentiation/process_banovich_qtls/covariates/"

eqtl_dir="/project2/gilad/bstrober/ipsc_differentiation/process_banovich_qtls/eqtl/"




##############################################################################
# Scripts
##############################################################################



######################
# Part 1: Convert from bams to normalized expression Data
if false; then
sh process_total_expression.sh $input_bam_dir $exon_file $processed_total_expression_dir $time_series_quantile_normalized $covariate_dir
fi




######################
# Part 2: Run independent time step e-qtl analysis

# Cis eqtl distance (EAGLE used 200 KB)
distance="50000"
# Minimum fraction of valid samples with less popular version of regulatory variant (homozygous reg variant vs heterozygous reg variant)
maf_cutoff=".1"
# The way in which we normalize the data
# Currently implemented for:
#####1. 'none'
#####2. 'standardize'
#####3. 'gaussian_projection'
normalization_method="none"
# The way in which the data was prepared
# Currently implemented for:
#####1. 'all_samples'  (quantile normalization and hidden factor correction was done across all cell lines)
#####2. 'cell_line_subset'  (quantile normalization and hidden factor correction was done for only cell lines found in our data)
data_prep_version="all_samples"
# The number of PCs to inlcude in the model 
num_pcs="3"
# Input expression data
quantile_normalized_expression=$processed_total_expression_dir"cell_line_subset_quantile_normalized.txt"
pca_loading_file=$covariate_dir"cell_line_subset_principal_components_3.txt"

sbatch eqtl_driver.sh $dosage_genotype_file $quantile_normalized_expression $gencode_gene_annotation_file $eqtl_dir $distance $maf_cutoff $normalization_method $data_prep_version $num_pcs $pca_loading_file

# Cis eqtl distance (EAGLE used 200 KB)
distance="50000"
# Minimum fraction of valid samples with less popular version of regulatory variant (homozygous reg variant vs heterozygous reg variant)
maf_cutoff=".1"
# The way in which we normalize the data
# Currently implemented for:
#####1. 'none'
#####2. 'standardize'
#####3. 'gaussian_projection'
normalization_method="none"
# The way in which the data was prepared
# Currently implemented for:
#####1. 'all_samples'  (quantile normalization and hidden factor correction was done across all cell lines)
#####2. 'cell_line_subset'  (quantile normalization and hidden factor correction was done for only cell lines found in our data)
data_prep_version="all_samples"
# The number of PCs to inlcude in the model 
num_pcs="10"
# Input expression data
quantile_normalized_expression=$processed_total_expression_dir"quantile_normalized.txt"
pca_loading_file=$covariate_dir"all_lines_principal_components_10.txt"

sbatch eqtl_driver.sh $dosage_genotype_file $quantile_normalized_expression $gencode_gene_annotation_file $eqtl_dir $distance $maf_cutoff $normalization_method $data_prep_version $num_pcs $pca_loading_file

