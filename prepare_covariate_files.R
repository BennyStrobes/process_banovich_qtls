args = commandArgs(trailingOnly=TRUE)
library(edgeR)
library(Rsubread)
library(Biobase)
library(preprocessCore)
library(ggplot2)
library(ggthemes)
library(glmnet)
library(reshape)
library(cowplot)
library(mvtnorm)
library(sva)
library(limma)

# Helper method to save DGE data structure to tab-deliminated text file
save_python_style_matrix <- function(counts, output_file, row_label_name) {
    #  Convert DGE data structure to matrix
    temp_mat <- as.matrix(counts)

    #  Edit colnames to include a header over the row-labels.
    revised_column_names <- colnames(temp_mat)
    revised_column_names[1] <- paste0(row_label_name,"\t",revised_column_names[1])

    write.table(temp_mat, output_file, quote=FALSE,col.names = revised_column_names, sep="\t")

}

# Write PC scores to output file
save_pcs <- function(sample_names, quant_expr, n, output_file) {
    svd1 <- svd(as.matrix(quant_expr))

    #  Scores of first n pc's across all samples
    pc <- svd1$v[,1:n]
    colnames(pc) <- paste0("PC",1:n)
    rownames(pc) <- sample_names
    save_python_style_matrix(pc, output_file, "Sample_id")
}



#########################################
# Command line args
#########################################

processed_total_expression_dir <- args[1]
covariate_dir <- args[2]

#########################################
# Load data
#########################################

#  Get sample information 
sample_info_file <- paste0(processed_total_expression_dir, "sample_info.txt")
sample_info <- read.table(sample_info_file, header=TRUE)

#  Get quantile normalized expression data
quantile_normalized_exp_file <- paste0(processed_total_expression_dir, "quantile_normalized.txt")
quant_expr <- read.csv(quantile_normalized_exp_file, header=TRUE, sep=" ")

#  Get quantile normalized expression data for cell line subset
cell_line_subset_quantile_normalized_exp_file <- paste0(processed_total_expression_dir, "cell_line_subset_quantile_normalized.txt")
cell_line_subset_quant_expr <- read.csv(cell_line_subset_quantile_normalized_exp_file, header=TRUE, sep=" ")

cell_line_subset <- colnames(cell_line_subset_quant_expr)

cell_line_subset <- substr(cell_line_subset, 2, nchar(cell_line_subset))
cell_line_subset[1] <-substr(cell_line_subset[1],8,12)

##############################################################################################################################
# Write PCs to output file
##############################################################################################################################
#  Number of pcs to save
n <- 10
#  Ouptut file to save PC loadings to 
pc_output_file <- paste0(covariate_dir, "all_lines_principal_components_", n, ".txt")
save_pcs(sample_info$Sample_name, quant_expr, n, pc_output_file)


#  Number of pcs to save
n <- 3
#  Ouptut file to save PC loadings to 
pc_output_file <- paste0(covariate_dir, "cell_line_subset_principal_components_", n, ".txt")
save_pcs(cell_line_subset, cell_line_subset_quant_expr, n, pc_output_file)
