import numpy as np
import os
import sys
import pdb
import gzip
import scipy.stats as ss
from scipy import stats
from sklearn import linear_model


# Create dictionary of all genes we are going to be testing (ie those that passed our filters)
def get_measured_genes(expression_mat):
    dicti = {}  # initialize
    head_count = 0  # For header
    f = open(expression_mat)
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:  # skip header
            head_count = head_count +1
            continue
        # Gene ids are found in the first column of each line
        gene_id = data[0]
        dicti[gene_id] = -1
    return dicti

# Create dictionary mapping all genes found in measured_genes, and are located on chromosome $chrom_num, to their tss
def get_mapping_from_gene_to_tss(gencode_gene_annotation_file, chrom_num, measured_genes):
    gene_to_tss_mapping = {}
    f = gzip.open(gencode_gene_annotation_file)
    count = 0
    for line in f:
        line = line.rstrip()
        data = line.split()
        if line.startswith('#'):  # ignore header lines
            continue
        gene_type = data[13].split('"')[1]  # ie protein_coding,pseudo_gene,etc
        gene_name = data[9].split('"')[1].split('.')[0]  # ensamble id
        line_chrom_num = data[0]
        start = int(data[3])  # Start  of gene
        end = int(data[4])  # End (downstream) of gene
        gene_part = data[2]  # gene,UTR,exon,etc
        strand = data[6]  # either positive or negative
        if line_chrom_num != 'chr' + chrom_num:  # limit to chromosome of interest
            continue
        if gene_name not in measured_genes:  # Only care about genes that we have measurements for
            continue
        # We now are limited to measured genes on the correct chromosome
        if gene_name not in gene_to_tss_mapping:  # We haven't seen this gene before
            if strand == '+':  # positive strand
                gene_to_tss_mapping[gene_name] = (start, strand)
            elif strand == '-':  # negative strand
                gene_to_tss_mapping[gene_name] = (end, strand)
        else:  # We've seen this gene before
            old_tuple = gene_to_tss_mapping[gene_name]
            old_tss = old_tuple[0]
            old_strand = old_tuple[1]
            tss = old_tss
            if old_strand != strand:  # Error checking
                print('ASSUMPTION ERROR')
                pdb.set_trace()
            if strand == '+' and start < old_tss: 
                tss = start
            elif strand == '-' and end > old_tss:
                tss = end
            gene_to_tss_mapping[gene_name] = (tss, strand)
    return gene_to_tss_mapping



# Prepare files for matrix eqtl related to gene expression:
def prepare_gene_expression_files(chrom_num, gene_expression_input_file, gencode_gene_annotation_file, expression_matrix_file, gene_location_file, normalization_method):
    # Create dictionary of all genes we are going to be testing (ie those that passed our filters)
    measured_genes = get_measured_genes(gene_expression_input_file)
    # Create dictionary mapping all genes found in measured_genes, and are located on chromosome $chrom_num, to their tss
    gene_to_tss = get_mapping_from_gene_to_tss(gencode_gene_annotation_file, chrom_num, measured_genes)

    # open file handle for gene expression matrix file
    t_expr = open(expression_matrix_file, 'w')
    # open file handle for gene location file
    t_loc = open(gene_location_file, 'w')
    t_loc.write('geneid\tchr\tsq\ts2\n')  # header

    # Loop through input quantile normalized expression file
    head_count = 0  # used to identify header of input file
    f = open(gene_expression_input_file)
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            t_expr.write('id')  # start writing header of output file
            head_count = head_count + 1
            # Keep track of which indices (samples) correspond to samples at the current time step
            indices = []
            # Keep ordered list of which cell lines are used
            cell_lines = []
            for i, sample_id in enumerate(data):
                if sample_id.endswith('id') == False:
                    indices.append(i)
                    cell_lines.append(sample_id.split('_')[0])
                    t_expr.write('\t' + sample_id.split('_')[0])
            t_expr.write('\n')
            continue
        gene_id = data[0]
        if gene_id not in gene_to_tss:  # Gene is not on this chromosome
            continue
        # Add line to gene location file
        t_loc.write(gene_id + '\tchr' + chrom_num + '\t' + str(gene_to_tss[gene_id][0]) + '\t' + str(gene_to_tss[gene_id][0]) + '\n')

        # normalize the expression data in various ways..:
        if normalization_method == 'none': # Don't do anything to the expression data
            expr_vec = np.asarray(data)[indices]
        elif normalization_method == 'standardize':  # standardize the expression data
            expr_vec_temp = np.asarray(data)[indices]
            expr_vec_float = expr_vec_temp.astype(float)  # convert to float
            expr_vec_standardized_float = (expr_vec_float - np.mean(expr_vec_float))/np.std(expr_vec_float)  # standardize
            expr_vec = expr_vec_standardized_float.astype(str)  # convert back to strings
        elif normalization_method == 'gaussian_projection':
            expr_vec_temp = np.asarray(data)[indices]
            expr_vec_float = expr_vec_temp.astype(float)  # convert to float
            expr_vec_ranked = ss.rankdata(expr_vec_float)/(len(expr_vec_float) + 1)
            expr_vec = stats.norm.ppf(expr_vec_ranked)
            expr_vec = expr_vec.astype(str)
            
        # Add line to gene expression matrix file
        t_expr.write(gene_id + '\t' + '\t'.join(expr_vec) + '\n')
    t_expr.close()
    t_loc.close()
    return cell_lines

# Compute maf of dosage based genotype
def get_maf(genotype_array):
    af = np.sum(genotype_array)/(2.0*len(genotype_array))
    return min(af, 1.0 - af)

# Extract an array of length number of time steps, where each element of the array is another array that contains indices of cell lines observed for this time step
def get_indices_for_each_time_step(data, cell_lines_all_time_steps):
    indices_all_time_steps = []
    for time_step, observed_cell_lines in enumerate(cell_lines_all_time_steps):
        indices = []
        for i, val in enumerate(data):
            if val in observed_cell_lines:
                indices.append(i)
        indices_all_time_steps.append(np.asarray(indices))
    return indices_all_time_steps

def pass_maf_cutoff_all_time_steps(data,indices_all_time_steps, maf_cutoff):
    pass_filter = True
    for time_step,index_vector in enumerate(indices_all_time_steps):
        genotype_for_one_time_step = data[index_vector].astype(float)
        if get_maf(genotype_for_one_time_step) < maf_cutoff:
            pass_filter = False
    return pass_filter

# Prepare files for matrix eqtl related to genotype
def prepare_genotype_files(chrom_num, ordered_cell_lines, dosage_genotype_file, genotype_matrix_file, variant_location_file, maf_cutoff):
    t_geno = open(genotype_matrix_file, 'w')
    t_loc = open(variant_location_file, 'w')
    # Print location file header
    t_loc.write('snp\tchr\tpos\n')

    f = gzip.open(dosage_genotype_file)
    chrom_num = 'chr' + chrom_num
    for line in f:
        line = line.rstrip()
        data = line.split()
        if line.startswith('#CHROM'):  # HEADER with sample names
            # Find indices that map genotype cell line order to that order found in the expression matrix
            indices = []
            for cell_line in ordered_cell_lines:
                for i,val in enumerate(data):
                    if val == cell_line:
                        indices.append(i)
            # Quick check to make sure this mapping is correct
            if np.array_equal(np.asarray(data)[indices], np.asarray(ordered_cell_lines)) == False:
                print('FATAL ERROR')
                pdb.set_trace()
            # Print genotype matrix header
            t_geno.write('id\t' + '\t'.join(np.asarray(data)[indices]) + '\n')
            continue
        if line.startswith('#'):  # headers we do not care about for now
            continue
        # normal line
        # Extract relavent features of line
        rs_id = data[2]
        chromer = data[0]
        pos = data[1]


        if chrom_num != data[0]:  # Not on correct chromsome
            continue

        if rs_id == '.':  #remove variants that do not have an rsID
            continue

        # Get genotype data in correct order
        genotype_data = np.asarray(data)[indices]
       
        # compute maf of our samples in this time step
        maf = get_maf(genotype_data.astype(float))

        # Ignore snps that have maf < cutoff in any time step
        if maf < maf_cutoff:
            continue

        # Print variant location
        t_loc.write(rs_id + '\t' + chromer + '\t' + pos + '\n')

        # Print genotype line
        t_geno.write(rs_id + '\t' + '\t'.join(genotype_data) + '\n')
    t_geno.close()
    t_loc.close()


# Extract list of cell lines that are found in all time steps (some time steps have different numbers of cell lines)
# When calling maf cutoff, we are to use only cell lines found in all time steps
def extract_list_of_cell_lines_in_all_time_steps(quantile_normalized_expression):
    # First, extract array of sample ids from quantile normalized expression matrix
    f = open(quantile_normalized_expression)
    head_count = 0
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            sample_ids = data[1:]
            continue
        continue
    f.close()
    return sample_ids

# Prepare files for matrix eqtl related to covariates when expression data is prepared independently at each time step
def prepare_covariate_files(ordered_cell_lines, pca_loading_file, num_pcs, covariate_matrix_file):
    t = open(covariate_matrix_file, 'w')  # open output file handle
    data = np.transpose(np.loadtxt(pca_loading_file, dtype=str))  # load in data
    # Make sure sample names are in correct order
    sample_names = data[0,1:]
    t.write('id')
    for i,ele in enumerate(sample_names):
        cell_line_id = ele
        t.write('\t' + cell_line_id)
        if cell_line_id != ordered_cell_lines[i]:
            print('FATAL ERROR IN COVARIATE PREP')
            pdb.set_trace()
    t.write('\n')
    # write Latent factors to covariate file
    for row_num in range(1, (num_pcs+1)):
        t.write('\t'.join(data[row_num,:]) + '\n')
    t.close()

# Prepare covariate (PCs) for when expression data is prepared in aggegrate
# SVA is used to estimate latent factors
def prepare_covariate_files_aggregrate(ordered_cell_lines, sva_loading_file, num_factors, covariate_matrix_file, time_step):
    t = open(covariate_matrix_file, 'w')  # open output file handle
    t.write('id')  # Write first element of header of output file

    aa = np.transpose(np.loadtxt(sva_loading_file,dtype=str))  # Load in SVA data
    
    # Learn indices of samples that correspond to samples at the current time step 
    # Also write header
    indices = [] 
    sample_ids = aa[0,1:]
    counter = 0
    for i, val in enumerate(sample_ids):
        i_cell_line = val.split('_')[0]
        i_time_step = val.split('_')[1]
        if i_time_step != time_step:
            continue
        if ordered_cell_lines[counter] == i_cell_line:
            indices.append(i)
            counter = counter + 1
            t.write('\t' + i_cell_line)
    # Simple check that what we did was correct
    if len(indices) != len(ordered_cell_lines):
        print('Fatal error')
        pdb.set_trace()
    t.write('\n')

    # Create 1 line for each additional latent factor
    for row_num in range(1, (num_factors +1)):
        row_name = aa[row_num, 0]
        data = np.asarray(aa[row_num,1:])[indices]  # limit to samples from this time step
        t.write(row_name + '\t' + '\t'.join(data) + '\n')
    t.close()

# Regress out the effects of latent variables, and save the corrected data to $gene_expression_input_file
def clean_expression_data(gene_expression_input_file, quantile_normalized_expression, num_pcs, sva_loading_file):
    #  Load in data
    expr_data = np.loadtxt(quantile_normalized_expression, dtype=str)
    sva_data = np.transpose(np.loadtxt(sva_loading_file, dtype=str))
    
    #  Parse loaded data
    expr_sample_ids = expr_data[0, 1:]
    gene_ids = expr_data[1:, 0]
    sva_sample_ids = sva_data[0, 1:]
    sva_factor_ids = sva_data[1:, 0]
    expr = expr_data[1:, 1:].astype(float)
    sva = np.transpose(sva_data[1:, 1:].astype(float))

    # Simple check to make sure assumptions are met
    if np.array_equal(expr_sample_ids, sva_sample_ids) == False:
        print('Fatal Error')
        pdb.set_trace()

    # Initialize corrected data matrix
    corrected_data = np.zeros(expr.shape)

    # Initialize output file
    t = open(gene_expression_input_file, 'w')
    # Write header to output file
    header = '\t'.join(expr_data[0, :])
    t.write(header + '\n')

    # Get number of genes
    num_genes = len(gene_ids)
    # Loop through all genes
    for gene_index in range(num_genes):
        gene_expression = expr[gene_index, :]  # Extract expression values for this gene
        # Fit linear model of latent factors onto the gene expression
        model = linear_model.LinearRegression(fit_intercept=True)
        modelfit = model.fit(sva, gene_expression)
        beta = modelfit.coef_
        # Regress out this model
        for i, val in enumerate(beta):
            gene_expression = gene_expression - sva[:, i]*val
        # Write corrected expression for this gene to the output file
        gene_name = gene_ids[gene_index]
        t.write(gene_name + '\t' + '\t'.join(gene_expression.astype(str)) + '\n')
    t.close()


########################
# Command line arguments
########################
chrom_num = sys.argv[1]  # Run each chromosome seperately (OK cause doing cis-eqtls)
dosage_genotype_file = sys.argv[2]  # File containing dosage-based genotype information
quantile_normalized_expression = sys.argv[3]  # Quantile normalized expression data (uncorrected data)
gencode_gene_annotation_file = sys.argv[4]  # hg19 gencodge gene annotation
output_root = sys.argv[5]  # Prefix of output files
maf_cutoff = float(sys.argv[6])  # Only use variants with maf >= maf_cutoff
normalization_method = sys.argv[7]  # String flag on how we want to normalize the expression data
data_prep_version = sys.argv[8]  # String flag to determine whether we use $quantile_normalized_expression or $quantile_normalized_time_step_independent_expression
num_pcs = (sys.argv[9])  # the number of PCs to use
pca_loading_file = sys.argv[10]


# Prepare files for matrix eqtl related to gene expression:
expression_matrix_file = output_root + 'expression.txt'
gene_location_file = output_root + 'gene_location.txt'
ordered_cell_lines = prepare_gene_expression_files(chrom_num, quantile_normalized_expression, gencode_gene_annotation_file, expression_matrix_file, gene_location_file, normalization_method)

# Prepare files for matrix eqtl related to genotype
genotype_matrix_file = output_root + 'genotype.txt'
variant_location_file = output_root + 'variant_location.txt'
prepare_genotype_files(chrom_num, ordered_cell_lines, dosage_genotype_file, genotype_matrix_file, variant_location_file, maf_cutoff)


# Prepare files for matrix eqtl related to covariates (not used because we are using corrected data)
covariate_matrix_file = output_root + 'covariate.txt'  # output file

# Prepare covariates seperately depending on whether using time_step_independent vs time_step_aggregrate
prepare_covariate_files(ordered_cell_lines, pca_loading_file, int(num_pcs), covariate_matrix_file)
