import numpy as np 
import os
import sys
import pdb



def get_mapping_from_gene_id_to_position(gene_location_file):
    # Initialize mapper
    dicti = {}
    head_count = 0  # Used to skip header
    f = open(gene_location_file)
    for line in f:
        if head_count == 0:  # Skip header
            head_count = head_count + 1
            continue
        line = line.rstrip()
        data = line.split()
        # Extract relevent fields
        gene_id = data[0]
        pos = data[2]
        dicti[gene_id] = pos
    return dicti


def get_mapping_from_variant_id_to_position(file_name):
    # Initialize mapper
    dicti = {}
    head_count = 0  # Used to skip header
    f = open(file_name)
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        # Extract relevent fields
        rs_id = data[0]
        pos = data[2]
        dicti[rs_id] = pos
    return dicti

# Compute maf of dosage based genotype
def get_maf(genotype_array):
    af = np.sum(genotype_array)/(2.0*len(genotype_array))
    return min(af, 1.0 - af)

def get_mapping_from_variant_id_to_allele_frequency(genotype_file):
    # Initialize mapper
    dicti = {}
    head_count = 0  # Used to skip header
    f = open(genotype_file)
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        # Extract relevent fields
        rs_id = data[0]
        geno = np.asarray(data[1:]).astype(float)
        dicti[rs_id] = str(get_maf(geno))
    return dicti

# Use dictionaries to add information to our matrix eqtl output file
def add_new_columns_to_eqtl_output(matrix_eqtl_output, new_output, chrom_num, gene_name_to_position, variant_id_to_position, variant_id_to_allele_frequency):
    f = open(matrix_eqtl_output)
    t = open(new_output, 'w')
    # Write header
    t.write('chrom_num\tgene_id\tgene_position\tregulatory_site_id\tregulatory_site_position\tregulatory_allele_frequency\tbeta\tpvalue\n')
    head_count = 0  # Will be used to skip header
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        # Extact relevent fields
        rs_id = data[0]
        gene_id = data[1]
        beta = data[2]
        pvalue = data[4]

        # Write to new_output
        t.write(chrom_num + '\t' + gene_id + '\t' + gene_name_to_position[gene_id] + '\t' + rs_id + '\t' + variant_id_to_position[rs_id] + '\t' + variant_id_to_allele_frequency[rs_id] + '\t' + beta + '\t' + pvalue + '\n')
    t.close()

matrix_eqtl_output = sys.argv[1]
new_output = sys.argv[2]
genotype_file = sys.argv[3]
variant_location_file = sys.argv[4]
gene_location_file = sys.argv[5]
chrom_num = sys.argv[6]

# Create dictionaries containing information we wish to add to our matrix eqtl output file
gene_name_to_position = get_mapping_from_gene_id_to_position(gene_location_file)
variant_id_to_position = get_mapping_from_variant_id_to_position(variant_location_file)
variant_id_to_allele_frequency = get_mapping_from_variant_id_to_allele_frequency(genotype_file)

# Use dictionaries to add information to our matrix eqtl output file
add_new_columns_to_eqtl_output(matrix_eqtl_output, new_output, chrom_num, gene_name_to_position, variant_id_to_position, variant_id_to_allele_frequency)