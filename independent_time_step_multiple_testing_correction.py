import numpy as np
import os
import sys
import pdb




# Run bonferonni correction treating every heterozygous site as a gene
def run_bonferonni_correction(testing_prefix, concatenated_output_file, bonferonni_corrected_output_file, testing_suffix):
    # Initialize dictionary that contains mapping from heterozygous site to a tuple that contains (pvalue,full_line, number of regulatory snps)
    sites = {}
    t = open(concatenated_output_file, 'w')
    for chrom_num in range(1,23):
        # Open up qtl output file for this chromosome
        chromosome_input_file = testing_prefix + 'chrom_' + str(chrom_num) + '_' + testing_suffix
        f = open(chromosome_input_file)
        head_count = 0 # used to find header
        # Loop through file
        for line in f:
            line = line.rstrip()
            data = line.split()
            if head_count == 0:
                head_count = head_count + 1
                if chrom_num == 1:
                    # Print header to concatenated output file. Only do for first chromosome b/c we only want header printed once
                    t.write(line + '\n')
                continue
            # print line to concatenated output file
            t.write(line + '\n')

            # Extract info from line of qtl file
            het_site_id = data[1]
            reg_site_id = data[3]
            pvalue = float(data[7])

            # If we've never seen this heterozygous site before
            if het_site_id not in sites:
                sites[het_site_id] = (pvalue, line,1)
            elif het_site_id in sites:  # We've seen this het site before
                # Check to see if the current pvalue is smaller than the previous best
                old_tuple = sites[het_site_id]
                old_pvalue = old_tuple[0]
                old_line = old_tuple[1]
                prev_count = old_tuple[2]
                if pvalue < old_pvalue: # new pvalue is better
                    sites[het_site_id] = (pvalue, line, prev_count + 1)
                else:  # old pvalue was better
                    sites[het_site_id] = (old_pvalue, old_line, prev_count + 1)          
        # Now that we have included this chromosome in the concatenate file, we can delete the original
        f.close()
        os.system('rm ' + chromosome_input_file)

    t.close()
    bf_tuples = []
    # Loop through sites and perform bonforonni correction
    for site_id in sites.keys():
        site_tuple = sites[site_id]
        bf_pvalue = site_tuple[0]*site_tuple[2]  # bf correction
        site_line = site_tuple[1]  # Full line corresponding to best site-regulatorySite pair
        new_line = site_line + '\t' + str(bf_pvalue)  # Add bf-corrected pvalue to line
        bf_tuples.append((bf_pvalue,new_line))  # Put back in tuple
    # Sort by bf-corrected pvalue
    sorted_bf_tuples = sorted(bf_tuples, key=lambda tup: tup[0])
    # print to output the best reg site per het. site
    t = open(bonferonni_corrected_output_file, 'w')
    t.write('chom_num\thet_site_id\thet_site_position\tregulatory_site_id\tregulatory_site_position\ttest_statistic\tregulatory_allele_frequency\tpvalue\tBF_pvalue\n')
    for tupler in sorted_bf_tuples:
        t.write(tupler[1] + '\n')
    t.close()

# Compute number of tests ran (for bh correction)
def get_number_of_tests_ran(file_name):
    f = open(file_name)
    count = 0
    for line in f:
        # 1 line for each test
        count = count + 1
    # header doesn't count
    count = count -1
    return count



# Run benjamin-hochberg on bf corrected data
def run_benjamini_hochberg_correction(bonferonni_corrected_output_file, bh_corrected_output_file, alpha):
    f = open(bonferonni_corrected_output_file)
    t = open(bh_corrected_output_file, 'w')
    head_count = 0  # to keep track of header
    number_of_tests = get_number_of_tests_ran(bonferonni_corrected_output_file)
    bf_count = 1  # parameter to keep track of for bh correction
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            # Print same header
            t.write(line + '\n')
            continue
        bf_corrected_pvalue = float(data[8])

        if bf_corrected_pvalue <= bf_count*alpha/number_of_tests:  # This test is significant
            t.write(line + '\n')
            bf_count = bf_count + 1
        else: # Test is not significant
            # break out of loop
            break
    t.close()
    t.close()






##########################################
# Input Data
##########################################

# There will be 22 files of form $testing_prefix + $chrom_num + '_aseqtl_results.txt'. They correspond to 22 seperate qtl runs for each chromosome
# This script merges those files. Runs bonferonni correction at the gene level. And then Benjamini-Hochberge at the genome level
testing_prefix = sys.argv[1]
testing_suffix = sys.argv[2]

alpha=.05 # Significance threshold for bh correction

concatenated_output_file = testing_prefix + testing_suffix
bonferonni_corrected_output_file = testing_prefix + 'bonferonni_correction.txt'
bh_corrected_output_file = testing_prefix + 'bh_correction_' + str(alpha) + '.txt'

# Run bonferonni correction treating every heterozygous site as a gene
run_bonferonni_correction(testing_prefix, concatenated_output_file, bonferonni_corrected_output_file, testing_suffix)

# Run benjamin-hochberg on bf corrected data
run_benjamini_hochberg_correction(bonferonni_corrected_output_file, bh_corrected_output_file, alpha)