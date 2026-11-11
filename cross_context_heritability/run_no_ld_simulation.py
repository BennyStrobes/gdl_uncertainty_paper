import numpy as np
import os
import sys
import pdb
from delta_beta_var_inference import delta_beta_var_inference





def simulate_data(sample_size, snps_per_gene):
    # First simulate genotype data
    # Each snp has allele frequency drawn uniformly between 0.05 and 0.5
    allele_freqs = np.random.uniform(0.05, 0.5, size=snps_per_gene)
    genotypes = np.random.binomial(2, allele_freqs, size=(sample_size, snps_per_gene))

    # Standardize each snp to mean 0 and variance 1 across individuals
    genotypes = (genotypes - np.mean(genotypes, axis=0)) / np.std(genotypes, axis=0)

    # Simulate causal effects for two contexts
    # 5 causal snps per context: 3 shared across contexts, 2 unique to each context
    per_snp_h2 = 0.015
    shared_corr = 0.8
    n_shared = 6
    n_unique = 1
    causal_indices = np.random.choice(snps_per_gene, size=n_shared + 2*n_unique, replace=False)
    shared_indices = causal_indices[:n_shared]
    unique_indices_1 = causal_indices[n_shared:(n_shared + n_unique)]
    unique_indices_2 = causal_indices[(n_shared + n_unique):]

    causal_effects_1 = np.zeros(snps_per_gene)
    causal_effects_2 = np.zeros(snps_per_gene)

    # Shared causal snps: effects drawn from bivariate normal with correlation shared_corr
    shared_cov = per_snp_h2*np.asarray([[1.0, shared_corr], [shared_corr, 1.0]])
    shared_effects = np.random.multivariate_normal(np.zeros(2), shared_cov, size=n_shared)
    causal_effects_1[shared_indices] = shared_effects[:, 0]
    causal_effects_2[shared_indices] = shared_effects[:, 1]

    # Unique causal snps
    causal_effects_1[unique_indices_1] = np.random.normal(0.0, np.sqrt(per_snp_h2), size=n_unique)
    causal_effects_2[unique_indices_2] = np.random.normal(0.0, np.sqrt(per_snp_h2), size=n_unique)
    

    # get genetic components
    genetic_comp_1 = genotypes.dot(causal_effects_1)
    genetic_comp_2 = genotypes.dot(causal_effects_2)

    # Add environmental noise
    env_noise_1 = np.random.normal(0.0, np.sqrt(1.0 - np.var(genetic_comp_1)), size=sample_size)
    env_noise_2 = np.random.normal(0.0, np.sqrt(1.0 - np.var(genetic_comp_2)), size=sample_size)
    expression_1 = genetic_comp_1 +env_noise_1
    expression_2 = genetic_comp_2 + env_noise_2

    return genotypes, genetic_comp_1, genetic_comp_2, expression_1, expression_2

def simulate_null_data(sample_size, snps_per_gene):
    # First simulate genotype data
    # Each snp has allele frequency drawn uniformly between 0.05 and 0.5
    allele_freqs = np.random.uniform(0.05, 0.5, size=snps_per_gene)
    genotypes = np.random.binomial(2, allele_freqs, size=(sample_size, snps_per_gene))

    # Standardize each snp to mean 0 and variance 1 across individuals
    genotypes = (genotypes - np.mean(genotypes, axis=0)) / np.std(genotypes, axis=0)

    # Simulate causal effects for two contexts
    # 5 causal snps per context: 3 shared across contexts, 2 unique to each context
    per_snp_h2 = 0.015
    shared_corr = 1.0
    n_shared = 7
    n_unique = 0
    causal_indices = np.random.choice(snps_per_gene, size=n_shared + 2*n_unique, replace=False)
    shared_indices = causal_indices[:n_shared]
    unique_indices_1 = causal_indices[n_shared:(n_shared + n_unique)]
    unique_indices_2 = causal_indices[(n_shared + n_unique):]

    causal_effects_1 = np.zeros(snps_per_gene)
    causal_effects_2 = np.zeros(snps_per_gene)

    # Shared causal snps: effects drawn from bivariate normal with correlation shared_corr
    shared_cov = per_snp_h2*np.asarray([[1.0, shared_corr], [shared_corr, 1.0]])
    shared_effects = np.random.multivariate_normal(np.zeros(2), shared_cov, size=n_shared)
    causal_effects_1[shared_indices] = shared_effects[:, 0]
    causal_effects_2[shared_indices] = shared_effects[:, 1]

    # Unique causal snps
    causal_effects_1[unique_indices_1] = np.random.normal(0.0, np.sqrt(per_snp_h2), size=n_unique)
    causal_effects_2[unique_indices_2] = np.random.normal(0.0, np.sqrt(per_snp_h2), size=n_unique)

    # Under the null, causal effects are identical across contexts (avoids numerical noise from singular covariance)
    causal_effects_2 = np.copy(causal_effects_1)
    

    # get genetic components
    genetic_comp_1 = genotypes.dot(causal_effects_1)
    genetic_comp_2 = genotypes.dot(causal_effects_2)

    # Add environmental noise
    env_noise_1 = np.random.normal(0.0, np.sqrt(1.0 - np.var(genetic_comp_1)), size=sample_size)
    env_noise_2 = np.random.normal(0.0, np.sqrt(1.0 - np.var(genetic_comp_2)), size=sample_size)
    expression_1 = genetic_comp_1 +env_noise_1
    expression_2 = genetic_comp_2 + env_noise_2

    return genotypes, genetic_comp_1, genetic_comp_2, expression_1, expression_2



####################
# Command line args
#####################
sample_size = int(sys.argv[1])
simulation_dir = sys.argv[2]



# Sim settings
n_genes = 10000
snps_per_gene = 500

###########################################################
# Run simulation with data with true changes between tissues
###########################################################
# Create output file
simulation_output_file = simulation_dir + 'no_ld_simulation_results_' + str(sample_size) + '.txt'
t = open(simulation_output_file,'w')
t.write("gene_iter\ttrue_delta_beta_var\ttrue_beta1_var\ttrue_beta2_var\test_delta_beta_var\test_delta_beta_var_se\test_delta_beta_var_ci_lb\test_delta_beta_var_ci_ub\test_delta_beta_var_p\test_delta_beta_var_jk_se\n")

for gene_iter in range(n_genes):
    # Generate simulated data
    genotypes, genetic_comp_1, genetic_comp_2, expression_1, expression_2 = simulate_data(sample_size, snps_per_gene)

    # Extract relevent simulated parameters
    true_delta_beta_var = np.var(genetic_comp_1 - genetic_comp_2)
    true_beta1_var = np.var(genetic_comp_1)
    true_beta2_var = np.var(genetic_comp_2)

    # Infer Var(\beta1 - \beta2)
    est_delta_beta_var, est_delta_beta_var_se, est_delta_beta_var_ci, est_delta_beta_var_p, est_delta_beta_var_jk_se = delta_beta_var_inference(expression_1, expression_2, genotypes)

    # print simulation results
    t.write(str(gene_iter) + "\t" + str(true_delta_beta_var) + "\t" + str(true_beta1_var) + "\t" + str(true_beta2_var) + "\t" + str(est_delta_beta_var) + "\t" + str(est_delta_beta_var_se) + "\t" + str(est_delta_beta_var_ci[0]) + "\t" + str(est_delta_beta_var_ci[1]) + "\t" + str(est_delta_beta_var_p) + "\t" + str(est_delta_beta_var_jk_se) + "\n")
t.close()


###########################################################
# Run Null simulation
###########################################################

# Create output file
simulation_output_file = simulation_dir + 'no_ld_null_simulation_results_' + str(sample_size) + '.txt'
t = open(simulation_output_file,'w')
t.write("gene_iter\ttrue_delta_beta_var\ttrue_beta1_var\ttrue_beta2_var\test_delta_beta_var\test_delta_beta_var_se\test_delta_beta_var_ci_lb\test_delta_beta_var_ci_ub\test_delta_beta_var_p\test_delta_beta_var_jk_se\n")

for gene_iter in range(n_genes):
    # Generate simulated data
    genotypes, genetic_comp_1, genetic_comp_2, expression_1, expression_2 = simulate_null_data(sample_size, snps_per_gene)

    # Extract relevent simulated parameters
    true_delta_beta_var = np.var(genetic_comp_1 - genetic_comp_2)
    true_beta1_var = np.var(genetic_comp_1)
    true_beta2_var = np.var(genetic_comp_2)

    # Infer Var(\beta1 - \beta2)
    est_delta_beta_var, est_delta_beta_var_se, est_delta_beta_var_ci, est_delta_beta_var_p, est_delta_beta_var_jk_se = delta_beta_var_inference(expression_1, expression_2, genotypes)

    # print simulation results
    t.write(str(gene_iter) + "\t" + str(true_delta_beta_var) + "\t" + str(true_beta1_var) + "\t" + str(true_beta2_var) + "\t" + str(est_delta_beta_var) + "\t" + str(est_delta_beta_var_se) + "\t" + str(est_delta_beta_var_ci[0]) + "\t" + str(est_delta_beta_var_ci[1]) + "\t" + str(est_delta_beta_var_p) + "\t" + str(est_delta_beta_var_jk_se) + "\n")
t.close()