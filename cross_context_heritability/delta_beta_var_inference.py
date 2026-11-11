import numpy as np
import scipy.stats


def he_regression_slope(n_pairs, sum_k, sum_kk, sum_y, sum_ky):
    # OLS slope (with intercept) of y_i*y_j on K_ij, computed from sufficient statistics
    return (sum_ky - sum_k*sum_y/n_pairs) / (sum_kk - np.square(sum_k)/n_pairs)


def delta_beta_var_inference(expression_1, expression_2, genotypes):
    delta_expr = expression_1 - expression_2
    delta_expr = delta_expr - np.mean(delta_expr)

    n_samples, n_snps = genotypes.shape

    # Haseman-Elston regression: regress delta_expr_i*delta_expr_j on GRM_ij over all pairs i != j
    # With GRM = GG^T/M, the slope estimates the genetic variance of delta_expr (= M*Var(beta1 - beta2))
    grm_full = genotypes.dot(genotypes.T) / n_snps
    grm = np.copy(grm_full)
    yy = np.outer(delta_expr, delta_expr)
    np.fill_diagonal(grm, 0.0)
    np.fill_diagonal(yy, 0.0)

    n_pairs = n_samples*(n_samples - 1)
    sum_k = np.sum(grm)
    sum_kk = np.sum(np.square(grm))
    sum_y = np.sum(yy)
    sum_ky = np.sum(grm*yy)
    est = he_regression_slope(n_pairs, sum_k, sum_kk, sum_y, sum_ky)

    # Leave-one-individual-out jackknife (dropping individual i removes row i and column i)
    jk_n_pairs = (n_samples - 1)*(n_samples - 2)
    jk_sum_k = sum_k - 2.0*np.sum(grm, axis=1)
    jk_sum_kk = sum_kk - 2.0*np.sum(np.square(grm), axis=1)
    jk_sum_y = sum_y - 2.0*np.sum(yy, axis=1)
    jk_sum_ky = sum_ky - 2.0*np.sum(grm*yy, axis=1)
    jk_est = he_regression_slope(jk_n_pairs, jk_sum_k, jk_sum_kk, jk_sum_y, jk_sum_ky)
    # Note: the delete-one jackknife overestimates the variance of this quadratic-form estimator (it roughly
    # double counts the pairwise noise term), so it is returned only for comparison
    est_jk_se = np.sqrt((n_samples - 1.0)/n_samples*np.sum(np.square(jk_est - np.mean(jk_est))))

    # Analytic SE: with intercept, est = y^T A y / sum(A^2), where A_ij = GRM_ij - mean off-diagonal GRM
    # (i != j) and A_ii = 0. For Gaussian y with covariance V, Var(y^T A y) = 2*tr(AVAV).
    # V is the plug-in fitted covariance: sigma_g^2 * GRM + sigma_e^2 * I
    a_mat = grm - sum_k/n_pairs
    np.fill_diagonal(a_mat, 0.0)
    # sigma_g^2 is plugged in unclamped below 0 (clamping at 0 biases the plug-in upward and inflates the SE);
    # it is only capped at var_y so that sigma_e^2 >= 0
    var_y = np.mean(np.square(delta_expr))
    sigma_g_sq = np.minimum(est, var_y)
    sigma_e_sq = var_y - sigma_g_sq
    v_mat = sigma_g_sq*grm_full + sigma_e_sq*np.eye(n_samples)
    av_mat = a_mat.dot(v_mat)
    var_quad = 2.0*np.sum(av_mat*av_mat.T)
    if var_quad <= 0.0:
        # Guard for extreme negative estimates: fall back to the null variance (sigma_g^2 = 0)
        var_quad = 2.0*np.square(var_y)*np.sum(np.square(a_mat))
    est_se = np.sqrt(var_quad) / np.sum(np.square(a_mat))

    est_ci = (est - 1.96*est_se, est + 1.96*est_se)
    # One-sided p-value (H1: genetic variance of delta_expr > 0)
    est_p = scipy.stats.norm.sf(est/est_se)

    return est, est_se, est_ci, est_p, est_jk_se
