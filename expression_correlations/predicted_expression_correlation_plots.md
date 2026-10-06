# Predicted vs observed expression correlation plots

Notation and definitions for the two "predicted expression correlation" calibration plots made by
`visualize_personalized_expression_correlations.R` from the per-gene output of
`personalized_expression_correlations_per_tissue.py`. Everything below is for one gene in one tissue
unless stated otherwise.

## Notation

| Symbol | Meaning |
|---|---|
| $n$ | number of individuals with expression in the tissue |
| $m$ | number of cis variants used for the gene |
| $X$ | $n \times m$ genotype matrix, each column standardized in sample to mean 0 and variance 1 |
| $R = X^\top X / n$ | in-sample LD matrix of the cis variants |
| $y$ | observed (residualized, renormalized) expression, length $n$, variance 1 |
| $b_j$ | Borzoi predicted effect of variant $j$, on the standardized-genotype scale |
| $s_{c(j)}$ | S-LDMC calibration slope of the Borzoi magnitude bin $c(j)$ that variant $j$ falls in |
| $\mu_j = s_{c(j)}\, b_j$ | rescaled Borzoi effect of variant $j$ |
| $\sigma_j^2$ | residual variance of the true effect around $\mu_j$ (two models, see below) |
| $\beta$ | true causal effects, length $m$ |
| $p = X\mu$ | rescaled predicted expression, length $n$ |

### Generative model

$$
\beta_j \sim N(\mu_j,\ \sigma_j^2) \ \text{independently},
\qquad
y = X\beta + \varepsilon,
\qquad
\operatorname{Var}(y_i) = 1 .
$$

$X\beta$ is the true genetic expression. $\varepsilon$ is independent of $X\beta$ and absorbs whatever is
needed for $y$ to have unit variance.

Two models for $\sigma_j^2$, giving the `_af_specific` twins of some columns:

- constant within Borzoi magnitude bin: $\sigma_j^2 = \text{resid\_var}_{c(j)}$
- allele-frequency specific: $\sigma_j^2 = 2 f_j (1-f_j)\, \tau^2_{c(j)}$, with $f_j$ the allele frequency

### Per-gene quantities written by the Python script

| Column | Definition |
|---|---|
| `rescaled_expression_correlation` | $r_{\text{obs}} = \operatorname{corr}(p, y)$, the realized correlation between the rescaled prediction and observed expression |
| `rescaled_predicted_expression_variance` | $V = \operatorname{Var}_i(p_i) = \mu^\top R \mu$, in-sample variance of the rescaled prediction |
| `predicted_cis_snp_h2` (and `_af_specific`) | $h^2_{\text{pred}} = \mathbb{E}[\operatorname{Var}_i(X\beta)] = V + \sum_j \sigma_j^2$ |
| `cis_snp_h2_he` | $\hat h^2$, Haseman-Elston estimate of the cis-SNP heritability of $y$ (unbiased, can be negative) |

## Two useful identities

**Expected covariance of the prediction with true genetic expression.**
Because $\mathbb{E}[\beta] = \mu$,

$$
\mathbb{E}\big[\operatorname{Cov}_i(p, X\beta)\big] = \mu^\top R\, \mathbb{E}[\beta] = \mu^\top R \mu = V .
$$

This does not involve $\sigma_j^2$. It is also the expected covariance with $y$, since $\varepsilon$ is independent of $p$.

**Observed-expression correlation is the genetic correlation shrunk by $\sqrt{h^2}$.**
With $y = X\beta + \varepsilon$ and $\operatorname{Var}(y) = 1$,

$$
\operatorname{corr}(p, y)
= \frac{\operatorname{Cov}(p, X\beta)}{\sqrt{\operatorname{Var}(p)}\cdot 1}
= \operatorname{corr}(p, X\beta)\, \sqrt{\operatorname{Var}(X\beta)}
= \operatorname{corr}(p, X\beta)\, \sqrt{h^2}.
$$

## Step by step: from Borzoi effects to $r_{\text{pred}}$ for one gene

These are the steps the Python script takes, in order. Steps 1 to 7 are shared by every quantity in the output.

**Step 1. Collect the cis variants.**
Take every variant with a Borzoi effect for this gene that is also in the genotype data. Call the Borzoi
per-allele effect of variant $j$ $b_j^{\text{raw}}$. Flip its sign if Borzoi's effect allele is not the allele the
genotype dosage counts, so that both refer to the same allele.

**Step 2. Build the genotype matrix.**
Take the dosages $g_{ij}$ (0, 1, 2) of variant $j$ in individual $i$ for the $n$ individuals with expression in this
tissue. Drop variants with any missing call or with no variance. Standardize each variant in sample:

$$
X_{ij} = \frac{g_{ij} - \bar g_j}{\operatorname{sd}(g_j)} ,
$$

so every column of $X$ has mean 0 and variance 1.

**Step 3. Put Borzoi effects on the standardized-genotype scale.**
An effect per allele copy becomes an effect per standard deviation of genotype:

$$
b_j = b_j^{\text{raw}} \cdot \operatorname{sd}(g_j) .
$$

This is the scale on which the S-LDMC slopes were estimated.

**Step 4. Assign each variant to a Borzoi magnitude bin.**
Bins are on $|b_j^{\text{raw}}|$ with edges $0, 0.001, 0.005, 0.01, 0.025, 0.05, 0.1, 0.2, 0.4, \infty$. Call the bin of
variant $j$ $c(j)$.

**Step 5. Rescale with the S-LDMC slope.**
Each bin $c$ has a slope $s_c$ from the S-LDMC results, the regression of true effects on Borzoi effects within that bin.
The rescaled effect is

$$
\mu_j = s_{c(j)} \, b_j .
$$

**Step 6. Form the prediction.**
One number per individual:

$$
p_i = \sum_{j=1}^{m} X_{ij}\, \mu_j, \qquad \text{i.e. } p = X\mu .
$$

**Step 7. Compute its variance across individuals.**

$$
V = \frac{1}{n}\sum_{i=1}^{n} (p_i - \bar p)^2 .
$$

This is the new column `rescaled_predicted_expression_variance`. Because the columns of $X$ are standardized,
$V = \mu^\top R \mu$ with $R$ the in-sample LD matrix.

**Step 8. Derive $r_{\text{pred}} = \sqrt{V}$.**

Setup. $p = X\mu$ is fixed (it depends only on genotypes and Borzoi). The data-generating model is

$$
y = X\beta + \varepsilon, \qquad \beta_j \sim N(\mu_j, \sigma_j^2) \text{ independently}, \qquad
\varepsilon \perp (X, \beta), \quad \mathbb{E}[\varepsilon] = 0, \qquad \operatorname{Var}_i(y_i) = 1 .
$$

All covariances and variances below are *across individuals* $i = 1, \dots, n$ (sample moments), and
$\mathbb{E}[\cdot]$ is the expectation over the random $\beta$ and $\varepsilon$.

*8.1 Write down what we want to predict.* The observed quantity is the sample correlation

$$
r_{\text{obs}} = \operatorname{corr}_i(p, y) = \frac{\operatorname{Cov}_i(p, y)}{\sqrt{\operatorname{Var}_i(p)\,\operatorname{Var}_i(y)}} .
$$

The predicted correlation is defined as the plug-in expectation, i.e. each moment replaced by its expectation under the model:

$$
r_{\text{pred}} := \frac{\mathbb{E}[\operatorname{Cov}_i(p, y)]}{\sqrt{\operatorname{Var}_i(p)\;\mathbb{E}[\operatorname{Var}_i(y)]}} .
$$

($\operatorname{Var}_i(p)$ needs no expectation because $p$ is fixed.)

*8.2 Centering.* Each column of $X$ has sample mean 0 (step 2). Hence $p = X\mu$ and $X\beta$ both have sample mean 0, and
sample covariances with them are plain averages of products:

$$
\operatorname{Cov}_i(p, X\beta) = \frac{1}{n}\sum_{i=1}^n p_i\,(X\beta)_i = \frac{1}{n}\,\mu^\top X^\top X\,\beta = \mu^\top R\,\beta,
\qquad R := \frac{X^\top X}{n} .
$$

*8.3 Split the covariance.* Substitute $y = X\beta + \varepsilon$ and use bilinearity:

$$
\operatorname{Cov}_i(p, y) = \operatorname{Cov}_i(p, X\beta) + \operatorname{Cov}_i(p, \varepsilon) = \mu^\top R\,\beta + \operatorname{Cov}_i(p, \varepsilon) .
$$

*8.4 Take the expectation over $\varepsilon$.* $\varepsilon$ is independent of $X$ with mean 0, so
$\mathbb{E}[\operatorname{Cov}_i(p, \varepsilon)] = 0$. Therefore

$$
\mathbb{E}[\operatorname{Cov}_i(p, y)] = \mathbb{E}[\mu^\top R\,\beta] .
$$

*8.5 Take the expectation over $\beta$.* $\mu^\top R\,\beta$ is linear in $\beta$ and $\mathbb{E}[\beta] = \mu$, so

$$
\mathbb{E}[\mu^\top R\,\beta] = \mu^\top R\,\mathbb{E}[\beta] = \mu^\top R\,\mu .
$$

Note that $\sigma_j^2$ never appears: the expectation of a linear function of $\beta$ depends only on the mean of $\beta$.

*8.6 Identify $\mu^\top R \mu$ as the variance of the prediction.* Again using that $p$ has mean 0,

$$
\operatorname{Var}_i(p) = \frac{1}{n}\sum_{i=1}^n p_i^2 = \frac{1}{n}\,\mu^\top X^\top X\,\mu = \mu^\top R\,\mu .
$$

Call this number $V$ (step 7). So steps 8.5 and 8.6 give the same number:

$$
\mathbb{E}[\operatorname{Cov}_i(p, y)] = V \qquad \text{and} \qquad \operatorname{Var}_i(p) = V .
$$

*8.7 Variance of expression.* The expression file is renormalized so that $\operatorname{Var}_i(y) = 1$, hence
$\mathbb{E}[\operatorname{Var}_i(y)] = 1$.

*8.8 Assemble.* Substitute 8.6 and 8.7 into the definition in 8.1:

$$
r_{\text{pred}} = \frac{V}{\sqrt{V \cdot 1}} = \frac{V}{\sqrt{V}} = \sqrt{V} .
$$

*Remarks.*
- $r_{\text{pred}}$ is the ratio of expected moments, not the expectation of the ratio. The two differ by a term of
  order $1/n$ that is ignored here.
- Only the mean of $\beta$ entered, so the constant-within-bin and allele-frequency-specific residual-variance models
  give the same $r_{\text{pred}}$.
- Interpretation: under the model, $p$ is a "calibrated" predictor of $y$ (regression slope of $y$ on $p$ equals 1 in
  expectation, since $\mathbb{E}[\operatorname{Cov}(p,y)] = \operatorname{Var}(p)$), and for a calibrated predictor the
  correlation with a unit-variance outcome is the predictor's standard deviation.

**Step 9. The observed counterpart.**
Replace the model with the data: $r_{\text{obs}} = \operatorname{corr}(p, y)$ computed from the same $p$ and the measured $y$.
This is `rescaled_expression_correlation`.

**Step 10 (Plot 2 only). Derive $r^{G}_{\text{pred}} = \sqrt{V / h^2_{\text{pred}}}$.**

Same setup as step 8, but the target is now the true genetic expression $X\beta$ instead of the observed expression
$y$. Sample moments are across individuals; $\mathbb{E}[\cdot]$ is over $\beta$.

*10.1 Write down what we want to predict.* The genetic-scale correlation is

$$
\operatorname{corr}_i(p, X\beta) = \frac{\operatorname{Cov}_i(p, X\beta)}{\sqrt{\operatorname{Var}_i(p)\,\operatorname{Var}_i(X\beta)}} ,
$$

and, as in 8.1, the predicted value is the plug-in expectation:

$$
r^{G}_{\text{pred}} := \frac{\mathbb{E}[\operatorname{Cov}_i(p, X\beta)]}{\sqrt{\operatorname{Var}_i(p)\;\mathbb{E}[\operatorname{Var}_i(X\beta)]}} .
$$

*10.2 Numerator.* By 8.2 and 8.5, $\mathbb{E}[\operatorname{Cov}_i(p, X\beta)] = \mu^\top R\,\mu = V$. (No $\varepsilon$ to
remove this time.)

*10.3 Variance of the prediction.* By 8.6, $\operatorname{Var}_i(p) = V$.

*10.4 Expected variance of true genetic expression.* Write $\beta = \mu + \delta$ with
$\delta_j \sim N(0, \sigma_j^2)$ independent across $j$, so $X\beta = X\mu + X\delta = p + X\delta$. Since $X\delta$ has
sample mean 0 (centered columns),

$$
\operatorname{Var}_i(X\beta) = \frac{1}{n}\,\beta^\top X^\top X\,\beta = \beta^\top R\,\beta
= \mu^\top R\,\mu + 2\,\mu^\top R\,\delta + \delta^\top R\,\delta .
$$

Take the expectation over $\delta$ term by term:

- $\mathbb{E}[\mu^\top R \mu] = V$ (constant);
- $\mathbb{E}[2\,\mu^\top R\,\delta] = 0$ since $\mathbb{E}[\delta] = 0$;
- $\mathbb{E}[\delta^\top R\,\delta] = \operatorname{tr}\!\big(R\,\mathbb{E}[\delta\delta^\top]\big)
  = \operatorname{tr}\!\big(R\,\operatorname{diag}(\sigma^2)\big) = \sum_j R_{jj}\,\sigma_j^2 = \sum_j \sigma_j^2$,
  because $R_{jj} = 1$ for standardized columns and the $\delta_j$ are independent, so the off-diagonal LD terms drop out.

Hence

$$
\mathbb{E}[\operatorname{Var}_i(X\beta)] = V + \sum_{j=1}^m \sigma_j^2 =: h^2_{\text{pred}} .
$$

This is `predicted_cis_snp_h2` (constant-within-bin $\sigma_j^2$) or `predicted_cis_snp_h2_af_specific`
(allele-frequency-specific $\sigma_j^2$). It is the expected cis-SNP heritability because $y$ has unit variance, so the
variance of its genetic component *is* its heritability.

*10.5 Assemble.* Substitute 10.2, 10.3 and 10.4 into 10.1:

$$
r^{G}_{\text{pred}} = \frac{V}{\sqrt{V \cdot h^2_{\text{pred}}}} = \sqrt{\frac{V}{h^2_{\text{pred}}}} .
$$

*Remarks.*
- Squaring gives $(r^G_{\text{pred}})^2 = V / h^2_{\text{pred}} = V / (V + \sum_j \sigma_j^2)$: the share of the true
  genetic variance that the point prediction $p$ accounts for. It lies in $[0, 1]$ because $\sum_j \sigma_j^2 \ge 0$.
- Unlike $r_{\text{pred}}$, this depends on the $\sigma_j^2$ through the denominator, so the two residual-variance
  models give different values and Plot 2 is made once per model.
- The two predictions are linked by $r_{\text{pred}} = r^{G}_{\text{pred}} \cdot \sqrt{h^2_{\text{pred}}}$, the model
  version of the identity $\operatorname{corr}(p, y) = \operatorname{corr}(p, X\beta)\,\sqrt{h^2}$.

**Step 11 (Plot 2 only). The observed counterpart of $r^G_{\text{pred}}$.**

$X\beta$ is not observed, so $\operatorname{corr}_i(p, X\beta)$ cannot be computed per gene. Instead use the identity
$\operatorname{corr}(p, y) = \operatorname{corr}(p, X\beta)\,\sqrt{h^2}$ (second identity above) to solve for it:

$$
\operatorname{corr}(p, X\beta) = \frac{r_{\text{obs}}}{\sqrt{h^2}} .
$$

Per gene, $\hat h^2$ (Haseman-Elston) is noisy and can be negative, so the ratio is formed at the bin level:
mean $r_{\text{obs}}$ over the bin, divided by the square root of the mean $\hat h^2$ over the bin. See Plot 2 below.

## Plot 1: observed-expression scale

File: `five_tissue_predicted_expression_r_calibration_per_tissue_panels*.pdf`

**Predicted correlation with observed expression.** Using the first identity, $\operatorname{Var}(p) = V$ and
$\operatorname{Var}(y) = 1$:

$$
r_{\text{pred}} = \frac{V}{\sqrt{V \cdot 1}} = \sqrt{V}
\qquad (\texttt{predicted\_expression\_r}).
$$

It does not depend on $\sigma_j^2$, so there is a single version.

**Observed counterpart.** $r_{\text{obs}}$ directly.

**Binning, per tissue.** Rank genes by $r_{\text{pred}}$ and split into $B = 10$ equal-count bins. For bin $b$ with $n_b$ genes:

$$
x_b = \frac{1}{n_b}\sum_{g \in b} r_{\text{pred}, g},
\qquad
y_b = \frac{1}{n_b}\sum_{g \in b} r_{\text{obs}, g},
$$

each with a 95% interval $\pm 1.96 \cdot \text{sd}/\sqrt{n_b}$ over genes in the bin. The dashed line is $y = x$.

## Plot 2: genetic-expression scale

Files: `five_tissue_predicted_genetic_r_calibration_per_tissue_panels*.pdf` and the `_af_specific` version.

**Predicted correlation with true genetic expression.** Using the first identity, $\operatorname{Var}(p) = V$ and
$\mathbb{E}[\operatorname{Var}(X\beta)] = h^2_{\text{pred}}$:

$$
r^{G}_{\text{pred}} = \frac{V}{\sqrt{V\, h^2_{\text{pred}}}} = \sqrt{\frac{V}{h^2_{\text{pred}}}}
\qquad (\texttt{predicted\_genetic\_r}).
$$

Equivalently $(r^G_{\text{pred}})^2 = V / h^2_{\text{pred}}$ is the fraction of true genetic expression variance
captured by the point prediction. $\sigma_j^2$ enters through $h^2_{\text{pred}}$, so there is one version per
residual-variance model.

**Observed counterpart.** By the second identity, $\operatorname{corr}(p, X\beta) = r_{\text{obs}} / \sqrt{h^2}$.
Per-gene $\hat h^2$ is noisy and can be negative, so the square root is taken once per bin, of the bin-mean heritability.

**Binning, per tissue.** Rank genes by $r^G_{\text{pred}}$ into $B = 10$ equal-count bins. For bin $b$:

$$
x_b = \frac{1}{n_b}\sum_{g \in b} r^{G}_{\text{pred}, g},
\qquad
y_b = \frac{\bar r_b}{\sqrt{\bar h^2_b}},
\quad
\bar r_b = \frac{1}{n_b}\sum_{g \in b} r_{\text{obs}, g},
\quad
\bar h^2_b = \frac{1}{n_b}\sum_{g \in b} \hat h^2_g .
$$

The interval on $x_b$ is $\pm 1.96 \cdot \text{sd}/\sqrt{n_b}$. The interval on $y_b$ is $\pm 1.96 \cdot \text{SE}(y_b)$ with the
delta method for $g(a, c) = a / \sqrt{c}$ at $a = \bar r_b$, $c = \bar h^2_b$:

$$
\operatorname{Var}(y_b) \approx
\frac{\operatorname{Var}(\bar r_b)}{\bar h^2_b}
\;-\; \frac{\bar r_b\, \operatorname{Cov}(\bar r_b, \bar h^2_b)}{(\bar h^2_b)^2}
\;+\; \frac{\bar r_b^2\, \operatorname{Var}(\bar h^2_b)}{4\, (\bar h^2_b)^3},
$$

where the variances and covariance of the bin means are the gene-level sample variances and covariance divided by $n_b$.
A bin with $\bar h^2_b \le 0$ gets no point.

## How the two plots relate

Dividing Plot 1's axes by the corresponding $\sqrt{h^2}$ gives Plot 2: the x-axis by $\sqrt{h^2_{\text{pred}}}$ and the
y-axis by $\sqrt{\hat h^2}$. So if Plot 1 is on $y = x$ but Plot 2 is not, the discrepancy is between predicted and
observed heritability, not between predicted and observed correlation.

## Gene sets

Each plot is made twice: all analyzed genes (`_all_genes` suffix) and heritable genes (cis-h2 LRT $p < 0.05$, no suffix).
Selecting on observed heritability inflates $\hat h^2$ relative to $h^2_{\text{pred}}$ (winner's curse), so the all-genes
versions are the unbiased ones.
