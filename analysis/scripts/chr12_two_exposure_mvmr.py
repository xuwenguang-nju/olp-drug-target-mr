"""Two-exposure multivariable MR at chromosome 12: can GDF11 and SUOX be separated?

The conditional analysis (chr12_conditional_analysis.py) showed that each SNP is a strong
cis-eQTL for BOTH genes and that the two SNPs are correlated at r^2 = 0.87. This script asks
the question the reviewer actually poses: with GDF11 and SUOX expression entered as two
separate exposures and both SNPs as instruments, is either gene's effect on OLP estimable
independently of the other?

Model: b_out = B_exp^T theta, with theta = (theta_GDF11, theta_SUOX)
       B_exp is 2 SNPs x 2 genes; the system is exactly identified and its conditioning
       determines whether the separation is meaningful.
Weighting uses the joint outcome covariance from the LD matrix (R^-1 Sigma R^-T).

Diagnostics reported: the condition number of B_exp and the correlation between the two
instruments' fitted exposure profiles, both of which quantify collinearity.

Output: results/tables/corrected/chr12_two_exposure_mvmr.csv
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import ndtr

OUT = Path("results/tables/corrected")
SNPS = ["rs7302200", "rs773107"]
GENES = ["GDF11", "SUOX"]

ld = pd.read_csv("data/processed/coloc/ld_1000g/SUOX_ld.csv.gz", index_col=0)
R = ld.loc[SNPS, SNPS].to_numpy()
r = R[0, 1]

eq = pd.read_csv("data/processed/coloc/coloc_eqtl_windows.csv")
eq = eq[eq.snp.isin(SNPS) & eq.gene.isin(GENES)]
fg = pd.read_csv("data/processed/coloc/coloc_fg_wide.csv").set_index("snp")


def harmonise(gene, snp):
    row = eq[(eq.gene == gene) & (eq.snp == snp)].iloc[0]
    ea, oa, z, n = row.effect_allele, row.other_allele, float(row.zscore), float(row.NrSamples)
    f = fg.loc[snp]
    af = float(f.af_alt) if ea == f.alt else 1.0 - float(f.af_alt)
    maf = min(af, 1 - af)
    se = 1.0 / np.sqrt(2 * maf * (1 - maf) * n)
    s = 1.0 if ea == f.alt else -1.0
    return z * se, se, s * float(f.beta), float(f.sebeta)


# ---- build the matrices -----------------------------------------------------
B_exp = np.zeros((len(SNPS), len(GENES)))
SE_exp = np.zeros_like(B_exp)
b_out = np.zeros(len(SNPS))
se_out = np.zeros(len(SNPS))

for i, snp in enumerate(SNPS):
    for k, gene in enumerate(GENES):
        be, se, bo, so = harmonise(gene, snp)
        B_exp[i, k] = be
        SE_exp[i, k] = se
        b_out[i] = bo
        se_out[i] = so

Rinv = np.linalg.inv(R)
Sigma_out = R * np.outer(se_out, se_out)
B_out = Rinv @ b_out
C_out = Rinv @ Sigma_out @ Rinv.T

# exposure side: joint effects and their covariance (used for the SE of theta)
B_exp_joint = np.zeros_like(B_exp)
for k in range(len(GENES)):
    Sigma_e = R * np.outer(SE_exp[:, k], SE_exp[:, k])
    B_exp_joint[:, k] = Rinv @ B_exp[:, k]
    if k == 0:
        C_exp = np.zeros((len(GENES), len(SNPS), len(SNPS)))
    C_exp[k] = Rinv @ Sigma_e @ Rinv.T

# ---- MVMR -------------------------------------------------------------------
W = np.linalg.inv(C_out)
A = B_exp_joint
theta = np.linalg.solve(A.T @ W @ A, A.T @ W @ B_out)

# variance of theta: A^T W A  plus the exposure-uncertainty term
# Var(beta_joint_exp) contributes through delta method: theta = M B_out with M = (A^T W A)^-1 A^T W
M = np.linalg.inv(A.T @ W @ A) @ A.T @ W
var_theta = M @ C_out @ M.T
for k in range(len(GENES)):
    G = np.zeros_like(M)          # d theta / d B_exp_joint[k, :]
    G[k, :] = -(M @ B_out)
    var_theta = var_theta + G @ C_exp[k] @ G.T
se_theta = np.sqrt(np.diag(var_theta))

print("=== two-exposure MVMR (GDF11 and SUOX entered jointly) ===")
rows = []
for k, gene in enumerate(GENES):
    lo, hi = theta[k] - 1.96 * se_theta[k], theta[k] + 1.96 * se_theta[k]
    p = 2 * ndtr(-abs(theta[k] / se_theta[k]))
    rows.append(dict(gene=gene, logOR=round(theta[k], 3), se=round(se_theta[k], 3),
                     OR=round(np.exp(theta[k]), 3),
                     OR_lo=round(np.exp(lo), 3), OR_hi=round(np.exp(hi), 3),
                     p=f"{p:.3g}"))
    print(f"  {gene:6s} OR {np.exp(theta[k]):6.3f} "
          f"({np.exp(lo):6.3f}-{np.exp(hi):7.3f})  SE {se_theta[k]:.3f}  P = {p:.3g}")

# ---- collinearity diagnostics ----------------------------------------------
# NOTE: with two instruments the correlation of the two *columns* of the exposure matrix is
# always +/-1, so it carries no information. What matters is whether the SUOX/GDF11 effect
# ratio differs between the instruments: if it were the same, one linear combination of the
# two exposures would be unidentifiable.
ratio = B_exp[:, 1] / B_exp[:, 0]
print("\n=== separability of the two exposures ===")
print(f"  SUOX/GDF11 marginal effect ratio, {SNPS[0]}: {ratio[0]:.3f}")
print(f"  SUOX/GDF11 marginal effect ratio, {SNPS[1]}: {ratio[1]:.3f}")
print(f"  spread: {abs(ratio[1] - ratio[0]) / abs(ratio.mean()) * 100:.1f}% of the mean")
print(f"  condition number of the exposure matrix: {np.linalg.cond(B_exp_joint):.2f}")
print("  With two instruments and two exposures the model is exactly identified, so the")
print("  conditional F statistic of Sanderson et al. is not defined (no instrument remains")
print("  once the other exposure is conditioned on). Separation rests only on the ratio spread.")

# ---- what a one-exposure model says, for contrast ---------------------------
print("\n=== single-exposure MVMR (each gene alone, both SNPs) ===")
for k, gene in enumerate(GENES):
    A1 = B_exp_joint[:, [k]]
    th = np.linalg.solve(A1.T @ W @ A1, A1.T @ W @ B_out)[0]
    M1 = np.linalg.inv(A1.T @ W @ A1) @ A1.T @ W
    v = (M1 @ C_out @ M1.T)[0, 0]
    G = np.zeros_like(M1)
    G[0, :] = -(M1 @ B_out)
    v += (G @ C_exp[k] @ G.T)[0, 0]
    s = np.sqrt(v)
    print(f"  {gene:6s} OR {np.exp(th):6.3f} "
          f"({np.exp(th - 1.96 * s):6.3f}-{np.exp(th + 1.96 * s):7.3f})  SE {s:.3f}")

pd.DataFrame(rows).to_csv(OUT / "chr12_two_exposure_mvmr.csv", index=False)
print("\nwritten:", OUT / "chr12_two_exposure_mvmr.csv")
