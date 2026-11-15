"""Conditional analysis at the chromosome-12 locus: can GDF11 and SUOX be separated?

Two instruments, rs7302200 (GDF11 top cis-eQTL) and rs773107 (SUOX top cis-eQTL), lie
79,929 bp apart with r = 0.934 (r^2 = 0.87) in 1000G EUR. Each SNP is a strong cis-eQTL for
BOTH genes, so a single-instrument Wald ratio cannot attribute the OLP association to one gene.

Three analyses, all two-SNP and exactly identified:

1. Conditional (COJO-style) statistics: the association of each SNP with each gene's
   expression and with OLP, after conditioning on the other SNP.
2. Multivariable MR (MVMR) per gene, using both SNPs as instruments, with the joint
   effect covariance derived from the marginal standard errors and the LD matrix.
3. The conditioning diagnostic: how much of each single-instrument estimate survives.

Marginal effects are modelled as b ~ N(B, Sigma) with Sigma_ij = r_ij * se_i * se_j, so
joint effects are R^-1 b and their covariance R^-1 Sigma R^-T.

Output: results/tables/corrected/chr12_conditional_analysis.csv
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import ndtr

OUT = Path("results/tables/corrected")
SNPS = ["rs7302200", "rs773107"]           # 0 = GDF11 IV, 1 = SUOX IV
GENES = ["GDF11", "SUOX"]

# ---------------------------------------------------------------- LD
ld = pd.read_csv("data/processed/coloc/ld_1000g/SUOX_ld.csv.gz", index_col=0)
R = ld.loc[SNPS, SNPS].to_numpy()
assert abs(R[0, 0] - 1) < 1e-6 and abs(R[1, 1] - 1) < 1e-6
r = R[0, 1]
print(f"r = {r:.4f}   r2 = {r**2:.4f}   1 - r2 = {1 - r**2:.4f}")

# ---------------------------------------------------------------- eQTL
eq = pd.read_csv("data/processed/coloc/coloc_eqtl_windows.csv")
eq = eq[eq.snp.isin(SNPS) & eq.gene.isin(GENES)]

# ---------------------------------------------------------------- outcome
fg = pd.read_csv("data/processed/coloc/coloc_fg_wide.csv")
fg = fg[fg.snp.isin(SNPS)][["snp", "ref", "alt", "beta", "sebeta", "af_alt"]].set_index("snp")
print("\noutcome (WIDE) at the two SNPs:")
print(fg.to_string())


def alleles_for(gene, snp):
    row = eq[(eq.gene == gene) & (eq.snp == snp)].iloc[0]
    return row.effect_allele, row.other_allele, float(row.zscore), float(row.NrSamples)


def af_for(snp, effect_allele):
    """Allele frequency of `effect_allele` from FinnGen af_alt (af_alt = ALT frequency)."""
    row = fg.loc[snp]
    if effect_allele == row.alt:
        return float(row.af_alt)
    if effect_allele == row.ref:
        return 1.0 - float(row.af_alt)
    raise ValueError(f"{snp}: effect allele {effect_allele} not in {row.ref}/{row.alt}")


def outcome_effect(snp, effect_allele):
    """FinnGen beta/se expressed per `effect_allele` (FinnGen beta is the ALT effect)."""
    row = fg.loc[snp]
    s = 1.0 if effect_allele == row.alt else -1.0
    return s * float(row.beta), float(row.sebeta)


print("\n=== 1. marginal effects, harmonised to each gene's own effect allele ===")
rows = []
for gene in GENES:
    z = {}
    se_exp = {}
    for j, snp in enumerate(SNPS):
        ea, oa, zz, n = alleles_for(gene, snp)
        af = af_for(snp, ea)
        maf = min(af, 1 - af)
        se_exp[snp] = 1.0 / np.sqrt(2 * maf * (1 - maf) * n)
        z[snp] = zz
        b_out, se_out = outcome_effect(snp, ea)
        rows.append(dict(gene=gene, snp=snp, effect_allele=ea, af=af, maf=maf,
                         z_exp=zz, b_exp=zz * se_exp[snp], se_exp=se_exp[snp],
                         b_out=b_out, se_out=se_out))
marg = pd.DataFrame(rows)
print(marg[["gene", "snp", "effect_allele", "maf", "z_exp", "b_out", "se_out"]]
      .to_string(index=False, float_format=lambda v: f"{v:.4g}"))


def joint(b, se):
    """Joint effects and their covariance for two correlated estimates."""
    Rinv = np.linalg.inv(R)
    Sigma = R * np.outer(se, se)
    return Rinv @ b, Rinv @ Sigma @ Rinv.T


print("\n=== 2. joint (conditional) effects and MVMR ===")
res = []
for gene in GENES:
    g = marg[marg.gene == gene].set_index("snp").loc[SNPS]
    b_exp, se_exp = g.b_exp.to_numpy(), g.se_exp.to_numpy()
    b_out, se_out = g.b_out.to_numpy(), g.se_out.to_numpy()

    B_exp, C_exp = joint(b_exp, se_exp)
    B_out, C_out = joint(b_out, se_out)

    # conditional standardised statistics
    z_exp = b_exp / se_exp
    z_out = b_out / se_out
    z_exp_cond = [(z_exp[j] - r * z_exp[1 - j]) / np.sqrt(1 - r ** 2) for j in (0, 1)]
    z_out_cond = [(z_out[j] - r * z_out[1 - j]) / np.sqrt(1 - r ** 2) for j in (0, 1)]

    # MVMR (exactly identified with two instruments)
    W = np.linalg.inv(C_out)
    A = B_exp.reshape(-1, 1)
    theta = np.linalg.inv(A.T @ W @ A) @ A.T @ W @ B_out
    var_theta = np.linalg.inv(A.T @ W @ A)[0, 0]
    se_theta = np.sqrt(var_theta)
    p_mvmr = 2 * ndtr(-abs(theta[0] / se_theta))

    # single-instrument estimates for reference
    single = [(b_out[j] / b_exp[j],
               abs(b_out[j] / b_exp[j]) * np.sqrt((se_out[j] / b_out[j]) ** 2 +
                                                  (se_exp[j] / b_exp[j]) ** 2))
              for j in (0, 1)]

    res.append(dict(
        gene=gene, r=round(r, 4), r2=round(r ** 2, 4),
        z_exp_rs7302200=round(z_exp[0], 2), z_exp_rs773107=round(z_exp[1], 2),
        z_exp_cond_rs7302200=round(z_exp_cond[0], 2),
        z_exp_cond_rs773107=round(z_exp_cond[1], 2),
        z_out_rs7302200=round(z_out[0], 2), z_out_rs773107=round(z_out[1], 2),
        z_out_cond_rs7302200=round(z_out_cond[0], 2),
        z_out_cond_rs773107=round(z_out_cond[1], 2),
        OR_single_GDF11IV=round(np.exp(single[0][0]), 3),
        OR_single_SUOXIV=round(np.exp(single[1][0]), 3),
        MVMR_logOR=round(theta[0], 3), MVMR_se=round(se_theta, 3),
        MVMR_OR=round(np.exp(theta[0]), 3),
        MVMR_OR_lo=round(np.exp(theta[0] - 1.96 * se_theta), 3),
        MVMR_OR_hi=round(np.exp(theta[0] + 1.96 * se_theta), 3),
        MVMR_p=f"{p_mvmr:.3g}",
    ))

out = pd.DataFrame(res)
out.to_csv(OUT / "chr12_conditional_analysis.csv", index=False)
print(out.T.to_string())
print("\nwritten:", OUT / "chr12_conditional_analysis.csv")

print("\n=== 3. conditioning diagnostic ===")
for gene in GENES:
    sub = marg[(marg.gene == gene)].set_index("snp").loc[SNPS]
    z = (sub.z_exp / sub.se_exp).to_numpy() if False else sub.z_exp.to_numpy()
    zc = [(z[j] - r * z[1 - j]) / np.sqrt(1 - r ** 2) for j in (0, 1)]
    print(f"{gene}: z(rs7302200) {z[0]:7.2f} -> conditional {zc[0]:7.2f} | "
          f"z(rs773107) {z[1]:7.2f} -> conditional {zc[1]:7.2f}")
