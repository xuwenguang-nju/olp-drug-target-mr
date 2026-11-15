"""Steiger directionality for the strict oral OLP endpoint.

Same estimator as run_steiger_test.py (correlations reconstructed from Z and
effective sample sizes, Fisher z comparison, two-sided erfc tail), applied to
K11_ORAL_LICHEN_PLANUS (920 cases / 499,266 controls).

Output: results/tables/corrected/steiger_oral_strict.csv
"""
import math
import pathlib

import numpy as np
import pandas as pd

ROOT = pathlib.Path(".")
CASES, CTRLS = 920, 499266
N_EFF = 4 * CASES * CTRLS / (CASES + CTRLS)

eq = pd.read_csv(ROOT / "data/processed/coloc/coloc_eqtl_windows.csv")
oral = pd.read_csv(ROOT / "results/tables/corrected/mr_olp_oral_results_corrected.csv")
IVS = {"GDF11": "rs7302200", "SUOX": "rs773107", "EGFL8": "rs41268896",
       "TAPBP": "rs3106190", "HLA-DPA1": "rs114390018"}

rows = []
for g, rsid in IVS.items():
    e = eq[(eq.gene == g) & (eq.snp == rsid)]
    o = oral[(oral.gene == g)]
    if e.empty or o.empty:
        print("skip", g)
        continue
    z_exp = float(e.iloc[0]["zscore"])
    n_exp = float(e.iloc[0]["NrSamples"])
    r_exp = math.sqrt(z_exp ** 2 / (z_exp ** 2 + n_exp))

    z_out = float(o.iloc[0]["Z_MR"])          # = Z_out * sign(Z_exp); magnitude is |Z_out|
    r_out = math.sqrt(z_out ** 2 / (z_out ** 2 + N_EFF))

    zf1 = 0.5 * math.log((1 + r_exp) / (1 - r_exp))
    zf2 = 0.5 * math.log((1 + r_out) / (1 - r_out))
    se_diff = math.sqrt(1 / (n_exp - 3) + 1 / (N_EFF - 3))
    z_diff = (zf1 - zf2) / se_diff
    p = math.erfc(abs(z_diff) / math.sqrt(2))
    rows.append(dict(gene=g, iv=rsid, outcome="K11_ORAL_LICHEN_PLANUS",
                     n_cases=CASES, z_exp=z_exp, z_out=z_out,
                     r_exp=r_exp, r_out=r_out,
                     direction="expression->outcome" if r_exp > r_out else "outcome->expression",
                     p_steiger=p))

res = pd.DataFrame(rows)
out = ROOT / "results/tables/corrected/steiger_oral_strict.csv"
res.to_csv(out, index=False)
print(res.to_string(index=False, float_format=lambda v: f"{v:.4g}"))
print("\nwritten:", out)
passed = (res.p_steiger < 0.05) & (res.direction == "expression->outcome")
print("Steiger-supported (P<0.05, expression->outcome):", int(passed.sum()), "/", len(res))
