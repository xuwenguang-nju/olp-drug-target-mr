"""Standard-error sensitivity analysis for the Wald-ratio MR estimates.

The manuscript's first-order standard error treats the exposure effect (b_exp) as
known, so the Wald Z equals Z_out x sign(Z_exp). That is what makes the reported
P values depend only on the outcome association. It also means the intervals do
not include uncertainty in b_exp.

This script recomputes the interval with the delta-method propagation of both
b_exp and b_out and reports how many Bonferroni hits survive.

    se(b) = sqrt( se_out^2 / b_exp^2 + b_out^2 * se_exp^2 / b_exp^4 )

No point estimate changes; only the standard error.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import ndtr

TAB = Path("results/tables/corrected")
OUT = Path("results/tables/corrected/mr_olp_se_propagation_sensitivity.csv")
BONF = 0.05 / 5040
CANDIDATES = ["GDF11", "SUOX", "EGFL8", "TAPBP", "HLA-DPA1"]


def analyse(tag):
    d = pd.read_csv(TAB / f"mr_olp_{tag}_results_corrected.csv")
    bx, bo = d.beta_exp, d.beta_out
    sx, so = d.se_exp, d.se_out

    se_known = so / bx.abs()
    se_full = np.sqrt(so ** 2 / bx ** 2 + bo ** 2 * sx ** 2 / bx ** 4)
    ratio = bo / bx

    d["se_known"] = se_known
    d["se_prop"] = se_full
    d["se_widen_pct"] = (se_full / se_known - 1) * 100
    d["OR_prop_lo"] = np.exp(ratio - 1.96 * se_full)
    d["OR_prop_hi"] = np.exp(ratio + 1.96 * se_full)
    d["p_prop"] = 2 * ndtr(-np.abs(ratio / se_full))
    d["bonf_known"] = d.p_MR < BONF
    d["bonf_prop"] = d.p_prop < BONF
    d["cohort"] = tag
    return d


frames = []
for tag in ("wide", "inclavo"):
    d = analyse(tag)
    frames.append(d)
    hits = d[d.bonf_known]
    kept = int(hits.bonf_prop.sum())
    print(f"{tag:8s} reported Bonferroni {len(hits):3d} -> retained {kept:3d} "
          f"(lost {len(hits) - kept})")
    lost = hits[~hits.bonf_prop]
    if len(lost):
        print("         lost:",
              ", ".join(f"{r.gene} (P {r.p_MR:.2g} -> {r.p_prop:.2g})"
                        for r in lost.itertuples()))

allr = pd.concat(frames)
keep = (["cohort", "gene", "snp", "OR", "OR_lo", "OR_hi", "p_MR", "se_known",
         "OR_prop_lo", "OR_prop_hi", "p_prop", "se_widen_pct", "bonf_known",
         "bonf_prop"])
allr[keep].to_csv(OUT, index=False)
print("\nwritten:", OUT, f"({len(allr)} rows)")

print("\n--- the five candidates ---")
sel = allr[allr.gene.isin(CANDIDATES)][
    ["cohort", "gene", "OR", "OR_lo", "OR_hi", "OR_prop_lo", "OR_prop_hi",
     "p_MR", "p_prop", "se_widen_pct", "bonf_prop"]]
print(sel.to_string(index=False, float_format=lambda v: f"{v:.4g}"))

w = allr[allr.cohort == "wide"]
print("\nWIDE se widening: min %.1f%% max %.1f%% median %.1f%%"
      % (w.se_widen_pct.min(), w.se_widen_pct.max(), w.se_widen_pct.median()))
