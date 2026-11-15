"""Leave-one-donor-out sensitivity analysis for the single-cell comparisons.

The full-atlas analyses compare three NEOLP with two EOLP libraries, so the sample-level
permutation test over C(5,3) splits cannot fall below P = 0.1 and every estimate rests on
five libraries. This script re-estimates the two headline single-cell claims after
dropping one donor at a time:

  (1) family-level composition shift, EOLP - NEOLP in percentage points
      (mean of per-sample fractions, as in the manuscript);
  (2) per-subtype expression of the five candidate genes,
      frac_expr difference and mean-expressing-level difference (EOLP - NEOLP),
      restricted to subtypes with cells in every remaining donor.

Dropping one of the three NEOLP donors leaves a balanced 2-versus-2 split, for which all
C(4,2) = 6 regroupings can be enumerated. Because a 2-versus-2 split and its complement give
the same |difference|, each value enters twice and the smallest attainable two-sided P is
2/6 = 0.333, not 1/6. Dropping one of the two EOLP donors leaves 3-versus-1: the difference
is still estimable but no balanced regrouping exists, so no P value is reported.

Inputs  (read only):
  results/tables/corrected/geo_gse211630_composition_per_sample.csv
  results/tables/corrected/geo_gse211630_erosion_per_sample.csv
Outputs:
  results/tables/corrected/scrna_leave_one_donor_out_composition.csv
  results/tables/corrected/scrna_leave_one_donor_out_expression.csv
"""
import itertools
from pathlib import Path

import numpy as np
import pandas as pd

TAB = Path("results/tables/corrected")
SAMPLES = ["NEOLP-1", "NEOLP-2", "NEOLP-3", "EOLP-1", "EOLP-2"]
GENES = ["GDF11", "SUOX", "EGFL8", "TAPBP", "HLA-DPA1"]


def perm_p(vals, lab):
    """Two-sided permutation p over every equal-sized split of a 4-donor set."""
    vals = np.asarray(vals, float)
    lab = np.asarray(lab, int)
    if len(vals) != 4 or not (vals.max() > vals.min()):
        return np.nan
    obs = vals[lab == 1].mean() - vals[lab == 0].mean()
    d = []
    for ii in itertools.combinations(range(4), 2):
        ll = np.zeros(4, int)
        ll[list(ii)] = 1
        d.append(vals[ll == 1].mean() - vals[ll == 0].mean())
    return float((np.abs(np.asarray(d)) >= abs(obs)).mean())


def composition_loo():
    comp = pd.read_csv(TAB / "geo_gse211630_composition_per_sample.csv",
                       encoding="utf-8-sig")
    fam = comp[comp.level == "family"].copy()
    fam["pct"] = fam.frac * 100

    rows = []
    # full set for reference (all five donors)
    for fam_name, sub in fam.groupby("group"):
        n = sub[sub.condition == "NEOLP"]["pct"]
        e = sub[sub.condition == "EOLP"]["pct"]
        rows.append(dict(dataset="all5", dropped="-", family=fam_name,
                         mean_NEOLP=round(n.mean(), 2), mean_EOLP=round(e.mean(), 2),
                         delta_pct=round(e.mean() - n.mean(), 2),
                         n_NEOLP=len(n), n_EOLP=len(e), perm_p_4donor=np.nan,
                         perm_estimable=False))

    for drop in SAMPLES:
        keep = [s for s in SAMPLES if s != drop]
        for fam_name, sub in fam.groupby("group"):
            sub = sub[sub["sample"].isin(keep)]
            if sub.empty:
                continue
            n = sub[sub.condition == "NEOLP"]["pct"]
            e = sub[sub.condition == "EOLP"]["pct"]
            if len(n) == 0 or len(e) == 0:
                continue
            vals = np.concatenate([n.to_numpy(), e.to_numpy()])
            lab = np.concatenate([np.zeros(len(n), int), np.ones(len(e), int)])
            # a permutation test needs a balanced split; dropping an EOLP donor gives 3v1
            balanced = len(n) == len(e) == 2
            rows.append(dict(dataset=f"drop_{drop}", dropped=drop, family=fam_name,
                             mean_NEOLP=round(n.mean(), 2), mean_EOLP=round(e.mean(), 2),
                             delta_pct=round(e.mean() - n.mean(), 2),
                             n_NEOLP=len(n), n_EOLP=len(e),
                             perm_p_4donor=round(perm_p(vals, lab), 4) if balanced else np.nan,
                             perm_estimable=balanced))
    return pd.DataFrame(rows).sort_values(["family", "dataset"])


def expression_loo():
    er = pd.read_csv(TAB / "geo_gse211630_erosion_per_sample.csv", encoding="utf-8-sig")
    er = er[er.gene.isin(GENES)]

    rows = []
    for (gene, cluster, celltype), sub in er.groupby(["gene", "cluster", "celltype"]):
        # full-set difference
        n = sub[sub.condition == "NEOLP"]
        e = sub[sub.condition == "EOLP"]
        if len(n) < 3 or len(e) < 2:
            continue
        base_frac = e.frac_expr.mean() - n.frac_expr.mean()
        base_mean = e.mean_exp.mean() - n.mean_exp.mean()
        rows.append(dict(dataset="all5", dropped="-", gene=gene, cluster=int(cluster),
                         celltype=celltype,
                         delta_frac=round(base_frac, 4), delta_mean_exp=round(base_mean, 4),
                         perm_p_4donor=np.nan, sign_preserved=np.nan))

        for drop in SAMPLES:
            s2 = sub[sub["sample"] != drop]
            n2 = s2[s2.condition == "NEOLP"]
            e2 = s2[s2.condition == "EOLP"]
            # keep any subset with at least one donor per condition: the difference is
            # estimable, even when the split is unbalanced (3 vs 1) and the permutation
            # test therefore is not
            if len(n2) < 1 or len(e2) < 1:
                continue
            d_frac = e2.frac_expr.mean() - n2.frac_expr.mean()
            # order rows by donor so the permutation labels line up with the values
            e2 = e2.sort_values("sample")
            n2 = n2.sort_values("sample")
            vals = np.concatenate([n2.frac_expr.to_numpy(), e2.frac_expr.to_numpy()])
            lab = np.concatenate([np.zeros(len(n2), int), np.ones(len(e2), int)])
            # the permutation test needs an equal-sized split; with an unbalanced
            # leave-one-out subset (3 NEOLP vs 1 EOLP) it is not defined
            perm = round(perm_p(vals, lab), 4) if abs(len(n2) - len(e2)) == 0 else np.nan
            rows.append(dict(dataset=f"drop_{drop}", dropped=drop, gene=gene,
                             cluster=int(cluster), celltype=celltype,
                             delta_frac=round(d_frac, 4),
                             delta_mean_exp=round(e2.mean_exp.mean() - n2.mean_exp.mean(), 4),
                             perm_p_4donor=perm,
                             perm_estimable=abs(len(n2) - len(e2)) == 0,
                             sign_preserved=np.sign(d_frac) == np.sign(base_frac)))
    return pd.DataFrame(rows)


def main():
    comp = composition_loo()
    comp.to_csv(TAB / "scrna_leave_one_donor_out_composition.csv", index=False)
    print("=== composition: EOLP - NEOLP (percentage points) ===")
    piv = comp.pivot_table(index="family", columns="dataset", values="delta_pct")
    print(piv.round(2).to_string())

    focus = ["Endothelial", "T/NK", "Fibroblast", "Myeloid", "Pericyte/SMC"]
    print("\n--- headline families: does the sign survive dropping any donor? ---")
    for fam_name in focus:
        r = comp[comp.family == fam_name].set_index("dataset")
        base = r.loc["all5", "delta_pct"]
        loo = r.drop(index="all5")
        same = (np.sign(loo.delta_pct) == np.sign(base)).sum()
        print(f"  {fam_name:14s} base {base:+7.2f}  LOO range {loo.delta_pct.min():+7.2f} to "
              f"{loo.delta_pct.max():+7.2f}  sign kept {same}/{len(loo)}")

    expr = expression_loo()
    expr.to_csv(TAB / "scrna_leave_one_donor_out_expression.csv", index=False)
    print("\n=== candidate expression: sign stability of the frac_expr difference ===")
    rows = []
    for gene in GENES:
        sub = expr[(expr.gene == gene) & (expr.dataset != "all5")]
        per = sub.groupby(["cluster", "celltype"])["sign_preserved"].all()
        eolp = sub[sub.dataset.str.startswith("drop_EOLP")]
        neolp = sub[sub.dataset.str.startswith("drop_NEOLP")]
        rows.append(dict(
            gene=gene,
            subtypes=int((expr[(expr.gene == gene) & (expr.dataset == "all5")]).shape[0]),
            sign_stable_all_loo=int(per.sum()),
            sign_stable_drop_EOLP=int(eolp.groupby(["cluster", "celltype"])["sign_preserved"].all().sum()),
            n_drop_EOLP=int(eolp.groupby(["cluster", "celltype"]).ngroups),
            sign_stable_drop_NEOLP=int(neolp.groupby(["cluster", "celltype"])["sign_preserved"].all().sum()),
            n_drop_NEOLP=int(neolp.groupby(["cluster", "celltype"]).ngroups),
        ))
    print(pd.DataFrame(rows).to_string(index=False))

    print("\n--- subtypes the manuscript relies on ---")
    hosts = {"GDF11": "Fibroblast", "SUOX": "Fibroblast", "TAPBP": "T/NK",
             "EGFL8": "Pericyte", "HLA-DPA1": "Myeloid"}
    for gene, family in hosts.items():
        sub = expr[(expr.gene == gene) &
                   (expr.celltype.str.contains(family, case=False, na=False))]
        if sub.empty:
            print(f"  {gene} x {family}: no rows")
            continue
        for _, rb in sub[sub.dataset == "all5"].iterrows():
            one = sub[(sub.dataset != "all5") & (sub.cluster == rb.cluster)]
            e = one[one.dataset.str.startswith("drop_EOLP")]
            n = one[one.dataset.str.startswith("drop_NEOLP")]
            def rng(x):
                return f"{x.delta_frac.min():+.4f} to {x.delta_frac.max():+.4f}" if len(x) else "n/a"
            print(f"  {gene} x {rb.celltype}: base {rb.delta_frac:+.4f}")
            print(f"      drop-EOLP ({len(e)} sets, 3v1, no test): {rng(e)} | "
                  f"sign kept {int((np.sign(e.delta_frac)==np.sign(rb.delta_frac)).sum())}/{len(e)}")
            print(f"      drop-NEOLP ({len(n)} sets, 2v2):            {rng(n)} | "
                  f"sign kept {int((np.sign(n.delta_frac)==np.sign(rb.delta_frac)).sum())}/{len(n)}")

    print("\nwritten: scrna_leave_one_donor_out_expression.csv")


if __name__ == "__main__":
    main()
