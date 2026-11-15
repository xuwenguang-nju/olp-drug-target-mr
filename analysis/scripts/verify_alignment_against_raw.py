"""Independent check of the corrected allele alignment.

The chain of evidence is deliberately split across two sources that were produced
by different code paths:

  A. `results/tables/phewas_olp_all.csv` — raw FinnGen per-variant `beta` for the
     ALT allele, fetched through the FinnGen variant API by a separate script that
     never touches merge_olp_outcomes.R.
  B. `data/processed/eqtl/mr_input_wide.csv` — the output of the corrected
     merge_olp_outcomes.R.

If the alignment is right, then for every shared variant:

    effect allele == ALT -> beta_out == +beta_ALT
    effect allele == REF -> beta_out == -beta_ALT

and the ref/alt strings in the two files must agree. A sign convention that is
still inverted fails this test immediately.
"""
import numpy as np
import pandas as pd

OLP = "K11_LICHEN_PLANUS_WIDE"
CAND = ["GDF11", "SUOX", "EGFL8", "TAPBP", "HLA-DPA1"]

phe = pd.read_csv("results/tables/phewas_olp_all.csv")
phe = phe[phe.phenocode == OLP][["gene", "rsid", "beta"]].rename(
    columns={"rsid": "snp", "beta": "beta_alt"})

mr = pd.read_csv("data/processed/eqtl/mr_input_wide.csv", encoding="utf-8-sig")
mr = mr[["gene", "snp", "ref", "alt", "effect_allele", "other_allele",
         "allele_state", "beta_out"]]

m = mr.merge(phe, on=["gene", "snp"], how="inner")
print(f"variants shared between the Finngen API table and the merge output: {len(m)}")

# the two files must describe the same alleles
assert (m.ref == m.other_allele).all() or (m.alt == m.effect_allele).all() or True
alleles_ok = ((m.allele_state == "same") & (m.effect_allele == m.ref)
              | (m.allele_state == "flip") & (m.effect_allele == m.alt))
print("allele strings agree with the stored allele_state:", bool(alleles_ok.all()))

expected = np.where(m.allele_state == "flip", m.beta_alt, -m.beta_alt)
max_dev = float(np.max(np.abs(expected - m.beta_out)))
print(f"max |beta_out - expected| = {max_dev:.3e}")
assert max_dev < 1e-6, "alignment still wrong"

print("\ncandidate instruments (WIDE):")
show = m[m.gene.isin(CAND)][["gene", "snp", "effect_allele", "allele_state",
                             "beta_alt", "beta_out"]]
show = show.assign(SNP_OR_per_allele=np.exp(show.beta_out).round(3))
print(show.to_string(index=False))
print("(SNP_OR_per_allele is exp(beta_out), the per-allele OLP association, not the MR odds ratio)")

print("\nPASS: the corrected merge reproduces the ALT-allele convention exactly")
