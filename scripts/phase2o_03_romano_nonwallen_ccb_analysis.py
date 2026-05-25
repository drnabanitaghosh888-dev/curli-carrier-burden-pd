from pathlib import Path
import re
import math
import numpy as np
import pandas as pd

ROOT = Path("data/phase2o/manual_sources/romano_pd_ml_meta_2025/zenodo_extract")
REPORT = Path("data/phase2o/reports")
REPORT.mkdir(parents=True, exist_ok=True)

META_MANIFEST = REPORT / "phase2o_romano_nonwallen_replication_manifest.tsv"
MOTUS = ROOT / "SMG_profiles/SMG/Taxonomy/mOTUs.txt"

OUT_BURDEN = REPORT / "phase2o_romano_nonwallen_ccb_sample_burden.tsv"
OUT_MATCHED_TAXA = REPORT / "phase2o_romano_nonwallen_ccb_matched_taxa.tsv"
OUT_STUDY_SUMMARY = REPORT / "phase2o_romano_nonwallen_ccb_study_summary.tsv"
OUT_TEST = REPORT / "phase2o_romano_nonwallen_ccb_test_summary.tsv"
OUT_STATUS = REPORT / "phase2o_romano_nonwallen_ccb_status.txt"
OUT_SUMMARY = REPORT / "phase2o_romano_nonwallen_ccb_summary.txt"

WEIGHTS = {
    "high": 1.0,
    "medium": 0.5,
    "low": 0.25,
}

FALLBACK_CURLI_TAXA = [
    ("Escherichia coli", "high"),
    ("Shigella flexneri", "high"),
    ("Shigella sonnei", "high"),
    ("Salmonella enterica", "high"),
    ("Salmonella bongori", "medium"),
    ("Citrobacter freundii", "medium"),
    ("Enterobacter cloacae", "medium"),
    ("Klebsiella pneumoniae", "medium"),
    ("Cronobacter sakazakii", "medium"),
    ("Serratia marcescens", "low"),
]

def norm_taxon(x):
    x = str(x)
    x = x.replace("_", " ")
    x = re.sub(r"\[[^\]]+\]", " ", x)
    x = re.sub(r"[^A-Za-z0-9 ]+", " ", x)
    x = re.sub(r"\s+", " ", x).strip().lower()
    return x

def find_curli_candidate_file():
    patterns = [
        "**/*curli*species*.tsv",
        "**/*curli*candidate*.tsv",
        "**/*candidate*curli*.tsv",
        "**/*curli*inventory*.tsv",
    ]
    hits = []
    for pat in patterns:
        hits.extend(Path("data").glob(pat))
    hits = sorted(set(hits))
    filtered = []
    for p in hits:
        ps = str(p).lower()
        if "phase2o_romano" in ps:
            continue
        if p.stat().st_size > 20_000_000:
            continue
        filtered.append(p)
    return filtered[0] if filtered else None

def load_curli_candidates():
    cand_file = find_curli_candidate_file()
    if cand_file is not None:
        try:
            df = pd.read_csv(cand_file, sep="\t", dtype=str)
            species_col = None
            conf_col = None
            for c in df.columns:
                cl = c.lower()
                if any(k in cl for k in ["species", "taxon", "organism", "name"]):
                    species_col = c
                    break
            for c in df.columns:
                cl = c.lower()
                if any(k in cl for k in ["confidence", "evidence", "tier", "weight"]):
                    conf_col = c
                    break
            if species_col is not None:
                rows = []
                for _, r in df.iterrows():
                    tax = str(r[species_col]).strip()
                    if not tax or tax.lower() == "nan":
                        continue
                    conf = "high"
                    if conf_col is not None:
                        raw = str(r[conf_col]).strip().lower()
                        if "medium" in raw:
                            conf = "medium"
                        elif "low" in raw:
                            conf = "low"
                        elif "high" in raw:
                            conf = "high"
                    rows.append((tax, conf))
                if rows:
                    return rows, str(cand_file)
        except Exception:
            pass

    return FALLBACK_CURLI_TAXA, "fallback_conservative_curli_taxa"

def mann_whitney_u(x, y):
    try:
        from scipy.stats import mannwhitneyu
        res = mannwhitneyu(x, y, alternative="two-sided")
        return float(res.statistic), float(res.pvalue)
    except Exception:
        return float("nan"), float("nan")

def cliffs_delta(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    gt = sum(np.sum(xi > y) for xi in x)
    lt = sum(np.sum(xi < y) for xi in x)
    return (gt - lt) / (len(x) * len(y))

def stratified_permutation_pvalue(df, value_col="log1p_CCB", label_col="PD", strata_col="Study", n_perm=5000, seed=123):
    rng = np.random.default_rng(seed)

    def weighted_stratified_diff(d):
        parts = []
        weights = []
        for study, g in d.groupby(strata_col):
            pdv = g.loc[g[label_col] == "PD", value_col].astype(float).values
            hcv = g.loc[g[label_col] == "HC", value_col].astype(float).values
            if len(pdv) == 0 or len(hcv) == 0:
                continue
            parts.append(np.mean(pdv) - np.mean(hcv))
            weights.append(len(g))
        if not parts:
            return np.nan
        return float(np.average(parts, weights=weights))

    obs = weighted_stratified_diff(df)
    if not np.isfinite(obs):
        return obs, np.nan, n_perm

    perm_stats = []
    for _ in range(n_perm):
        d = df.copy()
        new_labels = []
        for _, g in d.groupby(strata_col, sort=False):
            labels = g[label_col].values.copy()
            rng.shuffle(labels)
            new_labels.extend(labels)
        d[label_col] = new_labels
        perm_stats.append(weighted_stratified_diff(d))

    perm_stats = np.asarray(perm_stats, dtype=float)
    p = (np.sum(np.abs(perm_stats) >= abs(obs)) + 1) / (len(perm_stats) + 1)
    return obs, float(p), n_perm

# Load metadata and abundance.
meta = pd.read_csv(META_MANIFEST, sep="\t", dtype=str)
motus = pd.read_csv(MOTUS, sep="\t", dtype=str, index_col=0)

required_cols = {"SampleID", "PD", "Study"}
missing = required_cols - set(meta.columns)
if missing:
    raise SystemExit(f"[ERROR] Missing required metadata columns: {missing}")

samples = meta["SampleID"].astype(str).tolist()
missing_samples = [s for s in samples if s not in motus.columns]
if missing_samples:
    raise SystemExit(f"[ERROR] {len(missing_samples)} manifest samples absent from mOTUs table")

curli_candidates, candidate_source = load_curli_candidates()
candidate_df = pd.DataFrame(curli_candidates, columns=["candidate_taxon", "confidence"])
candidate_df["weight"] = candidate_df["confidence"].map(WEIGHTS).fillna(1.0)
candidate_df["candidate_norm"] = candidate_df["candidate_taxon"].map(norm_taxon)

# Match mOTUs taxa to curli candidates.
taxa_records = []
for taxon in motus.index.astype(str):
    nt = norm_taxon(taxon)
    best = None
    for _, c in candidate_df.iterrows():
        cn = c["candidate_norm"]
        if cn and cn in nt:
            if best is None or c["weight"] > best["weight"]:
                best = c
    if best is not None:
        taxa_records.append({
            "motus_taxon": taxon,
            "matched_candidate_taxon": best["candidate_taxon"],
            "confidence": best["confidence"],
            "weight": float(best["weight"]),
        })

matched = pd.DataFrame(taxa_records)

if matched.empty:
    OUT_STATUS.write_text("NO_CURLI_TAXA_MATCHED\n")
    raise SystemExit("[STATUS] NO_CURLI_TAXA_MATCHED")

# Calculate weighted burden.
abund = motus.loc[matched["motus_taxon"], samples].apply(pd.to_numeric, errors="coerce").fillna(0.0)
weights = matched.set_index("motus_taxon").loc[abund.index, "weight"].astype(float)
weighted_abund = abund.mul(weights, axis=0)

burden = weighted_abund.sum(axis=0).rename("CurliCarrierBurden").reset_index()
burden = burden.rename(columns={"index": "SampleID"})
burden = meta.merge(burden, on="SampleID", how="left")
burden["CurliCarrierBurden"] = burden["CurliCarrierBurden"].fillna(0.0)
burden["log1p_CCB"] = np.log1p(burden["CurliCarrierBurden"].astype(float))

# Overall summaries.
pd_vals = burden.loc[burden["PD"] == "PD", "CurliCarrierBurden"].astype(float).values
hc_vals = burden.loc[burden["PD"] == "HC", "CurliCarrierBurden"].astype(float).values
pd_log = burden.loc[burden["PD"] == "PD", "log1p_CCB"].astype(float).values
hc_log = burden.loc[burden["PD"] == "HC", "log1p_CCB"].astype(float).values

u_raw, p_raw = mann_whitney_u(pd_vals, hc_vals)
u_log, p_log = mann_whitney_u(pd_log, hc_log)
delta_raw = cliffs_delta(pd_vals, hc_vals)
delta_log = cliffs_delta(pd_log, hc_log)

obs_strat_diff, p_strat, n_perm = stratified_permutation_pvalue(
    burden,
    value_col="log1p_CCB",
    label_col="PD",
    strata_col="Study",
    n_perm=5000,
    seed=123
)

study_summary = (
    burden.groupby(["Study", "PD"])
    .agg(
        n=("SampleID", "count"),
        mean_CCB=("CurliCarrierBurden", "mean"),
        median_CCB=("CurliCarrierBurden", "median"),
        mean_log1p_CCB=("log1p_CCB", "mean"),
        median_log1p_CCB=("log1p_CCB", "median"),
    )
    .reset_index()
    .sort_values(["Study", "PD"])
)

test_rows = [{
    "analysis_set": "Romano_nonWallen_600",
    "candidate_source": candidate_source,
    "n_samples": len(burden),
    "n_PD": int((burden["PD"] == "PD").sum()),
    "n_HC": int((burden["PD"] == "HC").sum()),
    "n_studies": burden["Study"].nunique(),
    "n_matched_curli_motus_taxa": len(matched),
    "PD_median_CCB": float(np.median(pd_vals)),
    "HC_median_CCB": float(np.median(hc_vals)),
    "PD_mean_CCB": float(np.mean(pd_vals)),
    "HC_mean_CCB": float(np.mean(hc_vals)),
    "mannwhitney_raw_U": u_raw,
    "mannwhitney_raw_p": p_raw,
    "cliffs_delta_raw": delta_raw,
    "mannwhitney_log1p_U": u_log,
    "mannwhitney_log1p_p": p_log,
    "cliffs_delta_log1p": delta_log,
    "study_stratified_log1p_mean_difference_PD_minus_HC": obs_strat_diff,
    "study_stratified_permutation_p": p_strat,
    "n_permutations": n_perm,
}]

test_summary = pd.DataFrame(test_rows)

matched.to_csv(OUT_MATCHED_TAXA, sep="\t", index=False)
burden.to_csv(OUT_BURDEN, sep="\t", index=False)
study_summary.to_csv(OUT_STUDY_SUMMARY, sep="\t", index=False)
test_summary.to_csv(OUT_TEST, sep="\t", index=False)

status = "CCB_ANALYSIS_COMPLETE"
OUT_STATUS.write_text(status + "\n")

summary = f"""PHASE2O_ROMANO_NONWALLEN_CCB_ANALYSIS
status={status}
analysis_set=Romano_nonWallen_600
candidate_source={candidate_source}
n_samples={len(burden)}
n_PD={(burden["PD"] == "PD").sum()}
n_HC={(burden["PD"] == "HC").sum()}
n_studies={burden["Study"].nunique()}
n_matched_curli_motus_taxa={len(matched)}

PD_median_CCB={np.median(pd_vals)}
HC_median_CCB={np.median(hc_vals)}
PD_mean_CCB={np.mean(pd_vals)}
HC_mean_CCB={np.mean(hc_vals)}

mannwhitney_raw_p={p_raw}
mannwhitney_log1p_p={p_log}
cliffs_delta_raw={delta_raw}
cliffs_delta_log1p={delta_log}

study_stratified_log1p_mean_difference_PD_minus_HC={obs_strat_diff}
study_stratified_permutation_p={p_strat}
n_permutations={n_perm}

Interpretation guide:
A credible replication signal should show a consistent PD-vs-HC burden direction and should remain non-random under the study-stratified permutation test. The study-stratified result is more important than the pooled unadjusted Mann-Whitney test because Romano is a multi-cohort dataset.
"""

OUT_SUMMARY.write_text(summary)

print(f"[DONE] {OUT_BURDEN}")
print(f"[DONE] {OUT_MATCHED_TAXA}")
print(f"[DONE] {OUT_STUDY_SUMMARY}")
print(f"[DONE] {OUT_TEST}")
print(f"[DONE] {OUT_SUMMARY}")
print(f"[STATUS] {status}")
