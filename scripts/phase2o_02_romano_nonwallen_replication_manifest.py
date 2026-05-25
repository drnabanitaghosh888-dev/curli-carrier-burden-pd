from pathlib import Path
import pandas as pd

ROOT = Path("data/phase2o/manual_sources/romano_pd_ml_meta_2025/zenodo_extract")
OUT = Path("data/phase2o/reports")
OUT.mkdir(parents=True, exist_ok=True)

META = ROOT / "Metadata_All/Metadata/metadata.smg.txt"
MOTUS = ROOT / "SMG_profiles/SMG/Taxonomy/mOTUs.txt"

meta = pd.read_csv(META, sep="\t", dtype=str)
motus = pd.read_csv(MOTUS, sep="\t", dtype=str, index_col=0)

required = {"SampleID", "PD", "Study"}
missing = required - set(meta.columns)
if missing:
    raise SystemExit(f"[ERROR] Missing required metadata columns: {missing}")

# Keep only samples present in both metadata and mOTUs table.
motus_samples = set(motus.columns.astype(str))
meta["in_motus"] = meta["SampleID"].astype(str).isin(motus_samples)

# Exclude Wallen_2022 because it was already used as discovery cohort.
rep_meta = meta[(meta["Study"] != "Wallen_2022") & (meta["in_motus"])].copy()

# Basic summaries.
study_counts = (
    rep_meta.groupby(["Study", "PD"])
    .size()
    .reset_index(name="n_samples")
    .sort_values(["Study", "PD"])
)

overall_counts = (
    rep_meta["PD"]
    .value_counts()
    .rename_axis("PD")
    .reset_index(name="n_samples")
)

all_study_counts = (
    meta.groupby(["Study", "PD"])
    .size()
    .reset_index(name="n_samples")
    .sort_values(["Study", "PD"])
)

# Write outputs.
rep_manifest = OUT / "phase2o_romano_nonwallen_replication_manifest.tsv"
study_counts_path = OUT / "phase2o_romano_nonwallen_study_pd_counts.tsv"
overall_counts_path = OUT / "phase2o_romano_nonwallen_overall_pd_counts.tsv"
all_counts_path = OUT / "phase2o_romano_all_study_pd_counts.tsv"
status_path = OUT / "phase2o_romano_nonwallen_replication_status.txt"
summary_path = OUT / "phase2o_romano_nonwallen_replication_summary.txt"

rep_meta.to_csv(rep_manifest, sep="\t", index=False)
study_counts.to_csv(study_counts_path, sep="\t", index=False)
overall_counts.to_csv(overall_counts_path, sep="\t", index=False)
all_study_counts.to_csv(all_counts_path, sep="\t", index=False)

n_total = len(meta)
n_in_motus = int(meta["in_motus"].sum())
n_wallen = int(((meta["Study"] == "Wallen_2022") & meta["in_motus"]).sum())
n_rep = len(rep_meta)
n_pd = int((rep_meta["PD"] == "PD").sum())
n_hc = int((rep_meta["PD"] == "HC").sum())

ready = (
    n_total == 1324
    and n_in_motus == 1324
    and n_wallen == 724
    and n_rep == 600
    and n_pd > 0
    and n_hc > 0
)

status = "READY_FOR_NONWALLEN_CCB_ANALYSIS" if ready else "CHECK_REQUIRED_BEFORE_CCB"

summary = f"""PHASE2O_ROMANO_NONWALLEN_REPLICATION_MANIFEST
metadata_file={META}
motus_file={MOTUS}
n_total_metadata_samples={n_total}
n_samples_present_in_motus={n_in_motus}
n_wallen_samples_excluded={n_wallen}
n_nonwallen_replication_samples={n_rep}
n_nonwallen_PD={n_pd}
n_nonwallen_HC={n_hc}
status={status}

Interpretation:
Romano metadata and mOTUs are fully mapped. Wallen_2022 samples are excluded to avoid discovery-replication overlap. The remaining non-Wallen samples form the independent Romano replication subset for Curli Carrier Burden analysis.
"""

summary_path.write_text(summary)
status_path.write_text(status + "\n")

print(f"[DONE] {rep_manifest}")
print(f"[DONE] {study_counts_path}")
print(f"[DONE] {overall_counts_path}")
print(f"[DONE] {summary_path}")
print(f"[STATUS] {status}")
