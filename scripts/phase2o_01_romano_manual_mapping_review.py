from pathlib import Path
import re
import pandas as pd

ROOT = Path("data/phase2o/manual_sources/romano_pd_ml_meta_2025/zenodo_extract")
OUT = Path("data/phase2o/reports")
OUT.mkdir(parents=True, exist_ok=True)

META = ROOT / "Metadata_All/Metadata/metadata.smg.txt"
MOTUS = ROOT / "SMG_profiles/SMG/Taxonomy/mOTUs.txt"

summary_path = OUT / "phase2o_romano_mapping_summary.txt"
mapping_path = OUT / "phase2o_romano_sample_mapping.tsv"
meta_cols_path = OUT / "phase2o_romano_metadata_columns.tsv"
motus_cols_path = OUT / "phase2o_romano_motus_columns.tsv"
status_path = OUT / "phase2o_romano_mapping_status.txt"

def read_table(path):
    if not path.exists():
        raise FileNotFoundError(path)
    # Try tab first, then generic separator.
    try:
        return pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    except Exception:
        return pd.read_csv(path, sep=None, engine="python", dtype=str)

def norm(x):
    if pd.isna(x):
        return ""
    return str(x).strip()

def extract_accession(x):
    x = norm(x)
    m = re.search(r"\b[DES]RR\d+\b|\bSRR\d+\b|\bERR\d+\b|\bDRR\d+\b", x)
    return m.group(0) if m else ""

def infer_label(row):
    joined = " ".join(norm(v).lower() for v in row.values)
    # Conservative label inference.
    if re.search(r"\b(parkinson|pd|case|disease|patient)\b", joined):
        if not re.search(r"\b(control|healthy|hc|normal)\b", joined):
            return "PD_candidate"
    if re.search(r"\b(control|healthy|hc|normal)\b", joined):
        return "CONTROL_candidate"
    return "UNKNOWN"

def possible_id_columns(df):
    hits = []
    for c in df.columns:
        cl = c.lower()
        vals = df[c].dropna().astype(str).head(200).tolist()
        n_acc = sum(bool(extract_accession(v)) for v in vals)
        score = 0
        if any(k in cl for k in ["sample", "run", "accession", "id", "subject"]):
            score += 2
        if n_acc > 0:
            score += 3
        if score > 0:
            hits.append((c, score, n_acc))
    return sorted(hits, key=lambda x: (-x[1], x[0]))

def clean_motus_sample_columns(df):
    # mOTUs tables usually have feature/taxonomy columns followed by sample columns.
    non_sample_keywords = [
        "motu", "taxonomy", "tax", "consensus", "lineage", "species",
        "genus", "family", "order", "class", "phylum", "kingdom",
        "ncbi", "id", "name", "unnamed"
    ]
    cols = []
    for c in df.columns:
        cl = c.lower()
        if any(k in cl for k in non_sample_keywords):
            continue
        vals = pd.to_numeric(df[c], errors="coerce")
        numeric_fraction = vals.notna().mean()
        if numeric_fraction > 0.5:
            cols.append(c)
    return cols

def sample_key_variants(x):
    x = norm(x)
    variants = {x}
    acc = extract_accession(x)
    if acc:
        variants.add(acc)
    # Remove common table suffix/prefix punctuation.
    variants.add(re.sub(r"[^A-Za-z0-9]+", "", x))
    return {v for v in variants if v}

meta = read_table(META)
motus = read_table(MOTUS)

meta.to_csv(meta_cols_path, sep="\t", index=False)
pd.DataFrame({"motus_column": list(motus.columns)}).to_csv(motus_cols_path, sep="\t", index=False)

id_candidates = possible_id_columns(meta)
motus_sample_cols = clean_motus_sample_columns(motus)

# Choose the best metadata ID column.
best_id_col = id_candidates[0][0] if id_candidates else meta.columns[0]

meta2 = meta.copy()
meta2["__metadata_sample_id__"] = meta2[best_id_col].map(norm)
meta2["__metadata_accession__"] = meta2["__metadata_sample_id__"].map(extract_accession)
meta2["__disease_label_candidate__"] = meta2.apply(infer_label, axis=1)

# Build metadata lookup using exact IDs and accession IDs.
lookup = {}
for i, r in meta2.iterrows():
    for key in sample_key_variants(r["__metadata_sample_id__"]):
        lookup.setdefault(key, []).append(i)
    if r["__metadata_accession__"]:
        lookup.setdefault(r["__metadata_accession__"], []).append(i)

rows = []
for c in motus_sample_cols:
    matched_indices = set()
    for key in sample_key_variants(c):
        if key in lookup:
            matched_indices.update(lookup[key])
    if matched_indices:
        for idx in matched_indices:
            rows.append({
                "motus_sample_column": c,
                "metadata_row_index": idx,
                "metadata_id_column_used": best_id_col,
                "metadata_sample_id": meta2.loc[idx, "__metadata_sample_id__"],
                "metadata_accession": meta2.loc[idx, "__metadata_accession__"],
                "disease_label_candidate": meta2.loc[idx, "__disease_label_candidate__"],
                "mapping_status": "matched"
            })
    else:
        rows.append({
            "motus_sample_column": c,
            "metadata_row_index": "",
            "metadata_id_column_used": best_id_col,
            "metadata_sample_id": "",
            "metadata_accession": extract_accession(c),
            "disease_label_candidate": "",
            "mapping_status": "unmatched"
        })

mapping = pd.DataFrame(rows)
mapping.to_csv(mapping_path, sep="\t", index=False)

n_motus_samples = len(motus_sample_cols)
n_matched = int((mapping["mapping_status"] == "matched").sum())
n_unmatched = int((mapping["mapping_status"] == "unmatched").sum())
label_counts = mapping.loc[mapping["mapping_status"] == "matched", "disease_label_candidate"].value_counts().to_dict()

ready = (
    n_motus_samples > 0
    and n_matched > 0
    and n_unmatched == 0
    and label_counts.get("PD_candidate", 0) > 0
    and label_counts.get("CONTROL_candidate", 0) > 0
)

status = "READY_FOR_ANALYSIS" if ready else "NO_MAPPING_REQUIRED_FIRST"

summary = [
    "PHASE2O_ROMANO_MAPPING_REVIEW",
    f"metadata_file={META}",
    f"motus_file={MOTUS}",
    f"metadata_shape={meta.shape[0]}x{meta.shape[1]}",
    f"motus_shape={motus.shape[0]}x{motus.shape[1]}",
    f"metadata_id_column_used={best_id_col}",
    f"n_metadata_id_column_candidates={len(id_candidates)}",
    f"n_motus_sample_columns={n_motus_samples}",
    f"n_matched_motus_samples={n_matched}",
    f"n_unmatched_motus_samples={n_unmatched}",
    f"label_counts={label_counts}",
    f"READY_STATUS={status}",
    "",
    "TOP_METADATA_ID_COLUMN_CANDIDATES",
]
summary += [f"{c}\tscore={s}\tn_accession_like_values={n}" for c, s, n in id_candidates[:20]]

summary_path.write_text("\n".join(summary) + "\n")
status_path.write_text(status + "\n")

print(f"[DONE] {summary_path}")
print(f"[DONE] {mapping_path}")
print(f"[DONE] {meta_cols_path}")
print(f"[DONE] {motus_cols_path}")
print(f"[STATUS] {status}")
