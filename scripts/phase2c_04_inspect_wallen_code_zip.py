#!/usr/bin/env python3
from __future__ import annotations

import csv
from pathlib import Path
from zipfile import ZipFile


def flag_category(name: str) -> str:
    n = name.lower()
    tags = []
    if "meta" in n:
        tags.append("metadata")
    if "species" in n:
        tags.append("species")
    if "genus" in n:
        tags.append("genus")
    if "gene" in n or "family" in n:
        tags.append("gene_families")
    if "pathway" in n:
        tags.append("pathways")
    if "humann" in n:
        tags.append("HUMAnN")
    if "metaphlan" in n:
        tags.append("MetaPhlAn")
    if "stat" in n:
        tags.append("statistics")
    if "readme" in n:
        tags.append("README")
    return ";".join(tags) if tags else "none"


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    zpath = root / "data/phase2c/manual_sources/wallen_prjna834801/Supplementary_Code_24Oct2022.zip"
    out_dir = root / "data/phase2c/reports"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_tsv = out_dir / "phase2c_wallen_code_zip_inventory.tsv"
    out_status = out_dir / "phase2c_wallen_code_zip_status.txt"

    if not zpath.exists():
        with out_tsv.open("w", encoding="utf-8", newline="") as h:
            w = csv.DictWriter(
                h,
                fieldnames=["zip_entry", "file_size", "compressed_size", "category_flags"],
                delimiter="\t",
            )
            w.writeheader()
        out_status.write_text("MISSING_ZIP\n", encoding="utf-8")
        print(f"[DONE] wrote: {out_tsv}")
        print(f"[DONE] wrote: {out_status}")
        return 0

    rows = []
    with ZipFile(zpath, "r") as z:
        for zi in z.infolist():
            rows.append(
                {
                    "zip_entry": zi.filename,
                    "file_size": str(zi.file_size),
                    "compressed_size": str(zi.compress_size),
                    "category_flags": flag_category(zi.filename),
                }
            )

    with out_tsv.open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(
            h,
            fieldnames=["zip_entry", "file_size", "compressed_size", "category_flags"],
            delimiter="\t",
        )
        w.writeheader()
        for r in rows:
            w.writerow(r)

    out_status.write_text("ZIP_INVENTORY_READY\n", encoding="utf-8")
    print(f"[DONE] wrote: {out_tsv}")
    print(f"[DONE] wrote: {out_status}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
