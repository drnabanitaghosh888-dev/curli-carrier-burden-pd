from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXPECTED = {"csgA", "csgB", "csgC", "csgD", "csgE", "csgF", "csgG"}
CORE = {"csgA", "csgB", "csgD", "csgG"}
VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


def _parse_fasta(path: Path) -> list[dict[str, str]]:
    records = []
    header = None
    seq = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if header is not None:
                records.append({"header": header, "sequence": "".join(seq)})
            header = line[1:]
            seq = []
        else:
            seq.append(line.strip())
    if header is not None:
        records.append({"header": header, "sequence": "".join(seq)})
    return records


def _header_dict(h: str) -> dict[str, str]:
    out = {}
    for part in h.split("|"):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k] = v
    return out


def _read_meta(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8") as f:
        return [dict(r) for r in csv.DictReader(f, delimiter="\t")]


def test_mock_files_generated_and_have_seven_genes():
    subprocess.run(
        [sys.executable, "scripts/create_mock_curli_reference.py"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=True,
    )
    fasta = ROOT / "data/mock_reference/curli_mock_proteins.faa"
    meta = ROOT / "data/mock_reference/curli_mock_metadata.tsv"
    assert fasta.exists()
    assert meta.exists()
    recs = _parse_fasta(fasta)
    rows = _read_meta(meta)
    genes_fasta = {_header_dict(r["header"])["gene_symbol"] for r in recs}
    genes_meta = {r["gene_symbol"] for r in rows}
    assert genes_fasta == EXPECTED
    assert genes_meta == EXPECTED


def test_fasta_and_metadata_consistent_and_valid_symbols():
    fasta = ROOT / "data/mock_reference/curli_mock_proteins.faa"
    meta = ROOT / "data/mock_reference/curli_mock_metadata.tsv"
    recs = _parse_fasta(fasta)
    rows = _read_meta(meta)
    meta_by_gene = {r["gene_symbol"]: r for r in rows}
    accessions = set()
    for rec in recs:
        h = _header_dict(rec["header"])
        gene = h["gene_symbol"]
        acc = h["mock_accession"]
        seq = rec["sequence"]
        assert seq
        assert set(seq.upper()) <= VALID_AA
        assert gene in meta_by_gene
        m = meta_by_gene[gene]
        assert m["mock_accession"] == acc
        assert int(m["sequence_length"]) == len(seq)
        assert acc not in accessions
        accessions.add(acc)


def test_validator_and_scoring_work_for_complete_mock():
    cp = subprocess.run(
        [sys.executable, "scripts/validate_mock_curli_reference.py"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=True,
    )
    assert "validation passed" in cp.stdout.lower()

    out = ROOT / "data/mock_reference/.tmp_scores.json"
    try:
        subprocess.run(
            [sys.executable, "scripts/score_mock_curli_operon.py", "--out-json", str(out.relative_to(ROOT))],
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=True,
        )
        scores = json.loads(out.read_text(encoding="utf-8"))
        assert scores["Curli_Operon_Completeness"] == 1.0
        assert scores["Core_Curli_Score"] == 1.0
    finally:
        if out.exists():
            out.unlink()


def test_missing_one_gene_lowers_completeness_and_missing_core_lowers_core_score():
    fasta = ROOT / "data/mock_reference/curli_mock_proteins.faa"
    lines = fasta.read_text(encoding="utf-8").splitlines()
    # Drop csgC pair (non-core), then csgA pair (core), test via script output.
    def drop_gene(all_lines: list[str], gene: str) -> list[str]:
        out_lines = []
        i = 0
        while i < len(all_lines):
            line = all_lines[i]
            if line.startswith(">") and f"gene_symbol={gene}|" in line:
                i += 2
                continue
            out_lines.append(line)
            i += 1
        return out_lines

    tmp1 = ROOT / "data/mock_reference/.tmp_minus_csgC.faa"
    tmp2 = ROOT / "data/mock_reference/.tmp_minus_csgA.faa"
    out1 = ROOT / "data/mock_reference/.tmp_scores1.json"
    out2 = ROOT / "data/mock_reference/.tmp_scores2.json"
    try:
        tmp1.write_text("\n".join(drop_gene(lines, "csgC")) + "\n", encoding="utf-8")
        subprocess.run(
            [sys.executable, "scripts/score_mock_curli_operon.py", "--fasta", str(tmp1.relative_to(ROOT)), "--out-json", str(out1.relative_to(ROOT))],
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=True,
        )
        s1 = json.loads(out1.read_text(encoding="utf-8"))
        assert s1["Curli_Operon_Completeness"] < 1.0
        assert s1["Core_Curli_Score"] == 1.0

        tmp2.write_text("\n".join(drop_gene(lines, "csgA")) + "\n", encoding="utf-8")
        subprocess.run(
            [sys.executable, "scripts/score_mock_curli_operon.py", "--fasta", str(tmp2.relative_to(ROOT)), "--out-json", str(out2.relative_to(ROOT))],
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=True,
        )
        s2 = json.loads(out2.read_text(encoding="utf-8"))
        assert s2["Curli_Operon_Completeness"] < 1.0
        assert s2["Core_Curli_Score"] < 1.0
    finally:
        for p in [tmp1, tmp2, out1, out2]:
            if p.exists():
                p.unlink()


def test_no_heavy_external_command_invoked_via_dry_run():
    cp1 = subprocess.run(
        [sys.executable, "scripts/create_mock_curli_reference.py", "--dry-run"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=True,
    )
    cp2 = subprocess.run(
        [sys.executable, "scripts/validate_mock_curli_reference.py", "--dry-run"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=True,
    )
    cp3 = subprocess.run(
        [sys.executable, "scripts/score_mock_curli_operon.py", "--dry-run"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=True,
    )
    assert "[DRY-RUN]" in cp1.stdout
    assert "[DRY-RUN]" in cp2.stdout
    assert "[DRY-RUN]" in cp3.stdout
