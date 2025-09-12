# lib/kb_gtdbtk/merge_utils.py
import os
import io
import glob
import json
from typing import List, Tuple, Dict

import pandas as pd

def _safe_read_tsv(path: str) -> pd.DataFrame:
    """
    Read a TSV defensively:
      - treat everything as string
      - tolerate empty files
      - strip BOM / weird encodings
      - do not crash on bad lines
    """
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return pd.DataFrame()
    with open(path, "rb") as fh:
        raw = fh.read()
    # Strip UTF-8 BOM if present
    if raw.startswith(b"\xef\xbb\xbf"):
        raw = raw[3:]
    bio = io.BytesIO(raw)
    try:
        return pd.read_csv(
            bio,
            sep="\t",
            dtype=str,
            keep_default_na=False,   # keep blanks as ''
            na_values=[],
            engine="python",
            on_bad_lines="skip"      # never raise on malformed rows
        )
    except Exception:
        # last resort: return empty; log handled by caller
        return pd.DataFrame()

def merge_tsvs_safe(
    paths: List[str],
    out_path: str,
    min_nonempty: int = 1,
    fill_value: str = ""
) -> Dict:
    """
    Merge multiple TSVs by taking the UNION of all columns.
    Missing files or empty files are tolerated.
    Returns a 'manifest' with diagnostics for logging/reporting.
    """
    manifest = {
        "inputs_requested": sorted(paths),
        "inputs_found": [],
        "empty_files": [],
        "ignored_missing": [],
        "columns_per_file": {},
        "union_columns": [],
        "rows_per_file": {},
        "written_to": out_path,
        "notes": []
    }

    dfs = []
    for p in paths:
        if not os.path.exists(p):
            manifest["ignored_missing"].append(p)
            continue
        df = _safe_read_tsv(p)
        manifest["inputs_found"].append(p)
        if df.empty:
            manifest["empty_files"].append(p)
            manifest["columns_per_file"][p] = []
            manifest["rows_per_file"][p] = 0
        else:
            manifest["columns_per_file"][p] = list(df.columns)
            manifest["rows_per_file"][p] = int(df.shape[0])
            dfs.append(df)

    if len(dfs) < min_nonempty:
        # Write an empty sentinel file so downstream steps don't explode
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        pd.DataFrame().to_csv(out_path, sep="\t", index=False)
        manifest["union_columns"] = []
        manifest["notes"].append(
            f"Fewer than {min_nonempty} non-empty TSVs; wrote empty output."
        )
        _write_manifest(out_path + ".manifest.json", manifest)
        return manifest

    # Build union of columns while preserving a reasonable order:
    # 1) order from the widest frame, then
    # 2) append any remaining columns in sorted order.
    widest = max(dfs, key=lambda d: d.shape[1])
    ordered = list(widest.columns)
    all_cols = set(ordered)
    for df in dfs:
        for c in df.columns:
            if c not in all_cols:
                ordered.append(c)
                all_cols.add(c)

    manifest["union_columns"] = ordered

    # Reindex each frame to the union, filling gaps with fill_value
    normed = []
    for df in dfs:
        # ensure all union columns present
        for c in ordered:
            if c not in df.columns:
                df[c] = fill_value
        normed.append(df[ordered])

    merged = pd.concat(normed, ignore_index=True, sort=False)

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    merged.to_csv(out_path, sep="\t", index=False)

    _write_manifest(out_path + ".manifest.json", manifest)
    return manifest

def _write_manifest(path: str, data: Dict):
    try:
        with open(path, "w") as fh:
            json.dump(data, fh, indent=2, sort_keys=False)
    except Exception:
        pass
