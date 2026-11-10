#!/usr/bin/env python3
"""
GPMV proteomics QC + UniProt annotation pipeline

Purpose
-------
1. Read a Proteome Discoverer protein export (.xlsx or .csv).
2. Apply the thesis QC criteria consistently:
   - Protein FDR Confidence: Combined == High
   - Contaminant == False
   - # Unique Peptides >= 2
3. Restrict differential-abundance analysis to proteins with a numeric
   GPMV/Cell abundance ratio and adjusted p-value.
4. Define:
   - GPMV-enriched: adjusted p < 0.05 and ratio > 1
   - GPMV-depleted: adjusted p < 0.05 and ratio < 1
   - Not significant: adjusted p >= 0.05
5. Optionally retrieve UniProtKB annotations for enriched + depleted proteins.
6. Assign broad, NON-MUTUALLY-EXCLUSIVE localization categories using
   UniProt subcellular-location and GO Cellular Component annotations.
7. Save filtered tables, annotation tables, summaries and a comparison plot.

The XLSX reader uses only the Python standard library and reads the first
worksheet. This avoids changing the original workbook.

Example
-------
python gpmv_proteomics_uniprot_pipeline.py proteome.xlsx

By default, results are saved in a folder called
GPMV_proteomics_results next to the input file.

Filtering-only test:
python gpmv_proteomics_uniprot_pipeline.py proteome.xlsx --skip-uniprot

Dependencies
------------
Filtering: Python standard library only.
UniProt retrieval: internet connection only (standard library urllib).
Plotting: matplotlib (optional; CSV summaries are still produced without it).
"""

from __future__ import annotations

import argparse
import csv
import io
import math
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
import zipfile
from collections import Counter
from pathlib import Path


RATIO_CANDIDATES = [
    "Abundance Ratio: (GMPVs) / (Cells)",  # exact header in current export
    "Abundance Ratio: (GPMVs) / (Cells)",
]
PVAL_CANDIDATES = [
    "Abundance Ratio Adj. P-Value: (GMPVs) / (Cells)",
    "Abundance Ratio Adj. P-Value: (GPMVs) / (Cells)",
]

COL_FDR = "Protein FDR Confidence: Combined"
COL_CONTAMINANT = "Contaminant"
COL_UNIQUE = "# Unique Peptides"
COL_ACCESSION = "Accession"
COL_DESCRIPTION = "Description"

UNIPROT_FIELDS = [
    "accession",
    "id",
    "reviewed",
    "protein_name",
    "gene_names",
    "organism_name",
    "cc_subcellular_location",
    "go_c",
    "ft_transmem",
    "ft_signal",
]

MAIN_NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"


def is_number(value):
    if isinstance(value, bool):
        return False
    if isinstance(value, (int, float)):
        return math.isfinite(float(value))
    if isinstance(value, str):
        value = value.strip().replace(",", ".")
        if not value:
            return False
        try:
            return math.isfinite(float(value))
        except ValueError:
            return False
    return False


def to_float(value):
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return float(value)
    return float(str(value).strip().replace(",", "."))


def excel_col_to_index(cell_ref):
    """A1 -> 0, B1 -> 1, AA1 -> 26."""
    m = re.match(r"([A-Z]+)", cell_ref or "")
    if not m:
        return None
    n = 0
    for ch in m.group(1):
        n = n * 26 + (ord(ch) - ord("A") + 1)
    return n - 1


def xlsx_cell_value(cell, shared_strings):
    cell_type = cell.attrib.get("t")

    if cell_type == "inlineStr":
        inline = cell.find(MAIN_NS + "is")
        if inline is None:
            return ""
        return "".join((t.text or "") for t in inline.iter(MAIN_NS + "t"))

    value_node = cell.find(MAIN_NS + "v")
    if value_node is None:
        return ""

    raw = value_node.text or ""

    if cell_type == "s":
        return shared_strings[int(raw)]
    if cell_type == "b":
        return raw == "1"
    if cell_type in (None, "n"):
        try:
            number = float(raw)
            return int(number) if number.is_integer() else number
        except ValueError:
            return raw

    return raw


def read_first_sheet_xlsx(path):
    """
    Read the first worksheet from an .xlsx without modifying it.
    Supports both ordinary A1-referenced cells and sequential cells such as
    the current Proteome Discoverer export.
    """
    with zipfile.ZipFile(path) as zf:
        shared_strings = []
        if "xl/sharedStrings.xml" in zf.namelist():
            root = ET.fromstring(zf.read("xl/sharedStrings.xml"))
            for si in root.findall(MAIN_NS + "si"):
                shared_strings.append(
                    "".join((t.text or "") for t in si.iter(MAIN_NS + "t"))
                )

        workbook = ET.fromstring(zf.read("xl/workbook.xml"))
        sheets_node = workbook.find(MAIN_NS + "sheets")
        if sheets_node is None or len(sheets_node) == 0:
            raise ValueError("No worksheet found in the XLSX file.")

        first_sheet = sheets_node[0]
        rel_id = first_sheet.attrib[
            "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"
        ]

        rels_root = ET.fromstring(zf.read("xl/_rels/workbook.xml.rels"))
        target = None
        for rel in rels_root:
            if rel.attrib.get("Id") == rel_id:
                target = rel.attrib.get("Target")
                break
        if not target:
            raise ValueError("Could not resolve the first worksheet.")

        target = target.lstrip("/")
        if not target.startswith("xl/"):
            target = "xl/" + target

        sheet_root = ET.fromstring(zf.read(target))
        sheet_data = sheet_root.find(MAIN_NS + "sheetData")
        if sheet_data is None:
            raise ValueError("The first worksheet contains no sheetData.")

        matrix = []
        for row in sheet_data.findall(MAIN_NS + "row"):
            values = []
            sequential_col = 0
            for cell in row.findall(MAIN_NS + "c"):
                ref_index = excel_col_to_index(cell.attrib.get("r"))
                col_index = ref_index if ref_index is not None else sequential_col

                while len(values) < col_index:
                    values.append("")
                value = xlsx_cell_value(cell, shared_strings)

                if len(values) == col_index:
                    values.append(value)
                else:
                    values[col_index] = value
                sequential_col = col_index + 1
            matrix.append(values)

    if not matrix:
        raise ValueError("The worksheet is empty.")

    headers = [str(x).strip() for x in matrix[0]]
    ncols = len(headers)
    records = []
    for row in matrix[1:]:
        row = list(row) + [""] * max(0, ncols - len(row))
        records.append(dict(zip(headers, row[:ncols])))
    return headers, records


def read_csv_table(path):
    with open(path, "r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        return reader.fieldnames or [], list(reader)


def read_table(path):
    suffix = Path(path).suffix.lower()
    if suffix == ".xlsx":
        return read_first_sheet_xlsx(path)
    if suffix == ".csv":
        return read_csv_table(path)
    raise ValueError("Input must be .xlsx or .csv")


def find_column(headers, candidates):
    for c in candidates:
        if c in headers:
            return c
    raise KeyError(
        "None of the expected columns were found:\n  - "
        + "\n  - ".join(candidates)
    )


def bool_false(value):
    if value in (False, 0, "", None):
        return True
    text = str(value).strip().lower()
    return text in {"false", "0", "no", "n/a", "na", ""}


def write_csv(path, rows, fieldnames=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    rows = list(rows)
    if fieldnames is None:
        fieldnames = list(rows[0].keys()) if rows else []

    with open(path, "w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def filter_proteins(records, ratio_col, pval_col, min_unique, adj_p):
    qc_passed = []
    qc_failed = []

    for row in records:
        reasons = []

        fdr = str(row.get(COL_FDR, "")).strip().lower()
        if fdr and fdr != "high":
            reasons.append("FDR_not_High")

        if not bool_false(row.get(COL_CONTAMINANT, "")):
            reasons.append("Contaminant")

        unique = row.get(COL_UNIQUE, "")
        if not is_number(unique) or to_float(unique) < min_unique:
            reasons.append(f"Unique_peptides_lt_{min_unique}")

        if reasons:
            out = dict(row)
            out["QC_failure_reason"] = ";".join(reasons)
            qc_failed.append(out)
        else:
            qc_passed.append(dict(row))

    quantifiable = []
    no_quant = []

    for row in qc_passed:
        ratio = row.get(ratio_col, "")
        pval = row.get(pval_col, "")
        if is_number(ratio) and is_number(pval):
            out = dict(row)
            out["Ratio_numeric"] = to_float(ratio)
            out["Adj_P_numeric"] = to_float(pval)
            quantifiable.append(out)
        else:
            out = dict(row)
            out["Quantification_status"] = "Missing/non-numeric ratio or adjusted p-value"
            no_quant.append(out)

    enriched, depleted, not_significant = [], [], []

    for row in quantifiable:
        ratio = row["Ratio_numeric"]
        pval = row["Adj_P_numeric"]
        out = dict(row)

        if pval < adj_p and ratio > 1:
            out["Differential_group"] = "GPMV-enriched"
            enriched.append(out)
        elif pval < adj_p and ratio < 1:
            out["Differential_group"] = "GPMV-depleted"
            depleted.append(out)
        else:
            out["Differential_group"] = "Not significant"
            not_significant.append(out)

    return {
        "qc_passed": qc_passed,
        "qc_failed": qc_failed,
        "quantifiable": quantifiable,
        "no_quant": no_quant,
        "enriched": enriched,
        "depleted": depleted,
        "not_significant": not_significant,
    }


def canonical_accession(accession):
    """
    For localization annotation, map an isoform accession such as Q9UHI5-2
    to the canonical UniProtKB accession Q9UHI5.
    """
    return re.sub(r"-\d+$", "", str(accession).strip())


def find_header(row, contains):
    contains = [x.lower() for x in contains]
    for key in row:
        low = key.lower()
        if all(piece in low for piece in contains):
            return key
    return None


def load_uniprot_cache(path):
    path = Path(path)
    if not path.exists():
        return {}
    with open(path, "r", encoding="utf-8-sig", newline="") as f:
        rows = list(csv.DictReader(f))
    cache = {}
    for row in rows:
        key = row.get("Query_accession", "").strip()
        if key:
            cache[key] = row
    return cache


def save_uniprot_cache(path, cache):
    rows = list(cache.values())
    if not rows:
        return

    fieldnames = []
    seen = set()
    for row in rows:
        for key in row:
            if key not in seen:
                seen.add(key)
                fieldnames.append(key)

    write_csv(path, rows, fieldnames)


def fetch_uniprot_batch(accessions, timeout=60, attempts=4):
    query = "(" + " OR ".join(f"accession:{a}" for a in accessions) + ")"
    params = {
        "query": query,
        "format": "tsv",
        "fields": ",".join(UNIPROT_FIELDS),
        "size": "500",
    }
    url = "https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode(params)

    headers = {
        "User-Agent": "GPMV-proteomics-thesis-pipeline/1.0",
        "Accept": "text/tab-separated-values",
    }

    last_error = None
    for attempt in range(attempts):
        try:
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=timeout) as response:
                text = response.read().decode("utf-8")
            reader = csv.DictReader(io.StringIO(text), delimiter="\t")
            return list(reader)

        except urllib.error.HTTPError as e:
            last_error = e
            if e.code in {429, 500, 502, 503, 504} and attempt < attempts - 1:
                time.sleep(2 ** attempt)
                continue
            body = ""
            try:
                body = e.read().decode("utf-8", errors="replace")
            except Exception:
                pass
            raise RuntimeError(
                f"UniProt HTTP error {e.code}. Response:\n{body[:1000]}"
            ) from e

        except (urllib.error.URLError, TimeoutError) as e:
            last_error = e
            if attempt < attempts - 1:
                time.sleep(2 ** attempt)
                continue
            raise RuntimeError(f"UniProt request failed: {e}") from e

    raise RuntimeError(f"UniProt request failed: {last_error}")


def annotate_uniprot(accessions, cache_path, batch_size=50):
    cache = load_uniprot_cache(cache_path)

    canonical = []
    seen = set()
    for accession in accessions:
        c = canonical_accession(accession)
        if c and c not in seen:
            seen.add(c)
            canonical.append(c)

    missing = [a for a in canonical if a not in cache]
    print(f"UniProt: {len(canonical)} unique canonical accessions; {len(missing)} need retrieval.")

    for start in range(0, len(missing), batch_size):
        batch = missing[start:start + batch_size]
        print(f"  Fetching UniProt batch {start + 1}-{start + len(batch)} / {len(missing)}")
        result_rows = fetch_uniprot_batch(batch)

        found = set()
        for row in result_rows:
            entry_key = find_header(row, ["entry"])
            accession = row.get(entry_key, "").strip() if entry_key else ""
            if not accession:
                continue
            accession = canonical_accession(accession)
            found.add(accession)
            out = {"Query_accession": accession}
            out.update(row)
            cache[accession] = out

        for accession in batch:
            if accession not in found:
                cache[accession] = {
                    "Query_accession": accession,
                    "UniProt_retrieval_status": "Not found",
                }

        save_uniprot_cache(cache_path, cache)
        time.sleep(0.25)

    return cache


CATEGORY_TERMS = {
    "Plasma membrane / cell surface": [
        "plasma membrane",
        "cell membrane",
        "cell surface",
        "external side of plasma membrane",
        "apical membrane",
        "basal membrane",
        "basolateral membrane",
        "caveola",
    ],
    "Extracellular / secreted": [
        "secreted",
        "extracellular space",
        "extracellular region",
        "extracellular exosome",
        "extracellular vesicle",
    ],
    "Extracellular matrix": [
        "extracellular matrix",
        "basement membrane",
        "collagen-containing extracellular matrix",
        "matrix component",
    ],
    "Cytoskeleton / cell cortex": [
        "cytoskeleton",
        "cell cortex",
        "cortical cytoskeleton",
        "actin filament",
        "stress fiber",
        "focal adhesion",
        "microtubule",
        "intermediate filament",
    ],
    "Cytosol / cytoplasm": [
        "cytosol",
        "cytoplasm",
        "cytoplasmic",
    ],
    "Nucleus / chromatin": [
        "nucleus",
        "nuclear",
        "nucleoplasm",
        "nucleolus",
        "chromatin",
        "chromosome",
    ],
    "Mitochondrion": [
        "mitochondr",
    ],
    "ER / Golgi": [
        "endoplasmic reticulum",
        "golgi",
    ],
    "Endosome / lysosome / vesicle": [
        "endosome",
        "lysosome",
        "lysosomal",
        "phagosome",
        "autophagosome",
        "secretory vesicle",
        "transport vesicle",
        "cytoplasmic vesicle",
        "secretory granule",
    ],
    "Ribosome": [
        "ribosome",
        "ribosomal",
    ],
}


def text_from_uniprot_row(row):
    if not row:
        return "", "", "", "", ""

    subcell_key = find_header(row, ["subcellular", "location"])
    go_key = find_header(row, ["gene ontology", "cellular"])
    reviewed_key = find_header(row, ["reviewed"])
    tm_key = find_header(row, ["transmembrane"])
    signal_key = find_header(row, ["signal", "peptide"])

    subcell = row.get(subcell_key, "") if subcell_key else ""
    go_cc = row.get(go_key, "") if go_key else ""
    reviewed = row.get(reviewed_key, "") if reviewed_key else ""
    tm = row.get(tm_key, "") if tm_key else ""
    signal = row.get(signal_key, "") if signal_key else ""

    return str(subcell), str(go_cc), str(reviewed), str(tm), str(signal)


def classify_localization(subcell, go_cc):
    """
    Multi-label classification. A protein may belong to more than one category.
    This is deliberate because proteins can have multiple bona fide localizations.
    """
    text = f"{subcell}; {go_cc}".lower()
    labels = []

    for category, terms in CATEGORY_TERMS.items():
        if any(term in text for term in terms):
            labels.append(category)

    if not labels:
        labels = ["Other / no broad match"]

    return labels


def add_uniprot_annotations(rows, cache):
    annotated = []

    for row in rows:
        out = dict(row)
        original_acc = str(row.get(COL_ACCESSION, "")).strip()
        query_acc = canonical_accession(original_acc)
        uni = cache.get(query_acc, {})

        subcell, go_cc, reviewed, tm, signal = text_from_uniprot_row(uni)
        labels = classify_localization(subcell, go_cc)

        out["UniProt_query_accession"] = query_acc
        out["UniProt_reviewed"] = reviewed
        out["UniProt_subcellular_location"] = subcell
        out["UniProt_GO_cellular_component"] = go_cc
        out["UniProt_transmembrane"] = tm
        out["UniProt_signal_peptide"] = signal
        out["Localization_categories"] = "; ".join(labels)
        out["UniProt_retrieval_status"] = uni.get(
            "UniProt_retrieval_status",
            "Retrieved" if uni else "Missing from cache",
        )
        annotated.append(out)

    return annotated


def localization_summary(annotated_rows):
    groups = {
        "GPMV-enriched": [r for r in annotated_rows if r.get("Differential_group") == "GPMV-enriched"],
        "GPMV-depleted": [r for r in annotated_rows if r.get("Differential_group") == "GPMV-depleted"],
    }

    categories = list(CATEGORY_TERMS) + ["Other / no broad match"]
    summary = []

    for category in categories:
        row = {"Category": category}
        for group_name, group_rows in groups.items():
            count = 0
            for protein in group_rows:
                labels = {
                    x.strip()
                    for x in protein.get("Localization_categories", "").split(";")
                    if x.strip()
                }
                if category in labels:
                    count += 1
            denom = len(group_rows)
            row[f"{group_name}_count"] = count
            row[f"{group_name}_percent"] = (100.0 * count / denom) if denom else 0.0
        summary.append(row)

    return summary


def make_localization_plot(summary, path):
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not installed: skipping PNG plot.")
        return

    categories = [r["Category"] for r in summary]
    enriched = [r["GPMV-enriched_percent"] for r in summary]
    depleted = [r["GPMV-depleted_percent"] for r in summary]

    y = list(range(len(categories)))
    height = 0.38

    fig_height = max(6, 0.55 * len(categories))
    fig, ax = plt.subplots(figsize=(10, fig_height))
    ax.barh([v - height / 2 for v in y], enriched, height=height, label="GPMV-enriched")
    ax.barh([v + height / 2 for v in y], depleted, height=height, label="GPMV-depleted")

    ax.set_yticks(y)
    ax.set_yticklabels(categories)
    ax.invert_yaxis()
    ax.set_xlabel("Proteins in differential-abundance group (%)")
    ax.set_title("UniProt localization annotations of differentially abundant proteins")
    ax.legend()
    ax.text(
        0.0,
        -0.10,
        "Categories are non-mutually-exclusive; one protein may contribute to more than one category.",
        transform=ax.transAxes,
        fontsize=9,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_readme(path, counts, min_unique, adj_p, uniprot_done):
    text = f"""GPMV proteomics pipeline summary

FILTERING CRITERIA
- Protein FDR confidence: High
- Contaminant: False
- Minimum unique peptides: {min_unique}
- Quantitative comparison requires numeric GPMV/Cell ratio and adjusted p-value
- Differential significance: adjusted p < {adj_p}
- Enriched: significant and GPMV/Cell ratio > 1
- Depleted: significant and GPMV/Cell ratio < 1

COUNTS
- Total protein rows: {counts['total']}
- Passed QC: {counts['qc_passed']}
- Failed QC: {counts['qc_failed']}
- QC-passed and quantitatively comparable: {counts['quantifiable']}
- QC-passed but missing quantitative comparison: {counts['no_quant']}
- Significantly GPMV-enriched: {counts['enriched']}
- Significantly GPMV-depleted: {counts['depleted']}
- Not significant after QC: {counts['not_significant']}

UNIPROT
- Annotation performed: {'yes' if uniprot_done else 'no'}
- Localization categories are non-mutually-exclusive.
- Categories are generated from UniProt subcellular-location and
  GO Cellular Component text, not from protein names alone.
- UniProt reviewed/unreviewed status is retained as metadata, but it is
  not itself used as a biological localization category.
"""
    Path(path).write_text(text, encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(
        description="Filter GPMV proteomics data and annotate significant proteins with UniProt."
    )
    parser.add_argument("input_file", help="Proteome Discoverer .xlsx or .csv export")
    parser.add_argument(
        "--output-dir",
        default=None,
        help=(
            "Optional output folder. If omitted, results are saved in a "
            "'GPMV_proteomics_results' folder next to the input Excel/CSV file."
        ),
    )
    parser.add_argument(
        "--min-unique-peptides",
        type=int,
        default=2,
        help="Minimum unique peptides (default: 2)",
    )
    parser.add_argument(
        "--adj-p",
        type=float,
        default=0.05,
        help="Adjusted p-value threshold (default: 0.05)",
    )
    parser.add_argument(
        "--skip-uniprot",
        action="store_true",
        help="Run filtering only; do not access UniProt.",
    )
    parser.add_argument(
        "--uniprot-batch-size",
        type=int,
        default=50,
        help="Accessions per UniProt request (default: 50)",
    )
    args = parser.parse_args()

    input_path = Path(args.input_file).expanduser().resolve()

    if args.output_dir:
        output_dir = Path(args.output_dir).expanduser().resolve()
    else:
        output_dir = input_path.parent / "GPMV_proteomics_results"

    output_dir.mkdir(parents=True, exist_ok=True)

    headers, records = read_table(input_path)

    required = [COL_FDR, COL_CONTAMINANT, COL_UNIQUE, COL_ACCESSION]
    missing_required = [c for c in required if c not in headers]
    if missing_required:
        raise KeyError("Missing required columns: " + ", ".join(missing_required))

    ratio_col = find_column(headers, RATIO_CANDIDATES)
    pval_col = find_column(headers, PVAL_CANDIDATES)

    groups = filter_proteins(
        records,
        ratio_col=ratio_col,
        pval_col=pval_col,
        min_unique=args.min_unique_peptides,
        adj_p=args.adj_p,
    )

    counts = {
        "total": len(records),
        "qc_passed": len(groups["qc_passed"]),
        "qc_failed": len(groups["qc_failed"]),
        "quantifiable": len(groups["quantifiable"]),
        "no_quant": len(groups["no_quant"]),
        "enriched": len(groups["enriched"]),
        "depleted": len(groups["depleted"]),
        "not_significant": len(groups["not_significant"]),
    }

    print("\nFILTERING SUMMARY")
    print("-" * 60)
    for k, v in counts.items():
        print(f"{k:>18}: {v}")

    original_headers = headers[:]
    diff_headers = original_headers + ["Ratio_numeric", "Adj_P_numeric", "Differential_group"]

    write_csv(
        output_dir / "01_QC_failed.csv",
        groups["qc_failed"],
        original_headers + ["QC_failure_reason"],
    )
    write_csv(
        output_dir / "02_QC_passed_no_quantification.csv",
        groups["no_quant"],
        original_headers + ["Quantification_status"],
    )
    write_csv(
        output_dir / "03_GPMV_enriched_after_QC.csv",
        groups["enriched"],
        diff_headers,
    )
    write_csv(
        output_dir / "04_GPMV_depleted_after_QC.csv",
        groups["depleted"],
        diff_headers,
    )
    write_csv(
        output_dir / "05_not_significant_after_QC.csv",
        groups["not_significant"],
        diff_headers,
    )

    summary_rows = [
        {"Metric": "Total protein rows", "Count": counts["total"]},
        {"Metric": f"Passed QC (>= {args.min_unique_peptides} unique peptides)", "Count": counts["qc_passed"]},
        {"Metric": "Failed QC", "Count": counts["qc_failed"]},
        {"Metric": "QC-passed + quantitatively comparable", "Count": counts["quantifiable"]},
        {"Metric": "QC-passed but no quantitative comparison", "Count": counts["no_quant"]},
        {"Metric": "Significantly GPMV-enriched", "Count": counts["enriched"]},
        {"Metric": "Significantly GPMV-depleted", "Count": counts["depleted"]},
        {"Metric": "Not significant after QC", "Count": counts["not_significant"]},
    ]
    write_csv(output_dir / "00_filtering_summary.csv", summary_rows, ["Metric", "Count"])

    significant = groups["enriched"] + groups["depleted"]

    uniprot_done = False
    if not args.skip_uniprot:
        accessions = [r[COL_ACCESSION] for r in significant]
        cache_path = output_dir / "uniprot_annotation_cache.csv"

        cache = annotate_uniprot(
            accessions,
            cache_path=cache_path,
            batch_size=args.uniprot_batch_size,
        )

        annotated = add_uniprot_annotations(significant, cache)
        annotated_headers = diff_headers + [
            "UniProt_query_accession",
            "UniProt_reviewed",
            "UniProt_subcellular_location",
            "UniProt_GO_cellular_component",
            "UniProt_transmembrane",
            "UniProt_signal_peptide",
            "Localization_categories",
            "UniProt_retrieval_status",
        ]
        write_csv(
            output_dir / "06_significant_proteins_UniProt_annotated.csv",
            annotated,
            annotated_headers,
        )

        loc_summary = localization_summary(annotated)
        loc_fields = [
            "Category",
            "GPMV-enriched_count",
            "GPMV-enriched_percent",
            "GPMV-depleted_count",
            "GPMV-depleted_percent",
        ]
        write_csv(
            output_dir / "07_localization_summary.csv",
            loc_summary,
            loc_fields,
        )
        make_localization_plot(
            loc_summary,
            output_dir / "08_localization_enriched_vs_depleted.png",
        )

        # Reviewed/unreviewed metadata summary
        reviewed_counter = Counter()
        for r in annotated:
            reviewed_counter[
                (r.get("Differential_group", ""), r.get("UniProt_reviewed", "") or "Unknown")
            ] += 1

        reviewed_rows = [
            {"Differential_group": group, "UniProt_reviewed": status, "Count": count}
            for (group, status), count in sorted(reviewed_counter.items())
        ]
        write_csv(
            output_dir / "09_uniprot_reviewed_status_summary.csv",
            reviewed_rows,
            ["Differential_group", "UniProt_reviewed", "Count"],
        )

        uniprot_done = True

    write_readme(
        output_dir / "README_results.txt",
        counts,
        args.min_unique_peptides,
        args.adj_p,
        uniprot_done,
    )

    print(f"\nResults written to: {output_dir.resolve()}")
    if args.skip_uniprot:
        print("UniProt step skipped. Re-run without --skip-uniprot when internet is available.")


if __name__ == "__main__":
    main()