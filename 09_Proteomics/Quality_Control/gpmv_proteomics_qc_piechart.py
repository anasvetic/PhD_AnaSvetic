#!/usr/bin/env python3
"""
GPMV proteomics QC filtering + differential-abundance pie chart

This script starts from the COMPLETE Proteome Discoverer protein export and:
1. Retains proteins with:
   - Protein FDR Confidence: Combined == "High"
   - Contaminant == False
   - # Unique Peptides >= 2
2. Retains proteins with numeric GPMV/Cell abundance ratio and adjusted p-value.
3. Classifies proteins as:
   - GPMV-enriched: adjusted p < 0.05 and ratio > 1
   - GPMV-depleted: adjusted p < 0.05 and ratio < 1
   - Not significant: all remaining quantitatively comparable proteins
4. Saves the filtered protein table, a summary CSV, and a thesis-ready pie chart.

Expected result for the current dataset:
- 2,661 proteins retained for quantitative comparison
- 152 GPMV-enriched
- 470 GPMV-depleted
- 2,039 not significant

Usage in PyCharm terminal:
    python gpmv_proteomics_qc_piechart.py "C:\\path\\to\\protein_table.xlsx"

Results are automatically saved in a folder called:
    GPMV_QC_piechart_results
next to the original Excel file.

Dependencies:
    pandas
    matplotlib
    openpyxl

Install if necessary:
    pip install pandas matplotlib openpyxl
"""

from pathlib import Path
import argparse
import pandas as pd
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------
# COLUMN NAMES FROM THE PROTEOME DISCOVERER EXPORT
# ---------------------------------------------------------------------

COL_FDR = "Protein FDR Confidence: Combined"
COL_CONTAMINANT = "Contaminant"
COL_UNIQUE = "# Unique Peptides"
COL_RATIO = "Abundance Ratio: (GMPVs) / (Cells)"
COL_PADJ = "Abundance Ratio Adj. P-Value: (GMPVs) / (Cells)"


def contaminant_is_false(series):
    """
    Treat False, 0, blank and common textual false values as non-contaminants.
    """
    text = series.fillna("").astype(str).str.strip().str.lower()
    return text.isin(["false", "0", "no", "n", ""])


def autopct_with_count(values):
    """
    Show percentage and absolute protein count inside each pie slice.
    """
    total = sum(values)

    def _format(pct):
        count = int(round(pct * total / 100.0))
        return f"{pct:.1f}%\n(n={count})"

    return _format


def main():
    parser = argparse.ArgumentParser(
        description="Apply GPMV proteomics QC criteria and generate a differential-abundance pie chart."
    )
    parser.add_argument(
        "input_file",
        help="Path to the complete Proteome Discoverer .xlsx file",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Optional output folder. If omitted, results are saved next to the input Excel file.",
    )
    parser.add_argument(
        "--min-unique-peptides",
        type=int,
        default=2,
        help="Minimum number of unique peptides required per protein (default: 2)",
    )
    parser.add_argument(
        "--adj-p",
        type=float,
        default=0.05,
        help="Adjusted p-value significance threshold (default: 0.05)",
    )
    args = parser.parse_args()

    input_path = Path(args.input_file).expanduser().resolve()

    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found:\n{input_path}")

    if args.output_dir:
        output_dir = Path(args.output_dir).expanduser().resolve()
    else:
        output_dir = input_path.parent / "GPMV_QC_piechart_results"

    output_dir.mkdir(parents=True, exist_ok=True)

    # -----------------------------------------------------------------
    # READ DATA
    # -----------------------------------------------------------------

    df = pd.read_excel(input_path)

    required_columns = [
        COL_FDR,
        COL_CONTAMINANT,
        COL_UNIQUE,
        COL_RATIO,
        COL_PADJ,
    ]

    missing = [col for col in required_columns if col not in df.columns]
    if missing:
        raise KeyError(
            "The following required columns were not found:\n- "
            + "\n- ".join(missing)
        )

    total_rows = len(df)

    # -----------------------------------------------------------------
    # PROTEIN-LEVEL QC
    # -----------------------------------------------------------------

    fdr_mask = (
        df[COL_FDR]
        .fillna("")
        .astype(str)
        .str.strip()
        .str.lower()
        .eq("high")
    )

    contaminant_mask = contaminant_is_false(df[COL_CONTAMINANT])

    unique_peptides = pd.to_numeric(df[COL_UNIQUE], errors="coerce")
    peptide_mask = unique_peptides >= args.min_unique_peptides

    qc_mask = fdr_mask & contaminant_mask & peptide_mask
    qc_df = df.loc[qc_mask].copy()

    # -----------------------------------------------------------------
    # QUANTITATIVE FILTERING
    # -----------------------------------------------------------------

    qc_df["Ratio_numeric"] = pd.to_numeric(qc_df[COL_RATIO], errors="coerce")
    qc_df["Adj_P_numeric"] = pd.to_numeric(qc_df[COL_PADJ], errors="coerce")

    quant_df = qc_df.dropna(
        subset=["Ratio_numeric", "Adj_P_numeric"]
    ).copy()

    # -----------------------------------------------------------------
    # DIFFERENTIAL-ABUNDANCE CLASSIFICATION
    # -----------------------------------------------------------------

    quant_df["Differential_group"] = "Not significant"

    enriched_mask = (
        (quant_df["Adj_P_numeric"] < args.adj_p)
        & (quant_df["Ratio_numeric"] > 1)
    )

    depleted_mask = (
        (quant_df["Adj_P_numeric"] < args.adj_p)
        & (quant_df["Ratio_numeric"] < 1)
    )

    quant_df.loc[enriched_mask, "Differential_group"] = "GPMV-enriched"
    quant_df.loc[depleted_mask, "Differential_group"] = "GPMV-depleted"

    enriched_n = int(enriched_mask.sum())
    depleted_n = int(depleted_mask.sum())
    nonsig_n = int(
        (quant_df["Differential_group"] == "Not significant").sum()
    )

    retained_n = len(quant_df)
    qc_n = len(qc_df)
    missing_quant_n = qc_n - retained_n

    # -----------------------------------------------------------------
    # SAVE TABLES
    # -----------------------------------------------------------------

    quant_df.to_csv(
        output_dir / "01_QC_filtered_quantitative_proteins.csv",
        index=False,
    )

    quant_df.loc[
        quant_df["Differential_group"] == "GPMV-enriched"
    ].to_csv(
        output_dir / "02_GPMV_enriched_after_QC.csv",
        index=False,
    )

    quant_df.loc[
        quant_df["Differential_group"] == "GPMV-depleted"
    ].to_csv(
        output_dir / "03_GPMV_depleted_after_QC.csv",
        index=False,
    )

    quant_df.loc[
        quant_df["Differential_group"] == "Not significant"
    ].to_csv(
        output_dir / "04_not_significant_after_QC.csv",
        index=False,
    )

    summary = pd.DataFrame(
        {
            "Metric": [
                "Total proteins in original export",
                "Proteins passing protein-level QC",
                "QC-passed proteins retained for quantitative comparison",
                "QC-passed proteins lacking ratio and/or adjusted p-value",
                "GPMV-enriched",
                "GPMV-depleted",
                "Not significant",
            ],
            "Count": [
                total_rows,
                qc_n,
                retained_n,
                missing_quant_n,
                enriched_n,
                depleted_n,
                nonsig_n,
            ],
        }
    )

    summary.to_csv(
        output_dir / "00_filtering_summary.csv",
        index=False,
    )

    # -----------------------------------------------------------------
    # PIE CHART
    # -----------------------------------------------------------------

    labels = [
        "GPMV-enriched",
        "GPMV-depleted",
        "Not significant",
    ]

    values = [
        enriched_n,
        depleted_n,
        nonsig_n,
    ]

    fig, ax = plt.subplots(figsize=(8, 7))

    ax.pie(
        values,
        labels=labels,
        autopct=autopct_with_count(values),
        startangle=90,
        counterclock=False,
        wedgeprops={"edgecolor": "white", "linewidth": 1},
        textprops={"fontsize": 11},
    )

    ax.set_title(
        "Differential protein abundance after protein-level QC",
        fontsize=14,
        pad=18,
    )

    ax.axis("equal")

    fig.tight_layout()

    png_path = output_dir / "05_differential_abundance_pie_chart.png"
    pdf_path = output_dir / "05_differential_abundance_pie_chart.pdf"

    fig.savefig(
        png_path,
        dpi=300,
        bbox_inches="tight",
    )
    fig.savefig(
        pdf_path,
        bbox_inches="tight",
    )

    plt.close(fig)

    # -----------------------------------------------------------------
    # PRINT SUMMARY
    # -----------------------------------------------------------------

    print("\nGPMV PROTEOMICS QC + DIFFERENTIAL ABUNDANCE")
    print("-" * 60)
    print(f"Original protein rows:                  {total_rows}")
    print(f"Passed protein-level QC:                {qc_n}")
    print(f"Retained for quantitative comparison:   {retained_n}")
    print(f"Missing quantitative comparison:        {missing_quant_n}")
    print()
    print(f"GPMV-enriched:                          {enriched_n}")
    print(f"GPMV-depleted:                          {depleted_n}")
    print(f"Not significant:                        {nonsig_n}")
    print()
    print(f"Results saved to:\n{output_dir}")
    print(f"\nPie chart:\n{png_path}")


if __name__ == "__main__":
    main()
