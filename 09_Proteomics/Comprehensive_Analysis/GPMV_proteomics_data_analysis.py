# ============================================================
# GPMV PROTEOMICS ANALYSIS
# MDA-MB-231 cells vs GPMVs
#
# The script:
#   1. Lets you select the Excel files manually
#   2. Calculates proteome overlap
#   3. Calculates significant enrichment/depletion in GPMVs
#   4. Checks unique-peptide numbers
#   5. Checks fold-change reporting limits
#   6. Compares calculated enriched/depleted proteins with
#      collaborator-provided UP/DOWN files
#   7. Analyses grouped enriched-protein categories
#   8. Generates bar charts, pie charts and a volcano plot
#   9. Saves results automatically to the Desktop
# ============================================================


from pathlib import Path
from datetime import datetime
import os
import traceback

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import tkinter as tk
from tkinter import filedialog, messagebox


# ============================================================
# USER SETTINGS
# ============================================================

# Statistical significance threshold
ALPHA = 0.05

# Direction threshold:
# ratio > 1 = enriched in GPMVs
# ratio < 1 = depleted in GPMVs
FOLD_CHANGE_THRESHOLD = 1.0

# Proteome Discoverer apparent ratio reporting limits
UPPER_RATIO_LIMIT = 100
LOWER_RATIO_LIMIT = 0.01

# Resolution of saved figures
FIGURE_DPI = 300


# ============================================================
# FILE / FOLDER FUNCTIONS
# ============================================================

def find_desktop():
    """
    Try to find the user's Windows Desktop automatically.
    Also checks common OneDrive Desktop locations.
    """

    home = Path.home()

    candidates = []

    # Common OneDrive environment variables on Windows
    for variable in [
        "OneDrive",
        "OneDriveCommercial",
        "OneDriveConsumer"
    ]:
        location = os.environ.get(variable)

        if location:
            candidates.append(Path(location) / "Desktop")

    # OneDrive folders such as:
    # OneDrive - University Name
    try:
        for folder in home.glob("OneDrive*"):
            candidates.append(folder / "Desktop")
    except Exception:
        pass

    # Standard Windows Desktop
    candidates.append(home / "Desktop")

    # Return the first valid Desktop
    for candidate in candidates:
        if candidate.exists() and candidate.is_dir():
            return candidate

    return None


def choose_excel_file(title, required=True):
    """
    Opens a Windows file-selection dialog.
    """

    filename = filedialog.askopenfilename(
        title=title,
        filetypes=[
            ("Excel files", "*.xlsx *.xls"),
            ("All files", "*.*")
        ]
    )

    if not filename:
        if required:
            raise RuntimeError(
                f"No file was selected for:\n{title}"
            )
        return None

    return Path(filename)


def read_excel_file(path, preferred_sheet=None):
    """
    Reads an Excel file.

    If the preferred sheet exists, it is used.
    Otherwise, the first sheet is used.
    """

    excel_file = pd.ExcelFile(path)

    if preferred_sheet in excel_file.sheet_names:
        sheet = preferred_sheet
    else:
        sheet = excel_file.sheet_names[0]

    dataframe = pd.read_excel(
        path,
        sheet_name=sheet
    )

    return dataframe


# ============================================================
# COLUMN IDENTIFICATION
# ============================================================

def find_column(dataframe, exact_name=None, contains=None):
    """
    Finds a column either by exact name or by keywords.
    """

    columns = list(dataframe.columns)

    # Exact match first
    if exact_name is not None:
        for column in columns:
            if str(column).strip() == exact_name:
                return column

    # Keyword match
    if contains is not None:

        if isinstance(contains, str):
            contains = [contains]

        for column in columns:

            column_text = str(column).lower()

            if all(
                keyword.lower() in column_text
                for keyword in contains
            ):
                return column

    raise KeyError(
        f"Could not find the requested column.\n"
        f"Exact name: {exact_name}\n"
        f"Contains: {contains}"
    )


# ============================================================
# DETECTION / PRESENCE FUNCTIONS
# ============================================================

def value_means_detected(value):
    """
    Proteome Discoverer uses values such as 'High'
    when a protein was found and 'n/a' when it was not.
    """

    if pd.isna(value):
        return False

    text = str(value).strip().lower()

    missing_values = {
        "",
        "n/a",
        "na",
        "nan",
        "none"
    }

    return text not in missing_values


def calculate_presence_mask(dataframe, columns):
    """
    Returns True if the protein was found in at least
    one replicate belonging to the selected sample group.
    """

    detected = dataframe[columns].apply(
        lambda column: column.map(value_means_detected)
    )

    return detected.any(axis=1)


# ============================================================
# PROTEIN SET COMPARISON
# ============================================================

def compare_accession_sets(
        calculated_dataframe,
        provided_dataframe,
        accession_column_calculated,
        accession_column_provided
):
    """
    Compare a protein set calculated from TOTAL.xlsx
    with the collaborator-provided UP or DOWN file.
    """

    calculated = set(
        calculated_dataframe[
            accession_column_calculated
        ]
        .dropna()
        .astype(str)
        .str.strip()
    )

    provided = set(
        provided_dataframe[
            accession_column_provided
        ]
        .dropna()
        .astype(str)
        .str.strip()
    )

    all_accessions = sorted(
        calculated.union(provided)
    )

    rows = []

    for accession in all_accessions:

        in_calculated = accession in calculated
        in_provided = accession in provided

        if in_calculated and in_provided:
            status = "Match"

        elif in_calculated:
            status = "Calculated only"

        else:
            status = "Provided file only"

        rows.append(
            {
                "Accession": accession,
                "Calculated": in_calculated,
                "Provided file": in_provided,
                "Status": status
            }
        )

    comparison = pd.DataFrame(rows)

    summary = {
        "Calculated proteins": len(calculated),
        "Provided proteins": len(provided),
        "Matching proteins": len(
            calculated.intersection(provided)
        ),
        "Calculated only": len(
            calculated.difference(provided)
        ),
        "Provided only": len(
            provided.difference(calculated)
        )
    }

    return comparison, summary


# ============================================================
# PLOTTING FUNCTIONS
# ============================================================

def save_bar_chart(
        labels,
        values,
        title,
        ylabel,
        output_file,
        horizontal=False
):
    """
    Save a simple publication-quality bar chart.
    """

    plt.figure(figsize=(8, 5))

    if horizontal:

        bars = plt.barh(labels, values)

        plt.xlabel(ylabel)
        plt.gca().invert_yaxis()

        for bar, value in zip(bars, values):

            plt.text(
                value,
                bar.get_y() + bar.get_height() / 2,
                f" {value}",
                va="center"
            )

    else:

        bars = plt.bar(labels, values)

        plt.ylabel(ylabel)

        for bar, value in zip(bars, values):

            plt.text(
                bar.get_x() + bar.get_width() / 2,
                value,
                f"{value}",
                ha="center",
                va="bottom"
            )

    plt.title(title)
    plt.tight_layout()

    plt.savefig(
        output_file,
        dpi=FIGURE_DPI,
        bbox_inches="tight"
    )

    plt.close()


def save_pie_chart(
        labels,
        values,
        title,
        output_file
):
    """
    Save a pie chart.
    """

    plt.figure(figsize=(8, 7))

    plt.pie(
        values,
        labels=labels,
        autopct="%1.1f%%",
        startangle=90
    )

    plt.title(title)

    plt.tight_layout()

    plt.savefig(
        output_file,
        dpi=FIGURE_DPI,
        bbox_inches="tight"
    )

    plt.close()


# ============================================================
# MAIN ANALYSIS
# ============================================================

def main():

    # --------------------------------------------------------
    # Open invisible Tkinter window
    # --------------------------------------------------------

    root = tk.Tk()
    root.withdraw()

    try:

        # ====================================================
        # SELECT INPUT FILES
        # ====================================================

        messagebox.showinfo(
            "GPMV Proteomics Analysis",
            "You will now be asked to select your Excel files.\n\n"
            "The TOTAL file is required.\n"
            "The UP, DOWN and GROUPED files are optional."
        )

        total_path = choose_excel_file(
            "SELECT TOTAL PROTEOME FILE\n"
            "Example: GMPVs_TOTAL.xlsx",
            required=True
        )

        up_path = choose_excel_file(
            "SELECT GPMV-ENRICHED / UP FILE\n"
            "Press Cancel if you do not want to use it.",
            required=False
        )

        down_path = choose_excel_file(
            "SELECT GPMV-DEPLETED / DOWN FILE\n"
            "Press Cancel if you do not want to use it.",
            required=False
        )

        grouped_path = choose_excel_file(
            "SELECT GROUPED ENRICHED-PROTEIN FILE\n"
            "Press Cancel if you do not want to use it.",
            required=False
        )


        # ====================================================
        # CREATE OUTPUT FOLDER ON DESKTOP
        # ====================================================

        desktop = find_desktop()

        if desktop is None:

            messagebox.showwarning(
                "Desktop not found",
                "Python could not automatically find your Desktop.\n"
                "Please select the folder where you want the "
                "results to be saved."
            )

            selected_folder = filedialog.askdirectory(
                title="Select output folder"
            )

            if not selected_folder:
                raise RuntimeError(
                    "No output folder was selected."
                )

            desktop = Path(selected_folder)

        timestamp = datetime.now().strftime(
            "%Y-%m-%d_%H%M%S"
        )

        output_dir = desktop / (
            f"Proteomics_Analysis_{timestamp}"
        )

        output_dir.mkdir(
            parents=True,
            exist_ok=True
        )


        # ====================================================
        # READ TOTAL PROTEOME FILE
        # ====================================================

        total = read_excel_file(
            total_path,
            preferred_sheet="Proteins"
        )

        # Remove completely empty rows
        total = total.dropna(
            how="all"
        ).copy()


        # ====================================================
        # IDENTIFY IMPORTANT COLUMNS
        # ====================================================

        accession_col = find_column(
            total,
            exact_name="Accession"
        )

        description_col = find_column(
            total,
            exact_name="Description"
        )

        ratio_col = find_column(
            total,
            exact_name=(
                "Abundance Ratio: (GMPVs) / (Cells)"
            ),
            contains=[
                "Abundance Ratio",
                "GMPVs",
                "Cells"
            ]
        )

        pvalue_col = find_column(
            total,
            exact_name=(
                "Abundance Ratio Adj. P-Value: "
                "(GMPVs) / (Cells)"
            ),
            contains=[
                "Adj. P-Value",
                "GMPVs",
                "Cells"
            ]
        )

        unique_peptides_col = find_column(
            total,
            exact_name="# Unique Peptides",
            contains="Unique Peptides"
        )


        # ====================================================
        # IDENTIFY FOUND-IN-SAMPLE COLUMNS
        # ====================================================

        cell_detection_columns = [
            column
            for column in total.columns
            if (
                "found in sample" in str(column).lower()
                and
                "cells" in str(column).lower()
            )
        ]

        gpmv_detection_columns = [
            column
            for column in total.columns
            if (
                "found in sample" in str(column).lower()
                and
                "gmpvs" in str(column).lower()
            )
        ]

        if len(cell_detection_columns) == 0:
            raise RuntimeError(
                "No 'Found in Sample' columns for Cells "
                "could be identified."
            )

        if len(gpmv_detection_columns) == 0:
            raise RuntimeError(
                "No 'Found in Sample' columns for GPMVs "
                "could be identified."
            )


        # ====================================================
        # REMOVE ROWS WITHOUT ACCESSION
        # ====================================================

        total = total[
            total[accession_col].notna()
        ].copy()


        # ====================================================
        # PROTEOME OVERLAP
        # ====================================================

        cells_found = calculate_presence_mask(
            total,
            cell_detection_columns
        )

        gpmvs_found = calculate_presence_mask(
            total,
            gpmv_detection_columns
        )

        shared_mask = (
            cells_found
            &
            gpmvs_found
        )

        cells_only_mask = (
            cells_found
            &
            ~gpmvs_found
        )

        gpmvs_only_mask = (
            gpmvs_found
            &
            ~cells_found
        )

        neither_mask = (
            ~cells_found
            &
            ~gpmvs_found
        )


        # Dataframes for each overlap category

        cells_only = total[
            cells_only_mask
        ].copy()

        shared = total[
            shared_mask
        ].copy()

        gpmvs_only = total[
            gpmvs_only_mask
        ].copy()

        neither = total[
            neither_mask
        ].copy()


        # Counts

        number_cells = int(
            cells_found.sum()
        )

        number_gpmvs = int(
            gpmvs_found.sum()
        )

        number_shared = int(
            shared_mask.sum()
        )

        number_cells_only = int(
            cells_only_mask.sum()
        )

        number_gpmvs_only = int(
            gpmvs_only_mask.sum()
        )

        number_neither = int(
            neither_mask.sum()
        )

        union_detected = (
            number_cells_only
            +
            number_shared
            +
            number_gpmvs_only
        )


        # Percentages

        if number_gpmvs > 0:

            percent_gpmv_shared = (
                number_shared
                /
                number_gpmvs
                *
                100
            )

        else:
            percent_gpmv_shared = np.nan


        if number_cells > 0:

            percent_cells_in_gpmvs = (
                number_shared
                /
                number_cells
                *
                100
            )

        else:
            percent_cells_in_gpmvs = np.nan


        # ====================================================
        # QUANTITATIVE ANALYSIS
        # ====================================================

        total["_Ratio_numeric"] = pd.to_numeric(
            total[ratio_col],
            errors="coerce"
        )

        total["_AdjP_numeric"] = pd.to_numeric(
            total[pvalue_col],
            errors="coerce"
        )

        total["_UniquePeptides_numeric"] = (
            pd.to_numeric(
                total[unique_peptides_col],
                errors="coerce"
            )
        )


        quantified_mask = (
            total["_Ratio_numeric"].notna()
            &
            total["_AdjP_numeric"].notna()
        )


        enriched_mask = (
            quantified_mask
            &
            (
                total["_AdjP_numeric"]
                <
                ALPHA
            )
            &
            (
                total["_Ratio_numeric"]
                >
                FOLD_CHANGE_THRESHOLD
            )
        )


        depleted_mask = (
            quantified_mask
            &
            (
                total["_AdjP_numeric"]
                <
                ALPHA
            )
            &
            (
                total["_Ratio_numeric"]
                <
                FOLD_CHANGE_THRESHOLD
            )
        )


        nonsignificant_mask = (
            quantified_mask
            &
            ~(
                enriched_mask
                |
                depleted_mask
            )
        )


        quantified = total[
            quantified_mask
        ].copy()

        enriched = total[
            enriched_mask
        ].copy()

        depleted = total[
            depleted_mask
        ].copy()

        nonsignificant = total[
            nonsignificant_mask
        ].copy()


        number_quantified = len(
            quantified
        )

        number_enriched = len(
            enriched
        )

        number_depleted = len(
            depleted
        )

        number_nonsignificant = len(
            nonsignificant
        )


        if number_quantified > 0:

            percent_enriched = (
                number_enriched
                /
                number_quantified
                *
                100
            )

            percent_depleted = (
                number_depleted
                /
                number_quantified
                *
                100
            )

            percent_nonsignificant = (
                number_nonsignificant
                /
                number_quantified
                *
                100
            )

        else:

            percent_enriched = np.nan
            percent_depleted = np.nan
            percent_nonsignificant = np.nan


        # ====================================================
        # UNIQUE PEPTIDE QC
        # ====================================================

        unique_peptides = total[
            "_UniquePeptides_numeric"
        ]

        zero_unique = (
            unique_peptides == 0
        )

        one_unique = (
            unique_peptides == 1
        )

        two_or_more_unique = (
            unique_peptides >= 2
        )

        below_two_unique = (
            unique_peptides.notna()
            &
            (
                unique_peptides < 2
            )
        )


        n_zero_unique = int(
            zero_unique.sum()
        )

        n_one_unique = int(
            one_unique.sum()
        )

        n_two_or_more = int(
            two_or_more_unique.sum()
        )

        n_missing_unique = int(
            unique_peptides.isna().sum()
        )

        n_below_two = int(
            below_two_unique.sum()
        )


        enriched_below_two = int(
            (
                enriched[
                    "_UniquePeptides_numeric"
                ]
                <
                2
            ).sum()
        )


        depleted_below_two = int(
            (
                depleted[
                    "_UniquePeptides_numeric"
                ]
                <
                2
            ).sum()
        )


        below_two_dataframe = total[
            below_two_unique
        ].copy()


        # ====================================================
        # CHECK RATIO REPORTING LIMITS
        # ====================================================

        ratio_values = total[
            "_Ratio_numeric"
        ]

        n_ratio_100 = int(
            np.isclose(
                ratio_values,
                UPPER_RATIO_LIMIT,
                equal_nan=False
            ).sum()
        )

        n_ratio_001 = int(
            np.isclose(
                ratio_values,
                LOWER_RATIO_LIMIT,
                equal_nan=False
            ).sum()
        )


        # ====================================================
        # CREATE SUMMARY TABLE
        # ====================================================

        summary_rows = [

            {
                "Metric":
                    "Rows/protein entries in TOTAL file",
                "Value":
                    len(total)
            },

            {
                "Metric":
                    "Proteins detected in Cells",
                "Value":
                    number_cells
            },

            {
                "Metric":
                    "Proteins detected in GPMVs",
                "Value":
                    number_gpmvs
            },

            {
                "Metric":
                    "Proteins detected in both",
                "Value":
                    number_shared
            },

            {
                "Metric":
                    "Proteins detected only in Cells",
                "Value":
                    number_cells_only
            },

            {
                "Metric":
                    "Proteins detected only in GPMVs",
                "Value":
                    number_gpmvs_only
            },

            {
                "Metric":
                    "Proteins assigned to neither group",
                "Value":
                    number_neither
            },

            {
                "Metric":
                    "Union of detected proteins",
                "Value":
                    union_detected
            },

            {
                "Metric":
                    "% of GPMV proteome shared with Cells",
                "Value":
                    round(percent_gpmv_shared, 2)
            },

            {
                "Metric":
                    "% of Cell proteome detected in GPMVs",
                "Value":
                    round(percent_cells_in_gpmvs, 2)
            },

            {
                "Metric":
                    "Proteins with ratio + adjusted p-value",
                "Value":
                    number_quantified
            },

            {
                "Metric":
                    "Significantly enriched in GPMVs",
                "Value":
                    number_enriched
            },

            {
                "Metric":
                    "Significantly depleted in GPMVs",
                "Value":
                    number_depleted
            },

            {
                "Metric":
                    "Not significantly different",
                "Value":
                    number_nonsignificant
            },

            {
                "Metric":
                    "% significantly enriched",
                "Value":
                    round(percent_enriched, 2)
            },

            {
                "Metric":
                    "% significantly depleted",
                "Value":
                    round(percent_depleted, 2)
            },

            {
                "Metric":
                    "% not significantly different",
                "Value":
                    round(
                        percent_nonsignificant,
                        2
                    )
            },

            {
                "Metric":
                    "Proteins with 0 unique peptides",
                "Value":
                    n_zero_unique
            },

            {
                "Metric":
                    "Proteins with 1 unique peptide",
                "Value":
                    n_one_unique
            },

            {
                "Metric":
                    "Proteins with >=2 unique peptides",
                "Value":
                    n_two_or_more
            },

            {
                "Metric":
                    "Proteins with <2 unique peptides",
                "Value":
                    n_below_two
            },

            {
                "Metric":
                    "Enriched proteins with <2 unique peptides",
                "Value":
                    enriched_below_two
            },

            {
                "Metric":
                    "Depleted proteins with <2 unique peptides",
                "Value":
                    depleted_below_two
            },

            {
                "Metric":
                    "Missing unique-peptide values",
                "Value":
                    n_missing_unique
            },

            {
                "Metric":
                    "Proteins reported at ratio = 100",
                "Value":
                    n_ratio_100
            },

            {
                "Metric":
                    "Proteins reported at ratio = 0.01",
                "Value":
                    n_ratio_001
            }
        ]

        summary = pd.DataFrame(
            summary_rows
        )


        # ====================================================
        # OVERLAP SUMMARY TABLE
        # ====================================================

        overlap_summary = pd.DataFrame(
            {
                "Category": [
                    "Cells only",
                    "Shared",
                    "GPMVs only"
                ],
                "Protein count": [
                    number_cells_only,
                    number_shared,
                    number_gpmvs_only
                ]
            }
        )


        # ====================================================
        # DIFFERENTIAL ABUNDANCE SUMMARY
        # ====================================================

        differential_summary = pd.DataFrame(
            {
                "Category": [
                    "Enriched in GPMVs",
                    "Depleted in GPMVs",
                    "Not significant"
                ],
                "Protein count": [
                    number_enriched,
                    number_depleted,
                    number_nonsignificant
                ],
                "Percentage": [
                    percent_enriched,
                    percent_depleted,
                    percent_nonsignificant
                ]
            }
        )


        # ====================================================
        # OPTIONAL: COMPARE WITH PROVIDED UP FILE
        # ====================================================

        up_comparison = None
        up_summary = None

        if up_path is not None:

            provided_up = read_excel_file(
                up_path,
                preferred_sheet="Proteins"
            )

            provided_up_accession = find_column(
                provided_up,
                exact_name="Accession"
            )

            up_comparison, up_summary = (
                compare_accession_sets(
                    enriched,
                    provided_up,
                    accession_col,
                    provided_up_accession
                )
            )


        # ====================================================
        # OPTIONAL: COMPARE WITH PROVIDED DOWN FILE
        # ====================================================

        down_comparison = None
        down_summary = None

        if down_path is not None:

            provided_down = read_excel_file(
                down_path,
                preferred_sheet="Proteins"
            )

            provided_down_accession = find_column(
                provided_down,
                exact_name="Accession"
            )

            down_comparison, down_summary = (
                compare_accession_sets(
                    depleted,
                    provided_down,
                    accession_col,
                    provided_down_accession
                )
            )


        # ====================================================
        # OPTIONAL: GROUPED ENRICHED-PROTEIN FILE
        # ====================================================

        category_summary = None
        grouped_detail = None

        if grouped_path is not None:

            grouped_excel = pd.ExcelFile(
                grouped_path
            )

            # Prefer detailed grouped-protein sheet
            if (
                "Grouped proteins"
                in
                grouped_excel.sheet_names
            ):

                grouped_detail = pd.read_excel(
                    grouped_path,
                    sheet_name="Grouped proteins"
                )

                category_col = find_column(
                    grouped_detail,
                    exact_name="Category"
                )

                category_summary = (
                    grouped_detail
                    .groupby(category_col)
                    .size()
                    .reset_index(
                        name="Protein count"
                    )
                    .sort_values(
                        "Protein count",
                        ascending=False
                    )
                    .reset_index(drop=True)
                )

            # Otherwise use Summary sheet
            elif (
                "Summary"
                in
                grouped_excel.sheet_names
            ):

                category_summary = pd.read_excel(
                    grouped_path,
                    sheet_name="Summary"
                )

            else:

                grouped_detail = pd.read_excel(
                    grouped_path,
                    sheet_name=
                    grouped_excel.sheet_names[0]
                )

                category_col = find_column(
                    grouped_detail,
                    exact_name="Category"
                )

                category_summary = (
                    grouped_detail
                    .groupby(category_col)
                    .size()
                    .reset_index(
                        name="Protein count"
                    )
                    .sort_values(
                        "Protein count",
                        ascending=False
                    )
                    .reset_index(drop=True)
                )


        # ====================================================
        # SAVE PNG FIGURES
        # ====================================================

        # 1. Proteome overlap bar chart

        save_bar_chart(
            labels=[
                "Cells only",
                "Shared",
                "GPMVs only"
            ],
            values=[
                number_cells_only,
                number_shared,
                number_gpmvs_only
            ],
            title="Proteome overlap: Cells vs GPMVs",
            ylabel="Number of proteins",
            output_file=(
                output_dir
                /
                "01_proteome_overlap_bar.png"
            )
        )


        # 2. Proteome overlap pie chart

        save_pie_chart(
            labels=[
                "Cells only",
                "Shared",
                "GPMVs only"
            ],
            values=[
                number_cells_only,
                number_shared,
                number_gpmvs_only
            ],
            title="Proteome overlap: Cells vs GPMVs",
            output_file=(
                output_dir
                /
                "02_proteome_overlap_pie.png"
            )
        )


        # 3. Differential abundance bar chart

        save_bar_chart(
            labels=[
                "GPMV enriched",
                "GPMV depleted",
                "Not significant"
            ],
            values=[
                number_enriched,
                number_depleted,
                number_nonsignificant
            ],
            title=(
                "Differential protein abundance: "
                "GPMVs vs Cells"
            ),
            ylabel="Number of proteins",
            output_file=(
                output_dir
                /
                "03_differential_abundance_bar.png"
            )
        )


        # 4. Differential abundance pie chart

        save_pie_chart(
            labels=[
                "GPMV enriched",
                "GPMV depleted",
                "Not significant"
            ],
            values=[
                number_enriched,
                number_depleted,
                number_nonsignificant
            ],
            title=(
                "Differential protein abundance: "
                "GPMVs vs Cells"
            ),
            output_file=(
                output_dir
                /
                "04_differential_abundance_pie.png"
            )
        )


        # 5. Unique peptide QC chart

        peptide_labels = [
            "0 unique peptides",
            "1 unique peptide",
            ">=2 unique peptides"
        ]

        peptide_values = [
            n_zero_unique,
            n_one_unique,
            n_two_or_more
        ]

        save_bar_chart(
            labels=peptide_labels,
            values=peptide_values,
            title="Protein identification: unique peptides",
            ylabel="Number of proteins",
            output_file=(
                output_dir
                /
                "05_unique_peptide_QC.png"
            )
        )


        # ====================================================
        # 6. VOLCANO PLOT
        # ====================================================

        volcano = quantified.copy()

        # Ratios must be positive for log2
        volcano = volcano[
            volcano["_Ratio_numeric"] > 0
        ].copy()

        # Protect against p = 0
        positive_pvalues = volcano.loc[
            volcano["_AdjP_numeric"] > 0,
            "_AdjP_numeric"
        ]

        if len(positive_pvalues) > 0:

            minimum_p = positive_pvalues.min()

        else:

            minimum_p = 1e-300

        volcano[
            "_P_for_plot"
        ] = volcano[
            "_AdjP_numeric"
        ].clip(
            lower=minimum_p
        )

        volcano[
            "_log2_ratio"
        ] = np.log2(
            volcano[
                "_Ratio_numeric"
            ]
        )

        volcano[
            "_minus_log10_p"
        ] = -np.log10(
            volcano[
                "_P_for_plot"
            ]
        )


        volcano_enriched = volcano[
            (
                volcano["_AdjP_numeric"]
                <
                ALPHA
            )
            &
            (
                volcano["_Ratio_numeric"]
                >
                FOLD_CHANGE_THRESHOLD
            )
        ]


        volcano_depleted = volcano[
            (
                volcano["_AdjP_numeric"]
                <
                ALPHA
            )
            &
            (
                volcano["_Ratio_numeric"]
                <
                FOLD_CHANGE_THRESHOLD
            )
        ]


        volcano_other = volcano.drop(
            index=(
                volcano_enriched.index
                .union(
                    volcano_depleted.index
                )
            )
        )


        plt.figure(
            figsize=(8, 6)
        )

        plt.scatter(
            volcano_other[
                "_log2_ratio"
            ],
            volcano_other[
                "_minus_log10_p"
            ],
            s=15,
            alpha=0.5,
            label="Not significant"
        )

        plt.scatter(
            volcano_enriched[
                "_log2_ratio"
            ],
            volcano_enriched[
                "_minus_log10_p"
            ],
            s=18,
            alpha=0.7,
            label="Enriched in GPMVs"
        )

        plt.scatter(
            volcano_depleted[
                "_log2_ratio"
            ],
            volcano_depleted[
                "_minus_log10_p"
            ],
            s=18,
            alpha=0.7,
            label="Depleted in GPMVs"
        )

        plt.axhline(
            -np.log10(ALPHA),
            linestyle="--",
            linewidth=1
        )

        plt.axvline(
            0,
            linestyle="--",
            linewidth=1
        )

        plt.xlabel(
            "log2 abundance ratio (GPMVs / Cells)"
        )

        plt.ylabel(
            "-log10 adjusted p-value"
        )

        plt.title(
            "GPMV vs Cell proteome"
        )

        plt.legend()

        plt.tight_layout()

        plt.savefig(
            output_dir
            /
            "06_volcano_plot.png",
            dpi=FIGURE_DPI,
            bbox_inches="tight"
        )

        plt.close()


        # ====================================================
        # CATEGORY PLOTS
        # ====================================================

        if category_summary is not None:

            category_name_col = (
                category_summary.columns[0]
            )

            count_col = find_column(
                category_summary,
                exact_name="Protein count",
                contains="count"
            )

            category_summary[
                count_col
            ] = pd.to_numeric(
                category_summary[
                    count_col
                ],
                errors="coerce"
            )


            category_summary = (
                category_summary
                .dropna(
                    subset=[count_col]
                )
                .sort_values(
                    count_col,
                    ascending=False
                )
            )


            # Horizontal bar chart

            plt.figure(
                figsize=(10, 7)
            )

            bars = plt.barh(
                category_summary[
                    category_name_col
                ],
                category_summary[
                    count_col
                ]
            )

            plt.xlabel(
                "Number of proteins"
            )

            plt.title(
                "Functional categories of "
                "GPMV-enriched proteins"
            )

            plt.gca().invert_yaxis()

            for bar, value in zip(
                bars,
                category_summary[
                    count_col
                ]
            ):

                plt.text(
                    value,
                    (
                        bar.get_y()
                        +
                        bar.get_height() / 2
                    ),
                    f" {int(value)}",
                    va="center"
                )

            plt.tight_layout()

            plt.savefig(
                output_dir
                /
                "07_enriched_protein_categories_bar.png",
                dpi=FIGURE_DPI,
                bbox_inches="tight"
            )

            plt.close()


            # Category pie chart

            plt.figure(
                figsize=(10, 8)
            )

            wedges, texts, autotexts = plt.pie(
                category_summary[
                    count_col
                ],
                autopct="%1.1f%%",
                startangle=90,
                pctdistance=0.8
            )

            plt.legend(
                wedges,
                category_summary[
                    category_name_col
                ],
                title="Categories",
                loc="center left",
                bbox_to_anchor=(1, 0.5)
            )

            plt.title(
                "Functional categories of "
                "GPMV-enriched proteins"
            )

            plt.tight_layout()

            plt.savefig(
                output_dir
                /
                "08_enriched_protein_categories_pie.png",
                dpi=FIGURE_DPI,
                bbox_inches="tight"
            )

            plt.close()


        # ====================================================
        # SAVE CSV FILES
        # ====================================================

        summary.to_csv(
            output_dir
            /
            "summary_statistics.csv",
            index=False
        )

        enriched.to_csv(
            output_dir
            /
            "significantly_enriched_in_GPMVs.csv",
            index=False
        )

        depleted.to_csv(
            output_dir
            /
            "significantly_depleted_in_GPMVs.csv",
            index=False
        )

        below_two_dataframe.to_csv(
            output_dir
            /
            "proteins_below_2_unique_peptides.csv",
            index=False
        )

        overlap_summary.to_csv(
            output_dir
            /
            "proteome_overlap_summary.csv",
            index=False
        )

        differential_summary.to_csv(
            output_dir
            /
            "differential_abundance_summary.csv",
            index=False
        )

        if category_summary is not None:

            category_summary.to_csv(
                output_dir
                /
                "enriched_category_summary.csv",
                index=False
            )


        # ====================================================
        # SAVE COMPLETE EXCEL RESULTS WORKBOOK
        # ====================================================

        excel_output = (
            output_dir
            /
            "proteomics_analysis_results.xlsx"
        )

        with pd.ExcelWriter(
            excel_output,
            engine="openpyxl"
        ) as writer:

            # Input files used

            input_files = pd.DataFrame(
                {
                    "Input": [
                        "TOTAL",
                        "UP",
                        "DOWN",
                        "GROUPED"
                    ],
                    "File": [
                        str(total_path),
                        (
                            str(up_path)
                            if up_path
                            else "Not selected"
                        ),
                        (
                            str(down_path)
                            if down_path
                            else "Not selected"
                        ),
                        (
                            str(grouped_path)
                            if grouped_path
                            else "Not selected"
                        )
                    ]
                }
            )

            input_files.to_excel(
                writer,
                sheet_name="Input_files",
                index=False
            )

            summary.to_excel(
                writer,
                sheet_name="Summary",
                index=False
            )

            overlap_summary.to_excel(
                writer,
                sheet_name="Overlap_summary",
                index=False
            )

            differential_summary.to_excel(
                writer,
                sheet_name="Differential_summary",
                index=False
            )

            cells_only.to_excel(
                writer,
                sheet_name="Cells_only",
                index=False
            )

            shared.to_excel(
                writer,
                sheet_name="Shared",
                index=False
            )

            gpmvs_only.to_excel(
                writer,
                sheet_name="GPMVs_only",
                index=False
            )

            neither.to_excel(
                writer,
                sheet_name="Neither_group",
                index=False
            )

            enriched.to_excel(
                writer,
                sheet_name="GPMV_enriched",
                index=False
            )

            depleted.to_excel(
                writer,
                sheet_name="GPMV_depleted",
                index=False
            )

            nonsignificant.to_excel(
                writer,
                sheet_name="Not_significant",
                index=False
            )

            below_two_dataframe.to_excel(
                writer,
                sheet_name="Below_2_peptides",
                index=False
            )

            if category_summary is not None:

                category_summary.to_excel(
                    writer,
                    sheet_name="Category_summary",
                    index=False
                )

            if grouped_detail is not None:

                grouped_detail.to_excel(
                    writer,
                    sheet_name="Grouped_proteins",
                    index=False
                )

            if up_comparison is not None:

                up_comparison.to_excel(
                    writer,
                    sheet_name="UP_file_comparison",
                    index=False
                )

            if down_comparison is not None:

                down_comparison.to_excel(
                    writer,
                    sheet_name="DOWN_file_comparison",
                    index=False
                )


        # ====================================================
        # SAVE TEXT REPORT
        # ====================================================

        report_path = (
            output_dir
            /
            "proteomics_summary.txt"
        )

        with open(
            report_path,
            "w",
            encoding="utf-8"
        ) as report:

            report.write(
                "GPMV PROTEOMICS ANALYSIS\n"
            )

            report.write(
                "=" * 60 + "\n\n"
            )

            report.write(
                "INPUT FILES\n"
            )

            report.write(
                f"TOTAL: {total_path}\n"
            )

            report.write(
                f"UP: {up_path}\n"
            )

            report.write(
                f"DOWN: {down_path}\n"
            )

            report.write(
                f"GROUPED: {grouped_path}\n\n"
            )


            report.write(
                "PROTEOME OVERLAP\n"
            )

            report.write(
                "-" * 60 + "\n"
            )

            report.write(
                f"Cells: {number_cells}\n"
            )

            report.write(
                f"GPMVs: {number_gpmvs}\n"
            )

            report.write(
                f"Shared: {number_shared}\n"
            )

            report.write(
                f"Cells only: "
                f"{number_cells_only}\n"
            )

            report.write(
                f"GPMVs only: "
                f"{number_gpmvs_only}\n"
            )

            report.write(
                f"Neither: {number_neither}\n"
            )

            report.write(
                f"Union detected: "
                f"{union_detected}\n\n"
            )


            report.write(
                f"% of GPMV proteome shared "
                f"with Cells: "
                f"{percent_gpmv_shared:.2f}%\n"
            )

            report.write(
                f"% of Cell proteome represented "
                f"in GPMVs: "
                f"{percent_cells_in_gpmvs:.2f}%\n\n"
            )


            report.write(
                "DIFFERENTIAL ABUNDANCE\n"
            )

            report.write(
                "-" * 60 + "\n"
            )

            report.write(
                f"Quantified proteins: "
                f"{number_quantified}\n"
            )

            report.write(
                f"Enriched in GPMVs: "
                f"{number_enriched} "
                f"({percent_enriched:.2f}%)\n"
            )

            report.write(
                f"Depleted in GPMVs: "
                f"{number_depleted} "
                f"({percent_depleted:.2f}%)\n"
            )

            report.write(
                f"Not significant: "
                f"{number_nonsignificant} "
                f"({percent_nonsignificant:.2f}%)\n\n"
            )


            report.write(
                "UNIQUE-PEPTIDE QC\n"
            )

            report.write(
                "-" * 60 + "\n"
            )

            report.write(
                f"0 unique peptides: "
                f"{n_zero_unique}\n"
            )

            report.write(
                f"1 unique peptide: "
                f"{n_one_unique}\n"
            )

            report.write(
                f">=2 unique peptides: "
                f"{n_two_or_more}\n"
            )

            report.write(
                f"<2 unique peptides total: "
                f"{n_below_two}\n"
            )

            report.write(
                f"Enriched proteins with <2: "
                f"{enriched_below_two}\n"
            )

            report.write(
                f"Depleted proteins with <2: "
                f"{depleted_below_two}\n\n"
            )


            report.write(
                "RATIO REPORTING LIMITS\n"
            )

            report.write(
                "-" * 60 + "\n"
            )

            report.write(
                f"Ratio = 100: "
                f"{n_ratio_100}\n"
            )

            report.write(
                f"Ratio = 0.01: "
                f"{n_ratio_001}\n\n"
            )


            if up_summary is not None:

                report.write(
                    "UP FILE COMPARISON\n"
                )

                report.write(
                    "-" * 60 + "\n"
                )

                for key, value in (
                    up_summary.items()
                ):
                    report.write(
                        f"{key}: {value}\n"
                    )

                report.write("\n")


            if down_summary is not None:

                report.write(
                    "DOWN FILE COMPARISON\n"
                )

                report.write(
                    "-" * 60 + "\n"
                )

                for key, value in (
                    down_summary.items()
                ):
                    report.write(
                        f"{key}: {value}\n"
                    )

                report.write("\n")


            if category_summary is not None:

                report.write(
                    "GPMV-ENRICHED PROTEIN "
                    "CATEGORIES\n"
                )

                report.write(
                    "-" * 60 + "\n"
                )

                total_grouped = (
                    category_summary[
                        category_summary.columns[1]
                    ].sum()
                )

                for _, row in (
                    category_summary.iterrows()
                ):

                    category = row.iloc[0]
                    count = row.iloc[1]

                    percentage = (
                        count
                        /
                        total_grouped
                        *
                        100
                    )

                    report.write(
                        f"{category}: "
                        f"{int(count)} "
                        f"({percentage:.2f}%)\n"
                    )


        # ====================================================
        # PRINT RESULTS IN PYCHARM
        # ====================================================

        print("\n")
        print("=" * 70)
        print("GPMV PROTEOMICS ANALYSIS COMPLETE")
        print("=" * 70)

        print("\nPROTEOME OVERLAP")
        print(
            f"Cells:       {number_cells}"
        )
        print(
            f"GPMVs:       {number_gpmvs}"
        )
        print(
            f"Shared:      {number_shared}"
        )
        print(
            f"Cells only:  {number_cells_only}"
        )
        print(
            f"GPMVs only:  {number_gpmvs_only}"
        )
        print(
            f"Neither:     {number_neither}"
        )

        print(
            f"\nGPMV proteins shared with Cells: "
            f"{percent_gpmv_shared:.2f}%"
        )

        print(
            f"Cell proteins represented in GPMVs: "
            f"{percent_cells_in_gpmvs:.2f}%"
        )


        print("\nDIFFERENTIAL ABUNDANCE")

        print(
            f"Quantified:       "
            f"{number_quantified}"
        )

        print(
            f"GPMV enriched:    "
            f"{number_enriched} "
            f"({percent_enriched:.2f}%)"
        )

        print(
            f"GPMV depleted:    "
            f"{number_depleted} "
            f"({percent_depleted:.2f}%)"
        )

        print(
            f"Not significant:  "
            f"{number_nonsignificant} "
            f"({percent_nonsignificant:.2f}%)"
        )


        print("\nUNIQUE PEPTIDES")

        print(
            f"0 unique peptides: "
            f"{n_zero_unique}"
        )

        print(
            f"1 unique peptide:  "
            f"{n_one_unique}"
        )

        print(
            f">=2 peptides:      "
            f"{n_two_or_more}"
        )

        print(
            f"<2 total:          "
            f"{n_below_two}"
        )


        print(
            "\nResults saved to:"
        )

        print(
            output_dir
        )

        print("=" * 70)


        # ====================================================
        # FINISHED MESSAGE
        # ====================================================

        messagebox.showinfo(
            "Analysis complete",
            "Proteomics analysis completed successfully.\n\n"
            f"Results were saved here:\n\n"
            f"{output_dir}"
        )


        # Automatically open the result folder on Windows

        if os.name == "nt":

            try:
                os.startfile(output_dir)
            except Exception:
                pass


    except Exception as error:

        traceback.print_exc()

        messagebox.showerror(
            "Error",
            f"An error occurred:\n\n{error}\n\n"
            "Check the PyCharm console for more information."
        )

    finally:

        root.destroy()


# ============================================================
# RUN SCRIPT
# ============================================================

if __name__ == "__main__":
    main()