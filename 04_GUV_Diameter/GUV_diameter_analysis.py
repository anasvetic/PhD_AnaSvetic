import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path
from scipy.stats import shapiro, mannwhitneyu, ttest_ind


# ============================================================
# 1. FILE PATHS
# ============================================================

# CHANGE THESE TWO PATHS TO YOUR FILES

file_7 = Path(
    r"C:\Users\Korisnik\Desktop\GUVs\GUV_Electroformation_diameter\20250925_GUVs_7_diameter\20250925_GUVs_7_diameter.csv"
)

file_33 = Path(
    r"C:\Users\Korisnik\Desktop\GUVs\GUV_Electroformation_diameter\20251002_GUVs_33_diameter\20251002_GUVs_33_diameter.csv"
)


# ============================================================
# 2. OUTPUT FOLDER
# ============================================================

output_folder = (
    file_7.parent /
    "GUV_7_vs_33_diameter_analysis"
)

output_folder.mkdir(
    parents=True,
    exist_ok=True
)

print("\nResults will be saved in:")
print(output_folder.resolve())


# ============================================================
# 3. FUNCTION FOR READING TEXT-LIKE FILES
# ============================================================

def read_text_first_column(file_path):

    """
    Fallback reader for CSV/TXT files.

    The first row is treated as a header and ignored as data.
    The function attempts to detect whether the file is:
    - tab separated
    - semicolon separated
    - comma separated
    - or simply one value per line
    """

    with open(
        file_path,
        "r",
        encoding="utf-8-sig",
        errors="replace"
    ) as f:

        lines = [
            line.strip()
            for line in f
            if line.strip()
        ]


    if len(lines) < 2:

        raise ValueError(
            f"The file contains insufficient data:\n{file_path}"
        )


    header = lines[0]


    # --------------------------------------------------------
    # Detect separator based primarily on HEADER
    # --------------------------------------------------------

    if "\t" in header:

        df = pd.read_csv(
            file_path,
            sep="\t",
            header=0,
            dtype=str,
            encoding="utf-8-sig"
        )

        values = df.iloc[:, 0]


    elif ";" in header:

        df = pd.read_csv(
            file_path,
            sep=";",
            header=0,
            dtype=str,
            encoding="utf-8-sig"
        )

        values = df.iloc[:, 0]


    elif "," in header:

        # Comma-separated CSV
        df = pd.read_csv(
            file_path,
            sep=",",
            header=0,
            dtype=str,
            encoding="utf-8-sig"
        )

        values = df.iloc[:, 0]


    else:

        # ----------------------------------------------------
        # One value per line.
        #
        # This is particularly useful if values look like:
        #
        # Diameter
        # 18,56
        # 21,42
        # 14,70
        #
        # We DON'T treat the comma as a column separator.
        # ----------------------------------------------------

        values = pd.Series(
            lines[1:]
        )


    return values


# ============================================================
# 4. FUNCTION TO READ DIAMETER FILE
# ============================================================

def read_diameters(file_path):

    print("\n" + "=" * 65)
    print(f"READING FILE: {file_path.name}")
    print("=" * 65)


    # --------------------------------------------------------
    # Check that file exists
    # --------------------------------------------------------

    if not file_path.exists():

        raise FileNotFoundError(
            f"\nFile does not exist:\n{file_path}"
        )


    extension = file_path.suffix.lower()

    print(f"Detected extension: {extension}")


    # ========================================================
    # TRY READING THE FILE
    # ========================================================

    try:

        # ----------------------------------------------------
        # Modern Excel
        # ----------------------------------------------------

        if extension in [
            ".xlsx",
            ".xlsm",
            ".xltx",
            ".xltm"
        ]:

            df = pd.read_excel(
                file_path,
                header=0,
                usecols=[0],
                dtype=str,
                engine="openpyxl"
            )

            diameter_column = df.columns[0]

            print("File read as Excel workbook.")
            print(f"Column name: {diameter_column}")

            values = df.iloc[:, 0]


        # ----------------------------------------------------
        # Old Excel format
        # ----------------------------------------------------

        elif extension == ".xls":

            df = pd.read_excel(
                file_path,
                header=0,
                usecols=[0],
                dtype=str,
                engine="xlrd"
            )

            diameter_column = df.columns[0]

            print("File read as old Excel workbook.")
            print(f"Column name: {diameter_column}")

            values = df.iloc[:, 0]


        # ----------------------------------------------------
        # TXT / CSV
        # ----------------------------------------------------

        elif extension in [
            ".txt",
            ".csv",
            ".tsv"
        ]:

            values = read_text_first_column(
                file_path
            )

            print("File read as text/CSV file.")


        # ----------------------------------------------------
        # Unknown extension
        # ----------------------------------------------------

        else:

            print(
                "Unknown extension. "
                "Trying Excel format first..."
            )

            try:

                df = pd.read_excel(
                    file_path,
                    header=0,
                    usecols=[0],
                    dtype=str,
                    engine="openpyxl"
                )

                values = df.iloc[:, 0]

                print(
                    "Successfully read as Excel."
                )

            except Exception:

                print(
                    "Excel reading failed. "
                    "Trying text format..."
                )

                values = read_text_first_column(
                    file_path
                )


    # ========================================================
    # FALLBACK
    #
    # This is important if a file LOOKS like an Excel file
    # but is actually text/CSV internally.
    # ========================================================

    except Exception as excel_error:

        print(
            "\nStandard Excel reading failed."
        )

        print(
            f"Reason: {excel_error}"
        )

        print(
            "\nAttempting to read the file "
            "as text/CSV instead..."
        )


        try:

            values = read_text_first_column(
                file_path
            )

            print(
                "Fallback successful."
            )

        except Exception as text_error:

            raise ValueError(
                "\nThe file could not be read as either "
                "Excel or text/CSV.\n\n"
                f"Excel error:\n{excel_error}\n\n"
                f"Text error:\n{text_error}"
            )


    # ========================================================
    # CLEAN DIAMETER VALUES
    # ========================================================

    diameters = (
        values
        .astype(str)
        .str.strip()

        # Remove unusual spaces
        .str.replace(
            "\u00A0",
            "",
            regex=False
        )

        .str.replace(
            " ",
            "",
            regex=False
        )

        # European decimal comma -> decimal point
        .str.replace(
            ",",
            ".",
            regex=False
        )
    )


    # Convert to numeric
    diameters = pd.to_numeric(
        diameters,
        errors="coerce"
    )


    # Remove invalid / blank rows
    diameters = diameters.dropna()


    # Reset index
    diameters = diameters.reset_index(
        drop=True
    )


    # ========================================================
    # REPORT IMPORT
    # ========================================================

    print(
        f"\nValid diameter measurements found: "
        f"{len(diameters)}"
    )

    print("\nFirst 10 diameter measurements:")

    print(
        diameters
        .head(10)
        .to_string(index=False)
    )


    if len(diameters) == 0:

        raise ValueError(
            f"\nNo numeric diameter measurements "
            f"were found in:\n{file_path}"
        )


    return diameters


# ============================================================
# 5. LOAD BOTH GUV DATASETS
# ============================================================

diameter_7 = read_diameters(
    file_7
)

diameter_33 = read_diameters(
    file_33
)


# ============================================================
# 6. CHECK SAMPLE SIZES
# ============================================================

if len(diameter_7) < 3:

    raise ValueError(
        "7 mol% dataset contains fewer "
        "than 3 valid measurements."
    )


if len(diameter_33) < 3:

    raise ValueError(
        "33 mol% dataset contains fewer "
        "than 3 valid measurements."
    )


print("\n" + "=" * 65)
print("DATA LOADED SUCCESSFULLY")
print("=" * 65)

print(
    f"7 mol% cholesterol GUVs:  "
    f"n = {len(diameter_7)}"
)

print(
    f"33 mol% cholesterol GUVs: "
    f"n = {len(diameter_33)}"
)


# ============================================================
# 7. DESCRIPTIVE STATISTICS FUNCTION
# ============================================================

def descriptive_stats(data):

    n = len(data)

    mean = data.mean()

    median = data.median()

    sd = data.std(
        ddof=1
    )

    sem = sd / np.sqrt(n)

    q1 = data.quantile(
        0.25
    )

    q3 = data.quantile(
        0.75
    )

    iqr = q3 - q1

    minimum = data.min()

    maximum = data.max()


    return {

        "N": n,

        "Mean": mean,

        "Median": median,

        "SD": sd,

        "SEM": sem,

        "Q1": q1,

        "Q3": q3,

        "IQR": iqr,

        "Minimum": minimum,

        "Maximum": maximum
    }


# ============================================================
# 8. CALCULATE DESCRIPTIVE STATISTICS
# ============================================================

stats_7 = descriptive_stats(
    diameter_7
)

stats_33 = descriptive_stats(
    diameter_33
)


# ============================================================
# 9. PRINT DESCRIPTIVE STATISTICS
# ============================================================

print("\n" + "=" * 65)
print("7 mol% CHOLESTEROL GUVs")
print("=" * 65)


for key, value in stats_7.items():

    if key == "N":

        print(
            f"{key}: {value}"
        )

    else:

        print(
            f"{key}: {value:.3f} µm"
        )


print("\n" + "=" * 65)
print("33 mol% CHOLESTEROL GUVs")
print("=" * 65)


for key, value in stats_33.items():

    if key == "N":

        print(
            f"{key}: {value}"
        )

    else:

        print(
            f"{key}: {value:.3f} µm"
        )


# ============================================================
# 10. SHAPIRO-WILK NORMALITY TEST
# ============================================================

# Null hypothesis:
# Data are compatible with a normal distribution.
#
# p >= 0.05 -> normality is not rejected
# p < 0.05  -> significant deviation from normality


shapiro_7 = shapiro(
    diameter_7
)

shapiro_33 = shapiro(
    diameter_33
)


print("\n" + "=" * 65)
print("SHAPIRO-WILK NORMALITY TEST")
print("=" * 65)


print(
    f"7 mol%:\n"
    f"W = {shapiro_7.statistic:.4f}\n"
    f"p = {shapiro_7.pvalue:.6g}"
)


print(
    f"\n33 mol%:\n"
    f"W = {shapiro_33.statistic:.4f}\n"
    f"p = {shapiro_33.pvalue:.6g}"
)


normal_7 = (
    shapiro_7.pvalue >= 0.05
)

normal_33 = (
    shapiro_33.pvalue >= 0.05
)


# ============================================================
# 11. SELECT STATISTICAL TEST
# ============================================================

# If BOTH groups are compatible with normality:
#
# -> Welch's independent t-test
#
# Otherwise:
#
# -> Mann-Whitney U test


if normal_7 and normal_33:

    test_name = (
        "Welch's independent t-test"
    )


    result = ttest_ind(
        diameter_7,
        diameter_33,
        equal_var=False
    )


    test_statistic = (
        result.statistic
    )

    p_value = (
        result.pvalue
    )

    statistic_name = "t"


else:

    test_name = (
        "Mann-Whitney U test"
    )


    result = mannwhitneyu(
        diameter_7,
        diameter_33,
        alternative="two-sided",
        method="auto"
    )


    test_statistic = (
        result.statistic
    )

    p_value = (
        result.pvalue
    )

    statistic_name = "U"


# ============================================================
# 12. SIGNIFICANCE LEVEL
# ============================================================

if p_value < 0.0001:

    significance = "****"

    p_text = "p < 0.0001"


elif p_value < 0.001:

    significance = "***"

    p_text = "p < 0.001"


elif p_value < 0.01:

    significance = "**"

    p_text = "p < 0.01"


elif p_value < 0.05:

    significance = "*"

    p_text = "p < 0.05"


else:

    significance = "ns"

    p_text = (
        f"p = {p_value:.4f}"
    )


# ============================================================
# 13. PRINT SIGNIFICANCE RESULT
# ============================================================

print("\n" + "=" * 65)
print("STATISTICAL COMPARISON")
print("=" * 65)


print(
    f"Test used: {test_name}"
)


print(
    f"{statistic_name} = "
    f"{test_statistic:.4f}"
)


print(
    f"Exact p-value = "
    f"{p_value:.8g}"
)


print(
    f"Significance = "
    f"{significance}"
)


if p_value < 0.05:

    print(
        "\nCONCLUSION:"
        "\nThe GUV diameter distributions "
        "differ significantly between "
        "7 mol% and 33 mol% cholesterol."
    )

else:

    print(
        "\nCONCLUSION:"
        "\nNo statistically significant "
        "difference in GUV diameter was "
        "detected between 7 mol% and "
        "33 mol% cholesterol."
    )


# ============================================================
# 14. SAVE CLEANED DIAMETER VALUES
# ============================================================

diameter_7.to_csv(

    output_folder /
    "GUV_7mol_cleaned_diameters.txt",

    index=False,

    header=[
        "Diameter_um"
    ]
)


diameter_33.to_csv(

    output_folder /
    "GUV_33mol_cleaned_diameters.txt",

    index=False,

    header=[
        "Diameter_um"
    ]
)


# ============================================================
# 15. SAVE STATISTICAL RESULTS
# ============================================================

stats_file = (

    output_folder /

    "GUV_7_vs_33_statistical_analysis.txt"

)


with open(
    stats_file,
    "w",
    encoding="utf-8"
) as f:


    f.write(
        "GUV DIAMETER STATISTICAL ANALYSIS\n"
    )

    f.write(
        "7 mol% vs 33 mol% cholesterol\n"
    )

    f.write(
        "=" * 65 +
        "\n\n"
    )


    # ========================================================
    # 7 MOL%
    # ========================================================

    f.write(
        "7 mol% CHOLESTEROL GUVs\n"
    )

    f.write(
        "-" * 40 +
        "\n"
    )


    for key, value in stats_7.items():

        if key == "N":

            f.write(
                f"{key}: "
                f"{value}\n"
            )

        else:

            f.write(
                f"{key}: "
                f"{value:.3f} µm\n"
            )


    # ========================================================
    # 33 MOL%
    # ========================================================

    f.write(
        "\n33 mol% CHOLESTEROL GUVs\n"
    )

    f.write(
        "-" * 40 +
        "\n"
    )


    for key, value in stats_33.items():

        if key == "N":

            f.write(
                f"{key}: "
                f"{value}\n"
            )

        else:

            f.write(
                f"{key}: "
                f"{value:.3f} µm\n"
            )


    # ========================================================
    # NORMALITY
    # ========================================================

    f.write(
        "\nNORMALITY TEST\n"
    )

    f.write(
        "-" * 40 +
        "\n"
    )


    f.write(

        f"7 mol%: "
        f"Shapiro-Wilk W = "
        f"{shapiro_7.statistic:.4f}, "
        f"p = "
        f"{shapiro_7.pvalue:.8g}\n"

    )


    f.write(

        f"33 mol%: "
        f"Shapiro-Wilk W = "
        f"{shapiro_33.statistic:.4f}, "
        f"p = "
        f"{shapiro_33.pvalue:.8g}\n"

    )


    # ========================================================
    # COMPARISON
    # ========================================================

    f.write(
        "\nSTATISTICAL COMPARISON\n"
    )

    f.write(
        "-" * 40 +
        "\n"
    )


    f.write(
        f"Test used: "
        f"{test_name}\n"
    )


    f.write(
        f"{statistic_name} = "
        f"{test_statistic:.4f}\n"
    )


    f.write(
        f"Exact p-value = "
        f"{p_value:.8g}\n"
    )


    f.write(
        f"Significance = "
        f"{significance}\n"
    )


    # ========================================================
    # CONCLUSION
    # ========================================================

    if p_value < 0.05:

        f.write(

            "\nConclusion: "
            "A statistically significant difference "
            "in GUV diameter was detected between "
            "the 7 mol% and 33 mol% cholesterol groups.\n"

        )

    else:

        f.write(

            "\nConclusion: "
            "No statistically significant difference "
            "in GUV diameter was detected between "
            "the 7 mol% and 33 mol% cholesterol groups.\n"

        )


# ============================================================
# 16. CREATE COMPARATIVE BOXPLOT
# ============================================================

fig, ax = plt.subplots(
    figsize=(6, 7)
)


data = [
    diameter_7,
    diameter_33
]


# ============================================================
# 17. BOXPLOT
# ============================================================

ax.boxplot(

    data,

    widths=0.45,

    showmeans=True
)


# ============================================================
# 18. INDIVIDUAL DATA POINTS
# ============================================================

rng = np.random.default_rng(
    42
)


for position, values in enumerate(
    data,
    start=1
):


    jitter = rng.normal(

        loc=position,

        scale=0.035,

        size=len(values)

    )


    ax.scatter(

        jitter,

        values,

        alpha=0.5,

        s=22

    )


# ============================================================
# 19. AXIS LABELS
# ============================================================

ax.set_xticks(
    [1, 2]
)


ax.set_xticklabels(
    [
        "7 mol% Chol",
        "33 mol% Chol"
    ],
    fontsize=11
)


ax.set_ylabel(

    "GUV equivalent diameter (µm)",

    fontsize=12
)


ax.set_title(

    "GUV Diameter Distribution",

    fontsize=13
)


# Start y-axis at zero

ax.set_ylim(
    bottom=0
)


# ============================================================
# 20. SIGNIFICANCE BRACKET
# ============================================================

combined_max = max(

    diameter_7.max(),

    diameter_33.max()

)


combined_min = min(

    diameter_7.min(),

    diameter_33.min()

)


data_range = (

    combined_max -
    combined_min

)


# Prevent problems if all numbers happen to be identical

if data_range == 0:

    data_range = combined_max * 0.1


if data_range == 0:

    data_range = 1


bar_y = (

    combined_max +
    data_range * 0.10

)


bar_height = (

    data_range * 0.04

)


# Draw significance bracket

ax.plot(

    [1, 1, 2, 2],

    [
        bar_y,
        bar_y + bar_height,
        bar_y + bar_height,
        bar_y
    ],

    linewidth=1.2
)


# Significance stars / ns

ax.text(

    1.5,

    bar_y + bar_height * 1.5,

    significance,

    ha="center",

    va="bottom",

    fontsize=14
)


# p-value

ax.text(

    1.5,

    bar_y + bar_height * 0.15,

    p_text,

    ha="center",

    va="bottom",

    fontsize=9
)


# Add room at top

ax.set_ylim(

    0,

    bar_y + data_range * 0.20

)


# ============================================================
# 21. CLEAN FIGURE APPEARANCE
# ============================================================

ax.spines[
    "top"
].set_visible(
    False
)


ax.spines[
    "right"
].set_visible(
    False
)


plt.tight_layout()


# ============================================================
# 22. SAVE BOXPLOT
# ============================================================

graph_file = (

    output_folder /

    "GUV_7_vs_33_boxplot_significance.png"

)


plt.savefig(

    graph_file,

    dpi=600,

    bbox_inches="tight"

)


plt.show()


# ============================================================
# 23. FINAL SUMMARY
# ============================================================

print("\n" + "=" * 65)

print(
    "ANALYSIS COMPLETE"
)

print(
    "=" * 65
)


print(
    "\nResults saved in:"
)

print(
    output_folder.resolve()
)


print(
    "\nFiles created:"
)

print(
    "1. GUV_7mol_cleaned_diameters.txt"
)

print(
    "2. GUV_33mol_cleaned_diameters.txt"
)

print(
    "3. GUV_7_vs_33_statistical_analysis.txt"
)

print(
    "4. GUV_7_vs_33_boxplot_significance.png"
)