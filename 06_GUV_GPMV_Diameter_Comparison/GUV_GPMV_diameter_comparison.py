import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path
from scipy.stats import shapiro, mannwhitneyu, ttest_ind


# ============================================================
# 1. FILE PATHS
# ============================================================

# GPMV Nikon TXT file
gpmv_file = Path(
    r"C:\Users\Korisnik\Desktop\GPMVs\GPMVs_epifluo_diameter\MDA-MB-231\20241002_231_GPMV_diameter\gpmv 02102024 diameter and circularity.txt"
)

# GUV CSV file containing already extracted diameters
guv_file = Path(
    r"C:\Users\Korisnik\Desktop\GUVs\GUV_Electroformation_diameter\20251002_GUVs_33_diameter\20251002_GUVs_33_diameter.csv"
)


# ============================================================
# 2. GRAPH LABELS
# ============================================================

GPMV_LABEL = "GPMVs"

# Change this if appropriate, for example:
# "33 mol% Chol GUVs"
GUV_LABEL = "33 mol% Chol GUVs"


# ============================================================
# 3. GPMV DIAMETER INCLUSION RANGE
# ============================================================

# Same filtering logic as in the previous GPMV analysis

GPMV_MIN_DIAMETER = 1.0
GPMV_MAX_DIAMETER = 100.0


# ============================================================
# 4. OUTPUT FOLDER
# ============================================================

output_folder = (
    gpmv_file.parent /
    "GPMV_vs_GUV_diameter_analysis"
)

output_folder.mkdir(
    parents=True,
    exist_ok=True
)

print("\nResults will be saved in:")
print(output_folder.resolve())


# ============================================================
# 5. READ GPMV TXT FILE
# ============================================================

def read_gpmv_txt(file_path):

    print("\n" + "=" * 65)
    print("READING GPMV FILE")
    print("=" * 65)

    if not file_path.exists():

        raise FileNotFoundError(
            f"GPMV file not found:\n{file_path}"
        )


    # Nikon TXT files are tab-separated
    df = pd.read_csv(
        file_path,
        sep="\t",
        dtype=str,
        engine="python"
    )


    # Remove completely empty columns
    df = df.dropna(
        axis=1,
        how="all"
    )


    # Clean column names
    df.columns = (
        df.columns
        .str.replace(
            "\ufeff",
            "",
            regex=False
        )
        .str.strip()
    )


    print("\nColumns found:")
    print(df.columns.tolist())


    # --------------------------------------------------------
    # Make sure EqDiameter exists
    # --------------------------------------------------------

    if "EqDiameter" not in df.columns:

        raise ValueError(
            "\n'EqDiameter' column not found in GPMV TXT file.\n"
            f"Columns found: {df.columns.tolist()}"
        )


    # --------------------------------------------------------
    # Clean EqDiameter
    # --------------------------------------------------------

    df["EqDiameter_numeric"] = (
        df["EqDiameter"]
        .astype(str)
        .str.strip()
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
        .str.replace(
            ",",
            ".",
            regex=False
        )
    )


    df["EqDiameter_numeric"] = pd.to_numeric(
        df["EqDiameter_numeric"],
        errors="coerce"
    )


    # Remove non-numeric values
    numeric_df = df[
        df["EqDiameter_numeric"].notna()
    ].copy()


    # --------------------------------------------------------
    # Apply GPMV cutoff
    # --------------------------------------------------------

    excluded_df = numeric_df[
        (
            numeric_df["EqDiameter_numeric"]
            < GPMV_MIN_DIAMETER
        )
        |
        (
            numeric_df["EqDiameter_numeric"]
            > GPMV_MAX_DIAMETER
        )
    ].copy()


    filtered_df = numeric_df[
        (
            numeric_df["EqDiameter_numeric"]
            >= GPMV_MIN_DIAMETER
        )
        &
        (
            numeric_df["EqDiameter_numeric"]
            <= GPMV_MAX_DIAMETER
        )
    ].copy()


    diameters = (
        filtered_df["EqDiameter_numeric"]
        .reset_index(drop=True)
    )


    print(
        f"\nValid numeric measurements: "
        f"{len(numeric_df)}"
    )

    print(
        f"Measurements retained: "
        f"{len(diameters)}"
    )

    print(
        f"Measurements excluded: "
        f"{len(excluded_df)}"
    )


    print(
        f"\nGPMV inclusion range: "
        f"{GPMV_MIN_DIAMETER}–"
        f"{GPMV_MAX_DIAMETER} µm"
    )


    # Save excluded GPMV objects
    excluded_df.to_csv(
        output_folder /
        "excluded_GPMV_objects.txt",
        sep="\t",
        index=False
    )


    return diameters


# ============================================================
# 6. READ GUV CSV FILE
# ============================================================

def read_guv_csv(file_path):

    print("\n" + "=" * 65)
    print("READING GUV FILE")
    print("=" * 65)


    if not file_path.exists():

        raise FileNotFoundError(
            f"GUV file not found:\n{file_path}"
        )


    # --------------------------------------------------------
    # Read file as raw lines first.
    #
    # This is robust for a one-column CSV such as:
    #
    # Diameter
    # 15,34
    # 20,56
    # 18,91
    #
    # where comma is the DECIMAL separator rather than
    # the CSV column separator.
    # --------------------------------------------------------

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
            "GUV CSV file contains insufficient data."
        )


    # First row is the column name
    header = lines[0]

    print(
        f"\nColumn/header name: {header}"
    )


    # --------------------------------------------------------
    # If semicolon-separated
    # --------------------------------------------------------

    if ";" in header:

        df = pd.read_csv(
            file_path,
            sep=";",
            header=0,
            dtype=str,
            encoding="utf-8-sig"
        )

        values = df.iloc[:, 0]


    # --------------------------------------------------------
    # If tab-separated
    # --------------------------------------------------------

    elif "\t" in header:

        df = pd.read_csv(
            file_path,
            sep="\t",
            header=0,
            dtype=str,
            encoding="utf-8-sig"
        )

        values = df.iloc[:, 0]


    # --------------------------------------------------------
    # Otherwise assume one diameter per line.
    #
    # The first line is skipped because it is the header.
    # --------------------------------------------------------

    else:

        values = pd.Series(
            lines[1:]
        )


    # --------------------------------------------------------
    # Clean values
    # --------------------------------------------------------

    diameters = (
        values
        .astype(str)
        .str.strip()
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
        .str.replace(
            ",",
            ".",
            regex=False
        )
    )


    diameters = pd.to_numeric(
        diameters,
        errors="coerce"
    )


    diameters = (
        diameters
        .dropna()
        .reset_index(drop=True)
    )


    print(
        f"\nValid GUV diameter measurements: "
        f"{len(diameters)}"
    )


    print("\nFirst 10 GUV values:")

    print(
        diameters
        .head(10)
        .to_string(index=False)
    )


    return diameters


# ============================================================
# 7. LOAD BOTH DATASETS
# ============================================================

gpmv_diameters = read_gpmv_txt(
    gpmv_file
)

guv_diameters = read_guv_csv(
    guv_file
)


# ============================================================
# 8. CHECK SAMPLE SIZE
# ============================================================

if len(gpmv_diameters) < 3:

    raise ValueError(
        "GPMV dataset contains fewer than "
        "3 valid measurements."
    )


if len(guv_diameters) < 3:

    raise ValueError(
        "GUV dataset contains fewer than "
        "3 valid measurements."
    )


print("\n" + "=" * 65)
print("DATA LOADED")
print("=" * 65)

print(
    f"GPMVs: n = {len(gpmv_diameters)}"
)

print(
    f"GUVs:  n = {len(guv_diameters)}"
)


# ============================================================
# 9. DESCRIPTIVE STATISTICS
# ============================================================

def descriptive_stats(data):

    n = len(data)

    mean = data.mean()
    median = data.median()

    sd = data.std(
        ddof=1
    )

    sem = sd / np.sqrt(n)

    q1 = data.quantile(0.25)
    q3 = data.quantile(0.75)

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


gpmv_stats = descriptive_stats(
    gpmv_diameters
)

guv_stats = descriptive_stats(
    guv_diameters
)


# ============================================================
# 10. PRINT DESCRIPTIVE STATISTICS
# ============================================================

print("\n" + "=" * 65)
print("GPMV DIAMETERS")
print("=" * 65)

for key, value in gpmv_stats.items():

    if key == "N":

        print(
            f"{key}: {value}"
        )

    else:

        print(
            f"{key}: {value:.3f} µm"
        )


print("\n" + "=" * 65)
print("GUV DIAMETERS")
print("=" * 65)

for key, value in guv_stats.items():

    if key == "N":

        print(
            f"{key}: {value}"
        )

    else:

        print(
            f"{key}: {value:.3f} µm"
        )


# ============================================================
# 11. SHAPIRO-WILK NORMALITY TEST
# ============================================================

shapiro_gpmv = shapiro(
    gpmv_diameters
)

shapiro_guv = shapiro(
    guv_diameters
)


print("\n" + "=" * 65)
print("SHAPIRO-WILK NORMALITY TEST")
print("=" * 65)


print(
    f"GPMVs: "
    f"W = {shapiro_gpmv.statistic:.4f}, "
    f"p = {shapiro_gpmv.pvalue:.8g}"
)


print(
    f"GUVs: "
    f"W = {shapiro_guv.statistic:.4f}, "
    f"p = {shapiro_guv.pvalue:.8g}"
)


normal_gpmv = (
    shapiro_gpmv.pvalue >= 0.05
)

normal_guv = (
    shapiro_guv.pvalue >= 0.05
)


# ============================================================
# 12. SELECT STATISTICAL TEST
# ============================================================

if normal_gpmv and normal_guv:

    test_name = (
        "Welch's independent t-test"
    )


    test_result = ttest_ind(
        gpmv_diameters,
        guv_diameters,
        equal_var=False
    )


    statistic = (
        test_result.statistic
    )

    p_value = (
        test_result.pvalue
    )

    statistic_name = "t"


else:

    test_name = (
        "Mann-Whitney U test"
    )


    test_result = mannwhitneyu(
        gpmv_diameters,
        guv_diameters,
        alternative="two-sided",
        method="auto"
    )


    statistic = (
        test_result.statistic
    )

    p_value = (
        test_result.pvalue
    )

    statistic_name = "U"


# ============================================================
# 13. DETERMINE SIGNIFICANCE
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
    p_text = f"p = {p_value:.4f}"


# ============================================================
# 14. PRINT STATISTICAL COMPARISON
# ============================================================

print("\n" + "=" * 65)
print("STATISTICAL COMPARISON")
print("=" * 65)


print(
    f"Test used: {test_name}"
)

print(
    f"{statistic_name} = "
    f"{statistic:.4f}"
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
        "\nCONCLUSION:\n"
        "A statistically significant difference "
        "in diameter was detected between "
        "GPMVs and GUVs."
    )

else:

    print(
        "\nCONCLUSION:\n"
        "No statistically significant difference "
        "in diameter was detected between "
        "GPMVs and GUVs."
    )


# ============================================================
# 15. SAVE CLEANED DATA
# ============================================================

gpmv_diameters.to_csv(
    output_folder /
    "GPMV_cleaned_diameters.txt",
    index=False,
    header=[
        "Diameter_um"
    ]
)


guv_diameters.to_csv(
    output_folder /
    "GUV_cleaned_diameters.txt",
    index=False,
    header=[
        "Diameter_um"
    ]
)


# ============================================================
# 16. SAVE STATISTICAL REPORT
# ============================================================

stats_file = (
    output_folder /
    "GPMV_vs_GUV_statistical_analysis.txt"
)


with open(
    stats_file,
    "w",
    encoding="utf-8"
) as f:


    f.write(
        "GPMV vs GUV DIAMETER STATISTICAL ANALYSIS\n"
    )

    f.write(
        "=" * 65 +
        "\n\n"
    )


    # --------------------------------------------------------
    # GPMV
    # --------------------------------------------------------

    f.write(
        "GPMV DIAMETERS\n"
    )

    f.write(
        "-" * 40 +
        "\n"
    )


    for key, value in gpmv_stats.items():

        if key == "N":

            f.write(
                f"{key}: {value}\n"
            )

        else:

            f.write(
                f"{key}: "
                f"{value:.3f} µm\n"
            )


    # --------------------------------------------------------
    # GUV
    # --------------------------------------------------------

    f.write(
        "\nGUV DIAMETERS\n"
    )

    f.write(
        "-" * 40 +
        "\n"
    )


    for key, value in guv_stats.items():

        if key == "N":

            f.write(
                f"{key}: {value}\n"
            )

        else:

            f.write(
                f"{key}: "
                f"{value:.3f} µm\n"
            )


    # --------------------------------------------------------
    # Normality
    # --------------------------------------------------------

    f.write(
        "\nNORMALITY TEST\n"
    )

    f.write(
        "-" * 40 +
        "\n"
    )


    f.write(
        f"GPMVs: "
        f"Shapiro-Wilk W = "
        f"{shapiro_gpmv.statistic:.4f}, "
        f"p = "
        f"{shapiro_gpmv.pvalue:.8g}\n"
    )


    f.write(
        f"GUVs: "
        f"Shapiro-Wilk W = "
        f"{shapiro_guv.statistic:.4f}, "
        f"p = "
        f"{shapiro_guv.pvalue:.8g}\n"
    )


    # --------------------------------------------------------
    # Comparison
    # --------------------------------------------------------

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
        f"{statistic:.4f}\n"
    )


    f.write(
        f"Exact p-value = "
        f"{p_value:.8g}\n"
    )


    f.write(
        f"Significance = "
        f"{significance}\n"
    )


    if p_value < 0.05:

        f.write(
            "\nConclusion: "
            "A statistically significant difference "
            "in diameter was detected between "
            "GPMVs and GUVs.\n"
        )

    else:

        f.write(
            "\nConclusion: "
            "No statistically significant difference "
            "in diameter was detected between "
            "GPMVs and GUVs.\n"
        )


# ============================================================
# 17. CREATE COMPARATIVE BOXPLOT
# ============================================================

fig, ax = plt.subplots(
    figsize=(6, 7)
)


data = [
    gpmv_diameters,
    guv_diameters
]


# ============================================================
# 18. BOXPLOT
# ============================================================

ax.boxplot(
    data,
    widths=0.45,
    showmeans=True
)


# ============================================================
# 19. ADD INDIVIDUAL DATA POINTS
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
# 20. AXIS LABELS
# ============================================================

ax.set_xticks(
    [1, 2]
)


ax.set_xticklabels(
    [
        GPMV_LABEL,
        GUV_LABEL
    ],
    fontsize=11
)


ax.set_ylabel(
    "Equivalent diameter (µm)",
    fontsize=12
)


ax.set_title(
    "GPMV and GUV Diameter Distribution",
    fontsize=13
)


ax.set_ylim(
    bottom=0
)


# ============================================================
# 21. SIGNIFICANCE BRACKET
# ============================================================

combined_max = max(
    gpmv_diameters.max(),
    guv_diameters.max()
)


combined_min = min(
    gpmv_diameters.min(),
    guv_diameters.min()
)


data_range = (
    combined_max -
    combined_min
)


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


# Significance bracket

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


# Stars / ns

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


# Room for significance annotation

ax.set_ylim(
    0,
    bar_y + data_range * 0.20
)


# ============================================================
# 22. FIGURE APPEARANCE
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
# 23. SAVE FIGURE
# ============================================================

graph_file = (
    output_folder /
    "GPMV_vs_GUV_boxplot_significance.png"
)


plt.savefig(
    graph_file,
    dpi=600,
    bbox_inches="tight"
)


plt.show()


# ============================================================
# 24. FINAL SUMMARY
# ============================================================

print("\n" + "=" * 65)
print("ANALYSIS COMPLETE")
print("=" * 65)


print("\nResults saved in:")
print(output_folder.resolve())


print("\nFiles created:")

print(
    "1. GPMV_cleaned_diameters.txt"
)

print(
    "2. GUV_cleaned_diameters.txt"
)

print(
    "3. excluded_GPMV_objects.txt"
)

print(
    "4. GPMV_vs_GUV_statistical_analysis.txt"
)

print(
    "5. GPMV_vs_GUV_boxplot_significance.png"
)