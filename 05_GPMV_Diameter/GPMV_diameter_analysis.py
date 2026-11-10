import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


# ============================================================
# 1. SETTINGS
# ============================================================

# Change this to your TXT file
file_path = Path(
    r"C:\Users\Korisnik\Desktop\GPMVs\GPMVs_epifluo_diameter\MDA-MB-231\20241002_231_GPMV_diameter\gpmv 02102024 diameter and circularity.txt"
)

# ------------------------------------------------------------
# Inclusion criteria for GPMV equivalent diameter
# ------------------------------------------------------------

MIN_DIAMETER = 1.0       # µm
MAX_DIAMETER = 100.0     # µm

# IMPORTANT:
# Use the same inclusion criteria for all comparable datasets
# (for example D1, D2, D3, D4).


# ============================================================
# 2. CREATE OUTPUT FOLDER
# ============================================================

# Creates a folder next to the original TXT file
output_folder = file_path.parent / f"{file_path.stem}_analysis"
output_folder.mkdir(exist_ok=True)

print("\nResults will be saved in:")
print(output_folder.resolve())


# ============================================================
# 3. READ ORIGINAL TXT FILE
# ============================================================

# Read everything initially as text.
# The Nikon file is tab-separated.
df = pd.read_csv(
    file_path,
    sep="\t",
    dtype=str,
    engine="python"
)

# Remove completely empty columns
df = df.dropna(axis=1, how="all")

# Clean column names:
# - remove BOM if present
# - remove spaces
df.columns = (
    df.columns
    .str.replace("\ufeff", "", regex=False)
    .str.strip()
)


print("\nColumns found:")
print(df.columns.tolist())

print("\nFirst rows:")
print(df.head())


# ============================================================
# 4. CHECK THAT EqDiameter EXISTS
# ============================================================

if "EqDiameter" not in df.columns:
    raise ValueError(
        "\n'EqDiameter' column was not found.\n"
        f"Columns found: {df.columns.tolist()}"
    )


# ============================================================
# 5. CLEAN AND CONVERT EqDiameter
# ============================================================

# Keep original diameter column for checking later
df["EqDiameter_original"] = df["EqDiameter"]

# Clean EqDiameter
df["EqDiameter_numeric"] = (
    df["EqDiameter"]
    .astype(str)
    .str.strip()
    .str.replace("\u00A0", "", regex=False)  # remove non-breaking spaces
    .str.replace(" ", "", regex=False)       # remove normal spaces
    .str.replace(",", ".", regex=False)      # 18,89 -> 18.89
)

# Convert to real numeric values
df["EqDiameter_numeric"] = pd.to_numeric(
    df["EqDiameter_numeric"],
    errors="coerce"
)


# ============================================================
# 6. CHECK THE CONVERSION
# ============================================================

print("\nFirst converted diameter values:")

print(
    df[
        ["EqDiameter_original", "EqDiameter_numeric"]
    ].head(15)
)


# ============================================================
# 7. SHOW LARGEST VALUES BEFORE FILTERING
# ============================================================

print("\n" + "=" * 60)
print("15 LARGEST DIAMETERS BEFORE FILTERING")
print("=" * 60)

columns_to_show = []

for column in ["Item", "FieldID", "ObjID",
               "EqDiameter_original", "EqDiameter_numeric"]:

    if column in df.columns:
        columns_to_show.append(column)


print(
    df[columns_to_show]
    .sort_values(
        "EqDiameter_numeric",
        ascending=False
    )
    .head(15)
    .to_string(index=False)
)


# ============================================================
# 8. REMOVE NaN VALUES
# ============================================================

valid_numeric_df = df[
    df["EqDiameter_numeric"].notna()
].copy()


# ============================================================
# 9. APPLY DIAMETER CUTOFF
# ============================================================

# Objects outside the predefined range
excluded_df = valid_numeric_df[
    (valid_numeric_df["EqDiameter_numeric"] < MIN_DIAMETER) |
    (valid_numeric_df["EqDiameter_numeric"] > MAX_DIAMETER)
].copy()


# Objects retained for analysis
filtered_df = valid_numeric_df[
    (valid_numeric_df["EqDiameter_numeric"] >= MIN_DIAMETER) &
    (valid_numeric_df["EqDiameter_numeric"] <= MAX_DIAMETER)
].copy()


# Final diameter Series used for graph + statistics
diameters = filtered_df["EqDiameter_numeric"]


# ============================================================
# 10. REPORT FILTERING
# ============================================================

print("\n" + "=" * 60)
print("FILTERING SUMMARY")
print("=" * 60)

print(f"Rows in original file:          {len(df)}")
print(f"Valid numeric diameters:        {len(valid_numeric_df)}")
print(f"Measurements retained:          {len(filtered_df)}")
print(f"Measurements excluded:          {len(excluded_df)}")

print(
    f"\nDiameter inclusion range: "
    f"{MIN_DIAMETER:.1f}–{MAX_DIAMETER:.1f} µm"
)


if len(excluded_df) > 0:

    print("\nExcluded objects:")

    excluded_columns = []

    for column in [
        "Item",
        "Source",
        "FieldID",
        "ObjID",
        "EqDiameter_original",
        "EqDiameter_numeric"
    ]:

        if column in excluded_df.columns:
            excluded_columns.append(column)

    print(
        excluded_df[excluded_columns]
        .sort_values(
            "EqDiameter_numeric",
            ascending=False
        )
        .to_string(index=False)
    )


# ============================================================
# 11. SAVE EXCLUDED OBJECTS
# ============================================================

excluded_output = output_folder / "excluded_GPMV_objects.txt"

excluded_df.to_csv(
    excluded_output,
    sep="\t",
    index=False
)


# ============================================================
# 12. CHECK THAT DATA REMAIN AFTER FILTERING
# ============================================================

if len(diameters) == 0:

    raise ValueError(
        "\nNo GPMV diameters remain after filtering.\n"
        "Check MIN_DIAMETER and MAX_DIAMETER."
    )


# ============================================================
# 13. SAVE CLEAN DIAMETER DATA
# ============================================================

diameter_output = output_folder / "GPMV_diameters_filtered.txt"

diameters.to_csv(
    diameter_output,
    index=False,
    header=["EqDiameter_um"]
)


# Also save the complete retained object information

filtered_objects_output = (
    output_folder / "GPMV_objects_retained.txt"
)

filtered_df.to_csv(
    filtered_objects_output,
    sep="\t",
    index=False
)


# ============================================================
# 14. DESCRIPTIVE STATISTICS
# ============================================================

n = len(diameters)

mean = diameters.mean()
median = diameters.median()

# Sample standard deviation
std = diameters.std(ddof=1)

# Standard error of the mean
sem = std / np.sqrt(n)

minimum = diameters.min()
maximum = diameters.max()

q1 = diameters.quantile(0.25)
q3 = diameters.quantile(0.75)

iqr = q3 - q1

cv = (std / mean) * 100


statistics = {
    "N retained": n,
    "N excluded": len(excluded_df),

    "Mean (µm)": mean,
    "Median (µm)": median,

    "Standard deviation (µm)": std,
    "SEM (µm)": sem,

    "Minimum (µm)": minimum,
    "Q1 (µm)": q1,
    "Q3 (µm)": q3,
    "Maximum (µm)": maximum,

    "IQR (µm)": iqr,

    "Coefficient of variation (%)": cv
}


# ============================================================
# 15. PRINT STATISTICS
# ============================================================

print("\n" + "=" * 60)
print("GPMV DIAMETER STATISTICS")
print("=" * 60)

for name, value in statistics.items():

    if name.startswith("N"):
        print(f"{name}: {value}")

    else:
        print(f"{name}: {value:.3f}")


# ============================================================
# 16. SAVE STATISTICS AS TXT
# ============================================================

stats_txt = output_folder / "GPMV_diameter_statistics.txt"


with open(stats_txt, "w", encoding="utf-8") as f:

    f.write("GPMV DIAMETER STATISTICAL ANALYSIS\n")
    f.write("=" * 50 + "\n\n")

    f.write(
        f"Diameter inclusion range: "
        f"{MIN_DIAMETER:.1f}–{MAX_DIAMETER:.1f} µm\n\n"
    )

    for name, value in statistics.items():

        if name.startswith("N"):
            f.write(f"{name}: {value}\n")

        else:
            f.write(f"{name}: {value:.3f}\n")


# ============================================================
# 17. CREATE BOXPLOT
# ============================================================

fig, ax = plt.subplots(figsize=(5, 7))


# Boxplot
ax.boxplot(
    diameters,
    widths=0.40,
    showmeans=True
)


# ------------------------------------------------------------
# Add individual GPMV measurements
# ------------------------------------------------------------

# Fixed seed = same jitter arrangement each time
rng = np.random.default_rng(42)

x_jitter = rng.normal(
    loc=1,
    scale=0.025,
    size=n
)

ax.scatter(
    x_jitter,
    diameters,
    alpha=0.55,
    s=25
)


# ------------------------------------------------------------
# Labels
# ------------------------------------------------------------

ax.set_ylabel(
    "Equivalent diameter (µm)",
    fontsize=12
)

ax.set_xticks([1])
ax.set_xticklabels(["GPMVs"], fontsize=11)

ax.set_title(
    f"GPMV Diameter Distribution (n = {n})",
    fontsize=13
)


# ------------------------------------------------------------
# Start y-axis from 0 µm
#
# The upper limit is automatically determined from the
# FILTERED data, so the graph will no longer extend to 16,000.
# ------------------------------------------------------------

ax.set_ylim(bottom=0)


# Remove unnecessary top/right borders if desired
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)


plt.tight_layout()


# ============================================================
# 18. SAVE BOXPLOT
# ============================================================

boxplot_path = output_folder / "GPMV_diameter_boxplot.png"

plt.savefig(
    boxplot_path,
    dpi=600,
    bbox_inches="tight"
)

plt.show()


# ============================================================
# 19. STATISTICS TABLE FOR PNG
# ============================================================

table_parameters = [
    "N retained",
    "N excluded",
    "Mean (µm)",
    "Median (µm)",
    "Standard deviation (µm)",
    "SEM (µm)",
    "Minimum (µm)",
    "Q1 (µm)",
    "Q3 (µm)",
    "Maximum (µm)",
    "IQR (µm)",
    "Coefficient of variation (%)"
]


table_values = []

for parameter in table_parameters:

    value = statistics[parameter]

    if parameter.startswith("N"):
        table_values.append(str(value))

    else:
        table_values.append(f"{value:.3f}")


stats_df = pd.DataFrame({
    "Parameter": table_parameters,
    "Value": table_values
})


# ============================================================
# 20. CREATE STATISTICS PNG
# ============================================================

fig, ax = plt.subplots(figsize=(7, 6))

ax.axis("off")


table = ax.table(
    cellText=stats_df.values,
    colLabels=stats_df.columns,
    cellLoc="center",
    loc="center"
)

table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1, 1.4)


ax.set_title(
    "GPMV Diameter Descriptive Statistics\n"
    f"Inclusion range: "
    f"{MIN_DIAMETER:.1f}–{MAX_DIAMETER:.1f} µm",
    fontsize=13,
    pad=15
)


plt.tight_layout()


# ============================================================
# 21. SAVE STATISTICS PNG
# ============================================================

stats_png = output_folder / "GPMV_diameter_statistics.png"

plt.savefig(
    stats_png,
    dpi=600,
    bbox_inches="tight"
)

plt.show()


# ============================================================
# 22. FINAL OUTPUT SUMMARY
# ============================================================

print("\n" + "=" * 60)
print("ANALYSIS COMPLETE")
print("=" * 60)

print("\nResults saved here:")
print(output_folder.resolve())

print("\nFiles created:")

print(f"1. {diameter_output.name}")
print(f"2. {filtered_objects_output.name}")
print(f"3. {excluded_output.name}")
print(f"4. {stats_txt.name}")
print(f"5. {boxplot_path.name}")
print(f"6. {stats_png.name}")