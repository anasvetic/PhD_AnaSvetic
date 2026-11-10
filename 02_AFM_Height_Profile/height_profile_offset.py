import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import curve_fit
from scipy.stats import gaussian_kde
from scipy.signal import find_peaks
from pathlib import Path
import os


# ============================================================
# AFM HEIGHT PROFILE ANALYSIS
#
# Input:
#   Raw TXT height-profile data exported from Gwyddion.
#
# Processing:
#   - Skip the first 3 Gwyddion header rows
#   - Convert lateral distance from m to µm
#   - Convert height from m to nm
#   - Apply a vertical offset so that the minimum height = 0 nm
#   - Generate histogram and KDE representations
#   - Fit the height distribution using a two-component Gaussian model
#
# Output:
#   - Processed height-profile CSV
#   - Offset height-profile plot
#   - Height histogram
#   - KDE plot
#   - Bimodal Gaussian fit plot
# ============================================================


# ------------------------------------------------------------
# 1) SET INPUT AND OUTPUT PATHS
# ------------------------------------------------------------

file_path = r"copyandpasteyourheightprofilepath.txt"

# Folder where all output files will be saved
output_dir = r"copyandpasteyourpathhere"

# Ensure that the output folder exists
os.makedirs(output_dir, exist_ok=True)

# Derive the base filename from the input file
stem = Path(file_path).stem


# ------------------------------------------------------------
# 2) LOAD RAW HEIGHT-PROFILE DATA FROM GWYDDION
# ------------------------------------------------------------

# Gwyddion TXT exports contain 3 header rows
skip_lines = 3

# Expected columns:
# X_m = lateral position in metres
# Height_m = measured height in metres
col_names = ["X_m", "Height_m"]

df = pd.read_csv(
    file_path,
    skiprows=skip_lines,
    sep=r"\s+",
    header=None,
    names=col_names
)


# ------------------------------------------------------------
# 3) CONVERT UNITS
# ------------------------------------------------------------

# Lateral coordinate:
# metres -> micrometres
df["X_um"] = df["X_m"] * 1e6

# Height:
# metres -> nanometres
df["Height_nm"] = df["Height_m"] * 1e9


# ------------------------------------------------------------
# 4) APPLY VERTICAL HEIGHT OFFSET
# ------------------------------------------------------------

# The minimum measured height is used as a zero reference.
# This does NOT alter relative height differences within the
# profile; it simply shifts the complete profile vertically.

min_height = df["Height_nm"].min()

df["Height_nm_offset"] = (
    df["Height_nm"] - min_height
)


# ------------------------------------------------------------
# 5) SAVE PROCESSED DATA AS CSV
# ------------------------------------------------------------

csv_fname = (
    Path(output_dir)
    / f"{stem}_processed_height_profile.csv"
)

df.to_csv(
    csv_fname,
    index=False
)

print(
    f"Saved processed height-profile data to:\n{csv_fname}"
)


# ------------------------------------------------------------
# 6) PLOT OFFSET HEIGHT PROFILE
# ------------------------------------------------------------

plt.figure(figsize=(8, 5))

plt.plot(
    df["X_um"],
    df["Height_nm_offset"],
    color="darkblue",
    linewidth=2
)

# X corresponds to lateral distance along the selected profile.
# The offset is applied only to the height axis.
plt.xlabel("Distance (µm)")
plt.ylabel("Height (nm)")

# Change title if needed depending on sample type
# e.g. "SLB Height Profile" or "GUV Patch Height Profile"
plt.title("Height Profile")

plt.grid(False)


# Identify minimum and maximum points of the offset profile
min_idx = df["Height_nm_offset"].idxmin()
max_idx = df["Height_nm_offset"].idxmax()

min_y = df.loc[min_idx, "Height_nm_offset"]  # = 0 nm
max_y = df.loc[max_idx, "Height_nm_offset"]


# Horizontal reference lines
plt.axhline(
    y=min_y,
    color="red",
    linestyle="--",
    linewidth=1.2,
    label=f"Min: {min_y:.1f} nm"
)

plt.axhline(
    y=max_y,
    color="green",
    linestyle="--",
    linewidth=1.2,
    label=f"Max: {max_y:.1f} nm"
)

plt.legend()
plt.tight_layout()


# Save offset profile
profile_fname = (
    Path(output_dir)
    / f"{stem}_profile_offset.png"
)

plt.savefig(
    profile_fname,
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print(
    f"Saved offset height-profile plot to:\n{profile_fname}"
)


# ------------------------------------------------------------
# 7) HISTOGRAM OF ORIGINAL HEIGHT VALUES
# ------------------------------------------------------------

# The original, non-offset height values are used for the
# distribution analysis because applying a constant offset would
# only translate the distribution and would not change its shape.

plt.figure(figsize=(8, 5))

plt.hist(
    df["Height_nm"],
    bins=50,
    color="skyblue",
    edgecolor="black"
)

plt.xlabel("Height (nm)")
plt.ylabel("Frequency")
plt.title("Distribution of Height Profile Values")
plt.grid(False)
plt.tight_layout()


hist_fname = (
    Path(output_dir)
    / f"{stem}_histogram.png"
)

plt.savefig(
    hist_fname,
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print(
    f"Saved histogram to:\n{hist_fname}"
)


# ------------------------------------------------------------
# 8) KERNEL DENSITY ESTIMATION (KDE)
# ------------------------------------------------------------

plt.figure(figsize=(8, 5))

sns.kdeplot(
    df["Height_nm"],
    fill=True,
    color="teal",
    linewidth=2
)

plt.xlabel("Height (nm)")
plt.ylabel("Density")
plt.title("Smoothed Height Distribution")
plt.grid(False)
plt.tight_layout()


kde_fname = (
    Path(output_dir)
    / f"{stem}_kde.png"
)

plt.savefig(
    kde_fname,
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print(
    f"Saved KDE plot to:\n{kde_fname}"
)


# ------------------------------------------------------------
# 9) TWO-COMPONENT GAUSSIAN FIT
# ------------------------------------------------------------

# The Gaussian fit is applied to the KDE representation of the
# original height values.
#
# This analysis assumes that the selected profile contains two
# predominant height populations. Profiles containing more than
# two populations may not be adequately represented by this model.

data = df["Height_nm"].values

# Generate KDE
kde = gaussian_kde(data)

x_vals = np.linspace(
    np.min(data),
    np.max(data),
    1000
)

kde_vals = kde(x_vals)


# Define a two-component Gaussian function
def two_gaussian(
        x,
        amp1,
        mu1,
        sigma1,
        amp2,
        mu2,
        sigma2
):
    gaussian_1 = (
        amp1
        * np.exp(
            -((x - mu1) ** 2)
            / (2 * sigma1 ** 2)
        )
    )

    gaussian_2 = (
        amp2
        * np.exp(
            -((x - mu2) ** 2)
            / (2 * sigma2 ** 2)
        )
    )

    return gaussian_1 + gaussian_2


# ------------------------------------------------------------
# Initial parameter estimates
#
# These values may require adjustment for unusual datasets.
# ------------------------------------------------------------

initial_guess = [
    np.max(kde_vals),
    x_vals[np.argmax(kde_vals)],
    5,

    np.max(kde_vals) / 2,
    x_vals[len(x_vals) // 2],
    5
]


# Fit the two-component Gaussian model
popt, _ = curve_fit(
    two_gaussian,
    x_vals,
    kde_vals,
    p0=initial_guess,
    maxfev=7000
)


# Generate fitted curve
fitted_curve = two_gaussian(
    x_vals,
    *popt
)


# Identify local maxima in the fitted curve
peaks, _ = find_peaks(
    fitted_curve
)


# ------------------------------------------------------------
# 10) PLOT KDE AND BIMODAL GAUSSIAN FIT
# ------------------------------------------------------------

plt.figure(figsize=(8, 5))

plt.plot(
    x_vals,
    kde_vals,
    label="Smoothed KDE",
    color="blue"
)

plt.plot(
    x_vals,
    fitted_curve,
    label="Two-component Gaussian fit",
    color="red",
    linestyle="--"
)

plt.scatter(
    x_vals[peaks],
    fitted_curve[peaks],
    color="green",
    s=50,
    label="Detected peaks"
)

plt.xlabel("Height (nm)")
plt.ylabel("Density")
plt.title(f"{stem}: Two-component Gaussian Fit")

plt.legend()
plt.grid(False)
plt.tight_layout()


fit_fname = (
    Path(output_dir)
    / f"{stem}_bimodal_fit.png"
)

plt.savefig(
    fit_fname,
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print(
    f"Saved bimodal Gaussian fit plot to:\n{fit_fname}"
)


# ------------------------------------------------------------
# 11) PRINT COMPLETION MESSAGE
# ------------------------------------------------------------

print("\nAnalysis complete.")
print(f"All output files were saved in:\n{output_dir}")