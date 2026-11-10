import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.signal import find_peaks
from scipy.interpolate import make_interp_spline


# -----------------------------
# 1) Set your input/output
# -----------------------------
input_file = r"C:\Users\Korisnik\Desktop\LUVs\2mgml\20251127_Statistical_analysis\Ana_33chol_2mgml0003.txt"
output_folder = r"C:\Users\Korisnik\Desktop\LUVs\2mgml\20251127_Statistical_analysis\heightdistribution_corrected"

# ensure output folder exists
Path(output_folder).mkdir(parents=True, exist_ok=True)

# base name for output files
output_prefix = Path(output_folder) / "SLB33_2mgml_height_distribution_analysis"  # HERE CHANGE THE NAME FOR THE ONE YOU WANT


# -----------------------------
# 2) Load data, skip 3 header rows
# (as it is when one exports txt data from Gwyddion)
# -----------------------------
df = pd.read_csv(
    input_file,
    sep=r"\s+",
    skiprows=3,
    names=["z_m", "rho"]
)

# convert z from meters to nm
df["z_nm"] = df["z_m"] * 1e9


# -----------------------------
# 3) DESCRIPTIVE STATISTICS
# -----------------------------
z = df["z_nm"].values
rho = df["rho"].values

stats_dict = {
    "n": len(z),
    "mean (nm)": np.mean(z),
    "median (nm)": np.median(z),
    "std (nm)": np.std(z, ddof=1),
    "min (nm)": np.min(z),
    "max (nm)": np.max(z),
    "25% (nm)": np.percentile(z, 25),
    "75% (nm)": np.percentile(z, 75),
    "IQR (nm)": np.percentile(z, 75) - np.percentile(z, 25)
}


# -----------------------------
# 4) DETECT DISTRIBUTION PEAKS
# -----------------------------
peaks, properties = find_peaks(
    rho,
    prominence=np.max(rho) * 0.02
)

# First select the two strongest peaks based on rho height
sorted_peaks_by_rho = sorted(
    peaks,
    key=lambda p: rho[p],
    reverse=True
)

top_peaks = sorted_peaks_by_rho[:2]

# IMPORTANT:
# Reorder selected peaks according to z-height.
#
# Peak 1 = lower z-height
# Peak 2 = higher z-height
top_peaks = sorted(
    top_peaks,
    key=lambda p: z[p]
)


# -----------------------------
# 4a) ADD PEAK INFORMATION
# -----------------------------
for i, p in enumerate(top_peaks, start=1):
    stats_dict[f"Peak {i} z (nm)"] = z[p]
    stats_dict[f"Peak {i} ρ"] = rho[p]


# Calculate height difference:
# Peak 2 - Peak 1
if len(top_peaks) >= 2:

    peak1 = top_peaks[0]
    peak2 = top_peaks[1]

    peak_difference = z[peak2] - z[peak1]

    stats_dict["Peak-to-peak Δz (nm)"] = peak_difference

else:
    peak_difference = None


# -----------------------------
# 5) PRINT AND SAVE STATISTICS
# -----------------------------
print("\n--- Height Distribution Statistics ---")

for k, v in stats_dict.items():
    print(
        f"{k:25s}: {v:.4g}"
        if isinstance(v, (int, float, np.integer, np.floating))
        else f"{k:25s}: {v}"
    )

stats_file = Path(output_folder) / f"{output_prefix.stem}_stats.csv"

pd.DataFrame([stats_dict]).to_csv(
    stats_file,
    index=False
)

print(f"\nSaved statistics to: {stats_file}")


# -----------------------------
# 6) PLOT SMOOTH HEIGHT DISTRIBUTION
# -----------------------------
plt.figure(figsize=(8, 5))


# -----------------------------
# Smooth interpolation
# -----------------------------
z_smooth = np.linspace(
    z.min(),
    z.max(),
    300
)

rho_smooth = make_interp_spline(
    z,
    rho,
    k=3
)(z_smooth)


plt.plot(
    z_smooth,
    rho_smooth,
    color="blue",
    linewidth=2,
    label="Height Distribution"
)


# -----------------------------
# Mark detected peaks
# -----------------------------
if len(top_peaks) > 0:

    plt.scatter(
        z[top_peaks],
        rho[top_peaks],
        color="blue",
        zorder=5,
        label="Detected Peaks"
    )

    for i, p in enumerate(top_peaks, start=1):

        # Peak labels
        plt.text(
            z[p],
            rho[p] * 1.05,
            f"Peak {i}",
            color="blue",
            ha="center",
            va="bottom",
            fontsize=9
        )


# -----------------------------
# Peak-to-peak height difference
# -----------------------------
if len(top_peaks) >= 2:

    peak1 = top_peaks[0]
    peak2 = top_peaks[1]

    # Peak positions along the z-axis
    z1 = z[peak1]
    z2 = z[peak2]

    # Peak rho values
    rho1 = rho[peak1]
    rho2 = rho[peak2]

    # Difference is always:
    # higher z peak - lower z peak
    peak_difference = z2 - z1


    # -----------------------------
    # Position measurement bracket
    # -----------------------------

    # Put the horizontal line above both peaks
    y_measure = max(rho1, rho2) * 1.18


    # -----------------------------
    # Vertical line through Peak 1
    # -----------------------------
    plt.vlines(
        x=z1,
        ymin=0,
        ymax=y_measure,
        color="black",
        linestyle="--",
        linewidth=1
    )


    # -----------------------------
    # Vertical line through Peak 2
    # -----------------------------
    plt.vlines(
        x=z2,
        ymin=0,
        ymax=y_measure,
        color="black",
        linestyle="--",
        linewidth=1
    )


    # -----------------------------
    # Horizontal connecting line
    # -----------------------------
    plt.hlines(
        y=y_measure,
        xmin=z1,
        xmax=z2,
        color="black",
        linewidth=1.5
    )


    # -----------------------------
    # Δz label
    # -----------------------------
    plt.text(
        (z1 + z2) / 2,
        y_measure * 1.025,
        f"Δz = {peak_difference:.2f} nm",
        ha="center",
        va="bottom",
        fontsize=10
    )


# -----------------------------
# Plot formatting
# -----------------------------
plt.xlabel("Height z [nm]")
plt.ylabel("ρ")
plt.title("Height Distribution")

# Keep the horizontal scale identical across height-distribution plots.
# Change these values only if a wider height range needs to be shown.
X_AXIS_MIN_NM = 0
X_AXIS_MAX_NM = 5
X_AXIS_TICK_SPACING_NM = 1
plt.xlim(X_AXIS_MIN_NM, X_AXIS_MAX_NM)
plt.xticks(
    np.arange(
        X_AXIS_MIN_NM,
        X_AXIS_MAX_NM + X_AXIS_TICK_SPACING_NM,
        X_AXIS_TICK_SPACING_NM
    )
)

plt.grid(False)

plt.legend()

# Extra vertical space for peak labels
# and the Δz measurement bracket
plt.ylim(
    bottom=0,
    top=np.max(rho) * 1.40
)

plt.tight_layout()


# -----------------------------
# Save figure
# -----------------------------
line_plot_file = (
    Path(output_folder)
    / f"{output_prefix.stem}_line.png"
)

plt.savefig(
    line_plot_file,
    dpi=300,
    bbox_inches="tight"
)

plt.show()

print(f"Saved line graph to: {line_plot_file}")