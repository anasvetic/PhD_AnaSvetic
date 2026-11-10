import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from pathlib import Path
import tifffile
from tkinter import Tk, filedialog


# ============================================================
# SETTINGS
# ============================================================

# Percentiles used for the displayed fluorescence range.
# Using 1–99% prevents a few extreme pixels from dominating
# the intensity scale.
LOW_PERCENTILE = 1
HIGH_PERCENTILE = 99

# Number of labels shown on the fluorescence intensity bar
NUMBER_OF_TICKS = 5


# ============================================================
# 1. SELECT TIFF FILE
# ============================================================

root = Tk()
root.withdraw()
root.attributes("-topmost", True)

image_path = filedialog.askopenfilename(
    title="Select TIRF .TIF image",
    filetypes=[
        ("TIFF files", "*.tif *.tiff"),
        ("All files", "*.*")
    ]
)

root.destroy()

if not image_path:
    raise ValueError("No image was selected.")

image_path = Path(image_path)

# Folder containing the original image
output_dir = image_path.parent

# Original filename without .tif
base_name = image_path.stem

print("\nSelected image:")
print(image_path)

print("\nOutputs will be saved in:")
print(output_dir)


# ============================================================
# 2. LOAD TIFF IMAGE
# ============================================================

img = tifffile.imread(str(image_path))

print("\nIMAGE INFORMATION")
print("----------------------------------------")
print("Shape:", img.shape)
print("Data type:", img.dtype)
print("Dimensions:", img.ndim)


# ============================================================
# 3. EXTRACT FLUORESCENCE INTENSITY
# ============================================================

# ------------------------------------------------------------
# CASE 1: 2D grayscale image
# ------------------------------------------------------------

if img.ndim == 2:

    intensity = img.astype(np.float64)

    # The grayscale data will be displayed using a red LUT.
    display_image = intensity

    display_mode = "grayscale"


# ------------------------------------------------------------
# CASE 2: 3D image
#
# Could be:
#   RGB = (height, width, 3)
#   RGBA = (height, width, 4)
#   grayscale Z-stack = (z, height, width)
# ------------------------------------------------------------

elif img.ndim == 3:

    # RGB / RGBA image
    if img.shape[-1] in (3, 4):

        # Use the red channel as fluorescence intensity
        intensity = img[..., 0].astype(np.float64)

        # Keep the RGB image for display
        display_image = img[..., :3].astype(np.float64)

        display_mode = "rgb"

    else:

        # Assume grayscale Z-stack.
        # Average all Z slices.
        intensity = np.mean(
            img.astype(np.float64),
            axis=0
        )

        display_image = intensity

        display_mode = "grayscale"


# ------------------------------------------------------------
# CASE 3: 4D image
#
# Example:
#   RGB Z-stack = (z, height, width, 3)
# ------------------------------------------------------------

elif img.ndim == 4:

    if img.shape[-1] in (3, 4):

        # Average the red channel through the Z-stack
        intensity = np.mean(
            img[..., 0].astype(np.float64),
            axis=0
        )

        # Average RGB channels through the Z-stack
        display_image = np.mean(
            img[..., :3].astype(np.float64),
            axis=0
        )

        display_mode = "rgb"

    else:

        raise ValueError(
            f"Unsupported 4D TIFF format. Image shape: {img.shape}"
        )

else:

    raise ValueError(
        f"Unsupported TIFF dimensions: {img.ndim}"
    )


print("\nProcessed intensity image shape:")
print(intensity.shape)


# ============================================================
# 4. REMOVE NaN / INFINITE VALUES IF PRESENT
# ============================================================

valid_pixels = intensity[np.isfinite(intensity)]

if valid_pixels.size == 0:
    raise ValueError("Image contains no valid intensity pixels.")


# ============================================================
# 5. CALCULATE FLUORESCENCE INTENSITY STATISTICS
# ============================================================

minimum_intensity = np.min(valid_pixels)
maximum_intensity = np.max(valid_pixels)

mean_intensity = np.mean(valid_pixels)
median_intensity = np.median(valid_pixels)
std_intensity = np.std(valid_pixels)

low_intensity = np.percentile(
    valid_pixels,
    LOW_PERCENTILE
)

high_intensity = np.percentile(
    valid_pixels,
    HIGH_PERCENTILE
)


print("\nFLUORESCENCE INTENSITY STATISTICS")
print("----------------------------------------")
print(f"Minimum:              {minimum_intensity:.2f}")
print(f"Maximum:              {maximum_intensity:.2f}")
print(f"Mean:                 {mean_intensity:.2f}")
print(f"Median:               {median_intensity:.2f}")
print(f"Standard deviation:   {std_intensity:.2f}")
print(
    f"{LOW_PERCENTILE}th percentile:       "
    f"{low_intensity:.2f}"
)
print(
    f"{HIGH_PERCENTILE}th percentile:      "
    f"{high_intensity:.2f}"
)


# Avoid problems if the image has no meaningful intensity range
if high_intensity <= low_intensity:
    low_intensity = minimum_intensity
    high_intensity = maximum_intensity


# ============================================================
# 6. CREATE BLACK → RED FLUORESCENCE COLORMAP
# ============================================================

# This mimics the red LUT appearance:
#
# low fluorescence  = black
# high fluorescence = bright red

black_red_cmap = LinearSegmentedColormap.from_list(
    "black_to_red",
    [
        (0.0, 0.0, 0.0),   # black
        (1.0, 0.0, 0.0)    # red
    ]
)


# ============================================================
# 7. SAVE STATISTICS TO TXT FILE
# ============================================================

stats_path = output_dir / f"{base_name}_intensity_stats.txt"

with open(stats_path, "w", encoding="utf-8") as file:

    file.write("TIRF FLUORESCENCE INTENSITY ANALYSIS\n")
    file.write("========================================\n\n")

    file.write(f"Original file:\n{image_path}\n\n")

    file.write("IMAGE INFORMATION\n")
    file.write("----------------------------------------\n")
    file.write(f"Shape: {img.shape}\n")
    file.write(f"Data type: {img.dtype}\n")
    file.write(f"Display mode detected: {display_mode}\n\n")

    file.write("FLUORESCENCE INTENSITY STATISTICS\n")
    file.write("----------------------------------------\n")

    file.write(
        f"Minimum intensity: {minimum_intensity:.4f}\n"
    )

    file.write(
        f"Maximum intensity: {maximum_intensity:.4f}\n"
    )

    file.write(
        f"Mean intensity: {mean_intensity:.4f}\n"
    )

    file.write(
        f"Median intensity: {median_intensity:.4f}\n"
    )

    file.write(
        f"Standard deviation: {std_intensity:.4f}\n"
    )

    file.write(
        f"{LOW_PERCENTILE}th percentile: "
        f"{low_intensity:.4f}\n"
    )

    file.write(
        f"{HIGH_PERCENTILE}th percentile: "
        f"{high_intensity:.4f}\n"
    )

    file.write("\n")
    file.write(
        "Displayed intensity range: "
        f"{low_intensity:.4f} – "
        f"{high_intensity:.4f}\n"
    )


# ============================================================
# 8. CREATE THE FLUORESCENCE INTENSITY BAR
# ============================================================

# Creates a continuous gradient containing the intensity range
gradient = np.linspace(
    low_intensity,
    high_intensity,
    1000
).reshape(-1, 1)


# ============================================================
# 9. SAVE INTENSITY BAR AS A SEPARATE IMAGE
# ============================================================

fig_bar, ax_bar = plt.subplots(
    figsize=(2.2, 6)
)

ax_bar.imshow(
    gradient,
    cmap=black_red_cmap,
    origin="lower",
    aspect="auto",
    vmin=low_intensity,
    vmax=high_intensity
)

# Remove horizontal axis
ax_bar.set_xticks([])

# Create intensity labels
tick_positions = np.linspace(
    0,
    len(gradient) - 1,
    NUMBER_OF_TICKS
)

tick_values = np.linspace(
    low_intensity,
    high_intensity,
    NUMBER_OF_TICKS
)

ax_bar.set_yticks(tick_positions)

ax_bar.set_yticklabels(
    [f"{value:.1f}" for value in tick_values]
)

ax_bar.yaxis.tick_right()
ax_bar.yaxis.set_label_position("right")

ax_bar.set_ylabel(
    "Fluorescence intensity (a.u.)",
    fontsize=11,
    labelpad=12
)

# Remove box around bar
for spine in ax_bar.spines.values():
    spine.set_visible(False)

bar_path = output_dir / (
    f"{base_name}_intensity_bar.png"
)

plt.savefig(
    bar_path,
    dpi=600,
    bbox_inches="tight",
    transparent=False
)

plt.close(fig_bar)


# ============================================================
# 10. CREATE FINAL IMAGE + INTENSITY BAR
# ============================================================

fig = plt.figure(
    figsize=(8, 6)
)

# Image occupies most of the figure.
# Bar is narrow and positioned beside it.
grid = fig.add_gridspec(
    1,
    2,
    width_ratios=[20, 1],
    wspace=0.08
)

ax_image = fig.add_subplot(grid[0])
ax_bar = fig.add_subplot(grid[1])


# ============================================================
# 11. DISPLAY TIRF IMAGE
# ============================================================

if display_mode == "grayscale":

    # IMPORTANT:
    # This does NOT alter the original TIFF.
    #
    # It only displays the grayscale intensities using
    # black → red, reproducing the red LUT appearance.

    ax_image.imshow(
        display_image,
        cmap=black_red_cmap,
        vmin=low_intensity,
        vmax=high_intensity
    )


elif display_mode == "rgb":

    # The TIFF already contains RGB information.
    # Normalize only for matplotlib display.
    #
    # This does not alter the original file.

    rgb = display_image.copy()

    if np.issubdtype(img.dtype, np.integer):
        max_possible = np.iinfo(img.dtype).max
    else:
        max_possible = np.max(rgb)

    if max_possible > 0:
        rgb = rgb / max_possible

    rgb = np.clip(
        rgb,
        0,
        1
    )

    ax_image.imshow(rgb)


# Remove axes around fluorescence image
ax_image.axis("off")


# ============================================================
# 12. DISPLAY MATCHING FLUORESCENCE BAR
# ============================================================

ax_bar.imshow(
    gradient,
    cmap=black_red_cmap,
    origin="lower",
    aspect="auto",
    vmin=low_intensity,
    vmax=high_intensity
)

ax_bar.set_xticks([])

ax_bar.set_yticks(
    tick_positions
)

ax_bar.set_yticklabels(
    [f"{value:.1f}" for value in tick_values],
    fontsize=9
)

ax_bar.yaxis.tick_right()
ax_bar.yaxis.set_label_position("right")

ax_bar.set_ylabel(
    "Fluorescence intensity (a.u.)",
    fontsize=10,
    labelpad=10
)

for spine in ax_bar.spines.values():
    spine.set_visible(False)


# ============================================================
# 13. SAVE FINAL IMAGE
# ============================================================

combined_path = output_dir / (
    f"{base_name}_with_intensity_bar.png"
)

plt.savefig(
    combined_path,
    dpi=600,
    bbox_inches="tight",
    pad_inches=0.05,
    facecolor="white"
)

plt.close(fig)


# ============================================================
# 14. FINISHED
# ============================================================

print("\n========================================")
print("ANALYSIS COMPLETE")
print("========================================")

print("\nFiles saved:")

print(
    "\n1. Intensity statistics:\n"
    f"   {stats_path}"
)

print(
    "\n2. Fluorescence intensity bar:\n"
    f"   {bar_path}"
)

print(
    "\n3. TIRF image with intensity bar:\n"
    f"   {combined_path}"
)

print("\nOriginal TIFF was NOT modified.")