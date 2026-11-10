# Calculate % of Lo and Ld surface area, export CSV, plot bar + pie charts, and save them.


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os


# -----------------------------
# 1. INPUT SECTION
# -----------------------------
# Change these values as needed:
Lo_area = 0.85   # µm²
Ld_area = 13.67  # µm²
total_area = 14.52  # µm², or set to None to auto-sum Lo+Ld

# Output settings
output_dir = r"C:\Users\Korisnik\Desktop\temperature_cycling\GPMV_patch"
file_prefix = "20261023_GPMV_FINAL_phasesurfacearea"
dpi = 300
round_digits = 2

# -----------------------------
# 2. CALCULATE PERCENTAGES
# -----------------------------
if total_area is None:
    total_area = Lo_area + Ld_area

pct_Lo = (Lo_area / total_area) * 100
pct_Ld = (Ld_area / total_area) * 100

# Create dataframe
df = pd.DataFrame({
    "Phase": ["Lo", "Ld"],
    "Area_µm²": [Lo_area, Ld_area],
    "Percentage": [round(pct_Lo, round_digits), round(pct_Ld, round_digits)]
})

# -----------------------------
# 3. EXPORT TO CSV
# -----------------------------
os.makedirs(output_dir, exist_ok=True)
csv_path = os.path.join(output_dir, f"{file_prefix}.csv")
df.to_csv(csv_path, index=False)
print(f"✅ Data saved to: {csv_path}")

# -----------------------------
# 4. PLOTS
# -----------------------------

# Bar chart
fig_bar, ax_bar = plt.subplots(figsize=(5, 4))
bars = ax_bar.bar(df["Phase"], df["Percentage"])
ax_bar.set_ylabel("Percentage (%)")
ax_bar.set_title(f"Surface Area (%) (Total = {total_area} µm²)")
for bar in bars:
    height = bar.get_height()
    ax_bar.annotate(f"{height:.2f}%",
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3), textcoords="offset points",
                    ha="center", va="bottom")
fig_bar.tight_layout()
bar_path = os.path.join(output_dir, f"{file_prefix}_bar.png")
fig_bar.savefig(bar_path, dpi=dpi)
print(f"✅ Bar chart saved to: {bar_path}")

# Pie chart
fig_pie, ax_pie = plt.subplots(figsize=(5, 5))
ax_pie.pie(df["Percentage"], labels=df["Phase"], autopct="%1.1f%%", startangle=90)
ax_pie.set_title(f"Surface Area (%) (Total = {total_area} µm²)")
ax_pie.axis("equal")
fig_pie.tight_layout()
pie_path = os.path.join(output_dir, f"{file_prefix}_pie.png")
fig_pie.savefig(pie_path, dpi=dpi)
print(f"✅ Pie chart saved to: {pie_path}")

# -----------------------------
# 5. SHOW PLOTS
# -----------------------------
plt.show()

# -----------------------------
# 6. CONTINUE USING df IF NEEDED
# -----------------------------
# Now you can keep using df for further analysis in Python