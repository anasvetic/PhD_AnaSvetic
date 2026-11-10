#!/usr/bin/env python3
"""Compare Gwyddion grain exports before and after temperature cycling.

Requirements: Python 3.10+, numpy, pandas, matplotlib.
Install: python -m pip install numpy pandas matplotlib

In PyCharm, enter recorded projected ROI/masked areas in the settings below.
Missing file paths open file dialogs; missing ROI area/dimensions open ROI
dimension dialogs. Run once for GPMVs and
again for GUVs, changing SYSTEM_LABEL. Command-line options are also available.

Accepts the original Gwyddion export in SI units and the cleaned nm/nm^2
export created in this conversation. Units come from recognised headers,
never from the numerical magnitude. The original files are not modified.

Each row is a connected masked component, not necessarily a resolved domain.
No p-values, domain tracking, or mechanism assignments are made.
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime
import hashlib
import json
import math
from pathlib import Path
import re
import sys

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter


# ======================= SETTINGS TO EDIT ================================
SYSTEM_LABEL = "GPMV"  # Change to "33 mol% cholesterol GUV" for the GUV pair.

# Example Windows path: r"C:\Users\Ana\Documents\Ana_GPMVs_25C_initial.txt"
# None opens a file-selection dialog.
BEFORE_FILE = None
AFTER_FILE = None

# Enter (width, height) of the analysed ROI in micrometres.
# Example of syntax only: BEFORE_ROI_UM = (5.0, 5.0)
# (None, None) opens dialogs. No ROI dimensions are assumed.
BEFORE_ROI_UM = (None, None)
AFTER_ROI_UM = (None, None)

# Alternatively, enter your recorded ROI AREA directly in square micrometres.
# A supplied area takes precedence over width x height; no dimension dialogs.
BEFORE_ROI_AREA_UM2 = 30.66
AFTER_ROI_AREA_UM2 = 31.98

# Optional recorded TOTAL MASKED AREA, in square micrometres, for coverage.
# Use projected areas (footprints) for BOTH ROI and mask, not 3D surface area.
# These values are used as entered, independently of the grain-size cutoff.
# None uses the sum of exported grain areas passing the size cutoff.
BEFORE_MASKED_AREA_UM2 = 14.53
AFTER_MASKED_AREA_UM2 = 5.95

# Reject tiny mask remnants using equivalent DIAMETER, not radius or an axis.
# d_eq = 2 * sqrt(projected_area / pi). Keep d_eq >= 20 nm in BOTH scans.
# Rejected rows remain flagged in CSVs but do not enter domain counts/density,
# size/shape summaries or distributions. 0.0 disables this diameter cutoff.
MIN_DOMAIN_DIAMETER_NM = 20.0

# Optional additional area cutoff. 0.0 leaves the diameter cutoff in control.
# If both cutoffs are positive, a grain must satisfy both.
MIN_DOMAIN_AREA_NM2 = 0.0

# Optional: row/grain numbers of incomplete domains crossing the ROI boundary.
# They remain in exported area coverage and raw counts, but not size/shape
# summaries. These cannot be identified from the six measurements alone.
BEFORE_INCOMPLETE_GRAIN_IDS = []
AFTER_INCOMPLETE_GRAIN_IDS = []

# Set True ONLY if the export contains every masked component in the analysed
# membrane ROI, including small and boundary-crossing components. False or None
# means the calculated fraction refers only to the exported subset.
BEFORE_EXPORT_CONTAINS_FULL_MASK = None
AFTER_EXPORT_CONTAINS_FULL_MASK = None

OUTPUT_PARENT = None  # None saves next to the BEFORE input, in a new folder.
LOG_SIZE_AXES = True
# =========================================================================

COLORS = {"Before": "#285F83", "After": "#BD6745"}
MEASUREMENTS = (
    "projected_area_nm2",
    "surface_area_nm2",
    "equivalent_disc_radius_nm",
    "projected_boundary_length_nm",
    "major_semiaxis_nm",
    "minor_semiaxis_nm",
)


def normalise(text):
    return re.sub(r"[^a-z0-9]", "", text.lower().replace("²", "2").replace("µ", "u"))


# A unit-free Gwyddion symbol/full name has Gwyddion's base SI units.
# A cleaned header has an explicit nm or nm2 suffix.
ALIASES = {}
for name, symbols, titles, power in [
    (MEASUREMENTS[0], ["A_0"], ["Projected area"], 2),
    (MEASUREMENTS[1], ["A_s"], ["Surface area"], 2),
    (MEASUREMENTS[2], ["r_eq"], ["Equivalent disc radius", "Eq disc radius"], 1),
    (MEASUREMENTS[3], ["L_b0"], ["Projected boundary length", "Boundary length"], 1),
    (MEASUREMENTS[4], ["a_e1"], ["Major semiaxis of equivalent ellipse", "Major semiaxis"], 1),
    (MEASUREMENTS[5], ["a_e2"], ["Minor semiaxis of equivalent ellipse", "Minor semiaxis"], 1),
]:
    for alias in symbols + titles:
        ALIASES[normalise(alias)] = (name, 10.0 ** (9 * power), "SI")
        ALIASES[normalise(alias + (" (m^2)" if power == 2 else " (m)"))] = (
            name, 10.0 ** (9 * power), "SI"
        )
        ALIASES[normalise(alias + (" (nm^2)" if power == 2 else " (nm)"))] = (
            name, 1.0, "nm"
        )
    ALIASES[normalise(name)] = (name, 1.0, "nm")
for alias in ("Grain_ID", "Grain ID", "Grain number"):
    ALIASES[normalise(alias)] = ("grain_id", 1.0, "ID")


def tokens(line):
    text = line.strip().lstrip("#").strip()
    if "\t" in text:
        return [part.strip() for part in text.split("\t") if part.strip()]
    if "," in text:
        return [part.strip() for part in next(csv.reader([text]))]
    return text.split()


def read_export(path):
    """Validate a header, map its actual order, and convert its explicit units."""
    path = Path(path)
    if not path.is_file():
        raise ValueError(f"File not found: {path}")
    text = path.read_text(encoding="utf-8-sig")
    lines = text.splitlines()
    header = None
    header_line = None
    for i, line in enumerate(lines):
        parts = tokens(line)
        if not parts:
            continue
        mapped = [ALIASES.get(normalise(part)) for part in parts]
        if all(mapped) and any(item[0] == "projected_area_nm2" for item in mapped):
            names = [item[0] for item in mapped]
            if len(names) != len(set(names)):
                raise ValueError(f"Duplicate measurement columns in {path.name}.")
            header, header_line = mapped, i
            break
    if header is None:
        raise ValueError(
            f"{path.name}: no recognised header. Keep Gwyddion's header, e.g. "
            "'# A_0 A_s r_eq L_b0 a_e1 a_e2', or the cleaned nm/nm2 header. "
            "The program will not guess column order or units."
        )
    present = {item[0] for item in header}
    required = {"projected_area_nm2", "projected_boundary_length_nm",
                "major_semiaxis_nm", "minor_semiaxis_nm"}
    missing = sorted(required - present)
    if missing:
        raise ValueError(f"{path.name}: missing required columns: {', '.join(missing)}.")
    rows = []
    for line_num, line in enumerate(lines[header_line + 1:], header_line + 2):
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        parts = tokens(line)
        if len(parts) != len(header):
            raise ValueError(
                f"{path.name}, line {line_num}: expected {len(header)} values, "
                f"found {len(parts)}. Nothing was silently skipped."
            )
        try:
            rows.append([float(value) for value in parts])
        except ValueError as exc:
            raise ValueError(f"{path.name}, line {line_num}: non-numeric value.") from exc
    if not rows:
        raise ValueError(f"{path.name}: no grain rows. Check the mask/export.")
    df = pd.DataFrame(rows, columns=[item[0] for item in header])
    if not np.isfinite(df.to_numpy()).all():
        raise ValueError(f"{path.name}: missing, infinite, or non-finite values.")
    for name, factor, _ in header:
        df[name] *= factor
    if "grain_id" not in df:
        df.insert(0, "grain_id", np.arange(1, len(df) + 1))
    ids = df["grain_id"].to_numpy()
    if (ids < 1).any() or (ids != np.floor(ids)).any() or len(set(ids)) != len(ids):
        raise ValueError(f"{path.name}: grain IDs must be unique positive integers.")
    df["grain_id"] = df["grain_id"].astype(int)
    for name in MEASUREMENTS:
        if name in df and (df[name] <= 0).any():
            raise ValueError(f"{path.name}: {name} contains zero/negative measurements.")
    if (df["major_semiaxis_nm"] < df["minor_semiaxis_nm"]).any():
        raise ValueError(f"{path.name}: major semiaxis smaller than minor; check columns.")
    if "equivalent_disc_radius_nm" in df:
        expected_area = np.pi * df["equivalent_disc_radius_nm"] ** 2
        if not np.allclose(df["projected_area_nm2"], expected_area, rtol=5e-5, atol=0):
            raise ValueError(
                f"{path.name}: projected area does not agree with pi * radius^2. "
                "Check units and column mapping."
            )
    else:
        df["equivalent_disc_radius_nm"] = np.sqrt(df["projected_area_nm2"] / np.pi)
    if "surface_area_nm2" not in df:
        df["surface_area_nm2"] = np.nan
    return df[["grain_id", *MEASUREMENTS]], {
        "file": str(path.resolve()),
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "input_columns": [lines[header_line].strip()],
        "exported_components": len(df),
    }


def add_metrics(df, minimum_area, incomplete_ids, minimum_diameter=MIN_DOMAIN_DIAMETER_NM):
    df = df.copy()
    unknown = sorted(set(incomplete_ids) - set(df["grain_id"]))
    if unknown:
        raise ValueError(f"Incomplete-domain IDs not found in this export: {unknown}")
    area = df["projected_area_nm2"]
    a, b = df["major_semiaxis_nm"], df["minor_semiaxis_nm"]
    ratio = b / a
    df["equivalent_diameter_nm"] = 2 * np.sqrt(area / np.pi)
    df["ellipse_major_diameter_nm"] = 2 * a
    df["ellipse_minor_diameter_nm"] = 2 * b
    df["aspect_ratio"] = a / b
    df["ellipticity_1_minus_b_over_a"] = 1 - ratio
    # The maximum protects against numerical roundoff only, not bad axis ordering.
    df["eccentricity"] = np.sqrt(np.maximum(0.0, 1 - ratio ** 2))
    df["circularity_raw"] = 4 * np.pi * area / df["projected_boundary_length_nm"] ** 2
    df["circularity_above_one"] = df["circularity_raw"] > 1 + 1e-6
    df["below_minimum_area"] = area < minimum_area
    diameter = df["equivalent_diameter_nm"]
    df["below_minimum_diameter"] = (diameter < minimum_diameter) & ~np.isclose(
        diameter, minimum_diameter, rtol=1e-12, atol=0.0)
    df["below_size_cutoff"] = df["below_minimum_area"] | df["below_minimum_diameter"]
    df["included_in_domain_count"] = ~df["below_size_cutoff"]
    df["incomplete_domain"] = df["grain_id"].isin(incomplete_ids)
    df["included_in_size_shape"] = ~df["below_size_cutoff"] & ~df["incomplete_domain"]
    flags = []
    for row in df.itertuples():
        flags.append(";".join(
            flag for condition, flag in (
                (row.below_minimum_area, "below_minimum_area"),
                (row.below_minimum_diameter, "below_minimum_diameter"),
                (row.incomplete_domain, "incomplete_domain"),
                (row.circularity_above_one, "circularity_above_one"),
            ) if condition
        ))
    df["flags"] = flags
    return df


def describe(values):
    values = np.asarray(values, dtype=float)
    if len(values) == 0:
        return {key: None for key in ("mean", "sd", "median", "q1", "q3")}
    return {
        "mean": float(np.mean(values)),
        "sd": float(np.std(values, ddof=1)) if len(values) > 1 else None,
        "median": float(np.median(values)),
        "q1": float(np.percentile(values, 25)),
        "q3": float(np.percentile(values, 75)),
    }


def summarise(df, width_um, height_um, full_mask, roi_area_um2=None, masked_area_um2=None):
    if roi_area_um2 is None:
        if not all(v is not None and math.isfinite(v) and v > 0 for v in (width_um, height_um)):
            raise ValueError("Supply a positive ROI area or positive width and height in um.")
        roi_area = width_um * height_um
    else:
        roi_area = roi_area_um2
    if not math.isfinite(roi_area) or roi_area <= 0:
        raise ValueError("ROI area must be finite and positive, in um^2.")
    if masked_area_um2 is not None and (
        not math.isfinite(masked_area_um2) or not 0 <= masked_area_um2 <= roi_area
    ):
        raise ValueError("Recorded masked area must be between zero and ROI area, in um^2.")
    total_area = float(df["projected_area_nm2"].sum() / 1e6)
    if total_area > roi_area * (1 + 5e-5):
        raise ValueError(
            f"Exported projected area ({total_area:.6g} um^2) exceeds ROI area "
            f"({roi_area:.6g} um^2). Check ROI dimensions and units."
        )
    selected = df.loc[df["included_in_size_shape"]]
    retained = df.loc[df["included_in_domain_count"]]
    retained_area = float(retained["projected_area_nm2"].sum() / 1e6)
    coverage_area = masked_area_um2 if masked_area_um2 is not None else retained_area
    coverage_source = "manually_recorded_masked_area" if masked_area_um2 is not None else "size_filtered_export"
    result = {
        "roi_width_um": width_um,
        "roi_height_um": height_um,
        "roi_area_um2": roi_area,
        "roi_area_source": "manual_area" if roi_area_um2 is not None else "width_times_height",
        "recorded_masked_area_um2": masked_area_um2,
        "coverage_source": coverage_source,
        "coverage_projected_area_um2": coverage_area,
        "coverage_area_fraction_pct": 100 * coverage_area / roi_area,
        "retained_component_count": len(retained),
        "retained_component_density_per_um2": len(retained) / roi_area,
        "retained_total_projected_area_um2": retained_area,
        "size_shape_component_density_per_um2": len(selected) / roi_area,
        "export_contains_full_mask": full_mask,
        "exported_component_count": len(df),
        "exported_component_density_per_um2": len(df) / roi_area,
        "exported_total_projected_area_um2": total_area,
        "exported_area_fraction_pct": 100 * total_area / roi_area,
        "components_below_minimum_area": int(df["below_minimum_area"].sum()),
        "components_below_minimum_diameter": int(df["below_minimum_diameter"].sum()),
        "components_below_size_cutoff": int(df["below_size_cutoff"].sum()),
        "excluded_small_area_um2": total_area - retained_area,
        "incomplete_component_count": int(df["incomplete_domain"].sum()),
        "size_shape_component_count": len(selected),
        "circularity_above_one_all": int(df["circularity_above_one"].sum()),
        "circularity_above_one_selected": int(selected["circularity_above_one"].sum()),
        "circularity_requires_review": bool(selected.empty or selected["circularity_above_one"].any()),
    }
    for metric in ("projected_area_nm2", "equivalent_diameter_nm", "aspect_ratio",
                   "ellipticity_1_minus_b_over_a", "eccentricity", "circularity_raw"):
        result.update({f"{metric}_{stat}": value
                       for stat, value in describe(selected[metric]).items()})
    return result


COMPARISON_METRICS = [
    ("coverage_projected_area_um2", "um^2"),
    ("coverage_area_fraction_pct", "percentage points"),
    ("retained_component_count", "components"),
    ("retained_component_density_per_um2", "components/um^2"),
    ("retained_total_projected_area_um2", "um^2"),
    ("exported_total_projected_area_um2", "um^2"),
    ("exported_area_fraction_pct", "percentage points"),
    ("exported_component_count", "components"),
    ("exported_component_density_per_um2", "components/um^2"),
    ("size_shape_component_count", "components"),
    ("projected_area_nm2_mean", "nm^2"),
    ("projected_area_nm2_median", "nm^2"),
    ("equivalent_diameter_nm_mean", "nm"),
    ("equivalent_diameter_nm_median", "nm"),
    ("aspect_ratio_mean", "dimensionless"),
    ("aspect_ratio_median", "dimensionless"),
    ("ellipticity_1_minus_b_over_a_median", "dimensionless"),
    ("eccentricity_median", "dimensionless"),
    ("circularity_raw_mean", "dimensionless"),
    ("circularity_raw_median", "dimensionless"),
]


def make_comparison(before, after):
    rows = []
    for name, unit in COMPARISON_METRICS:
        b, a = before[name], after[name]
        available = b is not None and a is not None
        note = ""
        if name.startswith("coverage_"):
            note = f"Coverage source: before {before['coverage_source']}; after {after['coverage_source']}."
            if before["coverage_source"] != after["coverage_source"]:
                note += " REVIEW: coverage uses different definitions; align them before interpreting."
            elif before["coverage_source"] == "size_filtered_export" and not (
                before["export_contains_full_mask"] is True and after["export_contains_full_mask"] is True
            ):
                note += " Export completeness is unconfirmed; coverage describes the size-filtered exported subset."
        if name.startswith("retained_"):
            note = "Size-filtered exported components; boundary-crossing parts remain in count/area."
        if name.startswith("circularity"):
            if before["circularity_requires_review"] or after["circularity_requires_review"]:
                note = "REVIEW: raw circularity includes values >1 or an empty selected set; do not infer roundness."
        if name.startswith("exported_") and not (
            before["export_contains_full_mask"] is True and after["export_contains_full_mask"] is True
        ):
            note = "Describes the exported components; complete phase coverage/count is not confirmed."
        rows.append({
            "metric": name, "before": b, "after": a,
            "after_minus_before": a - b if available else None,
            "relative_change_pct": 100 * (a / b - 1) if available and b != 0 else None,
            "absolute_change_unit": unit, "note": note,
        })
    return pd.DataFrame(rows)


def export_plots(frames, summaries, label, outdir):
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10,
        "axes.spines.top": False, "axes.spines.right": False,
        "svg.fonttype": "none",
    })
    fig, axes = plt.subplots(2, 2, figsize=(10, 7.2), layout="constrained")
    panels = [
        ("projected_area_nm2", "Projected domain area (nm²)", True),
        ("equivalent_diameter_nm", "Equivalent disc diameter (nm)", True),
        ("aspect_ratio", "Aspect ratio (major / minor)", False),
        ("circularity_raw", "Circularity (raw, 4πA / P²)", False),
    ]
    for ax, (metric, xlabel, is_size) in zip(axes.flat, panels):
        arrays = {
            condition: np.sort(df.loc[df["included_in_size_shape"], metric].to_numpy())
            for condition, df in frames.items()
        }
        pooled = np.concatenate(list(arrays.values()))
        if is_size:
            lo = float(pooled.min()) * 0.85 if len(pooled) else 1.0
            hi = float(pooled.max()) * 1.15 if len(pooled) else 10.0
        elif metric == "aspect_ratio":
            lo, hi = 0.95, max(1.1, float(pooled.max()) * 1.05) if len(pooled) else 2.0
        else:
            lo, hi = 0.0, max(1.05, float(pooled.max()) * 1.05) if len(pooled) else 1.05
        for condition, df in frames.items():
            values = arrays[condition]
            if len(values):
                x = np.r_[lo, values, hi]
                y = np.r_[0, np.arange(1, len(values) + 1) / len(values), 1]
                ax.step(x, y, where="post", color=COLORS[condition], linewidth=1.9,
                        label=f"{condition} (n = {len(values)})")
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Cumulative fraction of selected domains")
        ax.set_ylim(0, 1.03)
        ax.set_xlim(lo, hi)
        ax.grid(alpha=0.18)
        if is_size and LOG_SIZE_AXES:
            ax.set_xscale("log")
            ax.xaxis.set_major_locator(LogLocator(base=10, subs=(1, 2, 5), numticks=8))
            ax.xaxis.set_major_formatter(FuncFormatter(lambda value, position: f"{value:g}"))
            ax.xaxis.set_minor_formatter(NullFormatter())
        if metric == "circularity_raw":
            ax.axvline(1, color="#777777", linestyle=":", linewidth=1)
            nb = summaries["Before"]["circularity_above_one_selected"]
            na = summaries["After"]["circularity_above_one_selected"]
            ax.set_title(f"Values > 1: before {nb}, after {na}", fontsize=10,
                         color="#A33D32" if nb + na else "#444444")
        if ax.get_legend_handles_labels()[0]:
            ax.legend(frameon=False, fontsize=9, loc="lower right")
        else:
            ax.text(.5, .5, "No domains pass the selection", ha="center", transform=ax.transAxes)
    fig.suptitle(f"{label} · Before and after temperature cycling", fontsize=14)
    for extension in ("png", "svg"):
        fig.savefig(outdir / f"domain_distributions.{extension}", dpi=300)
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.8), layout="constrained")
    overview = [
        ("coverage_area_fraction_pct", "Masked area (%)"),
        ("retained_component_density_per_um2", "Domains passing size cutoff per µm²"),
        ("equivalent_diameter_nm_median", "Median selected diameter (nm)"),
    ]
    for ax, (metric, ylabel) in zip(axes, overview):
        vals = [summaries[c][metric] for c in ("Before", "After")]
        ax.bar(["Before", "After"], [v if v is not None else 0 for v in vals],
               color=[COLORS["Before"], COLORS["After"]], width=.6)
        ax.set_ylabel(ylabel)
        ax.set_ylim(bottom=0)
        ax.margins(y=.22)
        ax.grid(axis="y", alpha=.18)
        ax.set_axisbelow(True)
        for i, value in enumerate(vals):
            ax.annotate("n/a" if value is None else f"{value:.3g}",
                        (i, value or 0), xytext=(0, 5), textcoords="offset points",
                        ha="center", fontsize=10)
    sources = [summaries[c]["coverage_source"] for c in ("Before", "After")]
    coverage_label = ("Recorded mask coverage" if all(s == "manually_recorded_masked_area" for s in sources)
                      else "Size-filtered export coverage" if all(s == "size_filtered_export" for s in sources)
                      else "REVIEW: different coverage sources")
    axes[0].set_title(coverage_label, fontsize=9)
    fig.suptitle(f"{label} · Descriptive comparison of one ROI pair", fontsize=13)
    for extension in ("png", "svg"):
        fig.savefig(outdir / f"comparison_overview.{extension}", dpi=300)
    plt.close(fig)


def write_report(summaries, comparison, settings, outdir):
    lines = [
        f"{settings['label']}: before/after temperature-cycling domain analysis",
        "",
        "DESCRIPTIVE RESULTS FOR THIS ROI PAIR",
        "No individual-domain matching, inferential p-values, or biological mechanism is assumed.",
        f"Minimum equivalent domain diameter: {settings['min_diameter_nm']:g} nm (inclusive).",
        f"Additional minimum domain area: {settings['min_area_nm2']:g} nm^2.",
        "",
    ]
    for condition, s in summaries.items():
        lines += [
            condition.upper(),
            f"ROI area: {s['roi_area_um2']:.6g} um^2 ({s['roi_area_source']}).",
            f"Coverage: {s['coverage_area_fraction_pct']:.6g}% from {s['coverage_source']} "
            f"({s['coverage_projected_area_um2']:.6g} um^2).",
            f"Excluded below size cutoff: {s['components_below_size_cutoff']} components "
            f"({s['excluded_small_area_um2']:.6g} um^2).",
            f"Retained domain count: {s['retained_component_count']}; "
            f"density: {s['retained_component_density_per_um2']:.6g} per um^2.",
            f"Exported components: {s['exported_component_count']}.",
            f"Total exported projected area: {s['exported_total_projected_area_um2']:.6g} um^2.",
            f"Exported area fraction: {s['exported_area_fraction_pct']:.6g}%.",
            f"Exported component density: {s['exported_component_density_per_um2']:.6g} per um^2.",
            f"Selected for size/shape: {s['size_shape_component_count']}.",
            f"Raw circularity >1: {s['circularity_above_one_all']} of all components; "
            f"{s['circularity_above_one_selected']} of selected components.",
        ]
        if s["export_contains_full_mask"] is not True:
            lines.append("Full-mask export is unconfirmed: export-derived quantities describe the exported subset; separately recorded coverage is unaffected.")
        if s["circularity_requires_review"]:
            lines.append("Circularity requires review. Values are retained, not clipped or selectively deleted.")
        if s["size_shape_component_count"] == 0:
            lines.append("No size/shape summaries are available for this condition.")
        lines.append("")
    lines += [
        "CHANGES (after minus before; relative change is 100*(after/before - 1))",
        comparison.to_string(index=False, float_format=lambda value: f"{value:.6g}"),
        "",
        "HOW TO INTERPRET",
        "Area fraction: uses recorded masked area when supplied; otherwise size-filtered exported area.",
        "Domain area/diameter: measures size. Smaller medians and a leftward size-distribution shift support smaller domains.",
        "Main count density: counts exported components passing the size cutoff per ROI area.",
        "Raw exported counts/areas are diagnostic records and include rejected tiny components.",
        "Aspect ratio near 1: the equivalent ellipse is less elongated; this alone does not establish a circular outline.",
        "Higher reliable circularity: closer to a circular outline, reflecting elongation and boundary irregularity.",
        "Ellipticity (1-b/a) and eccentricity are alternative expressions of the same axis-ratio information.",
        "",
        "A redistribution interpretation requires similar coverage, smaller domains, and compatible count-density changes.",
        "The program does not declare area conserved: inspect the difference and segmentation/repeatability uncertainty.",
        "Equal ROI dimensions alone do not establish that the same physical region was imaged.",
        "Use the same domain-phase definition, comparable physical resolution, and the same cutoff in both conditions.",
        "Manually recorded coverage is independent of the cutoff and is not altered by grain exclusions.",
        "If the recorded mask includes tiny remnants, they remain in that recorded coverage.",
        "Without a recorded mask area, coverage sums only size-filtered grains, including boundary-crossing parts.",
        "If grains were removed before export, their area/count cannot be recovered.",
        "Boundary-crossing components are incomplete; enter their IDs to exclude them from size/shape only.",
        "If these dominate an image, the remaining size/shape sample may be unrepresentative.",
        "Real holes affect interpretation: verify the boundary convention for perforated grains.",
        "Means/SD describe domain-to-domain spread, not uncertainty across independent patches.",
        "This ROI pair alone does not establish a general effect or prove fragmentation of tracked domains.",
        "",
        "See analysis_settings.json for input hashes, ROI dimensions, and selection settings.",
    ]
    (outdir / "analysis_report.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def id_list(text):
    if not text:
        return []
    try:
        ids = [int(x.strip()) for x in text.split(",")]
    except ValueError as exc:
        raise argparse.ArgumentTypeError("Use comma-separated integer grain IDs, e.g. 1,4,19.") from exc
    if len(ids) != len(set(ids)) or any(x < 1 for x in ids):
        raise argparse.ArgumentTypeError("Grain IDs must be unique positive integers.")
    return ids


def obtain_inputs(args):
    """Ask for dimensions only when neither a recorded area nor dimensions exist."""
    before = args.before or BEFORE_FILE
    after = args.after or AFTER_FILE
    rois = [args.before_roi or BEFORE_ROI_UM, args.after_roi or AFTER_ROI_UM]
    areas = [args.before_roi_area, args.after_roi_area]
    rois = [(None, None) if area is not None else roi for area, roi in zip(areas, rois)]
    if before and after and all(area is not None or all(v is not None for v in roi)
                               for area, roi in zip(areas, rois)):
        return Path(before), Path(after), rois
    try:
        import tkinter as tk
        from tkinter import filedialog, simpledialog
        root = tk.Tk()
        root.withdraw()
    except Exception as exc:
        raise ValueError(
            "Dialogs are unavailable. Enter paths and ROI values in the SETTINGS "
            "section, or use --before/--after and --before-roi-area/--after-roi-area."
        ) from exc
    try:
        paths = [before, after]
        for i, condition in enumerate(("BEFORE", "AFTER")):
            if not paths[i]:
                paths[i] = filedialog.askopenfilename(
                    title=f"{args.label}: select {condition} grain export",
                    filetypes=[("Text grain exports", "*.txt *.tsv *.csv"), ("All files", "*.*")],
                    parent=root,
                )
                if not paths[i]:
                    raise ValueError("Selection cancelled; no analysis was run.")
            if areas[i] is not None:
                continue
            dimensions = list(rois[i])
            for j, dimension in enumerate(("width", "height")):
                if dimensions[j] is None:
                    dimensions[j] = simpledialog.askfloat(
                        f"{condition} ROI", f"Enter analysed ROI {dimension} in micrometres (um):",
                        minvalue=1e-12, parent=root,
                    )
                    if dimensions[j] is None:
                        raise ValueError("ROI entry cancelled; no analysis was run.")
            rois[i] = dimensions
        return Path(paths[0]), Path(paths[1]), rois
    finally:
        root.destroy()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--before", type=Path)
    parser.add_argument("--after", type=Path)
    parser.add_argument("--before-roi", nargs=2, type=float, metavar=("WIDTH_UM", "HEIGHT_UM"))
    parser.add_argument("--after-roi", nargs=2, type=float, metavar=("WIDTH_UM", "HEIGHT_UM"))
    parser.add_argument("--before-roi-area", type=float, default=BEFORE_ROI_AREA_UM2,
                        help="Recorded BEFORE projected ROI area in um^2; overrides dimensions.")
    parser.add_argument("--after-roi-area", type=float, default=AFTER_ROI_AREA_UM2,
                        help="Recorded AFTER projected ROI area in um^2; overrides dimensions.")
    parser.add_argument("--before-masked-area", type=float, default=BEFORE_MASKED_AREA_UM2,
                        help="Recorded BEFORE projected masked area in um^2, for coverage.")
    parser.add_argument("--after-masked-area", type=float, default=AFTER_MASKED_AREA_UM2,
                        help="Recorded AFTER projected masked area in um^2, for coverage.")
    parser.add_argument("--label", default=SYSTEM_LABEL)
    parser.add_argument("--min-diameter-nm", type=float, default=MIN_DOMAIN_DIAMETER_NM,
                        help="Minimum equivalent diameter, inclusive; default 20 nm.")
    parser.add_argument("--min-area-nm2", type=float, default=MIN_DOMAIN_AREA_NM2)
    parser.add_argument("--before-edge-ids", type=id_list, default=BEFORE_INCOMPLETE_GRAIN_IDS)
    parser.add_argument("--after-edge-ids", type=id_list, default=AFTER_INCOMPLETE_GRAIN_IDS)
    parser.add_argument("--full-mask-export", action="store_true",
                        help="Confirm BOTH exports contain all masked components.")
    parser.add_argument("--output-parent", type=Path, default=OUTPUT_PARENT)
    args = parser.parse_args(argv)
    if not math.isfinite(args.min_area_nm2) or args.min_area_nm2 < 0:
        raise ValueError("Minimum domain area must be finite and non-negative.")
    if not math.isfinite(args.min_diameter_nm) or args.min_diameter_nm < 0:
        raise ValueError("Minimum domain diameter must be finite and non-negative.")
    before_path, after_path, rois = obtain_inputs(args)
    if before_path.resolve() == after_path.resolve():
        raise ValueError("The same file was selected twice. Choose distinct before and after exports.")
    frames, summaries, metadata = {}, {}, {}
    roi_areas = [args.before_roi_area, args.after_roi_area]
    masked_areas = [args.before_masked_area, args.after_masked_area]
    for i, (condition, path, edge_ids, complete) in enumerate((
        ("Before", before_path, args.before_edge_ids, BEFORE_EXPORT_CONTAINS_FULL_MASK),
        ("After", after_path, args.after_edge_ids, AFTER_EXPORT_CONTAINS_FULL_MASK),
    )):
        raw, source = read_export(path)
        df = add_metrics(raw, args.min_area_nm2, edge_ids, args.min_diameter_nm)
        width, height = rois[i]
        frames[condition] = df
        summaries[condition] = summarise(df, width, height,
                                        True if args.full_mask_export else complete,
                                        roi_area_um2=roi_areas[i], masked_area_um2=masked_areas[i])
        metadata[condition] = {**source, "roi_width_um": width, "roi_height_um": height,
                               "recorded_roi_area_um2": roi_areas[i],
                               "recorded_masked_area_um2": masked_areas[i],
                               "coverage_source": summaries[condition]["coverage_source"],
                               "incomplete_grain_ids": edge_ids,
                               "export_contains_full_mask": summaries[condition]["export_contains_full_mask"]}
    comparison = make_comparison(summaries["Before"], summaries["After"])
    parent = Path(args.output_parent) if args.output_parent else before_path.parent
    slug = re.sub(r"[^A-Za-z0-9_-]+", "_", args.label).strip("_") or "domains"
    outdir = parent / f"{slug}_temperature_cycling_{datetime.now():%Y%m%d_%H%M%S_%f}"
    outdir.mkdir(parents=True, exist_ok=False)
    settings = {
        "label": args.label, "min_area_nm2": args.min_area_nm2,
        "min_diameter_nm": args.min_diameter_nm,
        "effective_min_area_nm2": max(args.min_area_nm2, math.pi * (args.min_diameter_nm / 2) ** 2),
        "created_local_time": datetime.now().isoformat(timespec="seconds"),
        "inputs": metadata,
        "ellipticity_definition": "1 - minor_semiaxis / major_semiaxis",
        "circularity_definition": "4*pi*projected_area / projected_boundary_length**2",
        "circularity_policy": "Retain raw values; flag >1+1e-6; no clipping or automatic value-based deletion.",
        "coverage_policy": "Recorded projected masked area when supplied, unchanged by filtering; otherwise sum of size-filtered exported areas, including boundary-crossing parts.",
        "count_policy": "Main count/density excludes grains below size cutoff; boundary-crossing components remain. Raw counts are audit-only.",
        "log_size_axes": LOG_SIZE_AXES,
        "versions": {"python": sys.version.split()[0], "numpy": np.__version__,
                     "pandas": pd.__version__, "matplotlib": matplotlib.__version__},
    }
    (outdir / "analysis_settings.json").write_text(
        json.dumps(settings, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    for condition, df in frames.items():
        df.to_csv(outdir / f"{condition.lower()}_domain_measurements.csv", index=False,
                  float_format="%.12g", encoding="utf-8-sig")
    pd.DataFrame([{"condition": k, **v} for k, v in summaries.items()]).to_csv(
        outdir / "summary_by_condition.csv", index=False, float_format="%.12g",
        encoding="utf-8-sig")
    comparison.to_csv(outdir / "before_after_changes.csv", index=False,
                      float_format="%.12g", encoding="utf-8-sig")
    write_report(summaries, comparison, settings, outdir)
    export_plots(frames, summaries, args.label, outdir)
    print(f"\nAnalysis saved to:\n{outdir.resolve()}\n")
    for condition, s in summaries.items():
        print(f"{condition}: {s['exported_component_count']} exported components; "
              f"{s['components_below_size_cutoff']} below size cutoff; "
              f"{s['retained_component_count']} counted; "
              f"{s['size_shape_component_count']} selected for size/shape; "
              f"coverage {s['coverage_area_fraction_pct']:.4g}% ({s['coverage_source']}).")
        if s["circularity_requires_review"]:
            print(f"  REVIEW circularity: {s['circularity_above_one_selected']} selected values >1.")
    print("\nRead analysis_report.txt before interpreting the plots.")
    return outdir


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError) as exc:
        print(f"\nERROR: {exc}", file=sys.stderr)
        sys.exit(1)
