#!/usr/bin/env python3
"""Publication-ready diagram of the operational Gothic RGB decision trees."""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch

HERE = Path(__file__).resolve().parent
RULES = HERE / "output_11d_gothic/gothic_rgb_tempbin_decision_tree_rules.csv"
OUT = HERE / "output_11d_gothic"
COND = re.compile(r"(red|green|blue)\s*(<=|>)\s*([-+0-9.eE]+)")

FEATURE_LABEL = {"red": "C13 brightness temperature", "green": "C02 reflectance", "blue": "C05 reflectance"}
CLEAR = "#d9eaf7"
CLOUD = "#f4b183"
NODE = "#f4f4f4"
EDGE = "#333333"


def physical(feature: str, threshold: float) -> tuple[str, str]:
    if feature == "red":
        kelvin = 280.65 - 61.0 * threshold
        return "<" if False else "", f"{kelvin:.1f} K"
    if feature == "green":
        return "", f"{100 * 0.78 * threshold:.1f}%"
    return "", f"{100 * (0.01 + 0.58 * threshold):.1f}%"


def split_text(feature: str, threshold: float) -> str:
    _, units = physical(feature, threshold)
    return f"{FEATURE_LABEL[feature]}\n{feature.capitalize()} = {threshold:.3f}  ({units})"


def box(ax, xy, text, facecolor, width=0.32, height=0.14, fontsize=8.2, weight="normal"):
    x, y = xy
    patch = FancyBboxPatch(
        (x - width / 2, y - height / 2), width, height,
        boxstyle="round,pad=0.012,rounding_size=0.018",
        facecolor=facecolor, edgecolor=EDGE, linewidth=1.05,
        transform=ax.transAxes, clip_on=False,
    )
    ax.add_patch(patch)
    ax.text(x, y, text, ha="center", va="center", fontsize=fontsize,
            fontweight=weight, transform=ax.transAxes, linespacing=1.15)


def edge(ax, start, end, label):
    ax.annotate("", xy=end, xytext=start, xycoords=ax.transAxes,
                arrowprops={"arrowstyle": "-", "color": EDGE, "lw": 1.05})
    x = (start[0] + end[0]) / 2
    y = (start[1] + end[1]) / 2 + 0.012
    ax.text(x, y, label, ha="center", va="bottom", fontsize=7.8,
            fontweight="bold", transform=ax.transAxes,
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.6})


def draw_tree(ax, group: pd.DataFrame):
    first_conditions = [COND.findall(rule)[0] for rule in group.rule]
    root_feature = first_conditions[0][0]
    root_threshold = float(first_conditions[0][2])
    left = group[group.rule.str.startswith(f"{root_feature} <=")].copy()
    right = group[group.rule.str.startswith(f"{root_feature} >")].copy()

    def second_split(branch):
        parsed = [COND.findall(rule) for rule in branch.rule]
        return parsed[0][1][0], float(parsed[0][1][2])

    left_feature, left_threshold = second_split(left)
    right_feature, right_threshold = second_split(right)
    box(ax, (0.5, 0.79), split_text(root_feature, root_threshold), NODE, width=0.43)
    box(ax, (0.26, 0.49), split_text(left_feature, left_threshold), NODE, width=0.40)
    box(ax, (0.74, 0.49), split_text(right_feature, right_threshold), NODE, width=0.40)
    edge(ax, (0.45, 0.72), (0.29, 0.57), "≤")
    edge(ax, (0.55, 0.72), (0.71, 0.57), ">")

    leaves = [(left.iloc[0], 0.12), (left.iloc[1], 0.38), (right.iloc[0], 0.62), (right.iloc[1], 0.88)]
    for row, x in leaves:
        cloudy = int(row.prediction) == 1
        probability = float(row.leaf_cloudy_weight)
        label = f"{'CLOUD' if cloudy else 'CLEAR'}\nP(cloud) = {probability:.2f}"
        box(ax, (x, 0.17), label, CLOUD if cloudy else CLEAR, width=0.21, height=0.13,
            fontsize=8.0, weight="bold")
    edge(ax, (0.22, 0.42), (0.13, 0.24), "≤")
    edge(ax, (0.30, 0.42), (0.37, 0.24), ">")
    edge(ax, (0.70, 0.42), (0.63, 0.24), "≤")
    edge(ax, (0.78, 0.42), (0.87, 0.24), ">")

    first = group.iloc[0]
    ax.set_title(
        f"ERA5-Land $T_{{2m}}$: {first.temp_bin} °C\n"
        f"$n$ = {int(first.n_rows):,}   |   Macro-F1 = {first.macro_f1:.3f}",
        fontsize=11.3, fontweight="bold", pad=7,
    )
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")


def main() -> int:
    rules = pd.read_csv(RULES)
    rules = rules[(rules.status == "trained")].copy()
    order = ["<-10", "[-10, 0)", "[0, 10)", ">=10"]
    fig, axes = plt.subplots(2, 2, figsize=(12.2, 9.0))
    for ax, temp_bin in zip(axes.flat, order):
        draw_tree(ax, rules[rules.temp_bin == temp_bin].sort_values("leaf_id"))

    fig.suptitle(
        "Temperature-conditioned GOES Day Cloud Phase RGB decision trees",
        fontsize=16, fontweight="bold", y=0.985,
    )
    fig.text(
        0.5, 0.025,
        "At each pixel, select the tree using ERA5-Land 2-m air temperature. "
        "Conditions along a path are combined with AND; cloudy terminal paths are combined with OR.\n"
        "Red = inverted C13 brightness temperature normalized over 219.65–280.65 K; "
        "Green = C02 reflectance normalized over 0–0.78; Blue = C05 reflectance normalized over 0.01–0.59.",
        ha="center", va="bottom", fontsize=9.2, linespacing=1.35,
    )
    handles = [
        Line2D([0], [0], marker="s", linestyle="", markersize=13, markerfacecolor=CLEAR,
               markeredgecolor=EDGE, label="Clear terminal leaf"),
        Line2D([0], [0], marker="s", linestyle="", markersize=13, markerfacecolor=CLOUD,
               markeredgecolor=EDGE, label="Cloudy terminal leaf"),
    ]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.085),
               ncol=2, frameon=False, fontsize=10)
    fig.subplots_adjust(left=0.035, right=0.98, top=0.90, bottom=0.145, hspace=0.28, wspace=0.13)
    OUT.mkdir(parents=True, exist_ok=True)
    png = OUT / "gothic_rgb_tempbin_decision_tree_thresholds_paper.png"
    pdf = OUT / "gothic_rgb_tempbin_decision_tree_thresholds_paper.pdf"
    fig.savefig(png, dpi=400, facecolor="white")
    fig.savefig(pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(png)
    print(pdf)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
