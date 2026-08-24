#!/usr/bin/env python3
"""Append or replace the western-U.S. site map at the end of notebook 16."""

from pathlib import Path

import nbformat


NOTEBOOK = Path(__file__).with_name("16_plot_gothic_terrain_FIGURE.ipynb")
SECTION_TAG = "## Western U.S. evaluation sites"

markdown = f"""{SECTION_TAG}

Terrain context and locations of the four model-training/evaluation sites."""

code = r"""import geopandas as gpd
import rasterio
from rasterio.windows import from_bounds
from shapely.geometry import box

MAP_WEST, MAP_EAST = -125.0, -102.0
MAP_SOUTH, MAP_NORTH = 31.0, 49.0
NATURAL_EARTH_DIR = Path("analysis/data/natural_earth")
if not NATURAL_EARTH_DIR.exists():
    NATURAL_EARTH_DIR = Path("data/natural_earth")

RELIEF_TIF = NATURAL_EARTH_DIR / "HYP_50M_SR_W.tif"
STATES_SHP = NATURAL_EARTH_DIR / "ne_50m_admin_1_states_provinces.shp"
ANALYSIS_DIR = (
    Path.cwd()
    if (Path.cwd() / "data" / "natural_earth").exists()
    else Path.cwd() / "analysis"
)
MAP_OUTPUT = ANALYSIS_DIR / "output_16_western_us_site_map.png"

sites = {
    "Gothic, CO": (-106.9886, 38.9580),
    "Table Mountain, CO": (-105.2368, 40.12498),
    "Senator Beck, CO": (-107.7200, 37.9000),
    "CUES/Mammoth, CA": (-119.0290, 37.6430),
}
site_colors = {
    name: ("#111111" if name == "Gothic, CO" else "#D73027")
    for name in sites
}
label_offsets = {
    "Gothic, CO": (0.7, -0.15),
    "Table Mountain, CO": (-3.8, 0.65),
    "Senator Beck, CO": (-3.1, -1.0),
    "CUES/Mammoth, CA": (0.7, 0.55),
}

with rasterio.open(RELIEF_TIF) as src:
    window = from_bounds(
        MAP_WEST, MAP_SOUTH, MAP_EAST, MAP_NORTH, transform=src.transform
    ).round_offsets().round_lengths()
    relief = np.moveaxis(src.read(window=window), 0, -1)
    left, bottom, right, top = rasterio.windows.bounds(window, src.transform)
    relief_extent = (left, right, bottom, top)

states = gpd.read_file(STATES_SHP)
western_box = box(MAP_WEST, MAP_SOUTH, MAP_EAST, MAP_NORTH)
western_states = states.loc[
    (states.iso_a2 == "US") & states.geometry.intersects(western_box)
]

with plt.rc_context({
    "font.size": 14,
    "font.weight": "bold",
    "axes.labelweight": "bold",
}):
    fig, ax = plt.subplots(figsize=(10.5, 7.5))
    ax.imshow(
        relief,
        extent=relief_extent,
        origin="upper",
        interpolation="bilinear",
        alpha=0.92,
        zorder=0,
    )
    western_states.boundary.plot(
        ax=ax, color="white", linewidth=1.05, alpha=0.9, zorder=2
    )

    for name, (lon, lat) in sites.items():
        ax.scatter(
            lon, lat, s=145, marker="o", color=site_colors[name],
            edgecolor="white", linewidth=1.8, zorder=5,
        )
        dx, dy = label_offsets[name]
        ax.annotate(
            name,
            xy=(lon, lat),
            xytext=(lon + dx, lat + dy),
            textcoords="data",
            ha="left",
            va="center",
            fontsize=13,
            fontweight="bold",
            color="#111111",
            bbox=dict(
                boxstyle="round,pad=0.22",
                facecolor="white",
                edgecolor="none",
                alpha=0.84,
            ),
            arrowprops=dict(
                arrowstyle="-", color="#222222", linewidth=1.2,
                shrinkA=2, shrinkB=5,
            ),
            zorder=6,
        )

    ax.set_xlim(MAP_WEST, MAP_EAST)
    ax.set_ylim(MAP_SOUTH, MAP_NORTH)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_xticks(np.arange(-125, -101, 5))
    ax.set_yticks(np.arange(32, 50, 4))
    ax.grid(color="white", alpha=0.22, linewidth=0.8)
    for label in [*ax.get_xticklabels(), *ax.get_yticklabels()]:
        label.set_fontweight("bold")
    fig.tight_layout()
    fig.savefig(MAP_OUTPUT, dpi=300, bbox_inches="tight")
    plt.show()

print(f"Saved map to {MAP_OUTPUT.resolve()}")"""

nb = nbformat.read(NOTEBOOK, as_version=4)
for index, cell in enumerate(nb.cells):
    if cell.cell_type == "markdown" and SECTION_TAG in cell.source:
        nb.cells[index:index + 2] = [
            nbformat.v4.new_markdown_cell(markdown),
            nbformat.v4.new_code_cell(code),
        ]
        break
else:
    nb.cells.extend([
        nbformat.v4.new_markdown_cell(markdown),
        nbformat.v4.new_code_cell(code),
    ])

nbformat.write(nb, NOTEBOOK)
print(NOTEBOOK)
