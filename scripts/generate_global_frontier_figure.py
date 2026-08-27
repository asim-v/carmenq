"""Generate the publication graphic for the certified lambda=0.6 enclosure."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "data" / "global_frontier_l060_exact_assembly.json"
PDF_OUTPUT = ROOT / "figures" / "global_frontier_certificate.pdf"
PNG_OUTPUT = ROOT / "figures" / "global_frontier_certificate.png"

BLUE = "#174A7E"
GREEN = "#24735B"
ORANGE = "#C96A1B"
PALE_BLUE = "#EAF1F8"
PALE_GREEN = "#E8F3EE"
PALE_ORANGE = "#FFF0E3"
GREY = "#5F6B76"


def as_fraction(item: object) -> Fraction:
    if isinstance(item, list) and len(item) == 2:
        return Fraction(int(item[0]), int(item[1]))
    return Fraction(str(item))


def load_manifest(path: Path) -> dict[str, object]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("schema") != "carmenq.global-frontier-l060-exact-assembly.v1":
        raise RuntimeError("wrong global-frontier manifest schema")
    if not payload.get("complete"):
        raise RuntimeError("global-frontier manifest is incomplete")
    return payload


def labelled_box(
    axis: plt.Axes,
    x: float,
    y: float,
    width: float,
    height: float,
    text: str,
    *,
    facecolor: str,
    edgecolor: str,
) -> None:
    patch = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle="round,pad=0.012,rounding_size=0.018",
        facecolor=facecolor,
        edgecolor=edgecolor,
        linewidth=1.0,
    )
    axis.add_patch(patch)
    axis.text(
        x + width / 2,
        y + height / 2,
        text,
        ha="center",
        va="center",
        fontsize=7.7,
        linespacing=1.18,
    )


def draw_sector_map(axis: plt.Axes) -> None:
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.axis("off")
    axis.set_title("(a) Exhaustive terminal-readout partition", loc="left", pad=6)

    labelled_box(
        axis,
        0.25,
        0.84,
        0.50,
        0.105,
        "extreme qubit POVM\n(at most four active effects)",
        facecolor=PALE_BLUE,
        edgecolor=BLUE,
    )
    labelled_box(
        axis,
        0.02,
        0.57,
        0.29,
        0.13,
        "$w_{\\max}\\leq0.88325$\nexact dual replay\n$\\leq0.765893818$",
        facecolor=PALE_GREEN,
        edgecolor=GREEN,
    )
    labelled_box(
        axis,
        0.355,
        0.57,
        0.29,
        0.13,
        "$w_{\\max}>0.88325$\n2 or 3 active\n$\\leq0.76652$",
        facecolor=PALE_BLUE,
        edgecolor=BLUE,
    )
    labelled_box(
        axis,
        0.69,
        0.57,
        0.29,
        0.13,
        "$w_{\\max}>0.88325$\n4 active effects",
        facecolor=PALE_ORANGE,
        edgecolor=ORANGE,
    )
    labelled_box(
        axis,
        0.58,
        0.25,
        0.19,
        0.15,
        "$w_{\\min}\\geq0.0003$\nexact cover\n$\\leq0.76670$",
        facecolor=PALE_GREEN,
        edgecolor=GREEN,
    )
    labelled_box(
        axis,
        0.80,
        0.25,
        0.19,
        0.15,
        "$w_{\\min}<0.0003$\ndelete effect\n$\\leq0.76670$",
        facecolor=PALE_ORANGE,
        edgecolor=ORANGE,
    )

    arrow = dict(arrowstyle="-|>", color=GREY, linewidth=0.9, mutation_scale=8)
    axis.annotate("", xy=(0.165, 0.71), xytext=(0.43, 0.84), arrowprops=arrow)
    axis.annotate("", xy=(0.50, 0.71), xytext=(0.50, 0.84), arrowprops=arrow)
    axis.annotate("", xy=(0.835, 0.71), xytext=(0.57, 0.84), arrowprops=arrow)
    axis.annotate("", xy=(0.675, 0.41), xytext=(0.79, 0.57), arrowprops=arrow)
    axis.annotate("", xy=(0.895, 0.41), xytext=(0.87, 0.57), arrowprops=arrow)
    axis.text(
        0.50,
        0.08,
        r"maximum of all sectors: $\overline{\beta}_{2\mathrm{b}}(0.6)\leq0.76670$",
        ha="center",
        va="center",
        color=BLUE,
        fontsize=8.8,
        fontweight="bold",
    )


def draw_interval(axis: plt.Axes, payload: dict[str, object]) -> None:
    lower = float(as_fraction(payload["explicit_physical_lower_fraction"]))
    upper = float(as_fraction(payload["assembled_upper_fraction"]))
    sectors = payload["sector_bounds"]
    assert isinstance(sectors, dict)
    values = {
        "low-weight": float(
            as_fraction(sectors["maximum_effect_at_most_0p88325"]["fraction"])
        ),
        "projective": float(as_fraction(sectors["binary_projective"]["fraction"])),
        "ternary": float(as_fraction(sectors["ternary"]["fraction"])),
        "four-active / deletion": upper,
    }

    axis.set_title("(b) Certified support interval", loc="left", pad=6)
    xmin, xmax = 0.76584, 0.76674
    axis.axvspan(lower, upper, color=PALE_ORANGE, alpha=0.9, zorder=0)
    rows = list(values.items())
    for row, (label, value) in enumerate(rows[::-1], start=1):
        axis.hlines(row, xmin, value, color="#AEB8C2", linewidth=2.2)
        axis.scatter([value], [row], s=30, color=BLUE, zorder=3)
        axis.text(xmin + 0.000012, row + 0.15, label, fontsize=7.8, color=GREY)
        axis.text(
            value - 0.000006,
            row - 0.23,
            f"{value:.9f}",
            fontsize=7.2,
            ha="right",
            color=BLUE,
        )
    axis.axvline(lower, color=GREEN, linewidth=1.4)
    axis.axvline(upper, color=ORANGE, linewidth=1.4)
    axis.text(
        lower,
        0.26,
        "rational physical\nwitness",
        ha="left",
        va="bottom",
        fontsize=7.7,
        color=GREEN,
    )
    axis.text(
        upper,
        4.55,
        "global upper",
        ha="right",
        va="top",
        fontsize=7.7,
        color=ORANGE,
    )
    axis.text(
        (lower + upper) / 2,
        0.04,
        f"width = {upper - lower:.10f}",
        ha="center",
        va="bottom",
        fontsize=8.3,
        fontweight="bold",
        color="#7D4A13",
    )
    axis.set_xlim(xmin, xmax)
    axis.set_ylim(0, 4.75)
    axis.set_yticks([])
    axis.set_xlabel(r"support at $\lambda=0.6$")
    axis.ticklabel_format(axis="x", style="plain", useOffset=False)
    axis.tick_params(axis="x", labelsize=7.5)
    axis.grid(axis="x", color="#D7DEE5", linewidth=0.55)
    axis.spines[["left", "right", "top"]].set_visible(False)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    args = parser.parse_args()
    payload = load_manifest(args.manifest)

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 9.0,
            "axes.labelsize": 9.5,
            "axes.titlesize": 9.5,
        }
    )
    figure, (left, right) = plt.subplots(
        1,
        2,
        figsize=(7.35, 3.45),
        gridspec_kw={"width_ratios": (1.55, 1.0)},
        constrained_layout=True,
    )
    draw_sector_map(left)
    draw_interval(right, payload)

    PDF_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(PDF_OUTPUT, bbox_inches="tight")
    figure.savefig(PNG_OUTPUT, dpi=240, bbox_inches="tight")
    plt.close(figure)
    print(PDF_OUTPUT.relative_to(ROOT))
    print(PNG_OUTPUT.relative_to(ROOT))


if __name__ == "__main__":
    main()
