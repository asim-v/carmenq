"""Generate the publication figure for the order-sensitive support bounds."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from carmenq import (
    INTERLEAVED_BALANCED_COUNTEREXAMPLE,
    INTERLEAVED_ORDER_GAP_WEIGHT_THRESHOLD,
    interleaved_candidate_lower_bound,
    interleaved_support_upper_bound,
    rank_two_static_qubit_support,
)


ROOT = Path(__file__).resolve().parents[1]
PDF_OUTPUT = ROOT / "figures" / "order_support_gap.pdf"
PNG_OUTPUT = ROOT / "figures" / "order_support_gap.png"


def main() -> None:
    weights = np.linspace(0.0, 1.0, 241)
    grouped = np.array([rank_two_static_qubit_support(value) for value in weights])
    certificate = np.array(
        [interleaved_support_upper_bound(value) for value in weights]
    )
    candidate = np.array(
        [interleaved_candidate_lower_bound(value).support_value for value in weights]
    )
    no_record = 1.0 - weights / 2.0
    threshold = INTERLEAVED_ORDER_GAP_WEIGHT_THRESHOLD
    certified = weights >= threshold

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 9.5,
            "axes.labelsize": 10.5,
            "legend.fontsize": 8.5,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
        }
    )
    figure, axis = plt.subplots(figsize=(7.0, 3.65), constrained_layout=True)
    axis.fill_between(
        weights[certified],
        certificate[certified],
        grouped[certified],
        color="#F2A93B",
        alpha=0.20,
        linewidth=0,
        label="rigorously excluded by temporal order",
    )
    axis.plot(
        weights,
        grouped,
        color="#174A7E",
        linewidth=2.2,
        label=r"grouped exact / static ceiling $B_{4,2}$",
    )
    axis.plot(
        weights[certified],
        certificate[certified],
        color="#C96A1B",
        linewidth=2.2,
        label=r"interleaved rigorous certificate $U_{\rm I}$",
    )
    axis.plot(
        weights[~certified],
        certificate[~certified],
        color="#C96A1B",
        linewidth=1.2,
        linestyle=":",
        alpha=0.65,
    )
    axis.plot(
        weights,
        candidate,
        color="#188977",
        linewidth=1.8,
        linestyle="--",
        label="analytic achievable family",
    )
    axis.plot(
        weights,
        no_record,
        color="#6B7785",
        linewidth=1.1,
        linestyle=(0, (2, 2)),
        label="no-record strategy",
    )
    point = INTERLEAVED_BALANCED_COUNTEREXAMPLE
    axis.scatter(
        [point.audit_weight],
        [point.support_value],
        s=42,
        marker="D",
        color="#7A3E9D",
        edgecolor="white",
        linewidth=0.7,
        zorder=5,
        label="verified complete non-QND strategy",
    )
    axis.axvline(threshold, color="#C96A1B", linewidth=0.9, alpha=0.75)
    axis.annotate(
        r"$\lambda=3/7$",
        xy=(threshold, 0.733),
        xytext=(threshold + 0.025, 0.737),
        color="#9A4E11",
        fontsize=8.5,
    )
    axis.set_xlim(0.0, 1.0)
    axis.set_ylim(0.725, 1.008)
    axis.set_xlabel(r"AUDIT weight $\lambda$")
    axis.set_ylabel(r"support $\lambda P_{\rm A}+(1-\lambda)F_{\rm R}$")
    axis.set_xticks(np.linspace(0.0, 1.0, 6))
    axis.grid(True, color="#D7DEE5", linewidth=0.55, alpha=0.75)
    axis.spines[["top", "right"]].set_visible(False)
    axis.legend(loc="lower right", frameon=True, framealpha=0.96, ncol=1)

    PDF_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(PDF_OUTPUT, bbox_inches="tight")
    figure.savefig(PNG_OUTPUT, dpi=220, bbox_inches="tight")
    plt.close(figure)
    print(PDF_OUTPUT.relative_to(ROOT))
    print(PNG_OUTPUT.relative_to(ROOT))


if __name__ == "__main__":
    main()
