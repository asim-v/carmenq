"""Deterministic experiment suite and publication-figure generation."""

from __future__ import annotations

import csv
import json
import platform
from dataclasses import asdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy

from . import __version__
from .protocol import (
    NoiseModel,
    ProtocolConfig,
    conditional_record_information,
    environment_conditional_information,
    run_protocol,
)


PLOT_COLORS = {
    "visibility": "#0072B2",
    "reset_fidelity": "#D55E00",
    "predicate_fidelity": "#009E73",
    "transcript": "#CC79A7",
}


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError("Cannot write an empty table.")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _save_figure(figure: plt.Figure, base_path: Path) -> None:
    base_path.parent.mkdir(parents=True, exist_ok=True)
    # Suppress wall-clock metadata so two equivalent runs produce identical
    # PDF bytes.  This lets CI detect genuine figure drift.
    figure.savefig(
        base_path.with_suffix(".pdf"),
        bbox_inches="tight",
        metadata={"CreationDate": None, "ModDate": None},
    )
    figure.savefig(base_path.with_suffix(".png"), dpi=240, bbox_inches="tight")
    plt.close(figure)


def _plot_style() -> None:
    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.labelsize": 9,
            "axes.titlesize": 10,
            "legend.fontsize": 8,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.dpi": 120,
            "savefig.transparent": False,
        }
    )


def _sample_target_probability(
    probabilities: np.ndarray, target: int, shots: int, rng: np.random.Generator
) -> tuple[int, float]:
    counts = rng.multinomial(shots, np.clip(probabilities.real, 0.0, 1.0))
    return int(counts[target]), float(counts[target] / shots)


def benchmark_suite(
    data_dir: Path, figures_dir: Path, seed: int, shots: int
) -> list[dict[str, object]]:
    """Generate exact and finite-shot results for positive and negative controls."""
    cases: list[tuple[str, str, ProtocolConfig]] = [
        ("ideal", "Ideal coherent", ProtocolConfig()),
        ("classical_mixture", "Classical mixture", ProtocolConfig(coherent_input=False)),
        ("memory_retained", "Memory retained", ProtocolConfig(uncompute="leave_memory")),
        ("direct_phase", "Direct phase", ProtocolConfig(direct_phase=True)),
        ("action_bypassed", "Action bypassed", ProtocolConfig(enable_actions=False)),
        (
            "partial_environment_leakage",
            "Leakage, overlap 0.5",
            ProtocolConfig(noise=NoiseModel(environment_overlap=0.5)),
        ),
        (
            "imperfect_inversion",
            "Inverse failure 0.03",
            ProtocolConfig(noise=NoiseModel(inversion_error=0.03)),
        ),
    ]
    rng = np.random.default_rng(seed)
    rows: list[dict[str, object]] = []
    saved_states: dict[str, np.ndarray] = {}
    for case_id, label, config in cases:
        result = run_protocol(config)
        target_counts, sampled_fidelity = _sample_target_probability(
            np.diag(result.branch_state_after_readout).real, 3, shots, rng
        )
        record_information = conditional_record_information(config)
        environment_information = environment_conditional_information(
            config.noise.environment_overlap
        )
        row: dict[str, object] = {
            "case": case_id,
            "label": label,
            "visibility": result.metrics["visibility"],
            "reset_fidelity": result.metrics["reset_fidelity"],
            "predicate_fidelity": result.metrics["predicate_fidelity"],
            "target_contrast": result.metrics["target_contrast"],
            "record_chi_bits_given_parity": record_information,
            "environment_chi_bits_given_parity": environment_information,
            # This is a bookkeeping sum of separately evaluated marginals.  It
            # is not, without further assumptions, a bound on a correlated
            # joint record.
            "sum_of_marginal_residual_bits": record_information + environment_information,
            "target_counts": target_counts,
            "shots": shots,
            "sampled_predicate_fidelity": sampled_fidelity,
            "seed": seed,
        }
        rows.append(row)
        saved_states[f"{case_id}_branch_before"] = result.branch_state_before_readout
        saved_states[f"{case_id}_branch_after"] = result.branch_state_after_readout

    _write_csv(data_dir / "benchmark_summary.csv", rows)
    np.savez_compressed(data_dir / "reference_branch_states.npz", **saved_states)

    display_rows = rows[:6]
    x = np.arange(len(display_rows))
    width = 0.24
    figure, axis = plt.subplots(figsize=(7.2, 3.35))
    for offset, key, label in (
        (-width, "visibility", "Coherence visibility"),
        (0.0, "reset_fidelity", "Reset fidelity"),
        (width, "predicate_fidelity", "Predicate fidelity"),
    ):
        axis.bar(
            x + offset,
            [float(row[key]) for row in display_rows],
            width,
            label=label,
            color=PLOT_COLORS[key],
        )
    axis.set_ylim(0.0, 1.06)
    axis.set_ylabel("Exact metric")
    axis.set_xticks(x, [str(row["label"]) for row in display_rows], rotation=24, ha="right")
    axis.legend(
        frameon=False,
        ncol=3,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.035),
        borderaxespad=0.0,
    )
    axis.set_title("Four-history parity protocol and controls", pad=30)
    axis.grid(axis="y", alpha=0.2)
    _save_figure(figure, figures_dir / "protocol_controls")
    return rows


def noise_suite(data_dir: Path, figures_dir: Path) -> list[dict[str, object]]:
    """Sweep three explicit local channels at a per-operation probability."""
    levels = np.linspace(0.0, 0.12, 13)
    rows: list[dict[str, object]] = []
    for channel in ("dephasing", "depolarizing", "amplitude_damping"):
        for level in levels:
            noise = NoiseModel(**{channel: float(level)})
            result = run_protocol(ProtocolConfig(noise=noise))
            rows.append(
                {
                    "channel": channel,
                    "probability_per_operation": float(level),
                    "visibility": result.metrics["visibility"],
                    "reset_fidelity": result.metrics["reset_fidelity"],
                    "predicate_fidelity": result.metrics["predicate_fidelity"],
                    "target_contrast": result.metrics["target_contrast"],
                    "trace": result.metrics["trace"],
                }
            )
    _write_csv(data_dir / "noise_sweep.csv", rows)

    figure, axes = plt.subplots(1, 3, figsize=(8.3, 2.8), sharex=True, sharey=True)
    titles = {
        "dephasing": "Qudit dephasing",
        "depolarizing": "Qudit depolarization",
        "amplitude_damping": "Relaxation to zero",
    }
    for axis, channel in zip(axes, titles):
        selected = [row for row in rows if row["channel"] == channel]
        x = [float(row["probability_per_operation"]) for row in selected]
        for key, label in (
            ("visibility", "Visibility"),
            ("reset_fidelity", "Reset"),
            ("predicate_fidelity", "Predicate"),
        ):
            axis.plot(
                x,
                [float(row[key]) for row in selected],
                marker="o",
                markersize=2.6,
                linewidth=1.3,
                color=PLOT_COLORS[key],
                label=label,
            )
        axis.set_title(titles[channel])
        axis.set_xlabel("Error probability per operation")
        axis.grid(alpha=0.2)
    axes[0].set_ylabel("Exact metric")
    axes[0].set_ylim(-0.02, 1.03)
    axes[-1].legend(frameon=False, loc="best")
    figure.suptitle("Sensitivity to local Markovian channels", y=1.02)
    _save_figure(figure, figures_dir / "noise_sensitivity")
    return rows


def environment_suite(data_dir: Path, figures_dir: Path) -> list[dict[str, object]]:
    """Sweep the common overlap between branch-conditioned environment states."""
    rows: list[dict[str, object]] = []
    for overlap in np.linspace(0.0, 1.0, 21):
        result = run_protocol(
            ProtocolConfig(noise=NoiseModel(environment_overlap=float(overlap)))
        )
        rows.append(
            {
                "environment_overlap": float(overlap),
                "visibility": result.metrics["visibility"],
                "reset_fidelity": result.metrics["reset_fidelity"],
                "predicate_fidelity": result.metrics["predicate_fidelity"],
                "environment_chi_bits_given_parity": environment_conditional_information(
                    float(overlap)
                ),
            }
        )
    _write_csv(data_dir / "environment_sweep.csv", rows)

    figure, axis = plt.subplots(figsize=(4.8, 3.2))
    x = [float(row["environment_overlap"]) for row in rows]
    axis.plot(
        x,
        [float(row["visibility"]) for row in rows],
        color=PLOT_COLORS["visibility"],
        linewidth=2,
        label="Recovered visibility",
    )
    axis.plot(
        x,
        [float(row["environment_chi_bits_given_parity"]) for row in rows],
        color=PLOT_COLORS["transcript"],
        linewidth=2,
        label=r"Environment $\chi(H:E\mid P)$",
    )
    axis.set_xlabel(r"Environment-state overlap $\langle e_j|e_i\rangle$")
    axis.set_ylabel("Visibility or information (bits)")
    axis.set_xlim(0.0, 1.0)
    axis.set_ylim(0.0, 1.03)
    axis.grid(alpha=0.2)
    axis.legend(frameon=False)
    axis.set_title("Leakage–recoherence tradeoff for four histories")
    _save_figure(figure, figures_dir / "environment_tradeoff")
    return rows


def inversion_suite(data_dir: Path, figures_dir: Path) -> list[dict[str, object]]:
    """Sweep stochastic inverse-operation failures."""
    rows: list[dict[str, object]] = []
    for level in np.linspace(0.0, 0.12, 13):
        result = run_protocol(
            ProtocolConfig(noise=NoiseModel(inversion_error=float(level)))
        )
        rows.append(
            {
                "inverse_failure_probability": float(level),
                "visibility": result.metrics["visibility"],
                "reset_fidelity": result.metrics["reset_fidelity"],
                "predicate_fidelity": result.metrics["predicate_fidelity"],
            }
        )
    _write_csv(data_dir / "inversion_sweep.csv", rows)

    figure, axis = plt.subplots(figsize=(4.7, 3.15))
    x = [float(row["inverse_failure_probability"]) for row in rows]
    for key, label in (
        ("visibility", "Visibility"),
        ("reset_fidelity", "Reset fidelity"),
        ("predicate_fidelity", "Predicate fidelity"),
    ):
        axis.plot(
            x,
            [float(row[key]) for row in rows],
            marker="o",
            markersize=3,
            color=PLOT_COLORS[key],
            label=label,
        )
    axis.set_xlabel("Failure probability per inverse operation")
    axis.set_ylabel("Exact metric")
    axis.set_ylim(0.0, 1.03)
    axis.grid(alpha=0.2)
    axis.legend(frameon=False)
    axis.set_title("Sensitivity to imperfect inversion")
    _save_figure(figure, figures_dir / "inversion_sensitivity")
    return rows


def depth_suite(
    data_dir: Path, figures_dir: Path, seed: int
) -> list[dict[str, object]]:
    """Increase challenge-responsive history length at fixed local noise."""
    rng = np.random.default_rng(seed + 1)
    rows: list[dict[str, object]] = []
    for rounds in (1, 2, 3, 4, 6, 8):
        challenges = tuple(int(value) for value in rng.integers(0, 2, size=rounds))
        ideal = run_protocol(ProtocolConfig(challenges=challenges))
        noisy = run_protocol(
            ProtocolConfig(
                challenges=challenges,
                noise=NoiseModel(dephasing=0.003, depolarizing=0.001),
            )
        )
        rows.append(
            {
                "rounds": rounds,
                "logical_gate_count": 5 + 6 * rounds,
                "challenges": "".join(str(value) for value in challenges),
                "per_gate_dephasing": 0.003,
                "per_gate_depolarizing": 0.001,
                "ideal_visibility": ideal.metrics["visibility"],
                "ideal_reset_fidelity": ideal.metrics["reset_fidelity"],
                "ideal_predicate_fidelity": ideal.metrics["predicate_fidelity"],
                "noisy_visibility": noisy.metrics["visibility"],
                "noisy_reset_fidelity": noisy.metrics["reset_fidelity"],
                "noisy_predicate_fidelity": noisy.metrics["predicate_fidelity"],
                "seed": seed + 1,
            }
        )
    _write_csv(data_dir / "depth_sweep.csv", rows)

    figure, axis = plt.subplots(figsize=(4.8, 3.2))
    x = [int(row["logical_gate_count"]) for row in rows]
    for key, label, color_key in (
        ("noisy_visibility", "Visibility", "visibility"),
        ("noisy_reset_fidelity", "Reset fidelity", "reset_fidelity"),
        ("noisy_predicate_fidelity", "Predicate fidelity", "predicate_fidelity"),
    ):
        axis.plot(
            x,
            [float(row[key]) for row in rows],
            marker="o",
            color=PLOT_COLORS[color_key],
            label=label,
        )
    axis.set_xlabel("Logical operations")
    axis.set_ylabel("Exact metric")
    axis.set_ylim(0.7, 1.01)
    axis.grid(alpha=0.2)
    axis.legend(frameon=False)
    axis.set_title("Challenge-responsive depth at fixed local noise")
    _save_figure(figure, figures_dir / "depth_scaling")
    return rows


def generate_all(
    root: Path, seed: int = 20260812, shots: int = 8192
) -> dict[str, object]:
    """Regenerate every committed table, state file, figure, and metadata file."""
    root = Path(root).resolve()
    data_dir = root / "data"
    figures_dir = root / "figures"
    data_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)
    _plot_style()

    benchmark_rows = benchmark_suite(data_dir, figures_dir, seed, shots)
    noise_rows = noise_suite(data_dir, figures_dir)
    environment_rows = environment_suite(data_dir, figures_dir)
    inversion_rows = inversion_suite(data_dir, figures_dir)
    depth_rows = depth_suite(data_dir, figures_dir, seed)

    metadata: dict[str, object] = {
        "schema_version": 1,
        "simulator_version": __version__,
        "seed": seed,
        "shots": shots,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "matplotlib": matplotlib.__version__,
        "hilbert_space": "B(4) x W(4) x M(4) x G(2) x A(2)",
        "hilbert_dimension": 256,
        "predicate": "two-bit history-label parity",
        "default_challenges": list(ProtocolConfig().challenges),
        "command": f"python scripts/regenerate.py --seed {seed} --shots {shots}",
        "row_counts": {
            "benchmark_summary.csv": len(benchmark_rows),
            "noise_sweep.csv": len(noise_rows),
            "environment_sweep.csv": len(environment_rows),
            "inversion_sweep.csv": len(inversion_rows),
            "depth_sweep.csv": len(depth_rows),
        },
        "default_noise_model": asdict(NoiseModel()),
    }
    with (data_dir / "metadata.json").open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return metadata
