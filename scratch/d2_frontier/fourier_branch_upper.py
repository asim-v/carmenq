"""One spectral branch of the common-instrument Fourier upper bound.

The three nontrivial characters of Z_2^2 each have a qubit trace norm whose
active expression is either the positive scalar coefficient, the negative
scalar coefficient, or the Bloch-vector length.  This driver solves one member
of that finite cover for the fixed interior benchmark.
"""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import numpy as np

from choi_moment_reduced_upper import solve_povm
from two_block_choi_seesaw import canonical_three_effect_povm


BRANCH_NAMES = {
    "p": "scalar-positive",
    "n": "scalar-negative",
    "b": "bloch",
}
PRIOR_BOX = np.asarray(
    [
        [0.296875, 0.42596435546875],
        [0.224609375, 0.34832000732421875],
        [0.15234375, 0.258392333984375],
        [0.1083984375, 0.201324462890625],
    ],
    dtype=float,
)


def state_bounds_from_prior_box(prior_box: np.ndarray) -> np.ndarray:
    """Return coordinate bounds implied by subnormalised-state positivity."""

    bounds = np.zeros((4, 4, 2), dtype=float)
    bounds[:, 0] = prior_box
    bounds[:, 1:, 0] = -prior_box[:, 1, None]
    bounds[:, 1:, 1] = prior_box[:, 1, None]
    return bounds


def spectral_caps(
    code: str,
    plane_cap: int | None,
    sphere_cap: int | None,
    fine_sphere_cap: int | None = None,
) -> tuple[tuple[float, float, float, float] | None, ...]:
    """Caps after aligning the first and planarising the second Bloch mode."""

    caps: list[tuple[float, float, float, float] | None] = [None, None, None]
    bloch_indices = [index for index, letter in enumerate(code) if letter == "b"]
    if len(bloch_indices) >= 2 and plane_cap is not None:
        if not 0 <= plane_cap < 4:
            raise ValueError("plane cap must be 0,1,2,3")
        angle = (-3.0 + 2.0 * plane_cap) * np.pi / 8.0
        caps[bloch_indices[1]] = (
            float(np.cos(angle)),
            0.0,
            float(np.sin(angle)),
            float(np.cos(np.pi / 8.0)),
        )
    if sphere_cap is not None and fine_sphere_cap is not None:
        raise ValueError("choose either a coarse or a fine sphere cap")
    if len(bloch_indices) >= 3 and sphere_cap is not None:
        if not 0 <= sphere_cap < 6:
            raise ValueError("sphere cap must be 0,...,5")
        axis = np.zeros(3, dtype=float)
        axis[sphere_cap // 2] = 1.0 if sphere_cap % 2 == 0 else -1.0
        caps[bloch_indices[2]] = (
            float(axis[0]),
            float(axis[1]),
            float(axis[2]),
            float(1.0 / np.sqrt(3.0)),
        )
    if len(bloch_indices) >= 3 and fine_sphere_cap is not None:
        directions = [
            np.asarray(direction, dtype=float)
            for direction in itertools.product((-1.0, 0.0, 1.0), repeat=3)
            if direction != (0.0, 0.0, 0.0)
        ]
        if not 0 <= fine_sphere_cap < len(directions):
            raise ValueError("fine sphere cap must be 0,...,25")
        normal = directions[fine_sphere_cap]
        normal /= np.linalg.norm(normal)
        ratio = np.sqrt(3.0 / 2.0) - 1.0
        second = np.sqrt(2.0) - 1.0
        covering_cosine = 1.0 / np.sqrt(
            1.0 + second**2 + (ratio * (1.0 + second)) ** 2
        )
        caps[bloch_indices[2]] = (
            float(normal[0]),
            float(normal[1]),
            float(normal[2]),
            float(covering_cosine),
        )
    return tuple(caps)


def solve_branch(
    code: str,
    solver: str = "clarabel",
    plane_cap: int | None = None,
    sphere_cap: int | None = None,
    fine_sphere_cap: int | None = None,
) -> dict[str, object]:
    if len(code) != 3 or any(letter not in BRANCH_NAMES for letter in code):
        raise ValueError("branch code must contain exactly three letters from p,n,b")
    branch = tuple(BRANCH_NAMES[letter] for letter in code)
    payload = solve_povm(
        canonical_three_effect_povm(np.asarray([0.92, 0.64, 0.44])),
        0.55,
        (0, 1, 2, 3),
        "used",
        "quadratic",
        (0.25, 0.5, 1.0, 2.0, 4.0),
        state_bounds_from_prior_box(PRIOR_BOX),
        (),
        solver,
        False,
        fourier_trace_branches=branch,
        fourier_bloch_caps=spectral_caps(
            code, plane_cap, sphere_cap, fine_sphere_cap
        ),
    )
    payload["scope"] = "fixed ternary POVM, fixed prefix order, initial prior box"
    payload["branch_code"] = code
    payload["plane_cap"] = plane_cap
    payload["sphere_cap"] = sphere_cap
    payload["fine_sphere_cap"] = fine_sphere_cap
    payload["target"] = 0.758
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("branch")
    parser.add_argument("--plane-cap", type=int)
    parser.add_argument("--sphere-cap", type=int)
    parser.add_argument("--fine-sphere-cap", type=int)
    parser.add_argument("--solver", choices=("clarabel", "scs"), default="clarabel")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = solve_branch(
        args.branch.lower(),
        args.solver,
        args.plane_cap,
        args.sphere_cap,
        args.fine_sphere_cap,
    )
    rendered = json.dumps(payload, indent=2) + "\n"
    print(rendered)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()
