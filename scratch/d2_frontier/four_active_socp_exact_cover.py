"""Exact-dual SOCP cover for fully active four-effect terminal readouts.

The old global frontier used a 29-variable nonconvex spatial SCIP model for
this sector.  That geometry is unnecessary for an upper bound.  For terminal
effect traces ``w_i`` and syndrome priors ``p_i``, every physical strategy
satisfies three convex consequences:

* ``A <= sum_i w_i p_i`` (each rank-one effect has norm ``w_i``);
* the exact averaged aligned-projective comparisons
  ``A2_i = A-k_i(A-p_i)``, ``k_i=(1-w_i)/(2-w_i)``;
* the sixteen-path Hellinger RETURN cap.

The pulled-back prefix POVM also has every effect norm at most ``w_max``.
After fixing the order of the four prefix priors this is one linear audit
cap.  Every permutation of four binary labels is affine on ``F_2^2``;
translations canonicalise the prefix order and the six linear parts merely
permute the three nonzero syndrome labels.  Hence six SOCPs exhaust all 24
prefix-order cells.

Weight boxes are intersected with the exact sorted qubit-POVM polytope.  All
box-dependent coefficients are rounded in the relaxing direction.  Clarabel
only proposes canonical dual vectors.  Cone membership, stationarity, and
the resulting objective bounds are then recomputed with exact dyadic
``Fraction`` arithmetic, so no optimizer is trusted by the certificate.
"""

from __future__ import annotations

import argparse
import base64
from dataclasses import dataclass
from fractions import Fraction
import heapq
import itertools
import json
import math
from pathlib import Path
from typing import Any, Iterable
import zlib

import cvxpy as cp
from cvxpy.reductions.solvers.conic_solvers.clarabel_conif import CLARABEL
import numpy as np

from ternary_socp_exact_dual_probe import (
    canonical_hash,
    exact_dot,
    exact_sparse_stationarity,
    fraction_decimal,
    q,
    repair_dual_cones,
)


ROOT = Path(__file__).resolve().parent
SUPPORT_WEIGHT = Fraction(3, 5)
TARGET = Fraction(76652, 100000)
MAXIMUM_WEIGHT_FLOOR = Fraction(3533, 4000)
MINIMUM_ACTIVE_WEIGHT = Fraction(3, 10000)
PROJECTIVE_LINES = (
    (Fraction(11, 20), Fraction(7573, 10000)),
    (Fraction(3, 5), Fraction(766, 1000)),
)
NONZERO_PERMUTATIONS = tuple(itertools.permutations((1, 2, 3)))


def binary64_up(value: Fraction) -> float:
    candidate = float(value)
    while q(candidate) < value:
        candidate = float(np.nextafter(candidate, math.inf))
    return candidate


def binary64_down(value: Fraction) -> float:
    candidate = float(value)
    while q(candidate) > value:
        candidate = float(np.nextafter(candidate, -math.inf))
    return candidate


def fraction_pair(value: Fraction) -> list[int]:
    return [value.numerator, value.denominator]


@dataclass(frozen=True)
class WeightBox:
    """Box in ``(w0,w2,w3)``; ``w1=2-w0-w2-w3``."""

    bounds: tuple[tuple[Fraction, Fraction], ...]
    path: str = ""

    def serialise(self) -> dict[str, Any]:
        return {
            "path": self.path,
            "coordinates": [
                [fraction_pair(lower), fraction_pair(upper)]
                for lower, upper in self.bounds
            ],
        }

    @classmethod
    def deserialise(cls, payload: dict[str, Any]) -> "WeightBox":
        return cls(
            tuple(
                (
                    Fraction(*pair[0]),
                    Fraction(*pair[1]),
                )
                for pair in payload["coordinates"]
            ),
            str(payload["path"]),
        )


def determinant3(matrix: tuple[tuple[Fraction, ...], ...]) -> Fraction:
    a, b, c = matrix
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )


def solve3(
    matrix: tuple[tuple[Fraction, ...], ...],
    vector: tuple[Fraction, ...],
) -> tuple[Fraction, Fraction, Fraction] | None:
    determinant = determinant3(matrix)
    if determinant == 0:
        return None
    values = []
    for column in range(3):
        replaced = tuple(
            tuple(vector[row] if index == column else matrix[row][index] for index in range(3))
            for row in range(3)
        )
        values.append(determinant3(replaced) / determinant)
    return tuple(values)  # type: ignore[return-value]


def inequalities(
    box: WeightBox,
) -> tuple[tuple[tuple[Fraction, ...], Fraction], ...]:
    (xl, xu), (yl, yu), (zl, zu) = box.bounds
    return (
        ((-1, 0, 0), -xl),
        ((1, 0, 0), xu),
        ((0, -1, 0), -yl),
        ((0, 1, 0), yu),
        ((0, 0, -1), -zl),
        ((0, 0, 1), zu),
        # Sorted weights: w3 <= w2 <= w1 <= w0.
        ((0, -1, 1), Fraction(0)),
        ((1, 2, 1), Fraction(2)),
        ((-2, -1, -1), Fraction(-2)),
    )


def polytope_vertices(
    box: WeightBox,
) -> tuple[tuple[Fraction, Fraction, Fraction], ...]:
    rows = inequalities(box)
    vertices: set[tuple[Fraction, Fraction, Fraction]] = set()
    for selected in itertools.combinations(rows, 3):
        point = solve3(
            tuple(tuple(Fraction(item) for item in row[0]) for row in selected),
            tuple(row[1] for row in selected),
        )
        if point is None:
            continue
        if all(
            sum(Fraction(a) * value for a, value in zip(coefficients, point))
            <= right
            for coefficients, right in rows
        ):
            vertices.add(point)
    return tuple(sorted(vertices))


def weight_hull(
    box: WeightBox,
) -> tuple[tuple[Fraction, Fraction], ...] | None:
    vertices = polytope_vertices(box)
    if not vertices:
        return None
    weights = [
        (point[0], Fraction(2) - sum(point), point[1], point[2])
        for point in vertices
    ]
    result = tuple(
        (
            min(value[index] for value in weights),
            max(value[index] for value in weights),
        )
        for index in range(4)
    )
    if not all(
        MINIMUM_ACTIVE_WEIGHT <= lower <= upper <= 1
        for lower, upper in result
    ):
        raise ArithmeticError("exact weight hull escaped physical bounds")
    return result


def initial_box() -> WeightBox:
    floor = MAXIMUM_WEIGHT_FLOOR
    delta = MINIMUM_ACTIVE_WEIGHT
    return WeightBox(
        (
            (floor, Fraction(1)),
            (delta, (Fraction(2) - floor - delta) / 2),
            (delta, (Fraction(2) - floor) / 3),
        )
    )


def split_box(box: WeightBox, coordinate: int) -> tuple[WeightBox, WeightBox]:
    lower, upper = box.bounds[coordinate]
    middle = (lower + upper) / 2
    left = list(box.bounds)
    right = list(box.bounds)
    left[coordinate] = (lower, middle)
    right[coordinate] = (middle, upper)
    return WeightBox(tuple(left), box.path + "0"), WeightBox(tuple(right), box.path + "1")


def choose_split(box: WeightBox) -> int:
    root = initial_box()
    hull = weight_hull(box)
    if hull is None:
        raise ValueError("cannot split an empty weight box")
    # Coordinate hulls for w0,w2,w3.  Normalisation prevents the broad w2
    # chart from starving the more sensitive near-projective w0 direction.
    indices = (0, 2, 3)
    scores = []
    for coordinate, weight_index in enumerate(indices):
        width = hull[weight_index][1] - hull[weight_index][0]
        root_width = root.bounds[coordinate][1] - root.bounds[coordinate][0]
        scores.append(float(width / root_width))
    return max(range(3), key=scores.__getitem__)


class FourActiveOracle:
    """Persistent DPP SOCP for one canonical prefix-order cell."""

    def __init__(self) -> None:
        self.path = cp.Variable((4, 4), nonneg=True, name="path")
        self.audit = cp.Variable(nonneg=True, name="audit")
        self.returned = cp.Variable(nonneg=True, name="return")
        self.cross = cp.Variable(120, nonneg=True, name="hellinger_cross")
        self.syndrome_weight = cp.Parameter(4, nonneg=True)
        self.prefix_cap = cp.Parameter(4, nonneg=True)
        self.audit_coefficient = cp.Parameter((2, 4), nonneg=True)
        self.prior_coefficient = cp.Parameter((2, 4), nonneg=True)
        self.return_coefficient = cp.Parameter(2, nonneg=True)
        self.line_upper = cp.Parameter(2, nonneg=True)

        flat = cp.reshape(self.path, (16,), order="C")
        prefix = cp.sum(self.path, axis=1)
        self.syndrome = cp.hstack(
            [sum(self.path[z, z ^ s] for z in range(4)) for s in range(4)]
        )
        constraints: list[cp.Constraint] = [
            cp.sum(flat) == 1,
            self.audit <= 1,
            self.returned <= 1,
            self.syndrome <= self.audit,
            self.audit <= self.syndrome_weight @ self.syndrome,
            *(prefix[index] >= prefix[index + 1] for index in range(3)),
            self.audit <= self.prefix_cap @ prefix,
        ]
        cursor = 0
        for first in range(16):
            for second in range(first + 1, 16):
                constraints.append(
                    cp.SOC(
                        flat[first] + flat[second],
                        cp.hstack(
                            [
                                2 * self.cross[cursor],
                                flat[first] - flat[second],
                            ]
                        ),
                    )
                )
                cursor += 1
        if cursor != 120:
            raise AssertionError("wrong number of Hellinger cross terms")
        constraints.append(
            16 * self.returned <= cp.sum(flat) + 2 * cp.sum(self.cross)
        )
        for line in range(2):
            for syndrome in range(4):
                constraints.append(
                    self.audit_coefficient[line, syndrome] * self.audit
                    + self.prior_coefficient[line, syndrome]
                    * self.syndrome[syndrome]
                    + self.return_coefficient[line] * self.returned
                    <= self.line_upper[line]
                )

        objective_audit = binary64_up(SUPPORT_WEIGHT)
        objective_return = binary64_up(1 - SUPPORT_WEIGHT)
        objective = cp.Maximize(
            objective_audit * self.audit + objective_return * self.returned
        )
        base = cp.Problem(objective, constraints)
        variables = base.variables()
        if not all(bool(variable.attributes.get("nonneg")) for variable in variables):
            raise RuntimeError("every residual-controlled variable must be nonnegative")
        self.problem = cp.Problem(
            objective,
            [*constraints, *(variable <= 1 for variable in variables)],
        )
        if not self.problem.is_dpp():
            raise RuntimeError("four-active weight-box SOCP must be DPP")
        self.objective_encoding = (objective_audit, objective_return)

    def assign(
        self,
        hull: tuple[tuple[Fraction, Fraction], ...],
        permutation: tuple[int, int, int],
    ) -> dict[str, Any]:
        order = (0, *permutation)
        permuted = tuple(hull[index] for index in order)
        weight_upper = [binary64_up(upper) for _, upper in permuted]
        w0_lower, w0_upper = hull[0]
        prefix_exact = (
            w0_upper,
            w0_upper,
            2 * (1 - w0_lower),
            Fraction(0),
        )
        prefix_encoded = [binary64_up(value) for value in prefix_exact]
        audit_coefficients = np.zeros((2, 4))
        prior_coefficients = np.zeros((2, 4))
        return_coefficients = np.zeros(2)
        line_uppers = np.zeros(2)
        exact_coefficients: list[list[dict[str, Fraction]]] = []
        for line_index, (line_weight, line_upper) in enumerate(PROJECTIVE_LINES):
            return_coefficients[line_index] = binary64_down(1 - line_weight)
            line_uppers[line_index] = binary64_up(line_upper)
            rows = []
            for syndrome, (lower, _) in enumerate(permuted):
                loss = (1 - lower) / (2 - lower)
                audit_exact = line_weight * (1 - loss)
                prior_exact = line_weight * loss
                audit_coefficients[line_index, syndrome] = binary64_down(audit_exact)
                prior_coefficients[line_index, syndrome] = binary64_down(prior_exact)
                rows.append(
                    {
                        "loss_upper": loss,
                        "audit": audit_exact,
                        "prior": prior_exact,
                    }
                )
            exact_coefficients.append(rows)
        self.syndrome_weight.value = weight_upper
        self.prefix_cap.value = prefix_encoded
        self.audit_coefficient.value = audit_coefficients
        self.prior_coefficient.value = prior_coefficients
        self.return_coefficient.value = return_coefficients
        self.line_upper.value = line_uppers

        # Every comparison is weakened coefficientwise and every upper bound
        # is enlarged.  These exact checks make that direction auditable.
        for index, (_, upper) in enumerate(permuted):
            if q(weight_upper[index]) < upper:
                raise ArithmeticError("syndrome weight was not rounded upward")
        for index, exact in enumerate(prefix_exact):
            if q(prefix_encoded[index]) < exact:
                raise ArithmeticError("prefix cap was not rounded upward")
        for line_index, (line_weight, line_upper) in enumerate(PROJECTIVE_LINES):
            if q(return_coefficients[line_index]) > 1 - line_weight:
                raise ArithmeticError("RETURN coefficient was not rounded downward")
            if q(line_uppers[line_index]) < line_upper:
                raise ArithmeticError("line upper was not rounded upward")
            for syndrome in range(4):
                exact = exact_coefficients[line_index][syndrome]
                if q(audit_coefficients[line_index, syndrome]) > exact["audit"]:
                    raise ArithmeticError("audit coefficient was not relaxed")
                if q(prior_coefficients[line_index, syndrome]) > exact["prior"]:
                    raise ArithmeticError("prior coefficient was not relaxed")
        return {
            "syndrome_order": list(order),
            "weight_upper": [fraction_pair(q(value)) for value in weight_upper],
            "prefix_cap": [fraction_pair(q(value)) for value in prefix_encoded],
        }

    def canonical_data(self) -> dict[str, Any]:
        data, _, _ = self.problem.get_problem_data(cp.CLARABEL)
        return data


def exact_upper(data: dict[str, Any], dual: np.ndarray) -> tuple[Fraction, Fraction, Fraction]:
    residuals, correction = exact_sparse_stationarity(
        data["A"], dual, data["c"]
    )
    upper = exact_dot(np.asarray(data["b"]), dual) + correction
    maximum_residual = max(map(abs, residuals), default=Fraction(0))
    return upper, correction, maximum_residual


def encode_dual(vector: np.ndarray, dtype: str) -> str:
    array = np.asarray(vector, dtype="<f4" if dtype == "f32" else "<f8")
    return base64.b64encode(zlib.compress(array.tobytes(), 9)).decode("ascii")


def decode_dual(encoded: str, dtype: str) -> np.ndarray:
    raw = zlib.decompress(base64.b64decode(encoded))
    return np.frombuffer(raw, dtype="<f4" if dtype == "f32" else "<f8").astype(float)


class ExactCertifier:
    def __init__(self, target: Fraction = TARGET) -> None:
        self.oracle = FourActiveOracle()
        self.solver = CLARABEL()
        self.target = target

    def certify(
        self,
        hull: tuple[tuple[Fraction, Fraction], ...],
        permutation: tuple[int, int, int],
    ) -> dict[str, Any]:
        enclosure = self.oracle.assign(hull, permutation)
        data = self.oracle.canonical_data()
        result = self.solver.solve_via_data(
            data,
            warm_start=False,
            verbose=False,
            solver_opts={
                "tol_gap_abs": 1e-11,
                "tol_gap_rel": 1e-11,
                "tol_feas": 1e-11,
                "max_iter": 500,
            },
            solver_cache=None,
        )
        repaired, repaired_blocks = repair_dual_cones(
            np.asarray(result.z), data["dims"]
        )
        raw_upper, _, _ = exact_upper(data, repaired)
        storage_dtype = "f32"
        stored = decode_dual(encode_dual(repaired, storage_dtype), storage_dtype)
        stored, _ = repair_dual_cones(stored, data["dims"])
        upper, correction, residual = exact_upper(data, stored)
        if upper > self.target and raw_upper <= self.target:
            storage_dtype = "f64"
            stored = decode_dual(encode_dual(repaired, storage_dtype), storage_dtype)
            stored, _ = repair_dual_cones(stored, data["dims"])
            upper, correction, residual = exact_upper(data, stored)
        closed = upper <= self.target
        report: dict[str, Any] = {
            "syndrome_permutation": list(permutation),
            "coefficient_enclosure": enclosure,
            "canonical_shape": [int(data["A"].shape[0]), int(data["A"].shape[1])],
            "canonical_nonzeros": int(data["A"].nnz),
            "canonical_sha256": canonical_hash(data),
            "cone_dimensions": {
                "zero": int(data["dims"].zero),
                "nonnegative": int(data["dims"].nonneg),
                "soc": list(map(int, data["dims"].soc)),
            },
            "untrusted_solver_status": str(result.status),
            "untrusted_primal_objective": float(result.obj_val),
            "untrusted_dual_objective": float(result.obj_val_dual),
            "soc_heads_repaired": repaired_blocks,
            "certified_upper_fraction": fraction_pair(upper),
            "certified_upper_decimal": fraction_decimal(upper),
            "exact_residual_correction": fraction_pair(correction),
            "maximum_stationarity_residual_decimal": fraction_decimal(residual),
            "closed": closed,
            "trusted_optimizers": [],
            "untrusted_search_helpers": ["Clarabel dual-vector proposal"],
        }
        if closed:
            report.update(
                {
                    "dual_storage_dtype": storage_dtype,
                    "dual_zlib_base64": encode_dual(stored, storage_dtype),
                }
            )
        return report


def assess_box(certifier: ExactCertifier, box: WeightBox) -> dict[str, Any]:
    hull = weight_hull(box)
    if hull is None:
        return {"kind": "domain-empty", "box": box.serialise()}
    reports = [
        certifier.certify(hull, permutation)
        for permutation in NONZERO_PERMUTATIONS
    ]
    upper = max(
        Fraction(*report["certified_upper_fraction"]) for report in reports
    )
    return {
        "kind": "closed" if all(report["closed"] for report in reports) else "open",
        "box": box.serialise(),
        "exact_weight_hull": [
            [fraction_pair(lower), fraction_pair(upper)]
            for lower, upper in hull
        ],
        "maximum_certified_upper_fraction": fraction_pair(upper),
        "maximum_certified_upper_decimal": fraction_decimal(upper),
        "order_certificates": reports,
    }


def validate_leaf_tree(leaves: Iterable[dict[str, Any]]) -> None:
    paths = {str(leaf["box"]["path"]) for leaf in leaves}
    if len(paths) != len(list(leaves)):
        raise ArithmeticError("duplicate cover leaf path")
    pending = [""]
    while pending:
        path = pending.pop()
        if path in paths:
            continue
        if any(value.startswith(path) for value in paths):
            pending.extend((path + "0", path + "1"))
            continue
        raise ArithmeticError(f"weight-cover tree has a gap below {path!r}")


def write_payload(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", default="0.76652")
    parser.add_argument("--max-splits", type=int, default=5000)
    parser.add_argument("--checkpoint-every", type=int, default=10)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "four_active_socp_exact_cover_l060.json",
    )
    args = parser.parse_args()
    target = Fraction(args.target)
    certifier = ExactCertifier(target)
    root = initial_box()
    first = assess_box(certifier, root)
    leaves: list[dict[str, Any]] = []
    queue: list[tuple[float, int, WeightBox]] = []
    counter = 0
    if first["kind"] == "open":
        upper = Fraction(*first["maximum_certified_upper_fraction"])
        heapq.heappush(queue, (-float(upper), counter, root))
    else:
        leaves.append(first)
    splits = 0

    def payload() -> dict[str, Any]:
        complete = not queue
        closed = [leaf for leaf in leaves if leaf["kind"] == "closed"]
        maximum = max(
            (
                Fraction(*leaf["maximum_certified_upper_fraction"])
                for leaf in closed
            ),
            default=Fraction(0),
        )
        result = {
            "schema": "carmenq.four-active-weight-socp-exact-dual-cover.v1",
            "support_weight": fraction_pair(SUPPORT_WEIGHT),
            "target": fraction_pair(target),
            "maximum_weight_floor": fraction_pair(MAXIMUM_WEIGHT_FLOOR),
            "minimum_active_weight": fraction_pair(MINIMUM_ACTIVE_WEIGHT),
            "projective_lines": [
                [fraction_pair(weight), fraction_pair(upper)]
                for weight, upper in PROJECTIVE_LINES
            ],
            "prefix_order_reduction": (
                "AGL(2,2)=S4: translations canonicalise the prefix order; "
                "the six GL(2,2) parts permute the nonzero syndrome labels"
            ),
            "initial_box": root.serialise(),
            "boxes_split": splits,
            "leaf_count": len(leaves),
            "closed_leaf_count": len(closed),
            "domain_empty_leaf_count": sum(
                leaf["kind"] == "domain-empty" for leaf in leaves
            ),
            "boxes_remaining": len(queue),
            "maximum_open_upper": -queue[0][0] if queue else None,
            "maximum_certified_upper_decimal": fraction_decimal(maximum),
            "complete": complete,
            "all_cells_closed": complete,
            "leaves": sorted(leaves, key=lambda leaf: leaf["box"]["path"]),
            "trusted_optimizers": [],
            "untrusted_search_helpers": [
                "Clarabel dual-vector proposals",
                "best-first box ordering",
            ],
        }
        if complete:
            validate_leaf_tree(result["leaves"])
        return result

    while queue and splits < args.max_splits:
        _, _, box = heapq.heappop(queue)
        coordinate = choose_split(box)
        for child in split_box(box, coordinate):
            report = assess_box(certifier, child)
            if report["kind"] == "open":
                counter += 1
                upper = Fraction(*report["maximum_certified_upper_fraction"])
                heapq.heappush(queue, (-float(upper), counter, child))
            else:
                leaves.append(report)
        splits += 1
        if splits % args.checkpoint_every == 0:
            print(
                json.dumps(
                    {
                        "boxes_split": splits,
                        "closed_or_empty": len(leaves),
                        "open": len(queue),
                        "maximum_open_upper": -queue[0][0] if queue else None,
                    }
                ),
                flush=True,
            )
            write_payload(args.output, payload())
    result = payload()
    write_payload(args.output, result)
    print(json.dumps({key: value for key, value in result.items() if key != "leaves"}, indent=2))


if __name__ == "__main__":
    main()
