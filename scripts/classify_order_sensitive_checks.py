"""Enumerate the four-slot rank-two order classes used by CARMEN-Q.

The three nonzero columns of GF(2)^2 are encoded as 1, 2, and 3.  The
enumeration quotients full-rank length-four sequences by row-basis changes
(the six permutations induced by GL(2, 2)) and by time reversal.  It reports
the standard trellis-connectivity profile for one representative of every
class.  This is a finite structural check, not a proof of the AUDIT--RETURN
theorem.
"""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

from carmenq import trellis_connectivity_profile, trellis_connectivity_tau


NONZERO_COLUMNS = (1, 2, 3)


def _orbit(sequence: tuple[int, ...]) -> frozenset[tuple[int, ...]]:
    images: set[tuple[int, ...]] = set()
    for permutation in itertools.permutations(NONZERO_COLUMNS):
        relabel = dict(zip(NONZERO_COLUMNS, permutation))
        image = tuple(relabel[column] for column in sequence)
        images.add(image)
        images.add(image[::-1])
    return frozenset(images)


def _check_matrix(sequence: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    return (
        tuple((column >> 1) & 1 for column in sequence),
        tuple(column & 1 for column in sequence),
    )


def classify() -> dict[str, object]:
    """Return the complete quotient classification as JSON-ready data."""
    full_rank = tuple(
        sequence
        for sequence in itertools.product(NONZERO_COLUMNS, repeat=4)
        if len(set(sequence)) >= 2
    )
    seen: set[tuple[int, ...]] = set()
    classes: list[dict[str, object]] = []
    for sequence in full_rank:
        if sequence in seen:
            continue
        orbit = _orbit(sequence)
        seen.update(orbit)
        representative = min(orbit)
        matrix = _check_matrix(representative)
        classes.append(
            {
                "representative": representative,
                "orbit_size": len(orbit),
                "connectivity_profile": trellis_connectivity_profile(matrix),
                "tau": trellis_connectivity_tau(matrix),
            }
        )
    classes.sort(key=lambda row: row["representative"])
    return {
        "full_rank_sequences": len(full_rank),
        "equivalence_classes": len(classes),
        "tau_1_classes": sum(row["tau"] == 1 for row in classes),
        "tau_2_classes": sum(row["tau"] == 2 for row in classes),
        "classes": classes,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = classify()
    rendered = json.dumps(payload, indent=2) + "\n"
    print(rendered, end="")
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()
