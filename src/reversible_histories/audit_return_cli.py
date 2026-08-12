"""Command-line interface for the causal audit--return benchmark."""

from __future__ import annotations

import argparse
import json

from .audit_return import (
    BenchmarkCounts,
    PhenomenologicalNoise,
    certify_classical_memory,
    classical_memory_frontier,
    collective_classical_record_bound,
    plan_experiment,
)


def _common_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--steps", type=int, required=True, help="Number of stream slots.")
    parser.add_argument("--weight", type=float, default=0.5, help="AUDIT score weight.")
    parser.add_argument("--alpha", type=float, default=0.01, help="False-positive level.")
    parser.add_argument(
        "--audit-systematic",
        type=float,
        default=0.0,
        help="Preregistered upper bias bound for AUDIT probability.",
    )
    parser.add_argument(
        "--return-systematic",
        type=float,
        default=0.0,
        help="Preregistered upper bias bound for RETURN probability.",
    )
    parser.add_argument(
        "--null-slack",
        type=float,
        default=0.0,
        help="Preregistered enlargement of the classical null bound.",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Plan or analyse the causal audit--return quantum-memory benchmark."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    bound = subparsers.add_parser("bound", help="Print the exact null support point.")
    bound.add_argument("--steps", type=int, required=True)
    bound.add_argument("--weight", type=float, default=0.5)

    plan = subparsers.add_parser("plan", help="Plan a fixed-sample experiment.")
    _common_options(plan)
    plan.add_argument("--audit-probability", type=float)
    plan.add_argument("--return-fidelity", type=float)
    plan.add_argument(
        "--forecast-model",
        action="store_true",
        help="Use the documented default phenomenological forecast.",
    )
    plan.add_argument("--beta", type=float, default=0.1, help="False-negative target.")

    analyse = subparsers.add_parser("analyse", help="Analyse fixed observed counts.")
    _common_options(analyse)
    analyse.add_argument("--audit-successes", type=int, required=True)
    analyse.add_argument("--audit-trials", type=int, required=True)
    analyse.add_argument("--return-successes", type=int, required=True)
    analyse.add_argument("--return-trials", type=int, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "bound":
        point = classical_memory_frontier(args.steps, args.weight)
        payload = {
            "n_steps": args.steps,
            "audit_weight": args.weight,
            "streaming_classical_memory": point.__dict__,
            "collective_classical_record_bound": collective_classical_record_bound(
                args.weight
            ),
            "coherent_memory_algebraic_bound": 1.0,
        }
    elif args.command == "plan":
        if args.forecast_model:
            audit_probability, return_fidelity = PhenomenologicalNoise().point(
                args.steps
            )
        else:
            if args.audit_probability is None or args.return_fidelity is None:
                raise SystemExit(
                    "Provide both --audit-probability and --return-fidelity, "
                    "or select --forecast-model."
                )
            audit_probability = args.audit_probability
            return_fidelity = args.return_fidelity
        payload = plan_experiment(
            args.steps,
            audit_probability,
            return_fidelity,
            args.weight,
            args.alpha,
            args.beta,
            audit_systematic=args.audit_systematic,
            return_systematic=args.return_systematic,
            null_slack=args.null_slack,
        ).to_dict()
    else:
        counts = BenchmarkCounts(
            args.audit_successes,
            args.audit_trials,
            args.return_successes,
            args.return_trials,
        )
        payload = certify_classical_memory(
            counts,
            args.steps,
            args.weight,
            args.alpha,
            audit_systematic=args.audit_systematic,
            return_systematic=args.return_systematic,
            null_slack=args.null_slack,
        ).to_dict()
    print(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
