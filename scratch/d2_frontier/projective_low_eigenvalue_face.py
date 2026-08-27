"""Stable rank/rank interval bounds near a zero eigenvalue.

The generic generalized-eigenvalue jet is not Lipschitz when the low coarse
eigenvalue ``y`` vanishes: its off-diagonal entry is proportional to
``sqrt(y)``.  This module keeps the mature generic kernel intact and adds a
small subclass which evaluates the ``y=0`` face exactly, then controls motion
away from that face by a directed matrix perturbation of order ``sqrt(y)``.
The comparison applies whether or not the certified box itself touches zero.

The subclass changes neither the feasible set nor the certificate trust
boundary.  Tangents and branch choices remain untrusted suggestions; every
accepted upper bound uses outward-expanded binary64 intervals.
"""

from __future__ import annotations

from projective_tangent_interval_certificate import (
    BoundResult,
    Interval,
    Jet,
    ONE,
    ZERO,
    RankRankBox,
    RankRankCertifier,
    add,
    constant,
    divide,
    down,
    jadd,
    jdivide,
    jmaximum,
    jmul,
    jsub,
    jsquare,
    mul,
    root,
    scale,
    square,
    sub,
    up,
    variable,
)


class LowEigenvalueFaceCertifier(RankRankCertifier):
    """Rank/rank certifier with a regular ``y=0`` comparison bound."""

    def _rank_low_zero_generalized_jet(
        self,
        high: Jet,
        sine: Jet,
        label: int,
        tangent: float,
    ) -> Jet:
        tangent_j = constant(Interval.exact_float(tangent))
        inverse = jdivide(constant(ONE), tangent_j)
        intercept = jadd(constant(ONE), inverse)
        slope = jsub(tangent_j, inverse)
        n_high = jadd(intercept, jmul(slope, high))
        n_low = intercept
        z_value = jsquare(sine)
        selected = jsub(constant(ONE), z_value) if label == 0 else z_value
        d_high = jsub(
            constant(self.level),
            jmul(constant(self.weight), jmul(high, selected)),
        )
        return jmaximum(
            jdivide(n_high, d_high),
            jdivide(n_low, constant(self.level)),
        )

    def _low_eigenvalue_face_jet_expression(
        self, box: RankRankBox, tangents: tuple[float, ...]
    ) -> Jet:
        if self.first_kind != "rank" or self.second_kind != "rank":
            raise ValueError("low-eigenvalue face currently covers rank/rank only")
        if box.y.lo != 0.0 or box.y.hi != 0.0:
            raise ValueError("low-eigenvalue face requires y=0 exactly")
        residual = variable(box.residual, 1)
        first_sine = variable(box.first_sine, 2)
        second_sine = variable(box.second_sine, 3)
        first_high = jsub(constant(ONE), residual)
        terms = [
            self._rank_low_zero_generalized_jet(
                first_high, first_sine, label, tangents[label]
            )
            for label in (0, 1)
        ]
        for label in (0, 1):
            terms.append(
                self._generalized_jet(
                    constant(ONE),
                    residual,
                    second_sine,
                    label,
                    tangents[2 + label],
                )
            )
        total = constant(ZERO)
        for term in terms:
            total = jadd(total, term)
        return jmul(constant(self.scale), total)

    def _sqrt_y_jet_expression(
        self,
        coordinates: tuple[Interval, Interval, Interval, Interval],
        tangents: tuple[float, ...],
    ) -> Jet:
        """Full rank/rank sum in the regular coordinate ``u = sqrt(y)``."""

        if self.first_kind != "rank" or self.second_kind != "rank":
            raise ValueError("sqrt-y form currently covers rank/rank only")
        u_value = variable(coordinates[0], 0)
        residual = variable(coordinates[1], 1)
        first_sine = variable(coordinates[2], 2)
        second_sine = variable(coordinates[3], 3)
        y_value = jsquare(u_value)
        groups = (
            (
                jsub(jsub(constant(ONE), y_value), residual),
                y_value,
                first_sine,
            ),
            (
                jsub(constant(ONE), y_value),
                jadd(y_value, residual),
                second_sine,
            ),
        )
        total = constant(ZERO)
        cursor = 0
        for high, low, sine in groups:
            for label in (0, 1):
                total = jadd(
                    total,
                    self._generalized_jet(
                        high, low, sine, label, tangents[cursor]
                    ),
                )
                cursor += 1
        return jmul(constant(self.scale), total)

    def _sqrt_y_mean_value_upper(
        self, box: RankRankBox, tangents: tuple[float, ...]
    ) -> float:
        """Directed mean-value enclosure after substituting ``y = u^2``."""

        coordinates = (
            root(box.y), box.residual, box.first_sine, box.second_sine
        )
        enclosure = self._sqrt_y_jet_expression(coordinates, tangents)
        anchors: list[float] = []
        deltas: list[Interval] = []
        for coordinate, derivative in zip(coordinates, enclosure.gradient):
            if derivative.hi < 0.0:
                anchors.append(coordinate.lo)
                deltas.append(Interval(0.0, up(coordinate.width)))
            elif derivative.lo > 0.0:
                anchors.append(coordinate.hi)
                deltas.append(Interval(down(-coordinate.width), 0.0))
            else:
                anchor = coordinate.midpoint
                anchors.append(anchor)
                deltas.append(
                    Interval(
                        down(coordinate.lo - anchor),
                        up(coordinate.hi - anchor),
                    )
                )
        anchor_coordinates = tuple(
            Interval.exact_float(value) for value in anchors
        )
        mean_value = self._sqrt_y_jet_expression(
            anchor_coordinates, tangents
        ).value
        for derivative, delta in zip(enclosure.gradient, deltas):
            mean_value = add(mean_value, mul(derivative, delta))
        return mean_value.hi

    def _low_eigenvalue_face_lipschitz_upper(
        self,
        box: RankRankBox,
        tangents: tuple[float, ...] | None = None,
    ) -> float:
        """Compare any low-``y`` box with its exactly diagonal ``y=0`` face.

        If ``D,D0 >= d I``, the generalized Rayleigh quotient obeys

        ``lambda(N,D) <= lambda(N0,D0)
                         + (||N-N0|| + lambda(N0,D0)||D-D0||)/d``.

        The only nonsmooth change is an off-diagonal square root.  The bound
        ``|sqrt(a)-sqrt(b)| <= sqrt(|a-b|)`` and the Frobenius norm enclose it
        without differentiating at zero.
        """

        if self.first_kind != "rank" or self.second_kind != "rank":
            raise ValueError("low-eigenvalue face currently covers rank/rank only")
        face_box = RankRankBox(
            ZERO, box.residual, box.first_sine, box.second_sine
        )
        selected = tangents or self._select_tangents(face_box)
        enclosure = self._low_eigenvalue_face_jet_expression(face_box, selected)
        anchors: list[float] = []
        deltas: list[Interval] = []
        for coordinate, derivative in zip(face_box.coordinates, enclosure.gradient):
            if derivative.hi < 0.0:
                anchors.append(coordinate.lo)
                deltas.append(Interval(0.0, up(coordinate.width)))
            elif derivative.lo > 0.0:
                anchors.append(coordinate.hi)
                deltas.append(Interval(down(-coordinate.width), 0.0))
            else:
                anchor = coordinate.midpoint
                anchors.append(anchor)
                deltas.append(
                    Interval(
                        down(coordinate.lo - anchor),
                        up(coordinate.hi - anchor),
                    )
                )
        anchor_box = RankRankBox(
            *(Interval.exact_float(value) for value in anchors)
        )
        mean_value = self._low_eigenvalue_face_jet_expression(
            anchor_box, selected
        ).value
        for derivative, delta in zip(enclosure.gradient, deltas):
            mean_value = add(mean_value, mul(derivative, delta))

        y = box.y
        residual = box.residual
        groups = (
            (
                sub(sub(ONE, y), residual),
                y,
                sub(ONE, residual),
                ZERO,
                box.first_sine,
            ),
            (
                sub(ONE, y),
                add(y, residual),
                ONE,
                residual,
                box.second_sine,
            ),
        )
        perturbation = ZERO
        for group, (high, low, high_zero, low_zero, sine) in enumerate(groups):
            z_value = square(sine)
            one_minus_z = sub(ONE, z_value)
            for label in (0, 1):
                first_split, second_split = (
                    (one_minus_z, z_value)
                    if label == 0 else (z_value, one_minus_z)
                )
                diagonal_change_squared = mul(
                    square(y),
                    add(square(first_split), square(second_split)),
                )
                if group == 0:
                    off_change_squared = mul(
                        mul(high, y), mul(z_value, one_minus_z)
                    )
                else:
                    off_change_squared = mul(
                        mul(y, sub(sub(ONE, residual), y)),
                        mul(z_value, one_minus_z),
                    )
                delta_g = root(
                    add(
                        diagonal_change_squared,
                        scale(off_change_squared, 2.0),
                    )
                )
                trace_current = add(
                    mul(high, first_split), mul(low, second_split)
                )
                trace_face = add(
                    mul(high_zero, first_split),
                    mul(low_zero, second_split),
                )
                # The split matrix is rank one, hence its only nonzero
                # eigenvalue is exactly its trace.
                spectral_cap = Interval(
                    0.0,
                    up(max(trace_current.hi, trace_face.hi)),
                )
                denominator_floor = sub(
                    self.level, mul(self.weight, spectral_cap)
                )
                if denominator_floor.lo <= 0.0:
                    raise ZeroDivisionError("face denominator floor is not positive")
                tangent_i = Interval.exact_float(selected[2 * group + label])
                inverse = divide(ONE, tangent_i)
                slope = sub(tangent_i, inverse)
                slope_cap = Interval(
                    0.0, up(max(abs(slope.lo), abs(slope.hi)))
                )
                numerator_change = mul(slope_cap, y)
                if group == 0:
                    face_term = self._rank_low_zero_generalized_jet(
                        constant(high_zero), constant(sine), label,
                        selected[2 * group + label],
                    ).value
                else:
                    face_term = self._generalized_value(
                        high_zero,
                        low_zero,
                        sine,
                        label,
                        selected[2 * group + label],
                    )
                face_lambda_cap = Interval(
                    0.0, up(max(0.0, face_term.hi))
                )
                denominator_change = mul(self.weight, delta_g)
                perturbation = add(
                    perturbation,
                    add(
                        divide(numerator_change, denominator_floor),
                        divide(
                            mul(face_lambda_cap, denominator_change),
                            denominator_floor,
                        ),
                    ),
                )
        return add(mean_value, mul(self.scale, perturbation)).hi

    def bound(
        self, box: RankRankBox, tangents: tuple[float, ...] | None = None
    ) -> BoundResult:
        contracted = self.contract(box)
        if contracted is None:
            return BoundResult(-float("inf"), "domain-empty", (), None)
        selected = tangents or self._select_tangents(contracted)
        if contracted.y.lo == 0.0:
            # Tangent choice is untrusted but strongly affects sharpness.  A
            # tangent selected at the midpoint y>0 can be unnecessarily loose
            # on the comparison face, so anchor this singular branch at y=0.
            face_box = RankRankBox(
                ZERO, contracted.residual,
                contracted.first_sine, contracted.second_sine,
            )
            selected = self._select_tangents(face_box)
        try:
            upper = self._low_eigenvalue_face_lipschitz_upper(
                contracted, selected
            )
        except (ArithmeticError, ValueError, ZeroDivisionError):
            upper = float("inf")
        sqrt_y_upper = float("inf")
        if contracted.y.lo > 0.0:
            try:
                sqrt_y_upper = self._sqrt_y_mean_value_upper(
                    contracted, selected
                )
            except (ArithmeticError, ValueError, ZeroDivisionError):
                pass
        if sqrt_y_upper <= 1.0:
            return BoundResult(sqrt_y_upper, "sqrt-y-mean-value", selected, None)
        # The face comparison is often the cheapest closing method.  If it
        # remains open, however, the full secular suite can be substantially
        # sharper, so retain the minimum rather than forcing further face
        # subdivision.
        if upper <= 1.0:
            return BoundResult(
                upper,
                "zero-face-lipschitz",
                selected,
                None,
            )
        result = super().bound(contracted, selected)
        if sqrt_y_upper < result.upper:
            result = BoundResult(
                sqrt_y_upper,
                "sqrt-y-mean-value",
                selected,
                None,
            )
        if upper >= result.upper:
            return result
        # Reuse the established counter label: both methods are exact
        # perturbation bounds from a singular boundary face.
        return BoundResult(
            upper,
            "zero-face-lipschitz",
            result.tangents,
            result.gradient,
        )

    def _hessian_split_index(
        self,
        box: RankRankBox,
        scales: tuple[float, ...],
        tangents: tuple[float, ...] | None = None,
    ) -> int | None:
        """Choose the coordinate dominating the local Taylor remainder.

        This is an untrusted branching heuristic: the returned children are
        still certified independently by outward interval bounds.  A dyadic
        point Hessian is used only to rank coordinates, so failure at an
        eigenvalue crossing safely falls back to the established splitter.
        """

        coordinates = box.coordinates
        bisectable = [
            coordinate.lo < coordinate.midpoint < coordinate.hi
            for coordinate in coordinates
        ]
        if not any(bisectable):
            return None
        center_box = RankRankBox(
            *(
                Interval.exact_float(coordinate.midpoint)
                for coordinate in coordinates
            )
        )
        selected = tangents or self._select_tangents(center_box)
        center_jet, _ = self._jet2_expression(center_box, selected)
        widths = [coordinate.width for coordinate in coordinates]
        scores = []
        for index, width in enumerate(widths):
            score = width * sum(
                max(
                    abs(center_jet.hessian[index][other].lo),
                    abs(center_jet.hessian[index][other].hi),
                )
                * widths[other]
                for other in range(len(coordinates))
            )
            scores.append(score if bisectable[index] else 0.0)
        if max(scores) <= 0.0:
            return None
        index = max(range(len(coordinates)), key=lambda item: scores[item])
        normalized = [
            coordinate.width / scale
            for coordinate, scale in zip(coordinates, scales)
        ]
        widest = max(
            (item for item in range(len(coordinates)) if bisectable[item]),
            key=lambda item: normalized[item],
        )
        if normalized[index] < normalized[widest] / 16.0:
            return widest
        return index


    def split(
        self,
        box: RankRankBox,
        scales: tuple[float, ...],
        tangents: tuple[float, ...] | None = None,
    ) -> tuple[str, RankRankBox, RankRankBox]:
        contracted = self.contract(box)
        if contracted is None:
            raise ValueError("cannot split an empty box")
        y = contracted.y
        if y.lo == 0.0 and y.lo < y.midpoint < y.hi:
            try:
                face_box = RankRankBox(
                    ZERO,
                    contracted.residual,
                    contracted.first_sine,
                    contracted.second_sine,
                )
                selected = self._select_tangents(face_box)
                face_upper = self._low_eigenvalue_face_lipschitz_upper(
                    face_box, selected
                )
                near_upper = self._low_eigenvalue_face_lipschitz_upper(
                    contracted, selected
                )
                if face_upper <= 1.0 < near_upper:
                    left = list(contracted.coordinates)
                    right = list(contracted.coordinates)
                    left[0] = Interval(y.lo, y.midpoint)
                    right[0] = Interval(y.midpoint, y.hi)
                    children = (RankRankBox(*left), RankRankBox(*right))
                    children = tuple(
                        self.contract(child) or child for child in children
                    )
                    return "y", children[0], children[1]
            except (ArithmeticError, ValueError, ZeroDivisionError):
                pass
        try:
            index = self._hessian_split_index(contracted, scales, tangents)
            if index is not None:
                names = ("y", "residual", "first_sine", "second_sine")
                coordinates = list(contracted.coordinates)
                middle = coordinates[index].midpoint
                left = list(coordinates)
                right = list(coordinates)
                left[index] = Interval(coordinates[index].lo, middle)
                right[index] = Interval(middle, coordinates[index].hi)
                children = (RankRankBox(*left), RankRankBox(*right))
                children = tuple(
                    self.contract(child) or child for child in children
                )
                return names[index], children[0], children[1]
        except (ArithmeticError, ValueError, ZeroDivisionError):
            pass
        return super().split(contracted, scales, tangents)


__all__ = ["LowEigenvalueFaceCertifier"]
