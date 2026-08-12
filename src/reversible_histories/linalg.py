"""Small dense-linear-algebra utilities used by the reference simulator.

The largest default Hilbert space has dimension 256, so explicit density
matrices remain practical and make every approximation visible in the code.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence

import numpy as np

Array = np.ndarray


def basis(index: int, dimension: int) -> Array:
    """Return a computational-basis ket."""
    ket = np.zeros(dimension, dtype=complex)
    ket[index] = 1.0
    return ket


def density(ket: Array) -> Array:
    """Return ``|ket><ket|`` after normalizing numerical input."""
    ket = np.asarray(ket, dtype=complex)
    norm = np.linalg.norm(ket)
    if norm == 0:
        raise ValueError("A zero vector does not define a quantum state.")
    ket = ket / norm
    return np.outer(ket, ket.conj())


def partial_trace(rho: Array, dimensions: Sequence[int], keep: Sequence[int]) -> Array:
    """Trace out every subsystem not listed in ``keep``.

    Returned subsystems are ordered as in ``dimensions``.  The implementation
    intentionally rejects a reordered ``keep`` list to avoid silent basis
    conventions in scientific output.
    """
    dimensions = tuple(int(d) for d in dimensions)
    keep = tuple(int(k) for k in keep)
    if keep != tuple(sorted(keep)):
        raise ValueError("keep must be in ascending subsystem order")
    if len(set(keep)) != len(keep) or any(k < 0 or k >= len(dimensions) for k in keep):
        raise ValueError("keep contains an invalid subsystem index")
    total = int(np.prod(dimensions))
    if rho.shape != (total, total):
        raise ValueError("rho has a shape inconsistent with dimensions")

    tensor = np.asarray(rho, dtype=complex).reshape(dimensions + dimensions)
    current_dimensions = list(dimensions)
    trace_out = [axis for axis in range(len(dimensions)) if axis not in keep]
    for axis in sorted(trace_out, reverse=True):
        count = len(current_dimensions)
        tensor = np.trace(tensor, axis1=axis, axis2=axis + count)
        current_dimensions.pop(axis)
    retained_dimension = int(np.prod(current_dimensions, dtype=int)) if current_dimensions else 1
    return tensor.reshape(retained_dimension, retained_dimension)


def von_neumann_entropy(rho: Array, base: float = 2.0) -> float:
    """Compute ``-Tr(rho log rho)`` with small numerical eigenvalues removed."""
    hermitian = (rho + rho.conj().T) / 2
    eigenvalues = np.linalg.eigvalsh(hermitian).real
    eigenvalues = eigenvalues[eigenvalues > 1e-13]
    if not eigenvalues.size:
        return 0.0
    entropy = -np.sum(eigenvalues * np.log(eigenvalues)) / np.log(base)
    return float(max(0.0, entropy))


def pure_state_fidelity(rho: Array, ket: Array) -> float:
    """Return ``<ket|rho|ket>`` for a normalized pure target."""
    ket = np.asarray(ket, dtype=complex)
    ket = ket / np.linalg.norm(ket)
    value = np.vdot(ket, rho @ ket).real
    return float(np.clip(value, 0.0, 1.0))


def conditional_holevo_information(
    states: Sequence[Array],
    predicates: Sequence[int],
    probabilities: Sequence[float] | None = None,
) -> float:
    r"""Return the conditional Holevo quantity :math:`\chi(H:R\mid P)`.

    ``P`` is a deterministic coarse graining of history label ``H``.  This is
    zero when the retained record states depend only on the allowed predicate,
    and one bit for four equiprobable, orthogonal labels grouped into two
    parity classes.
    """
    states = [np.asarray(state, dtype=complex) for state in states]
    predicates = np.asarray(predicates)
    count = len(states)
    if len(predicates) != count or count == 0:
        raise ValueError("states and predicates must have equal nonzero length")
    if probabilities is None:
        probabilities_array = np.full(count, 1.0 / count)
    else:
        probabilities_array = np.asarray(probabilities, dtype=float)
        probabilities_array /= probabilities_array.sum()

    answer = 0.0
    for predicate in np.unique(predicates):
        members = np.flatnonzero(predicates == predicate)
        class_probability = float(probabilities_array[members].sum())
        conditional = probabilities_array[members] / class_probability
        average = sum(weight * states[index] for weight, index in zip(conditional, members))
        mean_entropy = sum(
            weight * von_neumann_entropy(states[index])
            for weight, index in zip(conditional, members)
        )
        answer += class_probability * (von_neumann_entropy(average) - mean_entropy)
    return float(max(0.0, answer))


def make_permutation(
    dimensions: Sequence[int],
    transform: Callable[[tuple[int, ...]], tuple[int, ...]],
) -> Array:
    """Compile a reversible basis-state transform into an output-index map."""
    dimensions = tuple(dimensions)
    total = int(np.prod(dimensions))
    mapping = np.empty(total, dtype=int)
    for flat_index in range(total):
        values = tuple(int(v) for v in np.unravel_index(flat_index, dimensions))
        output = transform(values)
        if len(output) != len(dimensions):
            raise ValueError("transform returned a tuple of the wrong length")
        mapping[flat_index] = np.ravel_multi_index(output, dimensions)
    if len(np.unique(mapping)) != total:
        raise ValueError("transform is not a permutation")
    return mapping


def apply_permutation(rho: Array, mapping: Array, phases: Array | None = None) -> Array:
    """Apply ``|i> -> phases[i] |mapping[i]>`` to a density matrix."""
    mapping = np.asarray(mapping, dtype=int)
    if phases is None:
        phases = np.ones(mapping.size, dtype=complex)
    phases = np.asarray(phases, dtype=complex)
    output = np.empty_like(rho, dtype=complex)
    output[np.ix_(mapping, mapping)] = phases[:, None] * rho * phases.conj()[None, :]
    return output


def apply_imperfect_permutation(
    rho: Array,
    mapping: Array,
    failure_probability: float,
) -> Array:
    """Model inversion failure as a stochastic skipped inverse operation."""
    failure_probability = float(failure_probability)
    if not 0.0 <= failure_probability <= 1.0:
        raise ValueError("failure_probability must lie in [0, 1]")
    ideal = apply_permutation(rho, mapping)
    return (1.0 - failure_probability) * ideal + failure_probability * rho


def apply_local_kraus(
    rho: Array,
    dimensions: Sequence[int],
    target: int,
    operators: Sequence[Array],
) -> Array:
    """Apply a local Kraus channel without constructing full Kronecker matrices."""
    dimensions = tuple(dimensions)
    count = len(dimensions)
    local_dimension = dimensions[target]
    tensor = rho.reshape(dimensions + dimensions)
    other_axes = [axis for axis in range(2 * count) if axis not in (target, count + target)]
    permutation = [target, count + target] + other_axes
    transposed = np.transpose(tensor, permutation)
    remaining_shape = transposed.shape[2:]
    blocks = transposed.reshape(local_dimension, local_dimension, -1)
    transformed = np.zeros_like(blocks, dtype=complex)
    for operator in operators:
        operator = np.asarray(operator, dtype=complex)
        if operator.shape != (local_dimension, local_dimension):
            raise ValueError("Kraus operator has the wrong local dimension")
        transformed += np.einsum(
            "ai,ijr,bj->abr", operator, blocks, operator.conj(), optimize=True
        )
    restored = transformed.reshape((local_dimension, local_dimension) + remaining_shape)
    inverse_permutation = np.argsort(permutation)
    return np.transpose(restored, inverse_permutation).reshape(rho.shape)


def apply_local_unitary(
    rho: Array, dimensions: Sequence[int], target: int, unitary: Array
) -> Array:
    """Apply a unitary to one subsystem."""
    return apply_local_kraus(rho, dimensions, target, [unitary])


def dephasing_kraus(dimension: int, probability: float) -> list[Array]:
    """Kraus operators whose off-diagonal multiplier is ``1-probability``."""
    probability = float(probability)
    operators = [np.sqrt(1.0 - probability) * np.eye(dimension, dtype=complex)]
    for index in range(dimension):
        projector = np.zeros((dimension, dimension), dtype=complex)
        projector[index, index] = np.sqrt(probability)
        operators.append(projector)
    return operators


def depolarizing_kraus(dimension: int, probability: float) -> list[Array]:
    r"""Kraus form of ``(1-p) rho + p I/d \otimes Tr_target(rho)``."""
    probability = float(probability)
    operators = [np.sqrt(1.0 - probability) * np.eye(dimension, dtype=complex)]
    coefficient = np.sqrt(probability / dimension)
    for row in range(dimension):
        for column in range(dimension):
            operator = np.zeros((dimension, dimension), dtype=complex)
            operator[row, column] = coefficient
            operators.append(operator)
    return operators


def amplitude_damping_to_ground_kraus(
    dimension: int, probability: float
) -> list[Array]:
    """A qudit relaxation channel taking every excited level toward level zero."""
    probability = float(probability)
    diagonal = np.full(dimension, np.sqrt(1.0 - probability), dtype=complex)
    diagonal[0] = 1.0
    operators = [np.diag(diagonal)]
    for level in range(1, dimension):
        operator = np.zeros((dimension, dimension), dtype=complex)
        operator[0, level] = np.sqrt(probability)
        operators.append(operator)
    return operators


def uniform_overlap_channel(
    rho: Array, dimensions: Sequence[int], target: int, overlap: float
) -> Array:
    r"""Trace a uniform-overlap environment out of ``target``.

    If environment records obey ``<e_j|e_i> = overlap`` for ``i != j``, this
    Schur channel multiplies every off-diagonal target block by that overlap.
    It is completely positive for the supported range ``0 <= overlap <= 1``.
    """
    overlap = float(overlap)
    if not 0.0 <= overlap <= 1.0:
        raise ValueError("overlap must lie in [0, 1]")
    dimensions = tuple(dimensions)
    count = len(dimensions)
    gram = np.full((dimensions[target], dimensions[target]), overlap, dtype=complex)
    np.fill_diagonal(gram, 1.0)
    shape = [1] * (2 * count)
    shape[target] = dimensions[target]
    shape[count + target] = dimensions[target]
    tensor = rho.reshape(dimensions + dimensions)
    return (tensor * gram.reshape(shape)).reshape(rho.shape)
