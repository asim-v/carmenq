"""Evaluate the certified temporal-order gap at balanced task weight."""

from carmenq import (
    INTERLEAVED_BALANCED_COUNTEREXAMPLE,
    INTERLEAVED_ORDER_GAP_WEIGHT_THRESHOLD,
    grouped_frontier,
    interleaved_support_upper_bound,
)


weight = 0.5
grouped = grouped_frontier(weight).support_value
interleaved_upper = interleaved_support_upper_bound(weight)
interleaved_lower = INTERLEAVED_BALANCED_COUNTEREXAMPLE.support_value

print(
    "Certified interval starts above "
    f"lambda={INTERLEAVED_ORDER_GAP_WEIGHT_THRESHOLD:.9f}"
)
print(f"Grouped exact support:          {grouped:.12f}")
print(f"Interleaved rigorous upper:     {interleaved_upper:.12f}")
print(f"Interleaved verified lower:     {interleaved_lower:.12f}")
print(f"Certified grouped-order gap:    {grouped - interleaved_upper:.12f}")
