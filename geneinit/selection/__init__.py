"""Deterministic evidence-based candidate selection."""

from .core import (
    SelectionError,
    build_loci_and_edges,
    evaluate_and_select,
    load_component_rows,
)

__all__ = [
    "SelectionError",
    "build_loci_and_edges",
    "evaluate_and_select",
    "load_component_rows",
]
