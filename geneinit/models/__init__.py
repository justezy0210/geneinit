"""Predictor-independent structural candidate models."""

from .candidate import (
    CandidateError,
    CandidateModel,
    CandidateUniverse,
    CdsPart,
    build_candidate_universe,
    parse_predictor_gtf,
)

__all__ = [
    "CandidateError",
    "CandidateModel",
    "CandidateUniverse",
    "CdsPart",
    "build_candidate_universe",
    "parse_predictor_gtf",
]
