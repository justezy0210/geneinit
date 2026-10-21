"""Same-run balanced-base selection used by the v2.2H supplement contract."""

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import combinations
from typing import Mapping, Sequence

from geneinit.selection.core import Edge, Record, _integer, _json_ints, _numeric


BASE_AUDIT_COLUMNS = (
    "candidate_id", "structure_id", "seqid", "strand", "start", "end",
    "locus_id", "eligible", "selected", "initial_status", "decision_status",
    "routes", "conflict_types", "reasons",
)


@dataclass(slots=True)
class BaseCall:
    record: Record
    eligible: bool = False
    selected: bool = False
    status: str = "DEFER_NO_EXTERNAL"
    initial_status: str = "DEFER_NO_EXTERNAL"
    routes: list[str] = field(default_factory=list)
    conflicts: set[str] = field(default_factory=set)
    reasons: set[str] = field(default_factory=set)
    protein_rank: int = 0
    rna_rank: int = 0
    supported_single_route: str = ""


def _core_metrics(record: Record, minimum: float) -> bool:
    if record.row.get("protein_chain_id", "NA") == "NA":
        return False
    values = (
        _numeric(record.row, "protein_model_precision", -1.0),
        _numeric(record.row, "protein_chain_recall", -1.0),
        _numeric(record.row, "protein_cds_f1", -1.0),
    )
    if min(values) < minimum:
        return False
    comparable = _integer(record.row, "protein_phase_comparable_parts")
    return comparable == 0 or _numeric(record.row, "protein_phase_fraction", -1.0) >= minimum


def _protein_strong(record: Record) -> bool:
    if not _core_metrics(record, 0.8):
        return False
    return not record.introns or _numeric(record.row, "protein_intron_precision", -1.0) >= 0.8


def _protein_supported_terminal(record: Record) -> bool:
    return (
        bool(record.introns)
        and _integer(record.row, "translated_protein_aa") >= 100
        and _core_metrics(record, 0.7)
        and _numeric(record.row, "protein_intron_precision", -1.0) >= 0.7
        and _numeric(record.row, "protein_intron_recall", -1.0) >= 0.7
        and _integer(record.row, "protein_start_match") + _integer(record.row, "protein_stop_match") >= 1
    )


def _protein_predictor_rescue(record: Record) -> bool:
    return (
        bool(record.introns)
        and _integer(record.row, "translated_protein_aa") >= 100
        and _core_metrics(record, 0.6)
        and _numeric(record.row, "protein_intron_precision", -1.0) >= 0.6
        and _numeric(record.row, "protein_intron_recall", -1.0) >= 0.6
        and _integer(record.row, "predictor_tier") >= 2
    )


def _strict_support(record: Record) -> bool:
    return bool(record.introns) and len(record.strong_introns) == len(record.introns)


def _full_support(record: Record) -> bool:
    if not record.introns:
        return False
    hq = _json_ints(record.row, "rna_hq_read_counts")
    samples = _json_ints(record.row, "rna_sample_hq_ge1_counts")
    if len(record.introns) == 1:
        return hq[0] >= 2 and samples[0] >= 2
    return all(value >= 1 for value in hq)


def _rna_predictor_complement(record: Record) -> bool:
    if len(record.introns) < 2 or _integer(record.row, "translated_protein_aa") < 100:
        return False
    hq = _json_ints(record.row, "rna_hq_read_counts")
    exact = sum(value >= 1 for value in hq)
    return (
        exact / len(record.introns) >= 0.8
        and sum(value == 0 for value in hq) <= 2
        and _integer(record.row, "predictor_tier") >= 2
    )


def _rna_protein_complement(record: Record) -> bool:
    if not record.introns:
        return False
    strong = len(record.strong_introns)
    phase_ok = (
        _integer(record.row, "protein_phase_comparable_parts") == 0
        or _numeric(record.row, "protein_phase_fraction", -1.0) >= 0.8
    )
    return (
        len(record.introns) - strong == 1
        and strong / len(record.introns) >= 0.8
        and _numeric(record.row, "protein_intron_precision", -1.0) >= 1.0
        and _numeric(record.row, "protein_intron_recall", -1.0) >= 1.0
        and _numeric(record.row, "protein_chain_recall", -1.0) >= 0.8
        and _numeric(record.row, "protein_cds_f1", -1.0) >= 0.8
        and phase_ok
    )


def _balanced_routes(record: Record) -> list[str]:
    if not record.introns or _integer(record.row, "translated_protein_aa") < 100:
        return []
    protein_tier = _integer(record.row, "protein_tier")
    predictor_tier = _integer(record.row, "predictor_tier")
    rna_tier = _integer(record.row, "rna_tier")
    recurrent = len(record.recurrent_introns)
    routes: list[str] = []
    if protein_tier >= 2 and predictor_tier >= 2:
        routes.append("PP_BALANCED")
    if protein_tier >= 2 and recurrent / len(record.introns) >= 0.8 and len(record.introns) - recurrent <= 1:
        routes.append("RP_BALANCED")
    if protein_tier >= 2 and rna_tier >= 3 and predictor_tier >= 1:
        routes.append("RPP_BALANCED")
    return routes


def _noncanonical_rescue(record: Record, strict: bool, full: bool) -> bool:
    noncanonical = _integer(record.row, "noncanonical_splice_count")
    if not noncanonical:
        return False
    motifs = tuple(__import__("json").loads(record.row["splice_motifs"]))
    hq = _json_ints(record.row, "rna_hq_read_counts")
    reps = _json_ints(record.row, "rna_condition_ge2rep_hq_ge1_counts")
    if len(motifs) != len(record.introns):
        return False
    noncanonical_supported = all(
        hq[index] >= 3 and reps[index] >= 1
        for index, motif in enumerate(motifs)
        if motif not in {"GT-AG", "GC-AG", "AT-AC"}
    )
    return noncanonical_supported and (strict or full)


from geneinit.evidence.te_classification import TE_COLUMNS, confirmed_only, repeat_neutral, partition_pools


def evaluate_base_record(record: Record) -> BaseCall:
    call = BaseCall(record)
    if (
        record.row["orf_status"] != "complete"
        or _integer(record.row, "common_qc_pass") != 1
        or _integer(record.row, "ambiguous_splice_count") != 0
    ):
        call.status = "REJECT_ORF"
        call.initial_status = call.status
        return call
    strict = _strict_support(record)
    full = _full_support(record)
    strong_protein = _protein_strong(record)
    terminals = _integer(record.row, "protein_start_match") + _integer(record.row, "protein_stop_match")
    balanced = _balanced_routes(record)
    if not record.introns:
        if strong_protein and terminals >= 1:
            call.routes.append("P_STRONG")
            call.status = "KEEP_P_STRONG"
            call.protein_rank = 3
        elif (
            _core_metrics(record, 0.9)
            and _integer(record.row, "predictor_tier") == 3
            and (repeat_neutral(record) or _numeric(record.row, "repeat_cds_fraction", 1.0) < 0.2)
        ):
            call.routes.append("P_SINGLE_RESCUE")
            call.status = "KEEP_P_SINGLE_RESCUE"
            call.protein_rank = 2
        else:
            quality = min(
                _numeric(record.row, key, -1.0)
                for key in (
                    "protein_model_precision", "protein_chain_recall",
                    "protein_cds_f1", "protein_phase_fraction",
                )
            )
            if (
                _integer(record.row, "translated_protein_aa") >= 100
                and (repeat_neutral(record) or _numeric(record.row, "repeat_cds_fraction", 1.0) < 0.1)
                and ((terminals >= 1 and quality >= 0.7) or (terminals == 0 and quality >= 0.85 and _integer(record.row, "predictor_tier") == 3))
            ):
                call.supported_single_route = "P_SINGLE_SUPPORTED"
            call.status = "REVIEW_SINGLE_EXON"
            call.initial_status = call.status
            return call
    else:
        if strict:
            call.routes.append("R_STRICT")
            call.status = "KEEP_R_STRICT"
            call.rna_rank = 3
        elif full:
            call.routes.append("R_FULL")
            call.status = "KEEP_R_FULL"
            call.rna_rank = 2
        if strong_protein:
            call.routes.append("P_STRONG")
            call.status = "KEEP_P_STRONG" if not call.routes[:-1] else call.status
            call.protein_rank = 3
        elif _protein_supported_terminal(record):
            call.routes.append("P_SUPPORTED_TERMINAL")
            call.status = "KEEP_P_SUPPORTED_TERMINAL" if not call.routes[:-1] else call.status
            call.protein_rank = 2
        elif _protein_predictor_rescue(record):
            call.routes.append("P_PREDICTOR_RESCUE")
            call.status = "KEEP_P_PREDICTOR_RESCUE" if not call.routes[:-1] else call.status
            call.protein_rank = 1
        if _rna_predictor_complement(record):
            call.routes.append("R_PREDICTOR_COMPLEMENT")
            call.status = "KEEP_R_PREDICTOR_COMPLEMENT" if len(call.routes) == 1 else call.status
            call.rna_rank = max(call.rna_rank, 1)
        if _rna_protein_complement(record):
            call.routes.append("RNA_PROTEIN_COMPLEMENT")
            call.status = "KEEP_RC" if len(call.routes) == 1 else call.status
            call.rna_rank = max(call.rna_rank, 1)
        for route in balanced:
            call.routes.append(route)
            call.rna_rank = max(call.rna_rank, 1 if route != "RPP_BALANCED" else 2)
            call.protein_rank = max(call.protein_rank, 1)
        if _integer(record.row, "noncanonical_splice_count"):
            if _noncanonical_rescue(record, strict, full):
                call.routes.append("NONCANONICAL_RNA")
                call.status = "KEEP_NONCANONICAL_RNA"
            else:
                call.status = "REVIEW_NONCANONICAL"
                call.initial_status = call.status
                return call
    if not call.routes:
        rna_moderate = bool(record.introns) and all(
            count >= 1 for count in _json_ints(record.row, "rna_hq_read_counts")
        )
        protein_moderate = _integer(record.row, "protein_tier") >= 2
        call.status = "DEFER_MODERATE" if rna_moderate or protein_moderate else "DEFER_NO_EXTERNAL"
        call.initial_status = call.status
        return call
    # Balanced routes are rescue routes: the historical v3.2 contract applied
    # them only to records that would otherwise be DEFER_MODERATE. Stronger
    # RNA/protein routes keep their own, higher-priority status.
    primary = [route for route in call.routes if route not in {"RP_BALANCED", "PP_BALANCED", "RPP_BALANCED"}]
    if not primary:
        rna_moderate = bool(record.introns) and all(
            count >= 1 for count in _json_ints(record.row, "rna_hq_read_counts")
        )
        protein_moderate = _integer(record.row, "protein_tier") >= 2
        if not (rna_moderate or protein_moderate):
            call.status = "DEFER_NO_EXTERNAL"
            call.initial_status = call.status
            return call
        route_priority = ("RP_BALANCED", "PP_BALANCED", "RPP_BALANCED")
        chosen = next(route for route in route_priority if route in balanced)
        call.status = f"KEEP_{chosen}"
    if not repeat_neutral(record) and _numeric(record.row, "repeat_cds_fraction") >= 0.8:
        call.initial_status = call.status
        call.status = "TE_RISK"
        return call
    call.initial_status = call.status
    call.eligible = True
    return call


def _protein_evidence_rank(call: BaseCall) -> int:
    routes = set(call.routes)
    if "P_STRONG" in routes:
        return 3
    if "P_SUPPORTED_TERMINAL" in routes:
        return 2
    if routes & {
        "P_PREDICTOR_RESCUE", "PP_BALANCED", "RP_BALANCED",
        "RPP_BALANCED", "P_SINGLE_RESCUE", "P_SINGLE_SUPPORTED",
    }:
        return 1
    return 0


def _rna_evidence_rank(call: BaseCall) -> int:
    routes = set(call.routes)
    if routes & {"R_STRICT", "NONCANONICAL_RNA"}:
        return 3
    if "R_FULL" in routes:
        return 2
    if routes & {"R_PREDICTOR_COMPLEMENT", "RP_BALANCED", "RPP_BALANCED"}:
        return 1
    return 0


def _protein_record_rank(call: BaseCall) -> tuple[object, ...]:
    record = call.record
    return (
        _protein_evidence_rank(call),
        _integer(record.row, "protein_start_match")
        + _integer(record.row, "protein_stop_match"),
        _integer(record.row, "protein_tier"),
        _numeric(record.row, "protein_model_precision", -1.0),
        _numeric(record.row, "protein_chain_recall", -1.0),
        _numeric(record.row, "protein_cds_f1", -1.0),
        _numeric(record.row, "protein_phase_fraction", -1.0),
    )


def _same_chain_rank(call: BaseCall) -> tuple[object, ...]:
    record = call.record
    terminal_count = (
        _integer(record.row, "protein_start_match")
        + _integer(record.row, "protein_stop_match")
    )
    return (
        int(terminal_count == 2),
        int(terminal_count >= 1),
        _protein_evidence_rank(call),
        int(_integer(record.row, "predictor_tier") == 3),
        _numeric(record.row, "protein_model_precision", -1.0),
        _numeric(record.row, "protein_chain_recall", -1.0),
        _numeric(record.row, "protein_cds_f1", -1.0),
    )


def _isoform_rank(call: BaseCall) -> tuple[object, ...]:
    return (
        _rna_evidence_rank(call),
        _protein_evidence_rank(call),
        len(call.record.recurrent_introns),
        _integer(call.record.row, "predictor_tier"),
        _protein_record_rank(call),
    )


def _single_rescue_sort_key(call: BaseCall) -> tuple[object, ...]:
    record = call.record
    quality = min(
        _numeric(record.row, key, -1.0)
        for key in (
            "protein_model_precision", "protein_chain_recall",
            "protein_cds_f1", "protein_phase_fraction",
        )
    )
    terminal_count = (
        _integer(record.row, "protein_start_match")
        + _integer(record.row, "protein_stop_match")
    )
    return (
        -terminal_count,
        -quality,
        -_integer(record.row, "predictor_tier"),
        0.0 if repeat_neutral(record) else _numeric(record.row, "repeat_cds_fraction"),
        -_integer(record.row, "translated_protein_aa"),
        record.key,
    )


def _unique_isoform_evidence(call: BaseCall, anchor: BaseCall) -> bool:
    unique_introns = set(call.record.introns) - set(anchor.record.introns)
    recurrent = bool(unique_introns & set(call.record.recurrent_introns))
    strong_protein = _protein_evidence_rank(call) == 3 and (
        (
            bool(unique_introns)
            and _numeric(call.record.row, "protein_intron_precision", -1.0) >= 1.0
            and _numeric(call.record.row, "protein_intron_recall", -1.0) >= 0.8
        )
        or (
            not unique_introns
            and _integer(call.record.row, "protein_start_match") == 1
            and _integer(call.record.row, "protein_stop_match") == 1
        )
    )
    return recurrent or strong_protein


def _spanning_introns(fusion: Record, left: Record, right: Record) -> tuple[tuple[int, int], ...]:
    if right.start < left.start:
        left, right = right, left
    return tuple(
        intron for intron in fusion.introns
        if left.end < intron[0] and intron[1] < right.start
    )


def _mark_rejected(call: BaseCall, status: str, reason: str) -> None:
    call.selected = False
    call.status = status
    call.reasons.add(reason)


def select_balanced_base(
    loci: Mapping[str, Sequence[Record]], edges: Sequence[Edge], policy=None
) -> tuple[set[str], dict[str, BaseCall]]:
    if policy is not None and confirmed_only(policy):
        for members in loci.values():
            for record in members:
                record.row["te_policy_id"] = policy["te_policy_id"]
        loci, edges = partition_pools(loci, edges)
    calls = {record.key: evaluate_base_record(record) for members in loci.values() for record in members}
    by_locus_edges: dict[str, list[Edge]] = __import__("collections").defaultdict(list)
    for edge in edges:
        by_locus_edges[edge.locus_id].append(edge)
    for locus_id, members in sorted(loci.items()):
        locus_calls = [calls[record.key] for record in members]
        edge_list = by_locus_edges.get(locus_id, ())
        relations = {
            tuple(sorted((edge.left, edge.right))): edge.relation
            for edge in edge_list
        }
        adjacency: dict[str, set[str]] = {call.record.key: set() for call in locus_calls}
        for edge in edge_list:
            adjacency[edge.left].add(edge.right)
            adjacency[edge.right].add(edge.left)

        # Structural flags are computed from every candidate, not only from
        # initially selected records.  A spanning model with two mutually
        # noncompeting neighbors is a fusion/split candidate.
        for call in locus_calls:
            neighbors = sorted(adjacency[call.record.key])
            if any(
                tuple(sorted((left, right))) not in relations
                for left, right in combinations(neighbors, 2)
            ):
                call.conflicts.add("FUSION_SPLIT_CANDIDATE")

        # New-locus single-exon rescue is permitted only when no main model is
        # already active and the candidate is not part of a fusion/split case.
        if not any(call.eligible for call in locus_calls):
            supported = [
                call for call in locus_calls
                if call.supported_single_route
                and "FUSION_SPLIT_CANDIDATE" not in call.conflicts
            ]
            if supported:
                winner = min(supported, key=_single_rescue_sort_key)
                winner.eligible = True
                winner.routes.append(winner.supported_single_route)
                winner.status = "KEEP_P_SINGLE_RESCUE"
                winner.initial_status = winner.status
                winner.protein_rank = 2
                for loser in supported:
                    if loser is not winner:
                        _mark_rejected(
                            loser, "CONFLICT_REVIEW",
                            "MULTIPLE_SINGLE_EXON_RESCUE_CANDIDATES",
                        )

        for call in locus_calls:
            call.selected = call.eligible

        # Same intron chain: RNA cannot determine alternative CDS boundaries.
        same_chain_groups: dict[tuple[tuple[int, int], ...], list[BaseCall]] = __import__("collections").defaultdict(list)
        for call in locus_calls:
            if call.record.introns:
                same_chain_groups[call.record.introns].append(call)
        for group in same_chain_groups.values():
            active = [call for call in group if call.selected]
            if len(active) < 2:
                continue
            for call in active:
                call.conflicts.add("SAME_INTRON_CHAIN_CDS_CONFLICT")
            best_rank = max(_same_chain_rank(call) for call in active)
            representative = min(
                (call for call in active if _same_chain_rank(call) == best_rank),
                key=lambda call: call.record.key,
            )
            representative.reasons.add("CDS_BOUNDARY_CONFLICT_REPRESENTATIVE")
            for call in active:
                if call is not representative:
                    _mark_rejected(
                        call, "CONFLICT_REVIEW", "RNA_CANNOT_RESOLVE_CDS_CONFLICT"
                    )

        # Fusion/split precedence is evaluated before ordinary isoform and hard
        # overlap resolution, matching the structural decision order.
        by_key = {call.record.key: call for call in locus_calls}
        for fusion in sorted(locus_calls, key=lambda call: (call.record.start, call.record.end, call.record.key)):
            neighbors = [
                by_key[key] for key in sorted(adjacency[fusion.record.key])
                if by_key[key].selected
            ]
            for left, right in combinations(neighbors, 2):
                if tuple(sorted((left.record.key, right.record.key))) in relations:
                    continue
                spanning = _spanning_introns(fusion.record, left.record, right.record)
                if not spanning:
                    continue
                for call in (fusion, left, right):
                    call.conflicts.add("FUSION_SPLIT_CANDIDATE")
                spanning_strict = bool(set(spanning) & set(fusion.record.strong_introns))
                split_strong = all(
                    _protein_evidence_rank(call) == 3
                    and _integer(call.record.row, "protein_start_match") == 1
                    and _integer(call.record.row, "protein_stop_match") == 1
                    for call in (left, right)
                )
                if spanning_strict and split_strong:
                    for call in (fusion, left, right):
                        _mark_rejected(call, "CONFLICT_REVIEW", "FUSION_SPLIT_BOTH_STRONG")
                elif spanning_strict:
                    fusion.reasons.add("STRICT_SPANNING_INTRON_FUSION_PRIORITY")
                    for call in (left, right):
                        if call.selected:
                            _mark_rejected(
                                call, "CONFLICT_REVIEW",
                                "STRICT_SPANNING_INTRON_FUSION_PRIORITY",
                            )
                elif left.selected and right.selected:
                    _mark_rejected(
                        fusion,
                        "FUSION_PROTEIN_REVIEW" if _protein_evidence_rank(fusion) == 3 else "CONFLICT_REVIEW",
                        "SPLIT_MODELS_WITHOUT_SPANNING_RNA",
                    )
                elif _protein_evidence_rank(fusion) == 3:
                    _mark_rejected(
                        fusion, "FUSION_PROTEIN_REVIEW", "PROTEIN_ONLY_FUSION"
                    )
                break

        # Alternative isoforms are resolved as connected components.  One
        # evidence-ranked anchor is kept; every additional isoform must carry a
        # recurrent unique intron or isoform-specific strong protein support.
        active = [call for call in locus_calls if call.selected]
        alt_adjacency: dict[str, set[str]] = {call.record.key: set() for call in active}
        for edge in edge_list:
            if edge.relation != "ALT_SPLICE":
                continue
            if by_key[edge.left].selected and by_key[edge.right].selected:
                alt_adjacency[edge.left].add(edge.right)
                alt_adjacency[edge.right].add(edge.left)
                by_key[edge.left].conflicts.add("ALTERNATIVE_SPLICE")
                by_key[edge.right].conflicts.add("ALTERNATIVE_SPLICE")
        remaining = set(alt_adjacency)
        while remaining:
            seed = min(remaining)
            stack = [seed]
            component: set[str] = set()
            while stack:
                key = stack.pop()
                if key in component:
                    continue
                component.add(key)
                stack.extend(sorted(alt_adjacency[key] - component, reverse=True))
            remaining -= component
            if len(component) < 2:
                continue
            isoforms = [by_key[key] for key in sorted(component) if by_key[key].selected]
            if len(isoforms) < 2:
                continue
            best_rank = max(_isoform_rank(call) for call in isoforms)
            anchor = min(
                (call for call in isoforms if _isoform_rank(call) == best_rank),
                key=lambda call: call.record.key,
            )
            anchor.reasons.add("ISOFORM_BASE_REPRESENTATIVE")
            for call in isoforms:
                if call is anchor:
                    continue
                if _unique_isoform_evidence(call, anchor):
                    call.reasons.add("ISOFORM_UNIQUE_EVIDENCE")
                else:
                    _mark_rejected(
                        call, "CONFLICT_REVIEW", "ISOFORM_WITHOUT_UNIQUE_EVIDENCE"
                    )

        # v3.2 narrow exception: a strong-protein single-exon model may beat
        # only moderate balanced spliced hard-overlap competitors.
        def strong_single(call: BaseCall) -> bool:
            return (
                call.selected and not call.record.introns
                and "FUSION_SPLIT_CANDIDATE" not in call.conflicts
                and call.initial_status == "KEEP_P_STRONG"
                and _protein_evidence_rank(call) == 3
                and _rna_evidence_rank(call) == 0
            )

        def moderate_spliced(call: BaseCall) -> bool:
            return (
                call.selected and bool(call.record.introns)
                and "FUSION_SPLIT_CANDIDATE" not in call.conflicts
                and call.initial_status in {"KEEP_RP_BALANCED", "KEEP_RPP_BALANCED"}
                and _rna_evidence_rank(call) == 1
                and _protein_evidence_rank(call) < 3
            )

        protected: set[str] = set()
        for single in locus_calls:
            if not strong_single(single):
                continue
            competitors = [
                by_key[key] for key in adjacency[single.record.key]
                if by_key[key].selected
                and relations[tuple(sorted((single.record.key, key)))] == "HARD_OVERLAP"
            ]
            if not any(moderate_spliced(call) for call in competitors):
                continue
            if any(
                not moderate_spliced(call)
                and (
                    _rna_evidence_rank(call) > 0
                    or _protein_evidence_rank(call) >= _protein_evidence_rank(single)
                )
                for call in competitors
            ):
                continue
            protected.add(single.record.key)

        ordered_calls = sorted(
            locus_calls,
            key=lambda call: (call.record.start, call.record.end, call.record.key),
        )
        for left, right in combinations(ordered_calls, 2):
            if not (left.selected and right.selected):
                continue
            # Empty intron chains are equal but do not constitute a
            # same-chain CDS-boundary decision in the v3.2 contract. Two
            # independently protein-supported single-exon models therefore
            # remain separate candidates even when their CDS overlaps.
            if left.record.introns == right.record.introns:
                continue
            if relations.get(tuple(sorted((left.record.key, right.record.key)))) != "HARD_OVERLAP":
                continue
            for call in (left, right):
                call.conflicts.add("OVERLAP_NO_SHARED_INTRON")
            priority = None
            for single, spliced in ((left, right), (right, left)):
                if single.record.key in protected and moderate_spliced(spliced):
                    priority = (single, spliced)
                    break
            if priority is not None:
                winner, loser = priority
                winner.reasons.add("STRONG_PROTEIN_SINGLE_PRIORITY_OVER_BALANCED_RNA")
                _mark_rejected(
                    loser, "CONFLICT_REVIEW",
                    "STRONG_PROTEIN_SINGLE_PRIORITY_OVER_BALANCED_RNA",
                )
                continue
            left_rna, right_rna = _rna_evidence_rank(left), _rna_evidence_rank(right)
            if left_rna != right_rna:
                loser = right if left_rna > right_rna else left
                _mark_rejected(loser, "CONFLICT_REVIEW", "RNA_PRIORITY_CONFLICT")
            elif left_rna > 0:
                for call in (left, right):
                    _mark_rejected(
                        call, "CONFLICT_REVIEW", "BOTH_RNA_SUPPORTED_CONFLICT"
                    )
            else:
                left_protein = _protein_evidence_rank(left)
                right_protein = _protein_evidence_rank(right)
                if left_protein != right_protein:
                    loser = right if left_protein > right_protein else left
                    _mark_rejected(
                        loser, "CONFLICT_REVIEW", "PROTEIN_PRIORITY_CONFLICT"
                    )
                else:
                    for call in (left, right):
                        _mark_rejected(
                            call, "CONFLICT_REVIEW", "UNRESOLVED_OVERLAP_CONFLICT"
                        )

    selected_keys = {key for key, call in calls.items() if call.selected}
    return selected_keys, calls


def base_audit_row(call: BaseCall) -> dict[str, object]:
    record = call.record
    return {
        "candidate_id": record.key,
        "structure_id": record.row["structure_id"],
        "seqid": record.seqid,
        "strand": record.strand,
        "start": record.start,
        "end": record.end,
        "locus_id": record.locus_id,
        "eligible": int(call.eligible),
        "selected": int(call.selected),
        "initial_status": call.initial_status,
        "decision_status": call.status,
        "routes": ",".join(call.routes) or "NA",
        "conflict_types": ",".join(sorted(call.conflicts)) or "NA",
        "reasons": ",".join(sorted(call.reasons)) or "NA",
    }
