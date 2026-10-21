"""Classified TE evidence contract shared by scoring and downstream gates."""

from urllib.parse import unquote


REPEAT_SCOPE = "classified_te_v1"

# A known TE order is sufficient: LTR/unknown is classified as LTR.
# Unknown repeats, RNA genes, simple repeats and tentative orders are excluded.
TE_ORDERS = frozenset({
    "ltr", "line", "sine", "dna", "rc", "mite", "tir", "helitron",
    "penelope", "dirs", "crypton",
})
TE_FEATURE_TYPES = {
    "ltr_retrotransposon": "LTR",
    "line_element": "LINE",
    "sine_element": "SINE",
    "terminal_inverted_repeat_element": "TIR",
    "tir_transposon": "TIR",
    "dna_transposon": "DNA",
    "helitron": "RC/Helitron",
    "mite": "MITE",
}


def is_classified_te(classification: str) -> bool:
    order = unquote(classification).strip().split("/", 1)[0].casefold()
    return order in TE_ORDERS


def repeat_classification(attributes: dict[str, str], feature_type: str) -> str:
    attributes = {key.casefold(): value for key, value in attributes.items()}
    for key in ("classification", "class"):
        if attributes.get(key):
            return unquote(attributes[key]).strip()
    # RepeatMasker library names can carry an explicit #class/family suffix.
    name = unquote(attributes.get("name", ""))
    if "#" in name:
        return name.rsplit("#", 1)[1].strip()
    return TE_FEATURE_TYPES.get(feature_type.casefold(), "Unknown")


def validate_repeat_scope(row: dict[str, str]) -> None:
    if row.get("repeat_scope") != REPEAT_SCOPE:
        raise ValueError(
            "component evidence must use repeat_scope=classified_te_v1; "
            "rerun score_component_evidence before selection/post-QC"
        )
