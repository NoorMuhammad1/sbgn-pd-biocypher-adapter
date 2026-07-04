"""Tests for the SBGN-PD to Biolink mapping tables."""

from __future__ import annotations

from sbgn_pd_adapter.biolink_mapper import (
    ARC_TO_BIOLINK,
    GLYPH_TO_BIOLINK,
    category_for_glyph,
    is_logical_operator,
    is_process_glyph,
    predicate_for_arc,
)


def test_every_glyph_mapping_has_rationale() -> None:
    # If someone adds a new mapping without documenting why, this test fails.
    for glyph_class, mapping in GLYPH_TO_BIOLINK.items():
        assert mapping.category, f"missing category for {glyph_class}"
        assert mapping.rationale.strip(), f"missing rationale for {glyph_class}"


def test_every_arc_mapping_has_rationale() -> None:
    for arc_class, mapping in ARC_TO_BIOLINK.items():
        assert mapping.predicate, f"missing predicate for {arc_class}"
        assert mapping.rationale.strip(), f"missing rationale for {arc_class}"


def test_biolink_categories_are_prefixed() -> None:
    for _, mapping in GLYPH_TO_BIOLINK.items():
        assert mapping.category is not None
        assert mapping.category.startswith("biolink:")


def test_biolink_predicates_are_prefixed() -> None:
    for _, mapping in ARC_TO_BIOLINK.items():
        assert mapping.predicate is not None
        assert mapping.predicate.startswith("biolink:")


def test_process_reification() -> None:
    for cls in ("process", "association", "dissociation", "omitted process"):
        assert is_process_glyph(cls) is True
    for cls in ("macromolecule", "simple chemical", "compartment"):
        assert is_process_glyph(cls) is False


def test_logical_operator_flags() -> None:
    assert is_logical_operator("and") is True
    assert is_logical_operator("or") is True
    assert is_logical_operator("not") is True
    assert is_logical_operator("process") is False


def test_case_insensitive_lookup() -> None:
    # Callers should not have to normalise case.
    assert category_for_glyph("Macromolecule") == category_for_glyph("macromolecule")
    assert predicate_for_arc("Consumption") == predicate_for_arc("consumption")


def test_unknown_class_returns_none() -> None:
    assert category_for_glyph("nonsense_class") is None
    assert predicate_for_arc("nonexistent arc") is None
