"""Tests for the similarity-threshold entity matcher."""

from __future__ import annotations

import pytest

from sbgn_pd_adapter.entity_matcher import (
    MatcherContracts,
    SimilarityThresholdMatcher,
)
from sbgn_pd_adapter.sbgn_ml_parser import SBGNGlyph


def _make_glyph(glyph_id: str, label: str, glyph_class: str = "simple chemical", **kwargs: object) -> SBGNGlyph:
    return SBGNGlyph(
        glyph_id=glyph_id,
        glyph_class=glyph_class,
        label=label,
        compartment=kwargs.pop("compartment", None),
        parent_id=kwargs.pop("parent_id", None),
        annotations=kwargs.pop("annotations", {}),  # type: ignore[arg-type]
    )


def test_weight_normalisation() -> None:
    m = SimilarityThresholdMatcher(
        label_weight=2.0, annotation_weight=2.0, compartment_weight=1.0
    )
    assert m.label_weight == pytest.approx(0.4)
    assert m.annotation_weight == pytest.approx(0.4)
    assert m.compartment_weight == pytest.approx(0.2)


def test_zero_weights_rejected() -> None:
    with pytest.raises(ValueError):
        SimilarityThresholdMatcher(
            label_weight=0.0, annotation_weight=0.0, compartment_weight=0.0
        )


def test_identical_labels_merge_above_threshold() -> None:
    m = SimilarityThresholdMatcher(threshold=0.7)
    a = _make_glyph("a", "ATP", annotations={"ref": "CHEBI:15422"})
    b = _make_glyph("b", "ATP", annotations={"ref": "CHEBI:15422"})
    decision = m.decide(a, b)
    assert decision.merged is True
    assert decision.score.composite > 0.7


def test_different_class_never_merges() -> None:
    m = SimilarityThresholdMatcher(threshold=0.0)  # would merge anything otherwise
    a = _make_glyph("a", "ATP", glyph_class="simple chemical")
    b = _make_glyph("b", "ATP", glyph_class="macromolecule")
    decision = m.decide(a, b)
    assert decision.merged is False
    assert "class" in decision.reason


def test_self_loop_forbidden() -> None:
    m = SimilarityThresholdMatcher()
    a = _make_glyph("a", "ATP")
    decision = m.decide(a, a)
    assert decision.merged is False
    assert "self-loop" in decision.reason


def test_annotation_overlap_dominates() -> None:
    # Different labels, identical annotations. High-weight-on-annotation matcher
    # should merge, low-weight matcher should not.
    a = _make_glyph("a", "ERK2", annotations={"ref": "UniProt:P28482"})
    b = _make_glyph("b", "MAPK1", annotations={"ref": "UniProt:P28482"})

    m_high = SimilarityThresholdMatcher(
        threshold=0.5, label_weight=0.1, annotation_weight=0.9, compartment_weight=0.0
    )
    m_low = SimilarityThresholdMatcher(
        threshold=0.5, label_weight=0.9, annotation_weight=0.1, compartment_weight=0.0
    )

    assert m_high.decide(a, b).merged is True
    assert m_low.decide(a, b).merged is False


def test_group_returns_canonical_ids() -> None:
    m = SimilarityThresholdMatcher(threshold=0.7)
    glyphs = [
        _make_glyph("g1", "ATP", annotations={"ref": "CHEBI:15422"}),
        _make_glyph("g2", "ATP", annotations={"ref": "CHEBI:15422"}),
        _make_glyph("g3", "ADP", annotations={"ref": "CHEBI:16761"}),
    ]
    id_map = m.group(glyphs)
    # g1 and g2 collapse; g3 stays alone.
    assert id_map["g1"] == id_map["g2"]
    assert id_map["g3"] == "g3"


def test_traversal_budget_contract() -> None:
    # With a small pair budget the matcher should not attempt every pair.
    m = SimilarityThresholdMatcher(
        threshold=0.7,
        contracts=MatcherContracts(max_pairs_per_glyph=1),
    )
    glyphs = [_make_glyph(f"g{i}", "ATP", annotations={"ref": "CHEBI:15422"}) for i in range(5)]
    m.group(glyphs)
    # Every glyph's pairs_per_glyph budget was 1; we should have made at most
    # ceil(5/2) = 3 pair comparisons if it were symmetric, but with the budget
    # cap, at most 5 total comparisons.
    assert len(m.decisions) <= 5


def test_decisions_are_recorded() -> None:
    m = SimilarityThresholdMatcher(threshold=0.7)
    m.group([
        _make_glyph("g1", "ATP", annotations={"ref": "CHEBI:15422"}),
        _make_glyph("g2", "ADP", annotations={"ref": "CHEBI:16761"}),
    ])
    assert len(m.decisions) == 1
    assert m.decisions[0].merged is False
