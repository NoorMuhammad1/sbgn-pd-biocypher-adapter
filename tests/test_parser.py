"""Tests for the SBGN-ML parser."""

from __future__ import annotations

from pathlib import Path

import pytest

from sbgn_pd_adapter.sbgn_ml_parser import (
    SBGNMLParseError,
    parse_sbgn_ml,
)


def test_parses_upper_glycolysis(upper_glycolysis_path: Path) -> None:
    doc = parse_sbgn_ml(upper_glycolysis_path)
    assert doc.source_path == str(upper_glycolysis_path)
    # 1 compartment + 4 small chemicals + 1 macromolecule + 1 process = 7 glyphs.
    assert len(doc.glyphs) == 7
    # 5 arcs in this file.
    assert len(doc.arcs) == 5


def test_parses_lower_glycolysis(lower_glycolysis_path: Path) -> None:
    doc = parse_sbgn_ml(lower_glycolysis_path)
    # 1 compartment + 4 small chemicals + 2 macromolecules + 2 processes = 9 glyphs.
    assert len(doc.glyphs) == 9
    assert len(doc.arcs) == 7


def test_annotations_are_extracted(upper_glycolysis_path: Path) -> None:
    doc = parse_sbgn_ml(upper_glycolysis_path)
    glucose = doc.glyph_by_id("glucose")
    assert glucose is not None
    # ChEBI reference should be captured under the annotation.resource key.
    assert any("CHEBI:17925" in v for v in glucose.annotations.values())


def test_compartment_reference_captured(upper_glycolysis_path: Path) -> None:
    doc = parse_sbgn_ml(upper_glycolysis_path)
    atp = doc.glyph_by_id("atp_upper")
    assert atp is not None
    assert atp.compartment == "cytosol"


def test_arc_class_lowercased(upper_glycolysis_path: Path) -> None:
    doc = parse_sbgn_ml(upper_glycolysis_path)
    classes = {arc.arc_class for arc in doc.arcs}
    # Every class should be lowercase after parsing so callers do not need to normalise.
    for cls in classes:
        assert cls == cls.lower()


def test_bad_xml_raises(tmp_path: Path) -> None:
    bad = tmp_path / "broken.sbgn"
    bad.write_text("<sbgn><not-closed", encoding="utf-8")
    with pytest.raises(SBGNMLParseError):
        parse_sbgn_ml(bad)


def test_non_pd_language_returns_empty(tmp_path: Path) -> None:
    # A minimal Activity Flow document. The parser should skip it rather than crash.
    af = tmp_path / "activity_flow.sbgn"
    af.write_text(
        '<?xml version="1.0" encoding="UTF-8"?>'
        '<sbgn xmlns="http://sbgn.org/libsbgn/0.2">'
        '<map language="activity flow"></map></sbgn>',
        encoding="utf-8",
    )
    doc = parse_sbgn_ml(af)
    assert doc.glyphs == []
    assert doc.arcs == []
