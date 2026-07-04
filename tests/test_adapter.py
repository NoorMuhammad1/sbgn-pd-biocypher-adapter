"""Tests for the SBGN-PD BioCypher adapter."""

from __future__ import annotations

from pathlib import Path

from sbgn_pd_adapter import SBGNPDAdapter, SBGNPDArcType, SBGNPDGlyphClass


def test_load_reads_all_documents(both_glycolysis_paths: list[Path]) -> None:
    adapter = SBGNPDAdapter(both_glycolysis_paths, matcher_threshold=0.7)
    adapter.load()
    assert adapter.stats.documents_parsed == 2


def test_atp_and_g6p_merge_across_files(both_glycolysis_paths: list[Path]) -> None:
    """Both files carry ATP (same ChEBI) and G6P (same ChEBI). Threshold 0.7
    should be enough for the matcher to collapse each pair."""
    adapter = SBGNPDAdapter(both_glycolysis_paths, matcher_threshold=0.7)
    adapter.load()
    # ATP appears in both files (atp_upper, atp_lower) and should merge.
    assert adapter._id_map["atp_lower"] == adapter._id_map["atp_upper"]
    # G6P appears as `g6p` and `g6p_b` and should merge.
    assert adapter._id_map["g6p"] == adapter._id_map["g6p_b"]


def test_high_threshold_prevents_merging(both_glycolysis_paths: list[Path]) -> None:
    """Threshold 0.99 should be too strict to merge the differently-labelled G6P."""
    adapter = SBGNPDAdapter(both_glycolysis_paths, matcher_threshold=0.99)
    adapter.load()
    # G6P has different label styles in the two files; at 0.99 they should not merge.
    assert adapter._id_map["g6p"] != adapter._id_map["g6p_b"]


def test_get_nodes_emits_biolink_categories(both_glycolysis_paths: list[Path]) -> None:
    adapter = SBGNPDAdapter(both_glycolysis_paths, matcher_threshold=0.7)
    nodes = list(adapter.get_nodes())
    labels = {label for _, label, _ in nodes}
    assert "biolink:SmallMolecule" in labels
    assert "biolink:MacromolecularMachineMixin" in labels
    assert "biolink:BiologicalProcess" in labels
    assert "biolink:CellularComponent" in labels


def test_get_edges_emits_biolink_predicates(both_glycolysis_paths: list[Path]) -> None:
    adapter = SBGNPDAdapter(both_glycolysis_paths, matcher_threshold=0.7)
    list(adapter.get_nodes())  # BioCypher expects nodes before edges
    edges = list(adapter.get_edges())
    predicates = {label for _, _, _, label, _ in edges}
    assert "biolink:has_input" in predicates
    assert "biolink:has_output" in predicates
    assert "biolink:catalyzes" in predicates


def test_edges_reference_canonical_ids(both_glycolysis_paths: list[Path]) -> None:
    adapter = SBGNPDAdapter(both_glycolysis_paths, matcher_threshold=0.7)
    list(adapter.get_nodes())
    edges = list(adapter.get_edges())
    # Every edge should point to a canonical id in the id map.
    valid_ids = set(adapter._id_map.values())
    for _, src, tgt, _, _ in edges:
        assert src in valid_ids
        assert tgt in valid_ids


def test_node_type_filter(both_glycolysis_paths: list[Path]) -> None:
    """Restricting to compartments should emit only compartments."""
    adapter = SBGNPDAdapter(
        both_glycolysis_paths,
        matcher_threshold=0.7,
        node_types=[SBGNPDGlyphClass.COMPARTMENT],
    )
    labels = {label for _, label, _ in adapter.get_nodes()}
    assert labels == {"biolink:CellularComponent"}


def test_edge_type_filter(both_glycolysis_paths: list[Path]) -> None:
    """Restricting to catalysis should emit only catalyzes edges."""
    adapter = SBGNPDAdapter(
        both_glycolysis_paths,
        matcher_threshold=0.7,
        edge_types=[SBGNPDArcType.CATALYSIS],
    )
    list(adapter.get_nodes())
    predicates = {label for _, _, _, label, _ in adapter.get_edges()}
    assert predicates == {"biolink:catalyzes"}


def test_no_dangling_arcs_after_merge(both_glycolysis_paths: list[Path]) -> None:
    """Arcs whose endpoints collapse to the same canonical id must be dropped."""
    adapter = SBGNPDAdapter(both_glycolysis_paths, matcher_threshold=0.7)
    list(adapter.get_nodes())
    edges = list(adapter.get_edges())
    for _, src, tgt, _, _ in edges:
        assert src != tgt, "self-loop after merge should have been dropped"


def test_report_is_stable(both_glycolysis_paths: list[Path]) -> None:
    adapter = SBGNPDAdapter(both_glycolysis_paths, matcher_threshold=0.7)
    list(adapter.get_nodes())
    list(adapter.get_edges())
    report = adapter.report()
    assert "documents parsed:        2" in report
    assert "glyphs merged:" in report
