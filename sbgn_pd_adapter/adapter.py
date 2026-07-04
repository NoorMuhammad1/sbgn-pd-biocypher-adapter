"""SBGN-PD BioCypher adapter.

Follows the BioCypher project-template Adapter conventions. Exposes
`get_nodes()` and `get_edges()` generators that yield the tuples the BioCypher
core expects:

    node tuple: (node_id, label, properties)
    edge tuple: (edge_id, source_id, target_id, label, properties)

The adapter can be pointed at one or many SBGN-ML files. When more than one is
given, the SimilarityThresholdMatcher decides which glyphs across files should
collapse to a single knowledge-graph node. The matcher's threshold is exposed
as a constructor argument so pipeline authors can tune reuse-vs-create without
touching the code.

The labels emitted here are the ones declared in `config/schema_config.yaml`.
Downstream users can override them by editing the schema config; BioCypher will
handle the relabelling and category expansion via the Biolink model.
"""

from __future__ import annotations

import logging
from collections.abc import Iterator
from dataclasses import dataclass
from enum import Enum, auto
from pathlib import Path

from sbgn_pd_adapter.biolink_mapper import (
    category_for_glyph,
    is_process_glyph,
    predicate_for_arc,
)
from sbgn_pd_adapter.entity_matcher import (
    MatcherContracts,
    SimilarityThresholdMatcher,
)
from sbgn_pd_adapter.sbgn_ml_parser import SBGNArc, SBGNDocument, SBGNGlyph, parse_sbgn_ml

logger = logging.getLogger(__name__)


class SBGNPDNodeField(Enum):
    """Node property fields emitted by this adapter.

    Values match the keys in schema_config.yaml. Add here first, then update
    the schema config to expose the new field to Neo4j.
    """

    ID = "id"
    LABEL = "label"
    GLYPH_CLASS = "sbgn_glyph_class"
    COMPARTMENT = "compartment"
    PARENT_ID = "parent_id"
    SOURCE_FILE = "source_file"
    ANNOTATIONS = "annotations"
    MERGED_FROM = "merged_from"


class SBGNPDEdgeField(Enum):
    """Edge property fields emitted by this adapter."""

    ID = "id"
    ARC_CLASS = "sbgn_arc_class"
    SOURCE_FILE = "source_file"
    NECESSARY = "necessary"


class SBGNPDGlyphClass(Enum):
    """The subset of SBGN-PD glyph classes we accept.

    Bounded so a malformed file listing `class=foo` does not create a bogus
    node type in the knowledge graph. Additions here should also be added to
    biolink_mapper.GLYPH_TO_BIOLINK.
    """

    MACROMOLECULE = auto()
    MACROMOLECULE_MULTIMER = auto()
    SIMPLE_CHEMICAL = auto()
    SIMPLE_CHEMICAL_MULTIMER = auto()
    NUCLEIC_ACID_FEATURE = auto()
    COMPLEX = auto()
    COMPLEX_MULTIMER = auto()
    UNSPECIFIED_ENTITY = auto()
    PERTURBING_AGENT = auto()
    SOURCE_AND_SINK = auto()
    PROCESS = auto()
    OMITTED_PROCESS = auto()
    UNCERTAIN_PROCESS = auto()
    ASSOCIATION = auto()
    DISSOCIATION = auto()
    PHENOTYPE = auto()
    COMPARTMENT = auto()
    AND = auto()
    OR = auto()
    NOT = auto()

    @classmethod
    def from_string(cls, name: str) -> "SBGNPDGlyphClass | None":
        canonical = name.strip().lower().replace(" ", "_")
        return cls.__members__.get(canonical.upper())


class SBGNPDArcType(Enum):
    """Accepted SBGN-PD arc classes."""

    CONSUMPTION = auto()
    PRODUCTION = auto()
    CATALYSIS = auto()
    MODULATION = auto()
    STIMULATION = auto()
    NECESSARY_STIMULATION = auto()
    INHIBITION = auto()
    LOGIC_ARC = auto()
    EQUIVALENCE_ARC = auto()

    @classmethod
    def from_string(cls, name: str) -> "SBGNPDArcType | None":
        canonical = name.strip().lower().replace(" ", "_")
        return cls.__members__.get(canonical.upper())


@dataclass
class AdapterStats:
    """Diagnostics for pipeline reporting."""

    documents_parsed: int = 0
    glyphs_read: int = 0
    arcs_read: int = 0
    glyphs_emitted: int = 0
    arcs_emitted: int = 0
    glyphs_merged: int = 0
    glyphs_skipped_unknown_class: int = 0
    arcs_skipped_unknown_class: int = 0
    arcs_skipped_dangling: int = 0


class SBGNPDAdapter:
    """BioCypher adapter over one or many SBGN-PD (SBGN-ML) files.

    Args:
        source_paths: iterable of paths to SBGN-ML files. May be a single path.
        matcher_threshold: composite-score cutoff for merging glyphs across
            files. See entity_matcher.SimilarityThresholdMatcher. Defaults to
            0.7 (high-precision).
        node_types: allow-list of SBGNPDGlyphClass values. If None, all classes
            declared in the enum are emitted.
        edge_types: allow-list of SBGNPDArcType values.
        node_fields: allow-list of SBGNPDNodeField values.
        edge_fields: allow-list of SBGNPDEdgeField values.

    The BioCypher core calls `get_nodes()` first, then `get_edges()`. We keep
    the parsed documents and the entity-matching decisions in memory between
    the two calls so the edge generator can rewrite source/target ids to the
    canonical merged ids.
    """

    def __init__(
        self,
        source_paths: str | Path | list[str | Path],
        *,
        matcher_threshold: float = 0.7,
        matcher_contracts: MatcherContracts | None = None,
        node_types: list[SBGNPDGlyphClass] | None = None,
        edge_types: list[SBGNPDArcType] | None = None,
        node_fields: list[SBGNPDNodeField] | None = None,
        edge_fields: list[SBGNPDEdgeField] | None = None,
    ) -> None:
        if isinstance(source_paths, (str, Path)):
            source_paths = [source_paths]
        self.source_paths: list[Path] = [Path(p) for p in source_paths]

        self.node_types = set(node_types) if node_types else set(SBGNPDGlyphClass)
        self.edge_types = set(edge_types) if edge_types else set(SBGNPDArcType)
        self.node_fields = set(node_fields) if node_fields else set(SBGNPDNodeField)
        self.edge_fields = set(edge_fields) if edge_fields else set(SBGNPDEdgeField)

        self.matcher = SimilarityThresholdMatcher(
            threshold=matcher_threshold,
            contracts=matcher_contracts or MatcherContracts(),
        )
        self.stats = AdapterStats()
        self._documents: list[SBGNDocument] = []
        self._all_glyphs: list[SBGNGlyph] = []
        self._id_map: dict[str, str] = {}
        self._loaded = False

    # ------------------------------------------------------------------ core
    def load(self) -> None:
        """Parse all source files and run entity matching.

        Called automatically on the first `get_nodes()` invocation. Exposed so
        pipeline authors can eager-load and inspect the matcher's decisions
        (via `adapter.matcher.decisions`) before iterating.
        """
        if self._loaded:
            return
        for path in self.source_paths:
            doc = parse_sbgn_ml(path)
            self._documents.append(doc)
            self._all_glyphs.extend(doc.glyphs)
            self.stats.documents_parsed += 1
            self.stats.glyphs_read += len(doc.glyphs)
            self.stats.arcs_read += len(doc.arcs)
        self._id_map = self.matcher.group(self._all_glyphs)
        # Count merges (a glyph is "merged" when its canonical id is not
        # itself).
        self.stats.glyphs_merged = sum(
            1 for gid, canon in self._id_map.items() if gid != canon
        )
        self._loaded = True

    def get_nodes(self) -> Iterator[tuple[str, str, dict]]:
        """Emit canonical nodes.

        Each equivalence class of glyphs produces one node. The representative
        glyph is the first one encountered (deterministic on file order).
        Merged glyph ids are attached as `merged_from` so downstream queries
        can inspect the collapse.
        """
        self.load()
        emitted: set[str] = set()
        # Build a reverse map so we can attach merged_from information.
        canonical_to_merged: dict[str, list[str]] = {}
        for gid, canon in self._id_map.items():
            canonical_to_merged.setdefault(canon, []).append(gid)

        for doc in self._documents:
            for glyph in doc.glyphs:
                canonical_id = self._id_map.get(glyph.glyph_id, glyph.glyph_id)
                if canonical_id in emitted:
                    continue
                # Only emit for the representative glyph, not the merged ones.
                if canonical_id != glyph.glyph_id:
                    continue

                glyph_enum = SBGNPDGlyphClass.from_string(glyph.glyph_class)
                if glyph_enum is None:
                    self.stats.glyphs_skipped_unknown_class += 1
                    logger.debug("skipping unknown glyph class %s", glyph.glyph_class)
                    continue
                if glyph_enum not in self.node_types:
                    continue

                mapping = category_for_glyph(glyph.glyph_class)
                if mapping is None or mapping.category is None:
                    continue
                label = mapping.category

                properties = self._node_properties(glyph, doc, canonical_to_merged.get(canonical_id, []))
                emitted.add(canonical_id)
                self.stats.glyphs_emitted += 1
                yield (canonical_id, label, properties)

    def get_edges(self) -> Iterator[tuple[str, str, str, str, dict]]:
        """Emit canonical edges.

        Source and target ids are rewritten to the canonical (post-merge) ids.
        Arcs whose endpoints did not survive filtering (unknown glyph class,
        excluded node type) are dropped and counted in `stats.arcs_skipped_*`.
        """
        self.load()
        # Precompute the set of glyph ids that survived node emission.
        surviving_canonical_ids: set[str] = set()
        for doc in self._documents:
            for glyph in doc.glyphs:
                canonical_id = self._id_map.get(glyph.glyph_id, glyph.glyph_id)
                if canonical_id != glyph.glyph_id:
                    continue
                glyph_enum = SBGNPDGlyphClass.from_string(glyph.glyph_class)
                if glyph_enum is None or glyph_enum not in self.node_types:
                    continue
                if category_for_glyph(glyph.glyph_class) is None:
                    continue
                surviving_canonical_ids.add(canonical_id)

        for doc in self._documents:
            for arc in doc.arcs:
                arc_enum = SBGNPDArcType.from_string(arc.arc_class)
                if arc_enum is None:
                    self.stats.arcs_skipped_unknown_class += 1
                    continue
                if arc_enum not in self.edge_types:
                    continue

                mapping = predicate_for_arc(arc.arc_class)
                if mapping is None or mapping.predicate is None:
                    continue

                canonical_source = self._id_map.get(arc.source_id, arc.source_id)
                canonical_target = self._id_map.get(arc.target_id, arc.target_id)

                if canonical_source == canonical_target:
                    # Self-loop after merging: skip. This is expected when two
                    # source files list the same glyph on both endpoints of an
                    # arc.
                    self.stats.arcs_skipped_dangling += 1
                    continue
                if (
                    canonical_source not in surviving_canonical_ids
                    or canonical_target not in surviving_canonical_ids
                ):
                    self.stats.arcs_skipped_dangling += 1
                    continue

                properties = self._edge_properties(arc, doc)
                self.stats.arcs_emitted += 1
                yield (
                    arc.arc_id,
                    canonical_source,
                    canonical_target,
                    mapping.predicate,
                    properties,
                )

    # -------------------------------------------------------------- helpers
    def _node_properties(
        self,
        glyph: SBGNGlyph,
        doc: SBGNDocument,
        merged_from: list[str],
    ) -> dict[str, object]:
        properties: dict[str, object] = {}
        if SBGNPDNodeField.LABEL in self.node_fields:
            properties["label"] = glyph.label
        if SBGNPDNodeField.GLYPH_CLASS in self.node_fields:
            properties["sbgn_glyph_class"] = glyph.glyph_class
        if SBGNPDNodeField.COMPARTMENT in self.node_fields and glyph.compartment:
            properties["compartment"] = glyph.compartment
        if SBGNPDNodeField.PARENT_ID in self.node_fields and glyph.parent_id:
            properties["parent_id"] = glyph.parent_id
        if SBGNPDNodeField.SOURCE_FILE in self.node_fields:
            properties["source_file"] = doc.source_path
        if SBGNPDNodeField.ANNOTATIONS in self.node_fields and glyph.annotations:
            # BioCypher CSV mode does not like nested structures; flatten.
            properties["annotations"] = ";".join(
                f"{k}={v}" for k, v in sorted(glyph.annotations.items())
            )
        if SBGNPDNodeField.MERGED_FROM in self.node_fields and len(merged_from) > 1:
            properties["merged_from"] = ";".join(sorted(gid for gid in merged_from if gid != glyph.glyph_id))
        # Flag process reification so downstream queries can filter on it.
        if is_process_glyph(glyph.glyph_class):
            properties["is_reified_process"] = True
        return properties

    def _edge_properties(self, arc: SBGNArc, doc: SBGNDocument) -> dict[str, object]:
        properties: dict[str, object] = {}
        if SBGNPDEdgeField.ARC_CLASS in self.edge_fields:
            properties["sbgn_arc_class"] = arc.arc_class
        if SBGNPDEdgeField.SOURCE_FILE in self.edge_fields:
            properties["source_file"] = doc.source_path
        if (
            SBGNPDEdgeField.NECESSARY in self.edge_fields
            and arc.arc_class.lower() == "necessary stimulation"
        ):
            properties["necessary"] = True
        return properties

    # ---------------------------------------------------------- diagnostics
    def get_node_count(self) -> int:
        """Return the number of canonical nodes the adapter will emit.

        Enumerates the generator so runs after this call re-parse the sources.
        Use once at the end for reporting, not in the hot path.
        """
        return sum(1 for _ in self.get_nodes())

    def report(self) -> str:
        """Human-readable summary of the last load."""
        s = self.stats
        return (
            f"SBGN-PD adapter report\n"
            f"  documents parsed:        {s.documents_parsed}\n"
            f"  glyphs read:             {s.glyphs_read}\n"
            f"  arcs read:               {s.arcs_read}\n"
            f"  glyphs emitted:          {s.glyphs_emitted}\n"
            f"  arcs emitted:            {s.arcs_emitted}\n"
            f"  glyphs merged:           {s.glyphs_merged}\n"
            f"  glyphs skipped (class):  {s.glyphs_skipped_unknown_class}\n"
            f"  arcs skipped (class):    {s.arcs_skipped_unknown_class}\n"
            f"  arcs skipped (dangling): {s.arcs_skipped_dangling}\n"
        )
