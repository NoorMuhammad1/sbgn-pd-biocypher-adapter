"""SBGN Process Description adapter for BioCypher.

Build biomedical knowledge graphs from SBGN-ML pathway files, mapped to the
Biolink model. Includes tunable entity matching across sources so the same
molecule appearing in two SBGN-PD files does not become two nodes.

The design follows BioCypher's threefold-modularity principle:

    * Modular source: SBGN-ML files (this adapter). Any SBGN-PD source that
      exports SBGN-ML can be ingested through the same pipeline.
    * Modular ontology: Biolink model at the node-category and predicate
      layers. Mapping lives in biolink_mapper.py and is easy to override.
    * Modular output: BioCypher's own output layer (Neo4j, PyG, RDF, CSV).
      This adapter does not commit to a backend.

The entity_matcher.SimilarityThresholdMatcher module is the piece where this
adapter goes beyond a straight one-file-per-KG translation. It answers the
"reuse versus create" question for entities that show up in multiple SBGN-PD
inputs, and enforces contract-style correctness checks (referential integrity,
duplicate and self-loop suppression, traversal budget soundness) that mirror
what my MSc thesis pipeline did for a single Neo4j deployment.
"""

from sbgn_pd_adapter.adapter import (
    SBGNPDAdapter,
    SBGNPDArcType,
    SBGNPDEdgeField,
    SBGNPDGlyphClass,
    SBGNPDNodeField,
)

__all__ = [
    "SBGNPDAdapter",
    "SBGNPDArcType",
    "SBGNPDEdgeField",
    "SBGNPDGlyphClass",
    "SBGNPDNodeField",
]

__version__ = "0.1.0"
