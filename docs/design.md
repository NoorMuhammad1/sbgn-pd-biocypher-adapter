# Design notes

The choices below are documented here rather than in code comments so they are visible without needing to open any single module.

## Process reification

SBGN-PD is inherently n-ary. A process glyph has multiple consumption inputs, multiple production outputs, and one or more modifiers (catalysis, stimulation, inhibition). The Biolink model is predominantly binary: subject-predicate-object.

Two options:

- **Reify the process as a node** (`biolink:BiologicalProcess`) and connect participants with directional predicates (`has_input`, `has_output`, `catalyzes`, `regulates`).
- **Expand each SBGN process into all pairwise substrate-product edges** and drop the process node.

The first option is what this adapter does. It preserves process identity, keeps SBGN's n-ary semantics intact, and lets downstream users attach process-level annotations (rate constants, kinetic laws, etc.) later without touching the schema. The second option loses the process identity and needs edge groupings that BioCypher's schema layer would have to unwind at query time.

The cost of reification is one extra hop in Cypher path queries. That is worth paying.

## Entity-matcher defaults

The default weights split (label 0.4, annotation 0.5, compartment 0.1) is deliberate:

- **Annotation overlap is weighted highest** because a shared UniProt or ChEBI reference is a much stronger equivalence signal than a matching label. Two glyphs labelled `ATP` in different files could plausibly be the same molecule or two different processes coincidentally sharing a display string. Two glyphs both annotating `CHEBI:15422` almost never are.
- **Label similarity is second** because in practice many SBGN-PD files omit cross-references entirely and fall back to display labels.
- **Compartment agreement is a tiebreaker** rather than a primary signal because SBGN compartment names are not standardised (`cytosol`, `cytoplasm`, `cyto`) and users often leave them off.

The default composite threshold of 0.7 is calibrated to a **high-precision** regime: the matcher would rather miss a valid merge than create a spurious one. For high-recall use cases (e.g. deduplicating a corpus of SBGN-ML files from a single well-curated source), 0.5 is a reasonable starting point.

## Relationship to BioCypher's built-in deduplication

BioCypher's schema-config-level `preferred_id` deduplication handles the common case where two source rows already carry the same normalised identifier. This adapter's matcher handles the harder case where two SBGN-ML files describe the same molecule but were exported by different tools and never harmonised their identifiers. In that setting BioCypher would have no way to know the two rows are the same, because their ids differ.

The two mechanisms are complementary. BioCypher runs after the adapter yields, so any canonicalisation the matcher performs is visible to the schema-level dedup layer.

## No non-standard-library parser

`sbgn_ml_parser.py` uses only `xml.etree.ElementTree` from the standard library rather than `libsbgn-python`. Reasons:

- One less dependency the adapter has to pin against.
- The parser only touches the SBGN-PD subset it needs; a full libsbgn dependency would bring in Activity Flow and Entity Relationship code the adapter never runs.
- The custom parser is easy to point at slightly-malformed SBGN-ML in the wild (missing labels, casing inconsistencies) without having to catch libsbgn-specific exceptions.

The trade-off is that this parser is not a general-purpose SBGN-ML validator. For strict validation, `libsbgn-python` on the same file is the right check to run in addition.

## Contract-style correctness checks

The invariants the matcher enforces (`MatcherContracts`) mirror the correctness contracts that were part of my MSc thesis pipeline:

- **Referential integrity** — after merging, every arc endpoint must resolve to a surviving canonical id. The adapter drops arcs whose endpoints did not survive filtering.
- **No self-loops** — a glyph cannot merge with itself; after merging, arcs whose endpoints collapse to the same canonical id are dropped.
- **Duplicate suppression** — the same pair is never scored twice.
- **Traversal budget soundness** — `max_pairs_per_glyph` caps the per-glyph fan-out so a pathological input (e.g. a corpus with 1000 files each carrying the same 100 ubiquitous cofactors) does not turn the matcher into an O(n^2) sinkhole.

These contracts are on by default and can be toggled off (for benchmarking) via the `MatcherContracts` dataclass.

## Node id strategy

The adapter uses the SBGN glyph id verbatim as the BioCypher node id. This is intentional:

- Glyph ids in real SBGN-ML files are typically already stable across exports from the same tool (Newt, Reactome, etc.).
- Where two files use conflicting glyph ids, they refer to different glyphs anyway; the matcher decides whether to collapse them.
- Using the glyph id preserves round-trip fidelity: a Cypher query can print an id that a Newt user can search for in the source pathway.

The trade-off is that glyph ids can be non-descriptive (`glyph1`, `sbgn12345`). Where this matters for downstream analysis, the adapter attaches the human-readable label as a node property, so queries can filter by `label` and print the id.
