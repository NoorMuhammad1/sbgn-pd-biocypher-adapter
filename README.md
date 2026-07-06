# SBGN-PD BioCypher adapter

A [BioCypher](https://biocypher.org) adapter for **SBGN Process Description** pathway data. Turns SBGN-ML files (as exported by [Newt](https://newteditor.org), Reactome, VANTED, or any other SBGN-PD tool) into a [Biolink](https://biolink.github.io/biolink-model/)-aligned biomedical knowledge graph, with tunable entity matching across sources.

Follows the [biocypher/project-template](https://github.com/biocypher/project-template) adapter conventions and the [Biolink model](https://biolink.github.io/biolink-model/) at the schema layer.

## What this adapter does

1. **Parses SBGN-ML** into a flat, adapter-friendly representation (`sbgn_ml_parser.py`).
2. **Maps SBGN-PD glyphs and arcs to Biolink categories and predicates** (`biolink_mapper.py`). The mapping is explicit and documented row-by-row so it is easy to extend and to argue about. The adapter accepts the seven SBGN-PD Level 1 flow-and-modulation arc classes (`consumption`, `production`, `catalysis`, `modulation`, `stimulation`, `necessary stimulation`, `inhibition`) plus `logic arc` and `equivalence arc`. These collapse onto six distinct Biolink predicates (`has_input`, `has_output`, `catalyzes`, `regulates`, `positively_regulates`, `negatively_regulates`), with `necessary stimulation` distinguished from plain `stimulation` by a `necessary=true` edge property rather than a separate predicate. See `docs/sbgn_to_biolink_mapping.md` for the row-by-row table.
3. **Reifies process glyphs as `biolink:BiologicalProcess` nodes** so SBGN's effectively n-ary process semantics survive translation to Biolink's predominantly binary predicate model. Participants attach via directional predicates (`has_input`, `has_output`, `catalyzes`, `regulates`, `positively_regulates`, `negatively_regulates`) set by the arc class in the source SBGN-ML. See `docs/design.md#process-reification` for the argument against the pairwise-expansion alternative, and `sketches/reification_vs_missingness/NOTE.md` for the carve-out of which SBGN-PD process-glyph subtypes are genuinely n-ary versus fixed-arity or unary.
4. **Matches entities across sources** with a similarity-threshold rule (`entity_matcher.py`) so the same molecule appearing under different labels in two SBGN-ML files does not become two knowledge-graph nodes. Every merge is recorded with a scored explanation for auditability.
5. **Yields BioCypher-shaped node and edge tuples** via `SBGNPDAdapter.get_nodes()` and `SBGNPDAdapter.get_edges()`. The BioCypher core then routes those tuples to whichever backend the runtime config selects. The default demo pipeline writes Neo4j admin-import CSVs, but the same adapter output feeds BioCypher's PyG, RDF, or plain-CSV writers without any adapter-side changes.

## Repository layout

```
sbgn_pd_biocypher_adapter/
├── sbgn_pd_adapter/
│   ├── adapter.py            main SBGNPDAdapter class
│   ├── sbgn_ml_parser.py     SBGN-ML XML parser, standard library only
│   ├── biolink_mapper.py     SBGN-PD -> Biolink mapping tables with rationale
│   └── entity_matcher.py     similarity-threshold matcher with contract checks
├── config/
│   ├── schema_config.yaml    Biolink schema, matched to what the adapter emits
│   └── biocypher_config.yaml runtime config (offline mode by default)
├── data/
│   ├── glycolysis_upper.sbgn hexokinase phosphorylation of glucose
│   └── glycolysis_lower.sbgn G6P->F6P->F1,6BP with GPI and PFK1
├── docs/
│   ├── design.md
│   ├── sbgn_to_biolink_mapping.md
│   └── example_cypher_queries.md
├── tests/                    34 pytest tests
├── create_knowledge_graph.py end-to-end pipeline
├── pyproject.toml
├── LICENSE                   MIT
└── README.md
```

## Quickstart

```bash
# Install (uv recommended; pip works too)
uv sync

# End-to-end: parse, match, write Neo4j admin-import CSVs
uv run python create_knowledge_graph.py

# Smoke test the adapter without BioCypher installed
uv run python create_knowledge_graph.py --no-biocypher

# Tune the entity-matcher threshold (default 0.7, high-precision)
uv run python create_knowledge_graph.py --threshold 0.85

# Run the test suite
uv run pytest
```

The BioCypher output lands in `biocypher-out/`, along with a Neo4j admin-import script (`neo4j-admin-import-call.ps1` or `.sh` depending on OS).

## Sample run

With the two glycolysis fragments in `data/`, the default threshold (0.7) collapses the two ATP glyphs (matched on shared ChEBI reference `CHEBI:15422`) and the two glucose-6-phosphate glyphs (`CHEBI:17665`) despite different label styles. Output:

```
SBGN-PD adapter report
  documents parsed:        2
  glyphs read:             16
  arcs read:               12
  glyphs emitted:          14
  arcs emitted:            12
  glyphs merged:           2
```

BioCypher then writes:

- `SmallMolecule-part001.csv` (6 rows)
- `MacromolecularMachineMixin-part001.csv` (3 rows for HK2, GPI, PFK1)
- `BiologicalProcess-part001.csv` (3 reified processes)
- `CellularComponent-part001.csv` (2 cytosol compartments. Note these do *not* merge because they carry no cross-refs.)
- `HasInput`, `HasOutput`, `Catalyzes`, and `Regulates` (when modulation arcs are present) edge CSVs

Swapping Neo4j for a different BioCypher output backend is a `biocypher_config.yaml` change, not an adapter change. The neo4j-admin-import script is regenerated on every run. Point it at your Neo4j instance to load.

## The entity-matcher, in detail

The matcher is the piece where this adapter goes beyond a straight one-file-per-KG translation. It is deliberately simple:

```python
Priority = w_label * label_similarity(l, r)
         + w_annotation * annotation_overlap(l, r)
         + w_compartment * compartment_agreement(l, r)
```

Weights sum to 1. Default threshold 0.7. Every candidate pair is scored and the composite is compared to the threshold. See `docs/design.md` for the reasoning behind the default weights and for how the matcher relates to BioCypher's own `preferred_id` deduplication (the two are complementary, not redundant).

The matcher enforces contract-style invariants:

- No self-loops (a glyph does not merge with itself).
- No cross-class merges (a `macromolecule` never collapses with a `simple chemical`, even if labels match).
- Traversal budget cap (`max_pairs_per_glyph`, default 32) so pathological inputs do not blow up runtime.
- Every merge decision is recorded on `matcher.decisions` with its score breakdown, so every merge is auditable end-to-end.

## Extending

- **New SBGN-PD glyph class**: add a row to `GLYPH_TO_BIOLINK` in `biolink_mapper.py` and a corresponding entry in `config/schema_config.yaml`. Add the enum value to `SBGNPDGlyphClass` in `adapter.py`.
- **New arc type**: same pattern with `ARC_TO_BIOLINK` and `SBGNPDArcType`.
- **Different Biolink target**. Change the `category` field in the mapping table. The schema config's `is_a` chain does the rest.
- **New data source (e.g. Reactome bulk import)**: drop SBGN-ML files under `data/` and re-run. The adapter recurses so subdirectories are fine.

## Design decisions worth reading

- Process reification vs pairwise expansion. See `docs/design.md#process-reification`.
- Threshold defaults and the label / annotation / compartment weight split. See `docs/design.md#entity-matcher-defaults`.
- Why the parser stays in the standard library. See `docs/design.md#no-non-standard-library-parser`.

## Experimental follow-ups

Three exploratory sketches live under `sketches/`. They are not required to run the adapter and are excluded from the core test suite. They exist to show what a downstream user might build on top.

- [`sketches/rgcn_link_prediction/`](sketches/rgcn_link_prediction/) contains an R-GCN link-prediction pipeline over the adapter's output. Includes a frozen [`PREREGISTRATION.md`](sketches/rgcn_link_prediction/PREREGISTRATION.md) for a typed-versus-flat relation comparison at Reactome scale, plus two pilot runs at synthetic scale ([`PILOT_RESULTS.md`](sketches/rgcn_link_prediction/PILOT_RESULTS.md), [`PILOT_V2_RESULTS.md`](sketches/rgcn_link_prediction/PILOT_V2_RESULTS.md)).
- [`sketches/topological_followup/DESIGN_NOTE.md`](sketches/topological_followup/DESIGN_NOTE.md) is a design memo on directed topological descriptors (persistent path homology after Chowdhury and Mémoli 2018, directed flag complexes after Lütgehetmann et al. 2020) as a complementary view to R-GCN. Working memo, not a preregistration.
- [`sketches/reification_vs_missingness/NOTE.md`](sketches/reification_vs_missingness/NOTE.md) disentangles the schema-level lossy-compression fix (reification, this adapter) from the embedding-space geometry fix (soft manifolds after Marinoni et al. 2026 TPAMI), and names four testable predictions about how the two axes compose.

## Provenance

The similarity-threshold entity-matching mechanism is adapted from my MSc thesis at Bilkent University (defended March 2026, advisor Prof. Ugur Dogrusoz, iVis Information Visualization group), where I built a Neo4j-backed graph database of biological pathways with the same reuse-versus-create decision structure but for a single application. This adapter generalises that mechanism into a reusable BioCypher-facing module.

## License

MIT. See `LICENSE`.
