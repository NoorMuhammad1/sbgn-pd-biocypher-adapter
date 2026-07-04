# SBGN-PD to Biolink mapping reference

The tables below are authoritative in `sbgn_pd_adapter/biolink_mapper.py`. This file is a human-readable extract for review.

## Glyph classes

| SBGN-PD class          | Biolink category                     | Notes |
|------------------------|--------------------------------------|-------|
| macromolecule          | biolink:MacromolecularMachineMixin   | Downstream cross-ref resolution refines to Protein / RNAProduct. |
| macromolecule multimer | biolink:MacromolecularComplex        | Multimer = homogeneous complex; add stoichiometry qualifier later. |
| simple chemical        | biolink:SmallMolecule                | Preferred_id = ChEBI. |
| simple chemical multimer | biolink:ChemicalMixture            | |
| nucleic acid feature   | biolink:NucleicAcidEntity            | Downstream refinement to Gene / Transcript. |
| complex                | biolink:MacromolecularComplex        | Direct analog. |
| complex multimer       | biolink:MacromolecularComplex        | With stoichiometry annotation. |
| unspecified entity     | biolink:BiologicalEntity             | SBGN wildcard. |
| perturbing agent       | biolink:ChemicalOrDrugOrTreatment    | Drug / stimulus superclass. |
| source and sink        | biolink:BiologicalEntity             | Open-boundary convention. |
| process                | biolink:BiologicalProcess            | Reified. See design.md. |
| omitted process        | biolink:BiologicalProcess            | `omitted` flag becomes a node property. |
| uncertain process      | biolink:BiologicalProcess            | Uncertainty as a property. |
| association            | biolink:BiologicalProcess            | With subclass qualifier. |
| dissociation           | biolink:BiologicalProcess            | With subclass qualifier. |
| phenotype              | biolink:PhenotypicFeature            | |
| compartment            | biolink:CellularComponent            | Container, not participant. |
| and / or / not         | biolink:BiologicalProcess            | Logical operator preserved as node so regulation logic stays queryable. |

## Arc classes

| SBGN-PD arc class      | Biolink predicate               | Direction                    |
|------------------------|---------------------------------|------------------------------|
| consumption            | biolink:has_input               | entity -> process            |
| production             | biolink:has_output              | process -> entity            |
| catalysis              | biolink:catalyzes               | enzyme -> process            |
| modulation             | biolink:regulates               | entity -> process            |
| stimulation            | biolink:positively_regulates    | entity -> process            |
| necessary stimulation  | biolink:positively_regulates    | with `necessary=true` prop   |
| inhibition             | biolink:negatively_regulates    | entity -> process            |
| logic arc              | biolink:regulates               | operator -> process          |
| equivalence arc        | biolink:same_as                 | entity -> entity             |

## Open questions

- **Reactome-specific SBO annotations.** Reactome SBGN-PD exports frequently include `sboTerm` refs on the process glyph (e.g. `SBO:0000176` for biochemical reaction). The adapter currently stores these under the `annotations` property but does not use them for refinement. A future extension could route each SBO subclass to a more specific Biolink subcategory (`MolecularActivity`, `PhysicalProcess`, etc.).
- **UniProt vs Ensembl for macromolecules.** The adapter maps `macromolecule` to `MacromolecularMachineMixin` and lets the downstream Biolink schema decide `Protein` vs `RNAProduct` at query time. An alternative is to inspect the annotation namespace and refine at ingestion. That is faster to query but couples the adapter to identifier namespaces, which reduces portability.
- **Complex participation.** SBGN nested complexes are represented via `parent_id` on child glyphs. There is no explicit Biolink predicate for `part_of_complex`; the adapter currently only records `parent_id` as a property, so queries have to join manually. Adding a `has_part` edge on load would make this searchable at the cost of doubling the edge count for complex-heavy pathways.
