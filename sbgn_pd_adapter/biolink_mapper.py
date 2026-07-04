"""SBGN-PD to Biolink model mapping tables.

Design decision: SBGN-PD is inherently n-ary (a process node has multiple
consumption inputs and multiple production outputs, plus modifiers). Biolink is
predominantly binary (subject-predicate-object). We reify process glyphs as
first-class nodes (Biolink `BiologicalProcess` by default) and connect
participants with directional predicates. This preserves the SBGN semantics
without inventing edge groupings that BioCypher's schema layer would then have
to unwind. The alternative -- expanding each SBGN process into all pairwise
substrate-product edges -- loses the process identity and makes it impossible
to attach process-level annotations later.

If the biological context calls for it, the process reification can be
overridden by pointing schema_config.yaml at a different Biolink category (for
example `MolecularActivity`, or a subclass like `Catalysis`).

The mapping is deliberately explicit rather than clever. Every SBGN class and
arc class is listed with a rationale. Extending it should be adding rows to
GLYPH_TO_BIOLINK and ARC_TO_BIOLINK, not writing new dispatch logic.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class BiolinkMapping:
    """A single row of the SBGN -> Biolink mapping.

    `category` is the Biolink node category (in `biolink:...` prefixed form as
    BioCypher expects). `predicate` is only used for arcs. `rationale` is a
    short human-readable note that ends up in the schema documentation.
    """

    category: str | None = None
    predicate: str | None = None
    rationale: str = ""


# SBGN-PD glyph class -> Biolink node category.
#
# For entity pool nodes we lean toward broad Biolink categories rather than
# specific ones like `Protein` because SBGN-PD does not distinguish protein
# from RNA at the glyph level (both are `macromolecule`). The Biolink refinement
# to `Protein` or `RNAProduct` happens downstream, when cross-refs to UniProt
# or Ensembl are resolved. Doing it here would need domain heuristics that a
# schema-layer decision should be responsible for.
GLYPH_TO_BIOLINK: dict[str, BiolinkMapping] = {
    "macromolecule": BiolinkMapping(
        category="biolink:MacromolecularMachineMixin",
        rationale="SBGN macromolecule covers proteins, RNAs, and multimers. Downstream "
        "cross-ref resolution refines to biolink:Protein or biolink:RNAProduct.",
    ),
    "macromolecule multimer": BiolinkMapping(
        category="biolink:MacromolecularComplex",
        rationale="Multimer glyphs are homogeneous complexes; Biolink models these as "
        "MacromolecularComplex with an integer count qualifier.",
    ),
    "simple chemical": BiolinkMapping(
        category="biolink:SmallMolecule",
        rationale="SBGN simple chemical is the small-molecule class; Biolink SmallMolecule "
        "aligns directly and carries ChEBI as preferred_id.",
    ),
    "simple chemical multimer": BiolinkMapping(
        category="biolink:ChemicalMixture",
        rationale="Multimers of small molecules are chemical mixtures in the Biolink model.",
    ),
    "nucleic acid feature": BiolinkMapping(
        category="biolink:NucleicAcidEntity",
        rationale="SBGN nucleic acid feature (gene, transcript, etc.) maps to "
        "NucleicAcidEntity, refined downstream.",
    ),
    "complex": BiolinkMapping(
        category="biolink:MacromolecularComplex",
        rationale="SBGN complex glyph groups multiple entity pool nodes; MacromolecularComplex "
        "is the direct Biolink analog.",
    ),
    "complex multimer": BiolinkMapping(
        category="biolink:MacromolecularComplex",
        rationale="Multimer of complexes; same Biolink category, with a stoichiometry "
        "annotation.",
    ),
    "unspecified entity": BiolinkMapping(
        category="biolink:BiologicalEntity",
        rationale="The SBGN wildcard entity type. Biolink BiologicalEntity is the safest broad "
        "supertype; downstream refinement is optional.",
    ),
    "perturbing agent": BiolinkMapping(
        category="biolink:ChemicalOrDrugOrTreatment",
        rationale="Perturbing agents in SBGN-PD are the drug/stimulus class; ChemicalOrDrugOr"
        "Treatment covers all three subtypes.",
    ),
    "source and sink": BiolinkMapping(
        category="biolink:BiologicalEntity",
        rationale="Source-and-sink is an SBGN convention for open-boundary reactants. Kept as "
        "BiologicalEntity so its role is preserved without needing a Biolink subtype.",
    ),
    # Process nodes -- reified to Biolink nodes so participants can attach.
    "process": BiolinkMapping(
        category="biolink:BiologicalProcess",
        rationale="SBGN generic process. Reified as a BiologicalProcess node so participants "
        "(substrate, product, modifier) can attach via directional predicates.",
    ),
    "omitted process": BiolinkMapping(
        category="biolink:BiologicalProcess",
        rationale="Same reification as process; the `omitted` flag becomes a node property.",
    ),
    "uncertain process": BiolinkMapping(
        category="biolink:BiologicalProcess",
        rationale="Same reification; uncertainty becomes a node property.",
    ),
    "association": BiolinkMapping(
        category="biolink:BiologicalProcess",
        rationale="SBGN association is a specialised process (component binding). Modelled as "
        "BiologicalProcess with a subclass qualifier so the KG stays queryable by process type.",
    ),
    "dissociation": BiolinkMapping(
        category="biolink:BiologicalProcess",
        rationale="Symmetric with association.",
    ),
    "phenotype": BiolinkMapping(
        category="biolink:PhenotypicFeature",
        rationale="SBGN phenotype is a downstream observable; Biolink PhenotypicFeature matches.",
    ),
    # Compartments -- containers, not participants.
    "compartment": BiolinkMapping(
        category="biolink:CellularComponent",
        rationale="SBGN compartments are subcellular containers; Biolink CellularComponent "
        "captures the semantics without conflating with GO:process.",
    ),
    # Logical operators. Retained as Biolink nodes so the modifier logic can be reconstructed
    # from the graph (rather than baked into edge properties, which would hide it from queries).
    "and": BiolinkMapping(
        category="biolink:BiologicalProcess",
        rationale="Logical AND connector. Preserved as a first-class node with subclass "
        "qualifier so a query can filter for AND-gated regulation.",
    ),
    "or": BiolinkMapping(
        category="biolink:BiologicalProcess",
        rationale="Logical OR connector; same treatment.",
    ),
    "not": BiolinkMapping(
        category="biolink:BiologicalProcess",
        rationale="Logical NOT connector; same treatment.",
    ),
}


# SBGN-PD arc class -> Biolink predicate.
#
# The direction convention in SBGN-PD:
#   consumption arc runs entity -> process (the entity is consumed)
#   production arc runs process -> entity (the entity is produced)
#   catalysis / modulation / stimulation / inhibition / necessary stimulation
#       all run entity -> process (the entity modifies the process)
#
# We map each to the closest Biolink predicate, keeping directionality explicit
# via the schema_config.yaml source/target order.
ARC_TO_BIOLINK: dict[str, BiolinkMapping] = {
    "consumption": BiolinkMapping(
        predicate="biolink:has_input",
        rationale="Substrate consumed by the process; biolink:has_input is the direct predicate.",
    ),
    "production": BiolinkMapping(
        predicate="biolink:has_output",
        rationale="Product generated by the process.",
    ),
    "catalysis": BiolinkMapping(
        predicate="biolink:catalyzes",
        rationale="Enzyme catalysing the process. Biolink has a first-class predicate.",
    ),
    "modulation": BiolinkMapping(
        predicate="biolink:regulates",
        rationale="Generic (up or down) modulation.",
    ),
    "stimulation": BiolinkMapping(
        predicate="biolink:positively_regulates",
        rationale="Positive modulation.",
    ),
    "necessary stimulation": BiolinkMapping(
        predicate="biolink:positively_regulates",
        rationale="Positive modulation that is required for the process. The `necessary` "
        "qualifier becomes an edge property so a query can filter on it.",
    ),
    "inhibition": BiolinkMapping(
        predicate="biolink:negatively_regulates",
        rationale="Negative modulation.",
    ),
    "logic arc": BiolinkMapping(
        predicate="biolink:regulates",
        rationale="Connects logical operators to their downstream target. Kept as regulates so "
        "the logical-operator subgraph is queryable as regulation.",
    ),
    "equivalence arc": BiolinkMapping(
        predicate="biolink:same_as",
        rationale="Marks equivalent entities across sub-maps.",
    ),
}


def category_for_glyph(glyph_class: str) -> BiolinkMapping | None:
    """Look up the Biolink node category for an SBGN-PD glyph class."""
    return GLYPH_TO_BIOLINK.get(glyph_class.lower())


def predicate_for_arc(arc_class: str) -> BiolinkMapping | None:
    """Look up the Biolink predicate for an SBGN-PD arc class."""
    return ARC_TO_BIOLINK.get(arc_class.lower())


def is_process_glyph(glyph_class: str) -> bool:
    """True if the glyph is a process-type node (needs reification)."""
    mapping = category_for_glyph(glyph_class)
    return mapping is not None and mapping.category == "biolink:BiologicalProcess"


def is_logical_operator(glyph_class: str) -> bool:
    """True if the glyph is a logical operator (AND / OR / NOT)."""
    return glyph_class.lower() in {"and", "or", "not"}
