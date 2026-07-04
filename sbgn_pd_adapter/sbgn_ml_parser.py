"""SBGN-ML parser.

SBGN-ML is the XML serialisation of the Systems Biology Graphical Notation. This
module reads Process Description (SBGN-PD) files and returns a flat, adapter-
friendly representation. We intentionally use only the standard-library
`xml.etree.ElementTree` so the adapter has no non-BioCypher runtime dependency
outside PyYAML.

Only the SBGN-PD subset of the spec is handled. Activity Flow and Entity
Relationship diagrams parse into empty node and arc lists (rather than
crashing) so the adapter can be pointed at a mixed corpus without pre-sorting.

The parser is written to survive real-world SBGN-ML: many exporters omit optional
elements, some produce inconsistent id casing, and some inline label text where
the spec places it in a child element. Where the file is malformed rather than
merely irregular, the parser raises `SBGNMLParseError` with the source line and
element so the caller can log and skip.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from xml.etree import ElementTree as ET

SBGN_NAMESPACES = {
    # SBGN-ML 0.2 and 0.3 both use this namespace.
    "sbgn": "http://sbgn.org/libsbgn/0.2",
    "sbgn03": "http://sbgn.org/libsbgn/0.3",
}


class SBGNMLParseError(ValueError):
    """Raised when the input file is not valid SBGN-ML."""


@dataclass(frozen=True)
class SBGNGlyph:
    """A node in the SBGN-PD graph.

    `glyph_class` is the SBGN-PD class (e.g. macromolecule, simple chemical,
    complex, process, association, dissociation, source and sink, and
    perturbing agent). `label` is the human-readable string that
    domain experts recognise, which is where entity matching starts. `parent_id`
    tracks nesting for complex containment.
    """

    glyph_id: str
    glyph_class: str
    label: str
    compartment: str | None = None
    parent_id: str | None = None
    # Free-form annotations copied through from SBGN-ML extensions (Newt writes
    # UniProt and ChEBI references here). Keys are lowercased.
    annotations: dict[str, str] = field(default_factory=dict)


@dataclass(frozen=True)
class SBGNArc:
    """A typed edge in the SBGN-PD graph.

    SBGN-PD arcs are directional between source and target glyphs. We store the
    raw SBGN arc class (consumption, production, catalysis, modulation,
    stimulation, inhibition, necessary stimulation, logic arc) so the Biolink
    mapper can decide the predicate.
    """

    arc_id: str
    arc_class: str
    source_id: str
    target_id: str


@dataclass(frozen=True)
class SBGNDocument:
    """A parsed SBGN-ML file.

    Kept flat on purpose. The adapter is easier to reason about when nested
    containment is expressed via `SBGNGlyph.parent_id` rather than nested lists.
    """

    source_path: str
    glyphs: list[SBGNGlyph]
    arcs: list[SBGNArc]

    def glyph_by_id(self, glyph_id: str) -> SBGNGlyph | None:
        for glyph in self.glyphs:
            if glyph.glyph_id == glyph_id:
                return glyph
        return None


def _local_name(tag: str) -> str:
    """Return the local tag name from a namespaced ElementTree tag."""
    return tag.rsplit("}", 1)[-1]


def _child_text(elem: ET.Element, local: str) -> str | None:
    """Return the text of a direct child with the given local name, if any."""
    for child in elem:
        if _local_name(child.tag) == local:
            return child.text
    return None


def _find_child(elem: ET.Element, local: str) -> ET.Element | None:
    for child in elem:
        if _local_name(child.tag) == local:
            return child
    return None


def _extract_annotations(glyph_elem: ET.Element) -> dict[str, str]:
    """Extract external references and notes as a flat dict.

    SBGN-ML files from Reactome / Newt / VANTED all put database cross-refs
    somewhere different. We collect any child element that has an `xlink:href`
    or a `ref` attribute, plus any `<extension>` block content, into a single
    flat namespace. Downstream Biolink mapping can then decide which keys it
    trusts as identifiers.
    """
    annotations: dict[str, str] = {}
    for descendant in glyph_elem.iter():
        for attr, value in descendant.attrib.items():
            attr_local = attr.rsplit("}", 1)[-1].lower()
            if attr_local in ("href", "ref", "reference", "resource"):
                key = f"{_local_name(descendant.tag).lower()}.{attr_local}"
                annotations[key] = value
        if _local_name(descendant.tag).lower() == "notes" and descendant.text:
            annotations["notes"] = descendant.text.strip()
    return annotations


def parse_sbgn_ml(path: str | Path) -> SBGNDocument:
    """Parse an SBGN-ML file into a flat SBGNDocument.

    Raises SBGNMLParseError on structural problems. Ignores glyphs and arcs
    that fall outside the SBGN-PD subset.
    """
    path = Path(path)
    try:
        tree = ET.parse(path)
    except ET.ParseError as exc:
        raise SBGNMLParseError(f"{path}: XML parse failed: {exc}") from exc

    root = tree.getroot()

    # SBGN-ML wraps everything in <sbgn><map>...</map></sbgn>. We accept both.
    map_elem = _find_child(root, "map")
    if map_elem is None and _local_name(root.tag) == "map":
        map_elem = root
    if map_elem is None:
        raise SBGNMLParseError(f"{path}: no <map> element found")

    # SBGN-PD maps declare `language="process description"`. If the language is
    # something else we return an empty document so the adapter can move on.
    language = map_elem.get("language", "process description")
    if language.lower() not in {"process description", "process_description", "pd"}:
        return SBGNDocument(source_path=str(path), glyphs=[], arcs=[])

    glyphs = list(_iter_glyphs(map_elem))
    arcs = list(_iter_arcs(map_elem))
    return SBGNDocument(source_path=str(path), glyphs=glyphs, arcs=arcs)


def _iter_glyphs(map_elem: ET.Element, parent_id: str | None = None):
    """Yield SBGNGlyph objects from a map or a container glyph.

    Recurses into containment glyphs (SBGN complexes and compartments) so the
    parent_id field is populated for every nested glyph.
    """
    for child in map_elem:
        if _local_name(child.tag) != "glyph":
            continue
        glyph_id = child.get("id")
        glyph_class = (child.get("class") or "").lower()
        if not glyph_id or not glyph_class:
            # SBGN-ML forbids these; skip rather than raise so a partially
            # broken file still yields the well-formed glyphs.
            continue
        label = _child_text(child, "label") or ""
        if not label:
            # `<label text="..."/>` variant.
            label_elem = _find_child(child, "label")
            if label_elem is not None:
                label = label_elem.get("text", "")
        compartment = child.get("compartmentRef")
        annotations = _extract_annotations(child)
        yield SBGNGlyph(
            glyph_id=glyph_id,
            glyph_class=glyph_class,
            label=label.strip(),
            compartment=compartment,
            parent_id=parent_id,
            annotations=annotations,
        )
        # Recurse for nested glyphs.
        yield from _iter_glyphs(child, parent_id=glyph_id)


def _iter_arcs(map_elem: ET.Element):
    for child in map_elem:
        if _local_name(child.tag) != "arc":
            continue
        arc_id = child.get("id")
        arc_class = (child.get("class") or "").lower()
        source_id = child.get("source")
        target_id = child.get("target")
        if not (arc_id and arc_class and source_id and target_id):
            continue
        yield SBGNArc(
            arc_id=arc_id,
            arc_class=arc_class,
            source_id=source_id,
            target_id=target_id,
        )
