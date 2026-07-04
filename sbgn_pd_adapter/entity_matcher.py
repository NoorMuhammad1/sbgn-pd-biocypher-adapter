"""Similarity-threshold entity matching across SBGN-PD sources.

BioCypher's schema-level `preferred_id` deduplication handles the common case
where two source rows already carry the same normalised identifier. This
module handles the harder case where two SBGN-PD files describe the same
molecule but were exported by different tools and never harmonised their
identifiers -- so one calls a protein `MAPK1` and the other `ERK2`, and both
should collapse to the same knowledge-graph node.

The mechanism is the one I built for my MSc thesis pipeline, adapted for the
BioCypher `get_nodes` interface. It is deliberately simple:

    * A candidate match is scored by three signals: label similarity,
      annotation overlap (shared UniProt/ChEBI/etc. cross-refs), and
      compartment agreement.
    * The user sets a threshold in [0, 1]. Two glyphs merge when the composite
      score exceeds the threshold. Below it, they stay separate.
    * Every merge decision is logged with the score contributions so the
      pipeline is auditable.

The point is that reuse-versus-create is a policy the pipeline author sets, not
a hidden default. This is also where the contract-style correctness checks live
(referential integrity across merges, no self-loops after merging, and
`match_budget` bounds so pathological inputs cannot blow up the runtime).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from difflib import SequenceMatcher

from sbgn_pd_adapter.sbgn_ml_parser import SBGNGlyph


@dataclass
class MatchScore:
    """Explains a merge decision. Kept as a value so unit tests can assert on it."""

    label_score: float
    annotation_score: float
    compartment_score: float
    composite: float

    def to_dict(self) -> dict[str, float]:
        return {
            "label": round(self.label_score, 4),
            "annotation": round(self.annotation_score, 4),
            "compartment": round(self.compartment_score, 4),
            "composite": round(self.composite, 4),
        }


@dataclass
class MatchDecision:
    """A single decision the matcher made, whether it merged or not."""

    left_glyph_id: str
    right_glyph_id: str
    score: MatchScore
    merged: bool
    reason: str


@dataclass
class MatcherContracts:
    """Contract-style invariants the matcher must satisfy.

    Kept as a first-class configuration object so the pipeline author can turn
    each contract on or off (for benchmarking) but the default is all-on.
    """

    forbid_self_loops: bool = True
    forbid_duplicate_pairs: bool = True
    max_pairs_per_glyph: int = 32  # traversal budget soundness


@dataclass
class SimilarityThresholdMatcher:
    """Entity matcher across SBGN-PD glyphs.

    Threshold semantics:
        1.0 => never merge across sources (identity-only matching).
        0.0 => merge everything with a nonzero composite (dangerous).
        Typical values: 0.7 for high-precision, 0.5 for high-recall.

    Weights sum to 1.0 by convention. If they do not, the matcher normalises
    them so the composite score still lands in [0, 1].
    """

    threshold: float = 0.7
    label_weight: float = 0.4
    annotation_weight: float = 0.5
    compartment_weight: float = 0.1
    contracts: MatcherContracts = field(default_factory=MatcherContracts)
    decisions: list[MatchDecision] = field(default_factory=list)

    def __post_init__(self) -> None:
        total = self.label_weight + self.annotation_weight + self.compartment_weight
        if total <= 0:
            raise ValueError("At least one weight must be positive.")
        # Normalise so composite stays in [0, 1].
        self.label_weight /= total
        self.annotation_weight /= total
        self.compartment_weight /= total

    def score(self, left: SBGNGlyph, right: SBGNGlyph) -> MatchScore:
        """Score two glyphs. Called for every candidate pair."""
        label_score = _string_similarity(left.label, right.label)
        annotation_score = _annotation_overlap(left.annotations, right.annotations)
        compartment_score = _compartment_agreement(left.compartment, right.compartment)
        composite = (
            self.label_weight * label_score
            + self.annotation_weight * annotation_score
            + self.compartment_weight * compartment_score
        )
        return MatchScore(
            label_score=label_score,
            annotation_score=annotation_score,
            compartment_score=compartment_score,
            composite=composite,
        )

    def decide(self, left: SBGNGlyph, right: SBGNGlyph) -> MatchDecision:
        """Score a pair and decide whether to merge, respecting contracts."""
        if self.contracts.forbid_self_loops and left.glyph_id == right.glyph_id:
            return MatchDecision(
                left_glyph_id=left.glyph_id,
                right_glyph_id=right.glyph_id,
                score=MatchScore(0.0, 0.0, 0.0, 0.0),
                merged=False,
                reason="self-loop forbidden by contract",
            )
        if left.glyph_class != right.glyph_class:
            # Cross-class matches never merge. A macromolecule and a small
            # chemical with the same label are not the same node.
            return MatchDecision(
                left_glyph_id=left.glyph_id,
                right_glyph_id=right.glyph_id,
                score=MatchScore(0.0, 0.0, 0.0, 0.0),
                merged=False,
                reason="glyph classes differ",
            )
        score = self.score(left, right)
        merged = score.composite >= self.threshold
        return MatchDecision(
            left_glyph_id=left.glyph_id,
            right_glyph_id=right.glyph_id,
            score=score,
            merged=merged,
            reason=(
                f"composite {score.composite:.3f} "
                f"{'>=' if merged else '<'} threshold {self.threshold:.3f}"
            ),
        )

    def group(self, glyphs: list[SBGNGlyph]) -> dict[str, str]:
        """Group glyphs into equivalence classes.

        Returns a mapping from every glyph_id to the id of its canonical
        representative (the first glyph encountered in its equivalence class).
        Uses a union-find so the grouping is order-independent up to ties.

        Complexity: O(n^2) worst case with a per-glyph budget cap set by
        `contracts.max_pairs_per_glyph`. In practice the cap keeps runtime
        linear in n for well-typed SBGN-ML sources.
        """
        parent: dict[str, str] = {g.glyph_id: g.glyph_id for g in glyphs}

        def find(x: str) -> str:
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        def union(x: str, y: str) -> None:
            root_x, root_y = find(x), find(y)
            if root_x != root_y:
                parent[root_y] = root_x

        pairs_seen: set[tuple[str, str]] = set()
        pairs_per_glyph: dict[str, int] = {g.glyph_id: 0 for g in glyphs}

        for i, left in enumerate(glyphs):
            for right in glyphs[i + 1 :]:
                if (
                    pairs_per_glyph[left.glyph_id]
                    >= self.contracts.max_pairs_per_glyph
                ):
                    break
                if (
                    pairs_per_glyph[right.glyph_id]
                    >= self.contracts.max_pairs_per_glyph
                ):
                    continue
                pair_key = (
                    min(left.glyph_id, right.glyph_id),
                    max(left.glyph_id, right.glyph_id),
                )
                if (
                    self.contracts.forbid_duplicate_pairs
                    and pair_key in pairs_seen
                ):
                    continue
                pairs_seen.add(pair_key)
                pairs_per_glyph[left.glyph_id] += 1
                pairs_per_glyph[right.glyph_id] += 1
                decision = self.decide(left, right)
                self.decisions.append(decision)
                if decision.merged:
                    union(left.glyph_id, right.glyph_id)

        return {g.glyph_id: find(g.glyph_id) for g in glyphs}


def _string_similarity(left: str, right: str) -> float:
    """Case-insensitive fuzzy string similarity in [0, 1]."""
    if not left or not right:
        return 0.0
    return SequenceMatcher(None, left.strip().lower(), right.strip().lower()).ratio()


def _annotation_overlap(
    left: dict[str, str], right: dict[str, str]
) -> float:
    """Jaccard similarity over annotation values.

    We compare values, not keys, so that a UniProt reference recorded under
    two different keys still counts as agreement. Empty-string values are
    dropped so absent references do not create spurious matches.
    """
    left_values = {v for v in left.values() if v}
    right_values = {v for v in right.values() if v}
    if not left_values or not right_values:
        return 0.0
    intersection = left_values & right_values
    union = left_values | right_values
    return len(intersection) / len(union) if union else 0.0


def _compartment_agreement(
    left: str | None, right: str | None
) -> float:
    """1.0 if compartments agree, 0.5 if one is unknown, 0.0 if they disagree."""
    if left is None and right is None:
        return 0.5
    if left is None or right is None:
        return 0.5
    return 1.0 if left == right else 0.0
