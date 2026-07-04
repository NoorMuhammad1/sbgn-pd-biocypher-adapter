"""End-to-end pipeline for the SBGN-PD BioCypher adapter.

Runs the adapter on every SBGN-ML file in `data/` and writes the resulting
knowledge graph to `biocypher-out/` in Neo4j admin-import CSV format.

Usage
-----
    uv run python create_knowledge_graph.py

    # or with a custom data directory:
    uv run python create_knowledge_graph.py --data-dir path/to/sbgn-files

    # or without BioCypher, for smoke-testing the adapter in isolation:
    uv run python create_knowledge_graph.py --no-biocypher

The pipeline also prints an adapter report so you can see how many glyphs
merged across sources, which arcs were dropped for dangling endpoints, and
which glyph classes fell outside the SBGN-PD subset.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from sbgn_pd_adapter import SBGNPDAdapter

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(name)s %(levelname)s: %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("sbgn_pd_pipeline")


def find_sbgn_files(data_dir: Path) -> list[Path]:
    """Return every SBGN-ML file in a directory tree."""
    if not data_dir.exists():
        raise FileNotFoundError(data_dir)
    candidates = sorted(
        list(data_dir.rglob("*.sbgn"))
        + list(data_dir.rglob("*.sbgnml"))
        + list(data_dir.rglob("*.xml"))
    )
    if not candidates:
        raise FileNotFoundError(f"No SBGN-ML files found under {data_dir}")
    return candidates


def run(
    data_dir: Path,
    threshold: float,
    use_biocypher: bool,
    out_dir: Path,
    biocypher_config: Path,
    schema_config: Path,
) -> int:
    sbgn_files = find_sbgn_files(data_dir)
    logger.info("found %d SBGN-ML files in %s", len(sbgn_files), data_dir)
    adapter = SBGNPDAdapter(sbgn_files, matcher_threshold=threshold)
    adapter.load()
    logger.info("\n%s", adapter.report())

    if not use_biocypher:
        # Smoke test: exhaust the generators to trigger stats without writing.
        for _ in adapter.get_nodes():
            pass
        for _ in adapter.get_edges():
            pass
        logger.info("smoke-test complete (no BioCypher output)\n%s", adapter.report())
        return 0

    try:
        from biocypher import BioCypher
    except ImportError:
        logger.error(
            "biocypher is not installed. install with `uv sync` or "
            "`pip install biocypher`, or re-run with --no-biocypher."
        )
        return 2

    bc = BioCypher(
        biocypher_config_path=str(biocypher_config),
        schema_config_path=str(schema_config),
        output_directory=str(out_dir),
    )
    bc.write_nodes(adapter.get_nodes())
    bc.write_edges(adapter.get_edges())
    bc.write_import_call()
    bc.summary()
    logger.info(
        "adapter final report after BioCypher write:\n%s", adapter.report()
    )
    return 0


def main(argv: list[str] | None = None) -> int:
    root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=root / "data",
        help="Directory of SBGN-ML files (recursed).",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.7,
        help="Composite matcher threshold in [0, 1]. 0.7 is high-precision.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=root / "biocypher-out",
        help="Directory to write BioCypher output files.",
    )
    parser.add_argument(
        "--biocypher-config",
        type=Path,
        default=root / "config" / "biocypher_config.yaml",
    )
    parser.add_argument(
        "--schema-config",
        type=Path,
        default=root / "config" / "schema_config.yaml",
    )
    parser.add_argument(
        "--no-biocypher",
        action="store_true",
        help="Skip BioCypher entirely and just run the adapter (useful when "
        "biocypher is not installed and you want to verify the adapter).",
    )
    args = parser.parse_args(argv)
    return run(
        data_dir=args.data_dir,
        threshold=args.threshold,
        use_biocypher=not args.no_biocypher,
        out_dir=args.out_dir,
        biocypher_config=args.biocypher_config,
        schema_config=args.schema_config,
    )


if __name__ == "__main__":
    sys.exit(main())
