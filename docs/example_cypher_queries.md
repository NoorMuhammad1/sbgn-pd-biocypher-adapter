# Example Cypher queries against the SBGN-PD knowledge graph

After running `python create_knowledge_graph.py`, load the generated CSVs into Neo4j:

```bash
# Windows PowerShell
biocypher-out/neo4j-admin-import-call.ps1

# Linux / macOS
bash biocypher-out/neo4j-admin-import-call.sh
```

Then connect to your Neo4j instance and try the queries below. They are ordered from simple to more interesting.

## 1. Count nodes by Biolink category

```cypher
MATCH (n)
RETURN labels(n)[0] AS category, count(*) AS n
ORDER BY n DESC;
```

Expected output for the shipped glycolysis sample:

| category                     | n |
|------------------------------|---|
| SmallMolecule                | 6 |
| MacromolecularMachineMixin   | 3 |
| BiologicalProcess            | 3 |
| CellularComponent            | 2 |

## 2. Every reaction and the enzyme that catalyses it

```cypher
MATCH (enzyme:MacromolecularMachineMixin)-[:CATALYZES]->(process:BiologicalProcess)
RETURN enzyme.label AS enzyme, process.label AS reaction;
```

## 3. What consumes ATP, and what produces ADP

Uses the reified process node to walk from substrate to product across the reaction.

```cypher
MATCH (atp {label: "ATP"})-[:HAS_INPUT]->(rxn:BiologicalProcess)-[:HAS_OUTPUT]->(product)
RETURN rxn.label AS reaction, product.label AS product;
```

## 4. Show entity-matcher merges

Entities collapsed across source files record their pre-merge ids in `merged_from`. Useful for auditing the matcher.

```cypher
MATCH (n)
WHERE n.merged_from IS NOT NULL
RETURN n.label AS canonical_label, n.merged_from AS merged_from, n.id AS canonical_id;
```

Expected for the shipped sample: ATP (`atp_upper` + `atp_lower`) and glucose-6-phosphate (`g6p` + `g6p_b`) should each show up.

## 5. Paths between two metabolites

Uses the SBGN process-reified graph to find any 2-4 hop path from glucose to fructose-1,6-bisphosphate.

```cypher
MATCH path = (start {label: "alpha-D-glucose"})-[*2..8]->(end {label: "fructose 1,6-bisphosphate"})
RETURN path
LIMIT 5;
```

Because processes are reified as nodes, the shortest path from a substrate to a downstream product runs `substrate -has_input-> process -has_output-> product` (two edges per reaction). Adjust the depth accordingly.

## 6. Reactions in a specific compartment

```cypher
MATCH (compartment:CellularComponent {label: "cytosol"})
MATCH (entity {compartment: compartment.id})-[:HAS_INPUT]->(process)
RETURN DISTINCT process.label AS reaction;
```

## 7. Traversal by SBGN arc class (preserved on the edge)

The BioCypher schema exposes the original SBGN arc class as an edge property, so a query can go both broad (via Biolink predicates) and specific (via SBGN class).

```cypher
MATCH ()-[e {sbgn_arc_class: "catalysis"}]->()
RETURN count(*) AS catalysis_edges;
```
