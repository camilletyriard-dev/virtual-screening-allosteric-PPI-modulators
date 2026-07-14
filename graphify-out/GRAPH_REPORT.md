# Graph Report - .  (2026-07-14)

## Corpus Check
- Large corpus: 2620 files · ~3,620,871 words. Semantic extraction will be expensive (many Claude tokens). Consider running on a subfolder.

## Summary
- 38 nodes · 37 edges · 9 communities (7 shown, 2 thin omitted)
- Extraction: 92% EXTRACTED · 5% INFERRED · 3% AMBIGUOUS · INFERRED: 2 edges (avg confidence: 0.68)
- Token cost: 0 input · 39,192 output

## Community Hubs (Navigation)
- Project Overview & Pipeline
- Supplier Availability Check
- Protein Interface Identification
- Candidate Filtering (Descriptors)
- Ligand Preparation
- Docking Results Ranking
- QSAR Model Training
- Receptor Preparation
- Docking Set Construction

## God Nodes (most connected - your core abstractions)
1. `Allosteric Modulator Discovery (project goal/rationale)` - 5 edges
2. `get_cid()` - 3 edges
3. `main()` - 3 edges
4. `get_heavy_atoms()` - 2 edges
5. `main()` - 2 edges
6. `descriptors_ok()` - 2 edges
7. `main()` - 2 edges
8. `embed_one()` - 2 edges
9. `main()` - 2 edges
10. `best_affinity()` - 2 edges

## Surprising Connections (you probably didn't know these)
- None detected - all connections are within the same source files.

## Import Cycles
- None detected.

## Hyperedges (group relationships)
- **AlphaFold + ZINC + DiffDock virtual screening pipeline** — readme_alphafold, readme_zinc, readme_diffdock, readme_scientific_agent_skills [EXTRACTED 0.90]

## Communities (9 total, 2 thin omitted)

### Community 0 - "Project Overview & Pipeline"
Cohesion: 0.33
Nodes (7): Allosteric Modulator Discovery (project goal/rationale), AlphaFold, DiffDock, Scientific Agent Skills (k-dense-ai/scientific-agent-skills), Target Protein (UniProt ID placeholder), Virtual Screening: Allosteric PPI Modulators, ZINC (compound database)

### Community 1 - "Supplier Availability Check"
Cohesion: 0.47
Nodes (5): get_cid(), get_vendor_sources(), main(), Check commercial supplier availability for the top-ranked candidates via the Pub, Returns (cid, status): status is 'found', 'not_in_pubchem' (CID=0     sentinel),

### Community 2 - "Protein Interface Identification"
Cohesion: 0.67
Nodes (3): get_heavy_atoms(), main(), Identify the MDM2-p53 interaction interface using the experimental complex 1YCR

### Community 3 - "Candidate Filtering (Descriptors)"
Cohesion: 0.67
Nodes (3): descriptors_ok(), main(), Filter the ZINC candidate pool (MW 300-500, logP 2-4 tranches) with RDKit/datamo

### Community 4 - "Ligand Preparation"
Cohesion: 0.67
Nodes (3): embed_one(), main(), Generate 3D conformers for the docking-set candidates and convert each to a PDBQ

### Community 5 - "Docking Results Ranking"
Cohesion: 0.67
Nodes (3): best_affinity(), main(), Parse AutoDock Vina batch docking output (one PDBQT per ligand, poses ranked bes

### Community 6 - "QSAR Model Training"
Cohesion: 0.67
Nodes (3): featurize(), main(), Build an MDM2-binding QSAR model from ChEMBL bioactivity data and use it to add

## Ambiguous Edges - Review These
- `Virtual Screening: Allosteric PPI Modulators` → `Target Protein (UniProt ID placeholder)`  [AMBIGUOUS]
  README.md · relation: references

## Knowledge Gaps
- **3 isolated node(s):** `Target Protein (UniProt ID placeholder)`, `ZINC (compound database)`, `Scientific Agent Skills (k-dense-ai/scientific-agent-skills)`
  These have ≤1 connection - possible missing edges or undocumented components.
- **2 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **What is the exact relationship between `Virtual Screening: Allosteric PPI Modulators` and `Target Protein (UniProt ID placeholder)`?**
  _Edge tagged AMBIGUOUS (relation: references) - confidence is low._
- **What connects `Target Protein (UniProt ID placeholder)`, `ZINC (compound database)`, `Scientific Agent Skills (k-dense-ai/scientific-agent-skills)` to the rest of the system?**
  _3 weakly-connected nodes found - possible documentation gaps or missing edges._