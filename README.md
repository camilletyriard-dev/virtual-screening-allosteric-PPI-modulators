# Virtual Screening: Allosteric PPI Modulators

**Target:** MDM2 (UniProt [Q00987](https://www.uniprot.org/uniprotkb/Q00987)), MDM2-p53 protein-protein interaction interface
**Goal:** Discover small molecules that disrupt the MDM2-p53 interaction via an AlphaFold + ZINC + AutoDock Vina + QSAR pipeline, screened for drug-likeness, commercial availability, and freedom-to-operate.

## Pipeline

| # | Script | Stage |
|---|--------|-------|
| 01 | `scripts/01_identify_interface.py` | Identify MDM2-p53 interface from experimental complex 1YCR, map onto AlphaFold MDM2 monomer (AF-Q00987) |
| 02 | `scripts/02_filter_candidates.py` | Filter ZINC pool for drug-likeness, structural alerts, similarity to known MDM2-p53 disruptors (Nutlin-3, idasanutlin) |
| 03 | `scripts/03_build_docking_set.py` | Select top-similarity + diverse subset for docking |
| 04 | `scripts/04_prepare_receptor.py` | Prepare AlphaFold MDM2 receptor for AutoDock Vina |
| 05 | `scripts/05_prepare_ligands.py` | Generate 3D conformers, convert to PDBQT (meeko) |
| 06 | `scripts/06_rank_docking_results.py` | Rank by best-pose Vina binding affinity |
| 07 | `scripts/07_train_qsar_model.py` | Train ChEMBL-based QSAR model, combine with docking score |
| 08 | `scripts/08_check_supplier_availability.py` | Check commercial availability via PubChem PUG |
| 09 | `scripts/09_check_patent_landscape.py` | Freedom-to-operate signal via PubChem/SureChEMBL patent cross-refs |
| 10 | `scripts/10_lead_optimization_profile.py` | Drug-likeness rules, NIBR filter, structural alerts, complexity, scaffold analysis |
| 11 | `scripts/11_final_shortlist.py` | Drop redundant enantiomer, produce final shortlist |

**Result:** 161,216 filtered ZINC candidates -> 50 docked -> 6 final shortlisted candidates across 6 distinct scaffolds. See [`results/reports/pipeline_report.md`](results/reports/pipeline_report.md) for the full writeup, funnel, and final candidate table, and [`results/rankings/final_shortlist.csv`](results/rankings/final_shortlist.csv) for the data.

Built using Scientific Agent Skills from https://github.com/k-dense-ai/scientific-agent-skills.

