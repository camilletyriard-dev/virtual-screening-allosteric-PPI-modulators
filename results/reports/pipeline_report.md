# Virtual Screening Pipeline Report: Allosteric MDM2-p53 PPI Modulators

**Target:** MDM2 (UniProt Q00987), MDM2-p53 protein-protein interaction interface
**Goal:** Identify commercially accessible, IP-clear, drug-like small molecules that disrupt the MDM2-p53 interaction, via an AlphaFold + ZINC + AutoDock Vina + QSAR pipeline.

## Pipeline stages

| # | Script | Stage | Output |
|---|--------|-------|--------|
| 01 | `01_identify_interface.py` | Identify MDM2-p53 interface from experimental complex 1YCR (AlphaFold monomer p53 TAD is disordered/low-confidence in isolation), map 22 contact residues onto AlphaFold MDM2 monomer (AF-Q00987) | `data/structures/mdm2_p53_interface.json` |
| 02 | `02_filter_candidates.py` | Filter ZINC pool (MW 300-500, logP 2-4) for drug-likeness (Lipinski/Veber), remove PAINS/Brenk structural alerts, score by ECFP4 similarity to Nutlin-3/idasanutlin | `data/compounds/zinc_filtered.csv` (161,216 candidates) |
| 03 | `03_build_docking_set.py` | Select top-30 by similarity + 20 diverse (from pool of 1,970) for a 50-compound docking set | `data/compounds/docking_set.csv` |
| 04 | `04_prepare_receptor.py` | Prepare AlphaFold MDM2 receptor for Vina; recompute pocket centroid in the AF model's own coordinate frame | receptor PDBQT + `pocket_centroid_xyz_af_frame` |
| 05 | `05_prepare_ligands.py` | Generate 3D conformers, convert to PDBQT via meeko | `data/docking/ligands_pdbqt/` |
| 06 | `06_rank_docking_results.py` | Parse Vina batch output, rank by best-pose binding affinity | `results/rankings/docking_ranked.csv` (50 candidates) |
| 07 | `07_train_qsar_model.py` | Train RandomForest QSAR model on ChEMBL MDM2 bioactivity (CHEMBL5023), predict pChEMBL, combine with docking score | `results/rankings/docking_ranked_with_qsar.csv` |
| 08 | `08_check_supplier_availability.py` | Check commercial availability for top 15 via PubChem PUG (CID lookup + vendor cross-refs) | `results/rankings/top_candidates_supplier_check.csv` |
| 09 | `09_check_patent_landscape.py` | Freedom-to-operate signal: PubChem PatentID cross-refs (SureChEMBL-mined) for candidates with resolved CID | `results/rankings/top_candidates_patent_check.csv` (7/9 checked candidates show zero patent xrefs) |
| 10 | `10_lead_optimization_profile.py` | Profile 7 IP-clear survivors: Ro5/Veber/leadlike-soft/CNS rules, NIBR screening-deck filter, structural alerts, SAS, Bertz complexity (p99), privileged scaffolds, Murcko scaffold grouping | `results/rankings/lead_optimization_profile.csv` |
| 11 | `11_final_shortlist.py` | Drop redundant enantiomer (ZINC5030666, weaker of an identical-property pair with ZINC5030661) | `results/rankings/final_shortlist.csv` (6 candidates) |

## Funnel

161,216 filtered ZINC candidates → 50 docked → 15 supplier-checked → 9 with resolved PubChem CID → 7 IP-clear (zero patent xrefs) → 7 pass lead-optimization filters → **6 final shortlist candidates** (6 distinct Murcko scaffolds)

## Final shortlist

| Rank | ZINC ID | Vina (kcal/mol) | QSAR pred. pChEMBL | Vendors | Patents | SAS | Scaffold class |
|------|---------|------------------|---------------------|---------|---------|-----|-----------------|
| 3 | 2382308675 | -6.712 | 6.48 | ChemBridge, Mcule | 0 | 3.15 | Triazole-piperidine |
| 5 | 4488995 | -8.426 | 5.44 | ZINC | 0 | 2.62 | Dihydroisoquinolinone |
| 6 | 55154650 | -7.699 | 5.60 | Ambinter, ChemBridge, Mcule, Molport, ZINC | 0 | 2.34 | Piperazinone-pyrazole (scaffold A) |
| 7 | 2325804728 | -8.126 | 5.28 | Enamine | 0 | 2.58 | Aminopyridinone |
| 12 | 14539886 | -7.751 | 5.22 | Ambinter, ChemBridge, Mcule, Molport | 0 | 2.50 | Piperazinone-pyrazole (scaffold B) |
| 13 | 5030661 | -7.176 | 5.45 | ZINC | 0 | 2.70 | Isoxazoline |

All 6 pass Rule of Five, Veber, NIBR screening-deck (severity 0), and structural-alert filters; excellent synthetic accessibility (SAS 2.3-3.2). ZINC55154650 and ZINC14539886 are independent hits on the same piperazinone-pyrazole pharmacophore, giving useful SAR redundancy without being structurally identical.

## Caveats

- **Interface/pocket:** derived from the experimental 1YCR complex, not from a co-folded AlphaFold-Multimer model; the AF monomer's pocket geometry may differ subtly from the bound conformation.
- **Docking:** single receptor conformation, no induced-fit or ensemble docking; Vina scores are a coarse ranking signal, not a binding-affinity prediction.
- **QSAR model:** trained on a relatively small/heterogeneous ChEMBL bioactivity set (mixed assay types, top- and mid-potency batches combined); predictions should be treated as a secondary prioritization signal alongside docking, not a standalone potency estimate.
- **IP screen:** absence of PubChem/SureChEMBL patent cross-references is a useful early signal, not formal FTO clearance -- Markush genus claims can cover a compound without naming it. A real FTO opinion requires attorney review.
- **Supplier check:** vendor cross-references indicate the exact structure appears in a catalog aggregator; stock/lead-time/price still need direct vendor confirmation before ordering.

## Suggested next steps

1. Confirm live supplier quotes and lead times for the 6 shortlisted ZINC IDs before committing to synthesis/purchase.
2. Prioritize in vitro binding assay (e.g. FP or SPR against MDM2) for the 2-3 top candidates by combined docking+QSAR score (ZINC4488995, ZINC2325804728, ZINC55154650).
3. If resources allow, run both piperazinone-pyrazole hits (ZINC55154650, ZINC14539886) in parallel as an early SAR pair on that pharmacophore.
4. For a real IP position, escalate the 6 shortlisted structures to attorney-level Markush/genus search before any development investment.
