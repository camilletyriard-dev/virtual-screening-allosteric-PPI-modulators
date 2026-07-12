"""
Filter the ZINC candidate pool (MW 300-500, logP 2-4 tranches) with RDKit/datamol:
  1. Parse + standardize
  2. Re-verify drug-like descriptor window (Lipinski/Veber)
  3. Remove structural alerts (PAINS, Brenk) via medchem
  4. Score by ECFP4 similarity to known MDM2-p53 PPI disruptors (Nutlin-3, Idasanutlin)
     -- validated chemotypes for exactly this interface, used only as a prioritization
        signal, not a hard filter
  5. Write annotated table and select a final docking set (top-similarity + diverse)
"""
import datamol as dm
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen
from rdkit.DataStructs import BulkTanimotoSimilarity
from rdkit.Chem import rdFingerprintGenerator
import medchem as mc
from medchem.catalogs import NamedCatalogs

IN_SMI = "data/compounds/zinc_candidates_dedup.smi"
OUT_CSV = "data/compounds/zinc_filtered.csv"
OUT_DOCK = "data/compounds/docking_set.csv"
N_JOBS = 4  # keep modest given constrained RAM

REFERENCE_LIGANDS = {
    "nutlin-3": "CC(C)OC1=C(C=CC(=C1)OC)C2=N[C@H]([C@H](N2C(=O)N3CCNC(=O)C3)C4=CC=C(C=C4)Cl)C5=CC=C(C=C5)Cl",
    "idasanutlin": "CC(C)(C)C[C@H]1[C@]([C@H]([C@@H](N1)C(=O)NC2=C(C=C(C=C2)C(=O)O)OC)C3=C(C(=CC=C3)Cl)F)(C#N)C4=C(C=C(C=C4)Cl)F",
}

FPGEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def descriptors_ok(mol):
    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)
    tpsa = Descriptors.TPSA(mol)
    rotb = Descriptors.NumRotatableBonds(mol)
    ok = (300 <= mw <= 500) and (2 <= logp <= 4) and hbd <= 5 and hba <= 10 and tpsa <= 140 and rotb <= 10
    return ok, mw, logp, hbd, hba, tpsa, rotb


def main():
    print(f"Loading {IN_SMI} ...")
    df = pd.read_csv(IN_SMI, sep=r"\s+", names=["smiles", "zinc_id"], skiprows=1)
    print(f"Loaded {len(df)} raw ZINC entries")

    mols = dm.utils.parallelized(dm.to_mol, df["smiles"].tolist(), n_jobs=N_JOBS, scheduler="threads", progress=True)
    df["mol"] = mols
    df = df[df["mol"].notna()].reset_index(drop=True)
    print(f"{len(df)} parsed successfully")

    df["mol"] = dm.utils.parallelized(
        lambda m: dm.standardize_mol(m, disconnect_metals=True, normalize=True, reionize=True),
        df["mol"].tolist(), n_jobs=N_JOBS, scheduler="threads", progress=True,
    )
    df = df[df["mol"].notna()].reset_index(drop=True)

    desc = [descriptors_ok(m) for m in df["mol"]]
    desc_df = pd.DataFrame(desc, columns=["druglike_ok", "mw", "logp", "hbd", "hba", "tpsa", "rotb"])
    df = pd.concat([df.reset_index(drop=True), desc_df], axis=1)
    print(f"Druglike window pass: {df['druglike_ok'].sum()} / {len(df)}")

    # structural alerts: PAINS + Brenk via medchem
    pains_catalog = NamedCatalogs.pains()
    brenk_catalog = NamedCatalogs.brenk()

    def has_alert(mol):
        return bool(pains_catalog.HasMatch(mol)) or bool(brenk_catalog.HasMatch(mol))

    df["has_alert"] = dm.utils.parallelized(has_alert, df["mol"].tolist(), n_jobs=N_JOBS, scheduler="threads", progress=True)
    print(f"Flagged by PAINS/Brenk: {df['has_alert'].sum()} / {len(df)}")

    passed = df[df["druglike_ok"] & ~df["has_alert"]].reset_index(drop=True)
    print(f"Passing drug-likeness + no structural alerts: {len(passed)}")

    # similarity to reference MDM2-p53 PPI disruptor chemotypes
    ref_mols = {name: Chem.MolFromSmiles(smi) for name, smi in REFERENCE_LIGANDS.items()}
    ref_fps = {name: FPGEN.GetFingerprint(m) for name, m in ref_mols.items()}

    passed_fps = dm.utils.parallelized(FPGEN.GetFingerprint, passed["mol"].tolist(), n_jobs=N_JOBS, scheduler="threads", progress=True)
    for name, ref_fp in ref_fps.items():
        passed[f"sim_{name}"] = BulkTanimotoSimilarity(ref_fp, passed_fps)
    passed["max_ref_similarity"] = passed[[f"sim_{n}" for n in REFERENCE_LIGANDS]].max(axis=1)

    passed = passed.sort_values("max_ref_similarity", ascending=False).reset_index(drop=True)
    passed_out = passed.drop(columns=["mol"])
    passed_out.to_csv(OUT_CSV, index=False)
    print(f"Saved {len(passed_out)} filtered candidates to {OUT_CSV}")

    # Build final docking set: top-N by reference similarity + a diverse subset
    top_similarity = passed.head(30)
    remaining_pool = passed.iloc[30:2000]  # cap diversity search to a manageable pool
    if len(remaining_pool) > 0:
        diverse_idx, _ = dm.pick_diverse(remaining_pool["mol"].tolist(), npick=min(20, len(remaining_pool)))
        diverse_subset = remaining_pool.iloc[diverse_idx]
    else:
        diverse_subset = remaining_pool

    docking_set = pd.concat([top_similarity, diverse_subset]).drop_duplicates(subset="zinc_id")
    docking_set = docking_set.drop(columns=["mol"])
    docking_set.to_csv(OUT_DOCK, index=False)
    print(f"Final docking set: {len(docking_set)} compounds -> {OUT_DOCK}")


if __name__ == "__main__":
    main()
