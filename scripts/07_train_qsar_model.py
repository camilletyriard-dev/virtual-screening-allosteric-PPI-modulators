"""
Build an MDM2-binding QSAR model from ChEMBL bioactivity data and use it to
add an ML-based potency prediction (predicted pChEMBL) alongside the Vina
docking score for the 50 candidates.

Training data: data/chembl/mdm2_bioactivity_raw.csv, pulled from ChEMBL target
CHEMBL5023 (E3 ubiquitin-protein ligase Mdm2) bioactivity records (IC50/Ki/Kd),
combining a top-potency batch and a mid-potency (1-100 uM) batch for a less
skewed training distribution.
"""
from pathlib import Path

import datamol as dm
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold, cross_val_predict
from sklearn.metrics import r2_score, mean_absolute_error

RAW_CSV = "data/chembl/mdm2_bioactivity_raw.csv"
TRAIN_CSV = "data/chembl/mdm2_qsar_train.csv"
DOCKING_SET_CSV = "data/compounds/docking_set.csv"
RANKED_CSV = "results/rankings/docking_ranked.csv"
OUT_CSV = "results/rankings/docking_ranked_with_qsar.csv"

FP_TYPE = "ecfp"
FP_RADIUS = 2
FP_NBITS = 2048


def featurize(smiles_list):
    fps = []
    valid_idx = []
    for i, smi in enumerate(smiles_list):
        mol = dm.to_mol(smi, sanitize=True)
        if mol is None:
            continue
        fp = dm.to_fp(mol, fp_type=FP_TYPE, radius=FP_RADIUS, fpSize=FP_NBITS)
        if fp is None:
            continue
        fps.append(fp)
        valid_idx.append(i)
    return np.array(fps), valid_idx


def main():
    raw = pd.read_csv(RAW_CSV)
    print(f"Loaded {len(raw)} raw bioactivity records")

    raw["canonical_smiles"] = dm.utils.parallelized(
        lambda s: dm.to_smiles(dm.to_mol(s)) if dm.to_mol(s) is not None else None,
        raw["smiles"].tolist(), n_jobs=4, scheduler="threads", progress=True,
    )
    raw = raw.dropna(subset=["canonical_smiles"])

    agg = raw.groupby("canonical_smiles")["pchembl_value"].mean().reset_index()
    agg.columns = ["smiles", "pchembl_value"]
    print(f"Aggregated to {len(agg)} unique compounds, pChEMBL range "
          f"[{agg.pchembl_value.min():.2f}, {agg.pchembl_value.max():.2f}]")
    agg.to_csv(TRAIN_CSV, index=False)

    X, valid_idx = featurize(agg["smiles"].tolist())
    y = agg["pchembl_value"].values[valid_idx]
    print(f"Featurized {len(X)} compounds with {FP_TYPE} fingerprints ({FP_NBITS} bits)")

    model = RandomForestRegressor(n_estimators=500, max_depth=None, n_jobs=4, random_state=0)
    cv = KFold(n_splits=5, shuffle=True, random_state=0)
    y_pred_cv = cross_val_predict(model, X, y, cv=cv, n_jobs=1)
    r2 = r2_score(y, y_pred_cv)
    mae = mean_absolute_error(y, y_pred_cv)
    print(f"5-fold CV: R2={r2:.3f}, MAE={mae:.3f} pChEMBL units")

    model.fit(X, y)

    dock_set = pd.read_csv(DOCKING_SET_CSV)
    X_cand, valid_cand_idx = featurize(dock_set["smiles"].tolist())
    pred = model.predict(X_cand)

    pred_df = pd.DataFrame({
        "zinc_id": dock_set["zinc_id"].values[valid_cand_idx].astype(str),
        "qsar_pred_pchembl": pred,
    })

    ranked = pd.read_csv(RANKED_CSV)
    ranked["zinc_id"] = ranked["zinc_id"].astype(str)
    merged = ranked.merge(pred_df, on="zinc_id", how="left")

    # Combined score: average of z-scored docking affinity (more negative = better,
    # so flip sign) and z-scored QSAR-predicted potency (higher = better).
    dock_z = -(merged["vina_affinity_kcal_mol"] - merged["vina_affinity_kcal_mol"].mean()) / merged["vina_affinity_kcal_mol"].std()
    qsar_z = (merged["qsar_pred_pchembl"] - merged["qsar_pred_pchembl"].mean()) / merged["qsar_pred_pchembl"].std()
    merged["combined_score"] = (dock_z + qsar_z) / 2
    merged = merged.sort_values("combined_score", ascending=False).reset_index(drop=True)
    merged.insert(0, "combined_rank", merged.index + 1)

    Path("results/rankings").mkdir(parents=True, exist_ok=True)
    merged.to_csv(OUT_CSV, index=False)
    print(f"\nWrote combined ranking to {OUT_CSV}")

    print("\nTop 10 by combined score (docking + QSAR):")
    cols = ["combined_rank", "zinc_id", "vina_affinity_kcal_mol", "qsar_pred_pchembl", "combined_score"]
    print(merged[cols].head(10).to_string(index=False))


if __name__ == "__main__":
    main()
