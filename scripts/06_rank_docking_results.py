"""
Parse AutoDock Vina batch docking output (one PDBQT per ligand, poses ranked
best-first) and merge best-pose binding affinities with candidate properties
from the docking set for a final ranked table.
"""
import re
from pathlib import Path

import pandas as pd

POSES_DIR = "data/docking/poses"
DOCKING_SET_CSV = "data/compounds/docking_set.csv"
OUT_CSV = "results/rankings/docking_ranked.csv"

VINA_RESULT_RE = re.compile(r"^REMARK VINA RESULT:\s*(-?\d+\.?\d*)")


def best_affinity(pdbqt_path):
    with open(pdbqt_path) as f:
        for line in f:
            m = VINA_RESULT_RE.match(line)
            if m:
                return float(m.group(1))
    return None


def main():
    rows = []
    for pdbqt in sorted(Path(POSES_DIR).glob("*_out.pdbqt")):
        zinc_id = pdbqt.stem.removesuffix("_out")
        affinity = best_affinity(pdbqt)
        rows.append({"zinc_id": zinc_id, "vina_affinity_kcal_mol": affinity})

    dock_df = pd.DataFrame(rows)
    dock_df["zinc_id"] = dock_df["zinc_id"].astype(str)
    print(f"Parsed docking results for {len(dock_df)} ligands")

    props_df = pd.read_csv(DOCKING_SET_CSV)
    props_df["zinc_id"] = props_df["zinc_id"].astype(str)

    merged = props_df.merge(dock_df, on="zinc_id", how="left")
    n_missing = merged["vina_affinity_kcal_mol"].isna().sum()
    if n_missing:
        print(f"WARNING: {n_missing} candidates missing docking results")

    merged = merged.sort_values("vina_affinity_kcal_mol", ascending=True).reset_index(drop=True)
    merged.insert(0, "docking_rank", merged.index + 1)

    Path("results/rankings").mkdir(parents=True, exist_ok=True)
    merged.to_csv(OUT_CSV, index=False)
    print(f"Wrote ranked results to {OUT_CSV}")

    print("\nTop 10 by predicted binding affinity:")
    cols = ["docking_rank", "zinc_id", "vina_affinity_kcal_mol", "mw", "logp", "max_ref_similarity"]
    print(merged[cols].head(10).to_string(index=False))


if __name__ == "__main__":
    main()
