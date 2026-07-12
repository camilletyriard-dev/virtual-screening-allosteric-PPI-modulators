"""
Resume point: zinc_filtered.csv was already produced by 02_filter_candidates.py
(sorted by max_ref_similarity, descending). This just rebuilds the final
docking set (top-similarity + diverse subset), which crashed previously due
to a dm.pick_diverse() return-value unpacking bug.
"""
import datamol as dm
import pandas as pd

IN_CSV = "data/compounds/zinc_filtered.csv"
OUT_DOCK = "data/compounds/docking_set.csv"

df = pd.read_csv(IN_CSV)
print(f"Loaded {len(df)} filtered candidates")

top_similarity = df.head(30).copy()
remaining_pool = df.iloc[30:2000].copy()
remaining_pool["mol"] = dm.utils.parallelized(dm.to_mol, remaining_pool["smiles"].tolist(), n_jobs=4, scheduler="threads", progress=True)
remaining_pool = remaining_pool[remaining_pool["mol"].notna()].reset_index(drop=True)

diverse_idx, _ = dm.pick_diverse(remaining_pool["mol"].tolist(), npick=min(20, len(remaining_pool)))
diverse_subset = remaining_pool.iloc[diverse_idx].drop(columns=["mol"])

docking_set = pd.concat([top_similarity, diverse_subset]).drop_duplicates(subset="zinc_id")
docking_set.to_csv(OUT_DOCK, index=False)
print(f"Final docking set: {len(docking_set)} compounds -> {OUT_DOCK}")
