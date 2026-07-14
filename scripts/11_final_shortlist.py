"""
Produce the final experimental shortlist from the 7 lead-optimization
survivors.

ZINC5030661 and ZINC5030666 are an enantiomeric pair (isoxazoline ring,
opposite stereocenter -- see murcko_scaffold/SMILES in stage 10) with nearly
identical properties (MW, logP, SAS, complexity all match). Keeping both adds
no chemotype diversity, so the weaker-ranked enantiomer (5030666, combined_rank
15, vina -7.109 kcal/mol) is dropped in favor of 5030661 (combined_rank 13,
vina -7.176 kcal/mol), leaving 6 chemically distinct candidates spanning 6
scaffolds for synthesis/assay prioritization.
"""
from pathlib import Path

import pandas as pd

IN_CSV = "results/rankings/lead_optimization_profile.csv"
OUT_CSV = "results/rankings/final_shortlist.csv"

DROP_ZINC_IDS = [5030666]  # redundant enantiomer of 5030661

SUMMARY_COLS = [
    "combined_rank", "zinc_id", "vina_affinity_kcal_mol", "qsar_pred_pchembl",
    "combined_score", "n_vendor_sources", "n_patents", "sas_score", "murcko_scaffold",
]


def main():
    df = pd.read_csv(IN_CSV)
    shortlist = df[~df["zinc_id"].isin(DROP_ZINC_IDS)].copy()
    shortlist = shortlist.sort_values("combined_rank").reset_index(drop=True)

    print(f"Loaded {len(df)} lead-optimization survivors, dropping {DROP_ZINC_IDS} (redundant enantiomer)")
    print(f"Final shortlist: {len(shortlist)} candidates\n")

    Path("results/rankings").mkdir(parents=True, exist_ok=True)
    shortlist.to_csv(OUT_CSV, index=False)
    print(f"Wrote final shortlist to {OUT_CSV}\n")

    print(shortlist[SUMMARY_COLS].to_string(index=False))

    print("\nScaffold diversity check:")
    for scaf, sub in shortlist.groupby("murcko_scaffold"):
        print(f"  {scaf}  ->  ZINC{sub['zinc_id'].tolist()}")


if __name__ == "__main__":
    main()
