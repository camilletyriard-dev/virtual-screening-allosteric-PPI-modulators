"""
Lead-optimization profiling of the 7 survivor candidates (docking + QSAR +
purchasable + clear IP). Applies medchem drug-likeness rules, NIBR
screening-deck severity, structural alerts, complexity/synthetic-accessibility
metrics, and privileged-scaffold detection.

Note: molfeat is incompatible with this environment's Python 3.11 (molfeat
requires-python caps at <3.11), so descriptor/fingerprint featurization below
uses datamol/RDKit directly instead, which covers the same ECFP/2D-descriptor
featurization molfeat would otherwise wrap.
"""
from pathlib import Path

import datamol as dm
import medchem as mc
import pandas as pd
from rdkit.Chem.Scaffolds import MurckoScaffold

IN_CSV = "results/rankings/top_candidates_patent_check.csv"
OUT_CSV = "results/rankings/lead_optimization_profile.csv"


def main():
    df = pd.read_csv(IN_CSV)
    survivors = df[df["likely_clear_ip"] == True].reset_index(drop=True)
    print(f"Profiling {len(survivors)} survivor candidates")

    mols = [dm.to_mol(s) for s in survivors["smiles"]]

    rules_df = mc.rules.RuleFilters(
        rule_list=["rule_of_five", "rule_of_veber", "rule_of_leadlike_soft", "rule_of_cns"]
    )(mols=mols, n_jobs=1, progress=False, keep_props=True)

    nibr_df = mc.structural.NIBRFilters()(mols=mols, n_jobs=1, progress=False)

    alerts_df = mc.structural.CommonAlertsFilters()(mols=mols, n_jobs=1, progress=False)

    sas = [dm.descriptors.sas(m) for m in mols]
    bertz_pass = mc.functional.complexity_filter(mols=mols, complexity_metric="bertz", limit="99", n_jobs=1)

    scaffold_group = mc.groups.ChemicalGroup(groups=["privileged_scaffolds"])
    has_priv_scaffold = [scaffold_group.has_match(m) for m in mols]

    murcko_scaffolds = [dm.to_smiles(MurckoScaffold.GetScaffoldForMol(m)) for m in mols]

    out = survivors.copy()
    out["rule_of_five"] = rules_df["rule_of_five"].values
    out["rule_of_veber"] = rules_df["rule_of_veber"].values
    out["rule_of_leadlike_soft"] = rules_df["rule_of_leadlike_soft"].values
    out["rule_of_cns"] = rules_df["rule_of_cns"].values
    out["nibr_pass"] = nibr_df["pass_filter"].values
    out["nibr_severity"] = nibr_df["severity"].values
    out["nibr_reasons"] = nibr_df["reasons"].values
    out["common_alerts_pass"] = alerts_df["pass_filter"].values
    out["common_alerts_reasons"] = alerts_df["reasons"].values
    out["sas_score"] = sas
    out["bertz_complexity_ok_p99"] = bertz_pass
    out["has_privileged_scaffold"] = has_priv_scaffold
    out["murcko_scaffold"] = murcko_scaffolds

    Path("results/rankings").mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_CSV, index=False)
    print(f"\nWrote lead optimization profile to {OUT_CSV}")

    print("\nScaffold grouping:")
    for scaf, sub in out.groupby("murcko_scaffold"):
        print(f"  {scaf}  ->  {sub['zinc_id'].tolist()}")

    print("\nFull profile:")
    cols = ["combined_rank", "zinc_id", "rule_of_five", "rule_of_veber", "rule_of_leadlike_soft",
            "nibr_pass", "nibr_severity", "common_alerts_pass", "sas_score", "bertz_complexity_ok_p99"]
    print(out[cols].to_string(index=False))


if __name__ == "__main__":
    main()
