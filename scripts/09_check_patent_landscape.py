"""
Screen top candidates for existing patent literature via PubChem's PatentID
cross-references (sourced from SureChEMBL patent-mining of full-text patents).
Absence of any patent hit is a useful (not definitive) freedom-to-operate
signal for early-stage virtual screening hits: it means no patent full-text
corpus mined by SureChEMBL mentions this exact structure. It does not
constitute formal FTO clearance (patents can claim Markush genus structures
that read on a compound without naming it explicitly) - a real FTO opinion
would need attorney review and genus/substructure search of active claims.

Only candidates with a resolved PubChem CID (i.e. that passed the supplier
check) are checked, since compounds absent from PubChem have no path to
purchase regardless of patent status.
"""
import time

import pandas as pd
import requests

SUPPLIER_CSV = "results/rankings/top_candidates_supplier_check.csv"
OUT_CSV = "results/rankings/top_candidates_patent_check.csv"
REQUEST_DELAY_S = 0.35

PUG_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"


def get_patent_ids(cid):
    url = f"{PUG_BASE}/compound/cid/{int(cid)}/xrefs/PatentID/JSON"
    r = requests.get(url, timeout=20)
    if r.status_code != 200:
        return []
    data = r.json()
    infos = data.get("InformationList", {}).get("Information", [])
    patents = []
    for info in infos:
        patents.extend(info.get("PatentID", []))
    return patents


def main():
    df = pd.read_csv(SUPPLIER_CSV)
    has_cid = df["pubchem_cid"].notna()
    print(f"Checking patent landscape for {has_cid.sum()} candidates with a resolved PubChem CID")

    results = []
    for _, row in df.iterrows():
        zinc_id = row["zinc_id"]
        cid = row["pubchem_cid"]
        if pd.isna(cid):
            results.append({"zinc_id": zinc_id, "n_patents": None, "sample_patent_ids": None,
                             "likely_clear_ip": None})
            continue
        patents = get_patent_ids(cid)
        time.sleep(REQUEST_DELAY_S)
        print(f"  ZINC{zinc_id} (CID {int(cid)}): {len(patents)} patent xrefs")
        results.append({
            "zinc_id": zinc_id,
            "n_patents": len(patents),
            "sample_patent_ids": "; ".join(patents[:5]),
            "likely_clear_ip": len(patents) == 0,
        })

    out = df.merge(pd.DataFrame(results), on="zinc_id", how="left")
    out.to_csv(OUT_CSV, index=False)
    print(f"\nWrote patent landscape check to {OUT_CSV}")
    cols = ["combined_rank", "zinc_id", "pubchem_cid", "n_vendor_sources", "n_patents", "likely_clear_ip"]
    print(out[cols].to_string(index=False))


if __name__ == "__main__":
    main()
