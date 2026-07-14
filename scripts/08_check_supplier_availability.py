"""
Check commercial supplier availability for the top-ranked candidates via the
PubChem PUG REST API. ZINC15's own site now gates programmatic access behind
a CAPTCHA, so instead we: (1) look up each candidate's PubChem CID by SMILES,
then (2) pull its cross-referenced source databases (xrefs/SourceName), which
for purchasable compounds includes commercial catalog vendors (Enamine,
MolPort, Mcule, ChemDiv, ChemBridge, etc.) alongside literature/patent
sources. Presence of known vendor names is used as a purchasability signal.

Respects PubChem's request-rate guidance (<=5 req/s) with a throttle.
"""
import time
import urllib.parse

import pandas as pd
import requests

RANKED_CSV = "results/rankings/docking_ranked_with_qsar.csv"
OUT_CSV = "results/rankings/top_candidates_supplier_check.csv"
TOP_N = 15
REQUEST_DELAY_S = 0.35

PUG_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

KNOWN_VENDOR_SOURCES = {
    "Enamine", "MolPort", "Mcule", "ChemDiv", "ChemBridge", "Vitas-M Laboratory",
    "Life Chemicals", "Otava Chemicals", "Aldrich Market Select", "eMolecules",
    "ZINC", "Ambinter", "AKos Consulting & Solutions", "Princeton BioMolecular Research",
    "ChemSpace", "Analyticon Discovery", "Apollo Scientific", "BLD Pharm",
    "Sigma-Aldrich", "TargetMol", "MedChemExpress",
}


def get_cid(smiles):
    """Returns (cid, status): status is 'found', 'not_in_pubchem' (CID=0
    sentinel), or 'lookup_error' (structure could not be parsed/standardized)."""
    url = f"{PUG_BASE}/compound/smiles/{urllib.parse.quote(smiles)}/cids/JSON"
    r = requests.get(url, timeout=20)
    if r.status_code != 200:
        return None, "lookup_error"
    data = r.json()
    cids = data.get("IdentifierList", {}).get("CID", [])
    if not cids or cids[0] == 0:
        return None, "not_in_pubchem"
    return cids[0], "found"


def get_vendor_sources(cid):
    url = f"{PUG_BASE}/compound/cid/{cid}/xrefs/SourceName/JSON"
    r = requests.get(url, timeout=20)
    if r.status_code != 200:
        return []
    data = r.json()
    infos = data.get("InformationList", {}).get("Information", [])
    sources = set()
    for info in infos:
        for s in info.get("SourceName", []):
            sources.add(s)
    return sorted(sources)


def main():
    df = pd.read_csv(RANKED_CSV).head(TOP_N).copy()
    print(f"Checking supplier availability for top {len(df)} candidates via PubChem")

    results = []
    for _, row in df.iterrows():
        zinc_id, smiles = row["zinc_id"], row["smiles"]
        cid, status = None, "lookup_error"
        vendors_found = []
        try:
            cid, status = get_cid(smiles)
            time.sleep(REQUEST_DELAY_S)
            if cid:
                sources = get_vendor_sources(cid)
                time.sleep(REQUEST_DELAY_S)
                vendors_found = [s for s in sources if any(v.lower() in s.lower() for v in KNOWN_VENDOR_SOURCES)]
        except requests.RequestException as e:
            print(f"  {zinc_id}: request failed ({e})")

        print(f"  ZINC{zinc_id}: status={status}, PubChem CID={cid}, vendor hits={len(vendors_found)}")
        results.append({
            "zinc_id": zinc_id,
            "pubchem_status": status,
            "pubchem_cid": cid,
            "n_vendor_sources": len(vendors_found),
            "vendor_sources": "; ".join(vendors_found),
            "likely_purchasable": bool(vendors_found),
        })

    out = df.merge(pd.DataFrame(results), on="zinc_id", how="left")
    out.to_csv(OUT_CSV, index=False)
    print(f"\nWrote supplier availability check to {OUT_CSV}")
    print(out[["combined_rank", "zinc_id", "pubchem_status", "pubchem_cid", "n_vendor_sources", "likely_purchasable"]].to_string(index=False))


if __name__ == "__main__":
    main()
