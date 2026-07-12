"""
Identify the MDM2-p53 interaction interface using the experimental complex
1YCR (Kussie et al. 1996), since isolated AlphaFold monomer predictions for
p53's transactivation domain (residues ~17-29) are low-confidence/disordered
(it folds into a helix only upon binding MDM2) and cannot be used to derive
a real MDM2-p53 interface on their own.

Interface residues found here are then mapped onto the AlphaFold MDM2
monomer model (AF-Q00987) for use as the docking pocket.
"""
import json
from Bio.PDB import PDBParser, NeighborSearch
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1

COMPLEX_PDB = "data/structures/1YCR.pdb"
AF_MDM2_PDB = "data/structures/AF-Q00987-MDM2.pdb"
OUT_JSON = "data/structures/mdm2_p53_interface.json"

CONTACT_CUTOFF = 4.5  # Angstrom, heavy-atom distance


def get_heavy_atoms(chain):
    atoms = []
    for res in chain:
        if not is_aa(res, standard=True):
            continue
        for atom in res:
            if atom.element != "H":
                atoms.append(atom)
    return atoms


def main():
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("1YCR", COMPLEX_PDB)
    model = structure[0]
    chain_mdm2 = model["A"]  # MDM2 hydrophobic cleft domain
    chain_p53 = model["B"]  # p53 transactivation domain peptide

    mdm2_atoms = get_heavy_atoms(chain_mdm2)
    p53_atoms = get_heavy_atoms(chain_p53)

    ns = NeighborSearch(mdm2_atoms)

    interface_residues = {}
    for atom in p53_atoms:
        close = ns.search(atom.coord, CONTACT_CUTOFF)
        for a in close:
            res = a.get_parent()
            resnum = res.id[1]
            resname = seq1(res.resname)
            interface_residues.setdefault(resnum, resname)

    interface_sorted = sorted(interface_residues.items())
    print(f"MDM2 interface residues (chain A, within {CONTACT_CUTOFF} A of p53 peptide):")
    for num, name in interface_sorted:
        print(f"  {name}{num}")

    # pocket centroid from interface residue CA atoms, for docking box definition
    import numpy as np
    centroid_atoms = []
    for res in chain_mdm2:
        if res.id[1] in interface_residues and "CA" in res:
            centroid_atoms.append(res["CA"].coord)
    centroid = np.mean(centroid_atoms, axis=0).tolist()

    # p53 peptide key anchor residues (well known from literature: F19, W23, L26)
    p53_seq = "".join(seq1(r.resname) for r in chain_p53 if is_aa(r, standard=True))
    p53_resnums = [r.id[1] for r in chain_p53 if is_aa(r, standard=True)]
    print(f"\np53 peptide (chain B) sequence: {p53_seq}")
    print(f"p53 peptide residue numbers: {p53_resnums[0]}-{p53_resnums[-1]}")

    result = {
        "source_complex": "1YCR",
        "contact_cutoff_angstrom": CONTACT_CUTOFF,
        "mdm2_interface_residues": [f"{name}{num}" for num, name in interface_sorted],
        "mdm2_interface_resnums": [num for num, _ in interface_sorted],
        "pocket_centroid_xyz": centroid,
        "p53_peptide_sequence": p53_seq,
        "p53_key_anchor_residues": ["F19", "W23", "L26"],
        "note": (
            "Interface derived from experimental complex 1YCR because AlphaFold "
            "monomer models of p53 show low/medium pLDDT (40-77) in this region "
            "(residues 15-30), consistent with it being an intrinsically "
            "disordered region that folds into a helix only upon MDM2 binding."
        ),
    }
    with open(OUT_JSON, "w") as f:
        json.dump(result, f, indent=2)
    print(f"\nSaved interface definition to {OUT_JSON}")

    # sanity check: confirm numbering matches the AlphaFold MDM2 model
    af_structure = parser.get_structure("AF-MDM2", AF_MDM2_PDB)
    af_chain = af_structure[0]["A"]
    mismatches = 0
    for num, name in interface_sorted:
        if num in af_chain:
            af_res = af_chain[num]
            af_name = seq1(af_res.resname)
            if af_name != name:
                mismatches += 1
                print(f"  MISMATCH at {num}: 1YCR={name} AF-model={af_name}")
        else:
            mismatches += 1
            print(f"  MISSING residue {num} in AlphaFold model")
    if mismatches == 0:
        print("Residue numbering confirmed consistent between 1YCR and AF-Q00987 model.")
    else:
        print(f"{mismatches} numbering mismatches found - re-check residue offset.")


if __name__ == "__main__":
    main()
