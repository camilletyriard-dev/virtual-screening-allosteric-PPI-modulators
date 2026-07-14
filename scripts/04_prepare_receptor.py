"""
Prepare the MDM2 receptor for AutoDock Vina docking.

Note: mdm2_p53_interface.json's pocket_centroid_xyz was computed from the
experimental 1YCR complex's own coordinate frame (see 01_identify_interface.py).
That frame is NOT compatible with AF-Q00987-MDM2.pdb -- the AlphaFold monomer
model has its own independent coordinate frame, even though residue numbering
matches. Docking against the AF model with the 1YCR-frame centroid would put
the search box in the wrong place in space. This script recomputes the pocket
centroid directly from the AF model's own interface-residue CA atoms and
records it in the interface JSON as pocket_centroid_xyz_af_frame.
"""
import json
import subprocess

import numpy as np
from Bio.PDB import PDBParser

AF_MDM2_PDB = "data/structures/AF-Q00987-MDM2.pdb"
INTERFACE_JSON = "data/structures/mdm2_p53_interface.json"
OUT_BASENAME = "data/docking/mdm2_receptor"
BOX_PADDING = 10.0  # Angstrom padding around interface-residue CA extent

MK_PREPARE_RECEPTOR = (
    r"C:\Users\PC\AppData\Local\Packages\PythonSoftwareFoundation.Python.3.11_qbz5n2kfra8p0"
    r"\LocalCache\local-packages\Python311\Scripts\mk_prepare_receptor.exe"
)


def main():
    with open(INTERFACE_JSON) as f:
        interface = json.load(f)
    resnums = set(interface["mdm2_interface_resnums"])

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("AF-MDM2", AF_MDM2_PDB)
    chain = structure[0]["A"]

    coords = []
    for res in chain:
        if res.id[1] in resnums and "CA" in res:
            coords.append(res["CA"].coord)
    coords = np.array(coords)
    print(f"Found {len(coords)} of {len(resnums)} interface-residue CA atoms in AF model")

    centroid = coords.mean(axis=0)
    extent = coords.max(axis=0) - coords.min(axis=0)
    box_size = extent + 2 * BOX_PADDING

    print(f"AF-frame pocket centroid: {centroid.tolist()}")
    print(f"Interface CA extent (A): {extent.tolist()}")
    print(f"Vina box size (A, with {BOX_PADDING} A padding): {box_size.tolist()}")

    interface["pocket_centroid_xyz_af_frame"] = centroid.tolist()
    interface["pocket_box_size_af_frame"] = box_size.tolist()
    interface["note_af_frame"] = (
        "pocket_centroid_xyz (original field) is in 1YCR's experimental coordinate "
        "frame and must NOT be used to define a docking box on AF-Q00987-MDM2.pdb. "
        "pocket_centroid_xyz_af_frame below is recomputed directly from the AF "
        "model's own interface-residue CA atoms and is the correct box center for "
        "docking against that structure."
    )
    with open(INTERFACE_JSON, "w") as f:
        json.dump(interface, f, indent=2)
    print(f"Updated {INTERFACE_JSON} with AF-frame centroid/box size")

    cx, cy, cz = centroid.tolist()
    sx, sy, sz = box_size.tolist()
    cmd = [
        MK_PREPARE_RECEPTOR,
        "--read_pdb", AF_MDM2_PDB,
        "-o", OUT_BASENAME,
        "-p",
        "-v",
        "--box_center", str(cx), str(cy), str(cz),
        "--box_size", str(sx), str(sy), str(sz),
        "--allow_bad_res",
    ]
    print("Running:", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    print(result.stdout)
    if result.returncode != 0:
        print(result.stderr)
        raise SystemExit(f"mk_prepare_receptor failed with code {result.returncode}")
    print(f"Receptor PDBQT + Vina box config written to {OUT_BASENAME}*")


if __name__ == "__main__":
    main()
