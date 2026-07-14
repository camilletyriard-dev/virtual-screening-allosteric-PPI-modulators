"""
Generate 3D conformers for the docking-set candidates and convert each to a
PDBQT ligand file for AutoDock Vina, using meeko (Gasteiger charges, rotatable
bonds detected automatically).
"""
import subprocess
from pathlib import Path

import datamol as dm
import pandas as pd
from rdkit import Chem

IN_CSV = "data/compounds/docking_set.csv"
OUT_SDF = "data/docking/ligands_3d.sdf"
OUT_LIGAND_DIR = "data/docking/ligands_pdbqt"

MK_PREPARE_LIGAND = (
    r"C:\Users\PC\AppData\Local\Packages\PythonSoftwareFoundation.Python.3.11_qbz5n2kfra8p0"
    r"\LocalCache\local-packages\Python311\Scripts\mk_prepare_ligand.exe"
)


def embed_one(row):
    mol = dm.to_mol(row.smiles, sanitize=True)
    if mol is None:
        return None
    mol = dm.conformers.generate(
        mol, n_confs=1, minimize_energy=True, forcefield="MMFF94s", align_conformers=False
    )
    mol = Chem.AddHs(mol, addCoords=True)
    mol.SetProp("_Name", str(row.zinc_id))
    return mol


def main():
    df = pd.read_csv(IN_CSV)
    print(f"Loaded {len(df)} candidates from {IN_CSV}")

    mols = dm.utils.parallelized(
        embed_one, list(df.itertuples(index=False)), n_jobs=4, scheduler="threads", progress=True
    )
    n_failed = sum(1 for m in mols if m is None)
    mols = [m for m in mols if m is not None]
    print(f"Embedded {len(mols)} conformers ({n_failed} failed)")

    Path("data/docking").mkdir(parents=True, exist_ok=True)
    writer = Chem.SDWriter(OUT_SDF)
    for m in mols:
        writer.write(m)
    writer.close()
    print(f"Wrote 3D conformers to {OUT_SDF}")

    Path(OUT_LIGAND_DIR).mkdir(parents=True, exist_ok=True)
    cmd = [
        MK_PREPARE_LIGAND,
        "-i", OUT_SDF,
        "--multimol_outdir", OUT_LIGAND_DIR,
    ]
    print("Running:", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    print(result.stdout)
    if result.stderr:
        print(result.stderr)

    pdbqt_files = list(Path(OUT_LIGAND_DIR).glob("*.pdbqt"))
    print(f"Produced {len(pdbqt_files)} ligand PDBQT files in {OUT_LIGAND_DIR}")
    if len(pdbqt_files) < len(mols):
        missing = len(mols) - len(pdbqt_files)
        print(f"WARNING: {missing} ligands failed PDBQT conversion (check meeko stderr above)")


if __name__ == "__main__":
    main()
