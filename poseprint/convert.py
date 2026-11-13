"""Turn a complex in PDB format into the MOL2 the ring perception needs.

The detectors read the protein from the PDB and the ligand from a MOL2 of the same
complex, because MOL2 carries the bond orders and the aromatic flags that ring
perception depends on and PDB does not.

Open Babel is used when it is on the path, which is what produced the files behind the
published labels. MDAnalysis is the fallback. It writes a usable MOL2 but it guesses the
bonds from distances, so a complex converted that way can differ on the aromatic
interactions. The hydrogen-bond and hydrophobic counts do not depend on this file.
"""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

__all__ = ["pdb_to_mol2", "have_openbabel"]


def have_openbabel() -> bool:
    return shutil.which("obabel") is not None


def pdb_to_mol2(pdb_path, mol2_path=None, prefer="obabel"):
    """Write a MOL2 beside the PDB, or at mol2_path, and return its path."""
    pdb_path = Path(pdb_path)
    if not pdb_path.is_file():
        raise FileNotFoundError(pdb_path)
    mol2_path = Path(mol2_path) if mol2_path else pdb_path.with_suffix(".mol2")
    mol2_path.parent.mkdir(parents=True, exist_ok=True)

    if prefer == "obabel" and have_openbabel():
        proc = subprocess.run(
            ["obabel", "-ipdb", str(pdb_path), "-omol2", "-O", str(mol2_path)],
            capture_output=True, text=True)
        if mol2_path.is_file() and mol2_path.stat().st_size > 0:
            return mol2_path
        raise RuntimeError("Open Babel could not convert %s\n%s"
                           % (pdb_path.name, proc.stderr.strip()))

    try:
        import MDAnalysis as mda
    except ImportError as exc:                       # pragma: no cover
        raise RuntimeError(
            "Converting a PDB to MOL2 needs either Open Babel on the path or "
            "MDAnalysis installed.") from exc
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        universe = mda.Universe(str(pdb_path), guess_bonds=True)
        universe.atoms.write(str(mol2_path))
    return mol2_path
