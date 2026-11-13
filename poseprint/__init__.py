"""POSEPRINT, protein-ligand interaction fingerprints from docked poses.

POSEPRINT reads a docked complex and reports, for every residue, which of seven
interaction types it makes with the ligand. It was written for the CLINCH study, where
it labelled 142,160 docked poses of 3,949 ligands against HIV-1 protease, but nothing in
it is specific to that target.

    from poseprint import analyse, fingerprint

    for result in analyse("complex.pdb"):
        print(result.ligand, result.counts())
        print(result.contact_table())

    for fp in fingerprint("complex.pdb"):
        print(fp["hbond_count"], fp["hbond_residues"])

    draw_diagram("complex.pdb", out="interactions.pdf")

The complex has to carry hydrogens. Hydrogen bonds are found through the donor hydrogen,
so a heavy-atom-only file returns nothing rather than failing.
"""
from __future__ import annotations

__version__ = "1.0.0"

from .core import (
    INTERACTION_COLUMNS,
    InteractionResult,
    analyse,
    analyse_many,
    fingerprint,
)
from .convert import have_openbabel, pdb_to_mol2
from .diagram import draw_diagram, interaction_records

__all__ = [
    "__version__",
    "analyse",
    "analyse_many",
    "fingerprint",
    "InteractionResult",
    "INTERACTION_COLUMNS",
    "pdb_to_mol2",
    "have_openbabel",
    "draw_diagram",
    "interaction_records",
]
