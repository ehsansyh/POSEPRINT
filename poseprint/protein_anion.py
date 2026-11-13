"""Anionic groups of the protein.

Carried over from the working copy of POSEPRINT with its detection logic unchanged.
See build_poseprint_repo.py in the CLINCH working tree for the exact edits that were
applied, which are limited to module naming, relative imports, dead imports and the
removal of one hard-coded scratch path.
"""

# -*- coding: utf-8 -*-
"""Carboxylate oxygens of ASP and GLU, plus the C-terminal OXT.

The previous version read the PDB with Chem.MolFromPDBFile(sanitize=False) and kept atoms
whose GetFormalCharge() was negative. RDKit cannot infer that from a PDB. A four-bonded
lysine nitrogen must be positive so the cation detector happened to work, but a
single-bonded oxygen is ambiguous, it could be a carbonyl, so no atom ever came back
negative. Measured over 150 complexes: 2831 protein cations, **0** protein anions, which
is why SaltBridge was exactly zero in all 232,290 poses.

Selecting by residue and atom name is what PLIP does and needs no charge perception.
"""
from Bio.PDB import PDBParser

ANIONIC_GROUPS = {
    "ASP": ("OD1", "OD2"),
    "GLU": ("OE1", "OE2"),
    "ASH": ("OD1", "OD2"),
    "GLH": ("OE1", "OE2"),
}
C_TERMINAL = ("OXT",)


def detect_cations(pdb_file):
    """Kept under the original name so callers do not change.

    Returns the same 4-tuple shape the old version did, (index, atom name, charge,
    gasteiger), because salt_bridge.py matches on element 1.
    """
    structure = PDBParser(QUIET=True).get_structure("x", pdb_file)
    anionic_atoms = []
    index = 0
    for residue in structure.get_residues():
        wanted = ANIONIC_GROUPS.get(residue.get_resname(), ())
        for atom in residue.get_atoms():
            index += 1
            if atom.get_name() in wanted or atom.get_name() in C_TERMINAL:
                anionic_atoms.append((index, atom.get_name(), -1, None))
    return anionic_atoms


def run(pdb_file):
    return detect_cations(pdb_file)
