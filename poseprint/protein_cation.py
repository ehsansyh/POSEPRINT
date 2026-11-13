"""Cationic groups of the protein.

Carried over from the working copy of POSEPRINT with its detection logic unchanged.
See build_poseprint_repo.py in the CLINCH working tree for the exact edits that were
applied, which are limited to module naming, relative imports, dead imports and the
removal of one hard-coded scratch path.
"""

# -*- coding: utf-8 -*-
"""
Created on Tue May  7 15:43:20 2024

@author: ehsan
"""

from rdkit import Chem
from rdkit.Chem import AllChem
import os

def detect_cations(pdb_file):
    # Load protein structure from Mol2 file
    # mol = Chem.MolFromMol2File(mol_file, sanitize=True)
    mol = Chem.MolFromPDBFile(pdb_file, sanitize=False)
    # with open(mol_file, 'r') as f:
    #     pdb_block = f.read()
    #     mol = Chem.MolFromPDBBlock(pdb_block, sanitize=True)
    if mol is None:
        raise ValueError("Failed to load protein structure from PDB file.")

    # Compute formal charges for the molecule
    AllChem.ComputeGasteigerCharges(mol)

    # Identify cationic atoms
    cationic_atoms = []
    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() > 0:
            cationic_atoms.append((atom.GetIdx(),atom.GetPDBResidueInfo().GetName(), atom.GetFormalCharge(), atom.GetPropsAsDict().get('_GasteigerHCharge')))

    return cationic_atoms
def run(pdb_file):
    cationic_atoms = detect_cations(pdb_file)
    
    return cationic_atoms

