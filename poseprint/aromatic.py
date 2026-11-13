"""Aromatic ring perception in the ligand.

Carried over from the working copy of POSEPRINT with its detection logic unchanged.
See build_poseprint_repo.py in the CLINCH working tree for the exact edits that were
applied, which are limited to module naming, relative imports, dead imports and the
removal of one hard-coded scratch path.
"""

# -*- coding: utf-8 -*-
"""
Created on Fri May  3 15:38:25 2024

@author: ehsan
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
import numpy as np
import networkx as nx
import MDAnalysis as md
import warnings
import subprocess
import os


# Planarity cut in Angstrom, RMS distance of the ring atoms from their own best-fit
# plane. Measured over 149 compounds against their SMILES: rings that should be aromatic
# have median 0.0035 and p90 0.0078, rings that should be saturated have median 0.1751.
# A cut of 0.02 sits above the aromatic p90 and far below the saturated median, and gave
# 100 percent agreement on "does this ligand have an aromatic ring" against 40.7 percent
# for the MOL2 route, whose bond blocks were flattened in place by the L5 fallback.
AROMATIC_PLANARITY_CUT = 0.02
# Planarity alone accepts any flat ring, and a flat ring is not necessarily aromatic.
# Measured over 200 complexes, the accepted rings included 6 three-membered and 2
# four-membered ones, which cannot be aromatic at all. Epoxides, cyclopropanes and
# azetidines are flat or nearly so. Aromatic rings in drug-like ligands are 5 or 6.
AROMATIC_RING_SIZES = (5, 6)


def ring_is_planar(atom_coords, cut=AROMATIC_PLANARITY_CUT):
    """Geometric aromaticity test, independent of MOL2 bond orders.

    The coordinates come from the PDB through detect_residue_ring_update and were never
    touched by the MOL2 rewriting, so this works on files whose bond block is already
    destroyed.
    """
    import numpy as _np
    pts = _np.asarray(atom_coords, dtype=float)
    if pts.ndim != 2 or pts.shape[0] < 3:
        return False
    if pts.shape[0] not in AROMATIC_RING_SIZES:
        return False
    centred = pts - pts.mean(axis=0)
    normal = _np.linalg.svd(centred)[2][-1]
    rmsd = float(_np.sqrt(_np.mean((centred @ normal) ** 2)))
    return rmsd <= cut


def extract_ligand_from_mol(mol_file, ligand_id, ligand_id_mol2, dirct):

    file = md.Universe(mol_file)
    try:
        fil1 = file.select_atoms(f'resname {ligand_id_mol2}')
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # Suppress the warning
            fil1.write(dirct+'ligand.mol2')
    except IndexError:
        fil1 = file.select_atoms(f'resname {ligand_id}')
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # Suppress the warning
            fil1.write(dirct+'ligand.mol2')


def _normalize_problematic_mol2_atom_types(mol2_path):
    """
    Normalize known-bad MOL2 atom types only in the temporary ligand file.
    RDKit cannot read some files where O.co2 atoms have degree > 1.
    """
    try:
        with open(mol2_path, 'r') as f:
            lines = f.readlines()

        atom_start = atom_end = bond_start = bond_end = None
        for idx, line in enumerate(lines):
            stripped = line.strip()
            if stripped == '@<TRIPOS>ATOM':
                atom_start = idx + 1
            elif stripped == '@<TRIPOS>BOND':
                atom_end = idx
                bond_start = idx + 1
            elif stripped.startswith('@<TRIPOS>') and bond_start is not None and bond_end is None:
                bond_end = idx
                break

        if atom_start is None or bond_start is None:
            return False
        if atom_end is None:
            atom_end = len(lines)
        if bond_end is None:
            bond_end = len(lines)

        atom_rows = {}
        degree_map = {}

        for idx in range(atom_start, atom_end):
            parts = lines[idx].split()
            if len(parts) >= 6:
                atom_id = int(parts[0])
                atom_rows[atom_id] = (idx, parts)
                degree_map.setdefault(atom_id, 0)

        for idx in range(bond_start, bond_end):
            parts = lines[idx].split()
            if len(parts) >= 3:
                a1 = int(parts[1])
                a2 = int(parts[2])
                degree_map[a1] = degree_map.get(a1, 0) + 1
                degree_map[a2] = degree_map.get(a2, 0) + 1

        modified = False
        for atom_id, (idx, parts) in atom_rows.items():
            if parts[5] == 'O.co2' and degree_map.get(atom_id, 0) > 1:
                parts[5] = 'O.2'
                lines[idx] = ' '.join(parts) + '\n'
                modified = True

        if modified:
            with open(mol2_path, 'w') as f:
                f.writelines(lines)
        return modified
    except Exception:
        return False


def _load_ligand_mol2(ligand_mol2_file):
    mol = Chem.MolFromMol2File(ligand_mol2_file)
    if mol is not None:
        return mol

    _normalize_problematic_mol2_atom_types(ligand_mol2_file)
    return Chem.MolFromMol2File(ligand_mol2_file)


def find_aromatic_rings(ligand_mol2_file):
    
    # mol = Chem.MolFromPDBFile(ligand_mol2_file)

    mol = _load_ligand_mol2(ligand_mol2_file)
    if mol is None:
        return []
    # AllChem.ComputeGasteigerCharges(mol)
    aromatic_atoms = set(atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic())
    # Detect aromatic atoms, including nitrogen
    # aromatic_atoms = set(atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic())

    # Get aromatic rings
    rings = Chem.GetSymmSSSR(mol)
    aromatic_rings = [ring for ring in rings if all(atom in aromatic_atoms for atom in ring)]

    return aromatic_rings

def get_atom_names(mol):
    """Get the names of atoms in the molecule."""
    
    # mol2 = Chem.MolFromPDBFile(mol)

    mol2 = _load_ligand_mol2(mol)
    if mol2 is None:
        return []
    atom_names = [atom.GetSymbol() for atom in mol2.GetAtoms()]
    return atom_names

def run(mol_file,pdb_file,ligand_id,ligand_id_mol2,dirct):
    """Aromatic rings of the ligand, or an empty list.

    This must never raise. ring_compares and pi_cation both fall back to ring_is_planar
    when this returns nothing, so an empty result costs nothing, whereas an exception
    kills the whole pose in run_parallel and loses every other interaction for it.
    A full PROLINT run over 232,290 poses died on exactly this.
    """
    try:
        extract_ligand_from_mol(mol_file, ligand_id, ligand_id_mol2, dirct)
    except Exception:
        return []
    if not os.path.isfile(dirct + 'ligand.mol2'):
        return []
    try:
        return _run_inner(mol_file, pdb_file, ligand_id, ligand_id_mol2, dirct)
    except Exception:
        return []


def _run_inner(mol_file,pdb_file,ligand_id,ligand_id_mol2,dirct):
    
    ligand_pdb_file = dirct+'ligand.pdb'
    ligand_mol2_file = dirct+'ligand.mol2'
    # mol2_file = dirct+'ligand.mol2'
    # sdf_file = dirct+'ligand.sdf'

    # pdb_to_mol2(ligand_pdb_file, ligand_mol2_file)
    # ligand_smiles = mol2_to_smiles(ligand_mol2_file)
    aromatic_rings = find_aromatic_rings(ligand_mol2_file)
    # print(aromatic_rings)
    
    ar_ring=[]
    atom_names = get_atom_names(ligand_mol2_file)
    for i, ring in enumerate(aromatic_rings, 1):
        ring_list = list(ring)  # Convert _vectint to list
        atom_names_in_ring = [atom_names[idx] for idx in ring_list if idx < len(atom_names)]
        # print(f"Aromatic ring {i}: {atom_names_in_ring}")
        # print(ring_list)
        ar_ring.append({'atom names' : atom_names_in_ring,
                   'atom numbers' : ring_list})
    return ar_ring


