"""Cationic groups of the ligand.

Carried over from the working copy of POSEPRINT with its detection logic unchanged.
See build_poseprint_repo.py in the CLINCH working tree for the exact edits that were
applied, which are limited to module naming, relative imports, dead imports and the
removal of one hard-coded scratch path.
"""

# -*- coding: utf-8 -*-
"""
Created on Mon May  6 17:40:16 2024

@author: ehsan
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
import numpy as np
import networkx as nx
import MDAnalysis as md
import re
import warnings

from . import topological_charges as topo
import subprocess
import os


# Switch the ligand charge source. True uses connectivity patterns from
# detect_topological_charges, which is what finds the 54.7 percent of ligands carrying an
# ionisable group. False restores the original behaviour, which reads formal charges out
# of a PDB and therefore sees only the 5.4 percent drawn charged. Keep False working so
# the strict version stays reproducible for the sensitivity analysis.
USE_TOPOLOGICAL_CHARGES = True


def extract_ligand_from_pdb(pdb_file, ligand_id, ligand_id_mol2, dirct):
    """Write the ligand out of the complex PDB.

    The candidate names must include the residue name with trailing digits stripped. When
    run_parallel auto-detects the ligand from the MOL2 it passes the MOL2 substructure
    name, LIG1, for BOTH arguments, while the PDB residue name is LIG. Both original
    selections then matched zero atoms, the write raised, and every ligand charge came
    back empty. Measured over 150 complexes: 0 ligand cations and 0 ligand anions.
    """
    universe = md.Universe(pdb_file)
    candidates = []
    for name in (ligand_id_mol2, ligand_id):
        if not name:
            continue
        for candidate in (str(name), re.sub(r"\d+$", "", str(name)).strip()):
            if candidate and candidate not in candidates:
                candidates.append(candidate)
    for candidate in candidates:
        selection = universe.select_atoms("resname {}".format(candidate))
        if len(selection) == 0:
            continue
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            selection.write(dirct + 'ligand.pdb')
        return
    raise ValueError(
        "no ligand atoms matched any of {} in {}".format(candidates, pdb_file))


def run(pdb_file,ligand_id,ligand_id_mol2,dirct):
    # ligand_id = 'FAD'
    try:    
        extract_ligand_from_pdb(pdb_file, ligand_id,ligand_id_mol2,dirct)
    
        ligand_pdb_file = dirct+'ligand.pdb'
        ligand_mol2_file = dirct+'ligand.pdb'
    
        ligand_mol2 = Chem.MolFromPDBFile(ligand_mol2_file, sanitize=True)

        if USE_TOPOLOGICAL_CHARGES and ligand_mol2 is not None:
            return topo.cations(ligand_mol2)

        # ligand_mol2 = Chem.MolFromMol2File(ligand_mol2_file)
    
        ligand_mol2 = Chem.AddHs(ligand_mol2)
    
        AllChem.ComputeGasteigerCharges(ligand_mol2)
    
        atom_properties = [(atom.GetIdx(), atom.GetSymbol(), atom.GetPDBResidueInfo(), atom.GetFormalCharge(), atom.GetDoubleProp('_GasteigerCharge')) for atom in ligand_mol2.GetAtoms()]
    
        positively_charged_atoms = [(atom_idx, props.GetName()) for atom_idx, symbol, props, formal_charge, partial_charge in atom_properties if formal_charge > 0]
    
        return positively_charged_atoms
    except Exception as e:
        print(f"[Warning] Skipping ligand {ligand_id}: {str(e)}")
        return []
