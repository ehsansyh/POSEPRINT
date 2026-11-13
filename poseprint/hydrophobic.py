"""Hydrophobic contacts between apolar ligand and protein atoms.

Carried over from the working copy of POSEPRINT with its detection logic unchanged.
See build_poseprint_repo.py in the CLINCH working tree for the exact edits that were
applied, which are limited to module naming, relative imports, dead imports and the
removal of one hard-coded scratch path.
"""

# -*- coding: utf-8 -*-
"""
Created on Wed May  8 12:54:47 2024

@author: ehsan
"""

import numpy as np
import networkx as nx
from Bio import PDB
from Bio.PDB import *
from . import rings as rg
import math
from Bio.PDB import NeighborSearch
from . import aromatic as far
import os


def is_aromatic_ring(rings_pdb, aromatic):
    for target in aromatic:
        if set(rings_pdb) == set(target['atom numbers']):
            #print(rings_pdb)
            return True

def calculate_distance(cent_p, cent_lig):
    """Calculate distance between two atoms."""
    distance = math.sqrt((cent_p[0]-(cent_lig[0]))**2+(cent_p[1]-(cent_lig[1]))**2+(cent_p[2]-(cent_lig[2]))**2)
    return distance

def hydrophobic(pdb_file,mol_file,ligand_id,ligand_id_mol2,dirct,structure):
    hph=[]
    data = rg.detect_rings_from_pdb(pdb_file,structure)
    prot_values = data.get('Amino Acid')
    lig_values = data.get('Ligand')
    aromatic= far.run(mol_file,pdb_file,ligand_id,ligand_id_mol2,dirct)
    
    if len(lig_values) == 1:
        for atom in lig_values[0]['pdb_atoms']:
            if atom.element.lower() in ['cl', 'c','s','f']:
                # print(atom.element)
                lig_atom = atom
                lig_coord = lig_atom.coord
    
                ns = NeighborSearch(list(structure.get_atoms()))
                neighbors = ns.search(atom.get_coord(), 4.0)
                for neighbor in neighbors:
                    if neighbor.get_parent().get_resname() in ['ALA','VAL','ILE','LEU','MET', 'PHE', 'TYR', 'TRP']:
                        if neighbor.element.lower() in ['c','s']:
                            res_coord = neighbor.coord
                            distance = calculate_distance(res_coord, lig_coord)
                            # print(distance)
                            if distance <= 4.0:
                                hph_data= (lig_atom.get_parent().get_resname(),lig_atom.get_parent().get_id()[1],neighbor.get_parent().get_resname(), neighbor.get_parent().get_id()[1],distance)
                                if hph_data not in hph:
                                    hph.append(hph_data)
                                
    else:
        for lig in lig_values:
            for atom in lig['pdb_atoms']:
                if atom.element.lower() in ['cl', 'c','s','f']:
                    # print(atom.element)
                    lig_atom = atom
                    lig_coord = lig_atom.coord
                    
                    ns = NeighborSearch(list(structure.get_atoms()))
                    neighbors = ns.search(atom.get_coord(), 4.0)
                    for neighbor in neighbors:
                        if neighbor.get_parent().get_parent().get_id() == lig['chain'].get_id():
                            if neighbor.get_parent().get_resname() in ['ALA','VAL','ILE','LEU','MET', 'PHE', 'TYR', 'TRP']:
                                if neighbor.element.lower() in ['c','s']:
                                    res_coord = neighbor.coord
                                    distance = calculate_distance(res_coord, lig_coord)
                                    # print(distance)
                                    if distance <= 4.0:
                                        hph_data= (lig_atom.get_parent().get_resname(),lig_atom.get_parent().get_id()[1],neighbor.get_parent().get_resname(), neighbor.get_parent().get_id()[1],distance)
                                        if hph_data not in hph:
                                            hph.append(hph_data)
    return hph
                
                
                
                
                
                
                
                