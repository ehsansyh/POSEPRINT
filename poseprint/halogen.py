"""Halogen bonds from ligand halogens to protein acceptors.

Carried over from the working copy of POSEPRINT with its detection logic unchanged.
See build_poseprint_repo.py in the CLINCH working tree for the exact edits that were
applied, which are limited to module naming, relative imports, dead imports and the
removal of one hard-coded scratch path.
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 14:23:21 2024

@author: ehsan
"""

from Bio.PDB import *
from Bio.PDB.vectors import calc_angle
import numpy as np
import math
from numpy import *
from numpy.linalg import norm
import os

def calculate_distance(acceptor_coords,halogen_coords):
    """Calculate distance between two atoms."""
    distance = math.sqrt((acceptor_coords[0]-(halogen_coords[0]))**2+(acceptor_coords[1]-(halogen_coords[1]))**2+(acceptor_coords[2]-(halogen_coords[2]))**2)
    return distance

def calculate_angle(atom1, atom2, atom3):
    """Calculate angle between three atoms."""
    f = atom2-atom1
    e = atom2-atom3
    abVec = norm(f)
    bcVec = norm(e)
    abNorm = f / abVec;
    bcNorm = e / bcVec;
    res = abNorm[0] * bcNorm[0] + abNorm[1] * bcNorm[1] + abNorm[2] * bcNorm[2];
    angle = arccos(res)*180.0/ pi
    return angle

def is_donor(atom):
    """Check if the atom can act as a halogen bond donor."""
    parent_residue = atom.get_parent().get_id()[0]
    if parent_residue != ' ':
        # Protein residue
        if atom.element in ("CL", "BR", "I"):
            return True
    return False


def prot_acceptor(atom):
    parent_residue = atom.get_parent().get_id()[0]
    if parent_residue == ' ':
        if atom.element in ('N', 'O', 'S'):

            return True
    return False


def detect_interactions(pdb_file,structure):
    """Detect interactions between ligand and protein atoms."""

    interactions = {}
    

    for res in structure.get_residues():
        if res.get_id()[0] != ' ':
            for atom in res.get_atoms():
                neighbor_atoms = []
                ns = NeighborSearch(list(structure.get_atoms()))
                neighbors = ns.search(atom.get_coord(), 3.0)
                
                for neighbor_atom in neighbors:
                    neighbor_residue = neighbor_atom.get_parent()
                    if neighbor_residue.get_id()[0] == ' ':
                        neighbor_atoms.append(neighbor_atom)
                
                for neighbor_atom in neighbor_atoms:
                    interaction_type = (res.get_resname(), neighbor_atom.get_parent().get_resname(), neighbor_atom.get_parent().get_id()[1])
                    interaction_atoms = (atom.get_name(), neighbor_atom.get_name())
                    interactions.setdefault(interaction_type, []).append(interaction_atoms)
    return interactions    
     
def halobond(pdb_file,structure, distance_max=4):

    dist=[]
    halocord=[]
    halobonds = {}
    for res in structure.get_residues():
        for atom in res.get_atoms():
            if is_donor(atom) == True:
                ns = NeighborSearch(list(structure.get_atoms()))
                neighbors = ns.search(atom.coord, distance_max)  
                for neighbor_atom in neighbors:
                    if prot_acceptor(neighbor_atom) == True:
                        donor_coords = atom.coord
                        acceptor_coords = neighbor_atom.coord
                        ns = NeighborSearch(neighbor_atom.parent.child_list)
                        neighbac = ns.search(neighbor_atom.get_coord(), 1.8)
                        ns = NeighborSearch(atom.parent.child_list)
                        neighbdn = ns.search(atom.get_coord(), 1.8)
                        for neighbs in neighbac:
                            covbond_coords = neighbs.coord
                            if list(covbond_coords) == list(acceptor_coords):
                                pass
                            else:
                                angle1 = calculate_angle(np.array(covbond_coords), np.array(acceptor_coords),np.array(donor_coords))
                                if angle1 >= 90:
                                    for neighbb in neighbdn:    
                                        cabond_coord = neighbb.coord                          
                                        if list(cabond_coord) == list(donor_coords):
                                            pass
                                        else:
                                            angle = calculate_angle(np.array(cabond_coord), np.array(donor_coords), np.array(acceptor_coords))
                                            if angle >= 140.0:
                                                distance= calculate_distance(acceptor_coords,donor_coords)
                                                if distance <= 4.0:  
                                                    interaction_type = (atom.get_parent().get_resname(), atom.get_parent().get_id()[1], neighbor_atom.get_parent().get_resname(), neighbor_atom.get_parent().get_id()[1])
                                                    interaction_atoms = (atom.get_name(),neighbs.get_name(), neighbor_atom.get_name())
                                                    dist.append(distance)
                                                    halobonds.setdefault(interaction_type, []).append((interaction_atoms, distance, angle))
                                            else:
                                                pass           
            else:
                continue
                            
    return halobonds



def main(pdb_file,structure):
    interactions = detect_interactions(pdb_file,structure)
    # Print detected interactions
    for interaction_type, atoms in interactions.items():
        ligand_residue, protein_residue, protein_id = interaction_type
    halobonds = halobond(pdb_file,structure)
    return halobonds
