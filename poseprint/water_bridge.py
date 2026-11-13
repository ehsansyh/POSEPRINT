"""Water-mediated hydrogen bonds.

Carried over from the working copy of POSEPRINT with its detection logic unchanged.
See build_poseprint_repo.py in the CLINCH working tree for the exact edits that were
applied, which are limited to module naming, relative imports, dead imports and the
removal of one hard-coded scratch path.
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 15:43:38 2024

@author: ehsan
"""

from Bio.PDB import *
from Bio.PDB.vectors import calc_angle
import numpy as np
import math
from numpy import *
from numpy.linalg import norm
import os



def calculate_distance(acceptor_coords,hydrogen_coords):
    """Calculate distance between two atoms."""
    distance = math.sqrt((acceptor_coords[0]-(hydrogen_coords[0]))**2+(acceptor_coords[1]-(hydrogen_coords[1]))**2+(acceptor_coords[2]-(hydrogen_coords[2]))**2)
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
    angle = (arccos(res)*180.0)/pi
    return angle

def is_donor(atom):
    """Check if the atom can act as a hydrogen bond donor."""
    parent_residue = atom.get_parent().get_id()[0]
    if parent_residue != ' ' and parent_residue != 'W':
        # Protein residue
        if atom.element in ('N', 'O'):
            for donor_neighbor in atom.parent.child_list:
                ns = NeighborSearch(atom.parent.child_list)
                neighbors = ns.search(atom.get_coord(), 1.1)
                for neighbor in neighbors:
                    if neighbor.element == 'H':
                        return True
    return False

def water_donor(atom):
    """Check if the atom can act as a hydrogen bond donor."""
    parent_residue = atom.get_parent().get_id()[0]
    if parent_residue == 'W':
        if atom.element in ('H'):
            return True
        else:
            return False
    return False

def prot_donor(atom):
    parent_residue = atom.get_parent().get_id()[0]
    if parent_residue == ' ':
        if atom.element in ('N', 'O'):
            for donor_neighbor in atom.parent.child_list:
                ns = NeighborSearch(atom.parent.child_list)
                neighbors = ns.search(atom.get_coord(), 1.1)
                for neighbor in neighbors:
                    if neighbor.element == 'H':
                        return True
    return False

def is_acceptor(atom):
    """Check if the atom can act as a hydrogen bond acceptor."""
    parent_residue = atom.get_parent().get_id()[0]
    if parent_residue != ' ' and parent_residue != 'W':
        # Non-protein residue
        if atom.element in ("O", "N", "F"):
            return True
    return False

def water_acceptor(atom):
    """Check if the atom can act as a hydrogen bond donor."""
    parent_residue = atom.get_parent().get_id()[0]
    if parent_residue == 'W':
        # Protein residue
        if atom.element in ('O'):
            return True
        else:
            return False
    return False

def prot_acceptor(atom):
    parent_residue = atom.get_parent().get_id()[0]
    if parent_residue == ' ':
        if atom.element in ("O", "N", "F"):

            return True
    return False

def detect_interactions(pdb_file,structure):
    """Detect interactions between ligand and protein atoms."""

    interactions = {}

    for res in structure.get_residues():
        if res.get_id()[0] == 'W':
            for atom in res.get_atoms():
                neighbor_atoms = []
                ns = NeighborSearch(list(structure.get_atoms()))
                neighbors = ns.search(atom.get_coord(), 3.0)
                
                for neighbor_atom in neighbors:
                    neighbor_residue = neighbor_atom.get_parent()
                    if neighbor_residue.get_id()[0] == ' ' or neighbor_residue.get_id()[0] != ' ' and neighbor_residue.get_id()[0] != 'W':
                        neighbor_atoms.append(neighbor_atom)
                for neighbor_atom in neighbor_atoms:
                    interaction_type = (res.get_resname(), neighbor_atom.get_parent().get_resname(), neighbor_atom.get_parent().get_id()[1])
                    interaction_atoms = (atom.get_name(), neighbor_atom.get_name())
                    interactions.setdefault(interaction_type, []).append(interaction_atoms)
    return interactions    
     
def wbbond(pdb_file,structure, distance_max=4):

    dist_wa=[]
    dist_wd=[]
    hcord=[]
    wbbonds = {}
    for res in structure.get_residues():
        for atom in res.get_atoms():
            if atom.get_parent().get_id()[0] != ' ' and atom.get_parent().get_id()[0] != 'W':
                if is_donor(atom) == True:
                    ns = NeighborSearch(list(structure.get_atoms()))
                    neighbors = ns.search(atom.coord, distance_max)  
                    for neighbor_atom in neighbors:
                        if water_acceptor(neighbor_atom) == True:

                            ns = NeighborSearch(list(structure.get_atoms()))
                            neighbors = ns.search(neighbor_atom.coord, 1.1)
                            for neighbor_atom1 in neighbors:
                                if water_donor(neighbor_atom1) == True:
                                    # print(neighbor_atom1)
                                    ns = NeighborSearch(list(structure.get_atoms()))
                                    neighbors = ns.search(neighbor_atom1.coord, distance_max)
                                    for neighbor_atom_p in neighbors:
                                        if prot_acceptor(neighbor_atom_p) == True:
                                            # print(neighbor_atom_p)
                                            ns = NeighborSearch(atom.parent.child_list)
                                            neighb = ns.search(atom.get_coord(), 1.1)
                                            for neighbs in neighb:
                                                if neighbs.element == 'H':
                                                    # print(neighbs)
                                                    hydrogen_coords = neighbs.coord
                                                    donor_coords = atom.coord
                                                    acceptor_coords = neighbor_atom_p.coord
                                                    water_donor_coord= neighbor_atom1.coord
                                                    water_acc_coord= neighbor_atom.coord

                                                    angle_dw = calculate_angle(np.array(donor_coords), np.array(hydrogen_coords), np.array(water_acc_coord))

                                                    if angle_dw >= 140.0:
                                                        # print(angle_dw,'***')
                                                        angle_aw = calculate_angle(np.array(acceptor_coords), np.array(water_acc_coord), np.array(hydrogen_coords))
                                                        if angle_aw >= 71.0 and angle_aw <= 140.0:
                                                            # print(angle_aw)
                                                            distance_wd= calculate_distance(donor_coords,water_acc_coord)
                                                            distance_wa= calculate_distance(acceptor_coords,water_acc_coord)
                                                            if distance_wd <= 3.0 and distance_wd >= 1.9:
                                                                if distance_wa <= 3.0 and distance_wa >= 1.9:

                                                                    interaction_type = (atom.get_parent().get_resname(), atom.get_parent().get_id()[1], neighbor_atom.get_parent().get_resname(), neighbor_atom.get_parent().get_id()[1],neighbor_atom_p.get_parent().get_resname(), neighbor_atom_p.get_parent().get_id()[1])
                                                                    interaction_atoms = (atom.get_name(), neighbor_atom.get_name(),neighbor_atom_p.get_name())
                                                                    data = (interaction_atoms,distance_wd,distance_wa,angle_dw, angle_aw)
                                                                    dist_wa.append(distance_wa)
                                                                    dist_wd.append(distance_wd)
                                                                    wbbonds.setdefault(interaction_type, [])
                                                                    if data not in wbbonds[interaction_type]:
                                                                        wbbonds[interaction_type].append(data)
                                                    else:
                                                        pass
                                                    
            elif atom.get_parent().get_id()[0] == ' ':
                if prot_donor(atom) == True:
                    ns = NeighborSearch(list(structure.get_atoms()))
                    neighbors = ns.search(atom.coord, distance_max)  # Convert generator to list
                    for neighbor_atom in neighbors:
                        if water_acceptor(neighbor_atom) == True:
                            ns = NeighborSearch(list(structure.get_atoms()))
                            neighbors = ns.search(neighbor_atom.coord, 1.1)
                            for neighbor_atom1 in neighbors:
                                if water_donor(neighbor_atom1) == True:
                                    ns = NeighborSearch(list(structure.get_atoms()))
                                    neighbors = ns.search(neighbor_atom1.coord, distance_max)
                                    for neighbor_atom_l in neighbors:
                                        if is_acceptor(neighbor_atom_l):
                                            ns = NeighborSearch(atom.parent.child_list)
                                            neighb = ns.search(atom.get_coord(), 1.1)
                                            for neighbs in neighb:
                                                if neighbs.element == 'H':
                                                    # print(neighbs)
                                                    hydrogen_coords = neighbs.coord
                                                    donor_coords = atom.coord
                                                    acceptor_coords = neighbor_atom_l.coord
                                                    water_donor_coord= neighbor_atom1.coord
                                                    water_acc_coord= neighbor_atom.coord
                                                    angle_dw = calculate_angle(np.array(donor_coords), np.array(hydrogen_coords), np.array(water_acc_coord))
                                                    if angle_dw >= 140.0:
                                                        angle_aw = calculate_angle(np.array(acceptor_coords), np.array(water_acc_coord), np.array(hydrogen_coords))
                                                        if angle_aw >= 71.0 and angle_aw <= 140.0:
                                                            distance_wd= calculate_distance(donor_coords,water_acc_coord)
                                                            distance_wa= calculate_distance(acceptor_coords,water_acc_coord)
                                                            if distance_wd <= 3.0 and distance_wd >= 1.9:
                                                                if distance_wa <= 3.0 and distance_wa >= 1.9:
                                                                    interaction_type = (atom.get_parent().get_resname(), atom.get_parent().get_id()[1], neighbor_atom.get_parent().get_resname(), neighbor_atom.get_parent().get_id()[1],neighbor_atom_l.get_parent().get_resname(), neighbor_atom_l.get_parent().get_id()[1])
                                                                    interaction_atoms = (atom.get_name(), neighbor_atom.get_name(),neighbor_atom_l.get_name())
                                                                    data = (interaction_atoms,distance_wd,distance_wa,angle_dw, angle_aw)
                                                                    dist_wa.append(distance_wa)
                                                                    dist_wd.append(distance_wd)
                                                                    wbbonds.setdefault(interaction_type, [])
                                                                    if data not in wbbonds[interaction_type]:
                                                                        wbbonds[interaction_type].append(data)
                                                    else:
                                                        pass
            else:
                continue
                            
    return wbbonds

def main(pdb_file,structure):
    interactions = detect_interactions(pdb_file,structure)
    wbbonds = wbbond(pdb_file,structure)
    return wbbonds

