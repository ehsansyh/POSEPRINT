"""Aromatic stacking, parallel, offset parallel and T-shaped.

Carried over from the working copy of POSEPRINT with its detection logic unchanged.
See build_poseprint_repo.py in the CLINCH working tree for the exact edits that were
applied, which are limited to module naming, relative imports, dead imports and the
removal of one hard-coded scratch path.
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 19:33:35 2024

@author: ehsan
"""

import numpy as np
import networkx as nx
from Bio import PDB
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

def generate_plane(vertices):
    normal = np.cross(vertices[1] - vertices[0], vertices[2] - vertices[0])
    normal = normal / np.linalg.norm(normal)
    return normal


def calculate_angle_between_planes(normal1, normal2):
    # Convert normal vectors to numpy arrays
    normal_vector_plane1 = np.array(normal1)
    normal_vector_plane2 = np.array(normal2)
    
    # Check if the normal vector of plane 1 is on the right side of plane 2
    if np.dot(normal_vector_plane1, normal_vector_plane2) < 0:
        # Swap the planes if plane 1 is on the right side
        normal_vector_plane1, normal_vector_plane2 = normal_vector_plane2, normal_vector_plane1
    
    # Calculate the dot product of the two normal vectors
    dot_product = np.dot(normal_vector_plane1, normal_vector_plane2)
    
    # Calculate the magnitudes of the normal vectors
    magnitude_plane1 = np.linalg.norm(normal_vector_plane1)
    magnitude_plane2 = np.linalg.norm(normal_vector_plane2)
    
    # Calculate the cosine of the angle between the planes
    cos_angle = np.abs(dot_product) / (np.abs(magnitude_plane1) * np.abs(magnitude_plane2))
    
    # Calculate the angle in radians
    angle_radians = np.arccos(cos_angle)
    
    # Convert angle to degrees
    angle_degrees = np.degrees(angle_radians)
    
    return angle_degrees


# PLIP accepts a stacking pair on three conditions, not two. Centroid distance and
# inter-plane angle are not enough, because two rings can sit 4.9 A apart at 62 degrees
# and still be side by side rather than stacked. The third is the offset, how far the
# second centre lies from the first ring's axis once projected onto its plane. Without it
# 55 percent of the pairs this module recorded were pairs PLIP rejects.
#
# The offset limit is deliberately wider than PLIP's 2.0 A. Three reasons, and all of
# them are properties of THIS dataset rather than preferences:
#   1. the labels are soft frequencies over a 45 pose ensemble, not a binary call on one
#      experimental structure, so a near-stacked geometry still carries real signal
#   2. the poses are Boltz-2 predictions and carry coordinate error
#   3. the parallel displaced benzene dimer stays attractive out to roughly this lateral
#      displacement, so 3.0 A is near the physical ceiling rather than an arbitrary number
# Measured on 200 complexes, the offsets of the pairs this module accepted had median
# 2.08 A and p75 3.00 A, so a 3.0 A cut keeps about three quarters of them while still
# discarding the side by side pairs that have no ring overlap at all.
# Set these back to 2.0 (and PICATION_DIST_MAX to 6.0) to reproduce PLIP exactly, which
# is the sensitivity analysis the manuscript should report alongside the main result.
PISTACK_OFFSET_MAX = 3.0


def ring_centre_offset(cent_a, normal_a, cent_b, normal_b):
    """PLIP's offset: project each centre onto the other ring's plane, take the smaller."""
    a = np.asarray(cent_a, dtype=float)
    b = np.asarray(cent_b, dtype=float)
    na = np.asarray(normal_a, dtype=float)
    nb = np.asarray(normal_b, dtype=float)
    d = b - a
    off_a = np.linalg.norm(d - np.dot(d, na) * na)
    off_b = np.linalg.norm(-d - np.dot(-d, nb) * nb)
    return float(min(off_a, off_b))


def chain_guard_is_meaningful(prot_values, lig_values):
    """Should the protein/ligand same-chain test be enforced at all?

    Boltz-2 writes the ligand as its own chain, protein A and ligand B, so requiring
    res['chain'] == lig['chain'] is never true and every ring pair is discarded before
    any geometry is computed. That is why PiStacking came out as exactly zero in all
    232,290 poses. The other PROLINT detectors never hit this because pi_cation,
    hydrophobic_interactions and salt_bridge all have a len(lig_values) == 1 fast path
    that skips the guard, and ring_compares has no such fast path.

    The guard still matters for a complex where the ligand really does sit inside a
    protein chain, so it is kept whenever the two share any chain id and dropped only
    when they cannot possibly match.
    """
    protein_chains = {r['chain'].get_id() for r in prot_values}
    ligand_chains = {l['chain'].get_id() for l in lig_values}
    return bool(protein_chains & ligand_chains)


def compare_rings(ring1, ring2):
    centroid1 = np.mean(ring1, axis=0)
    centroid2 = np.mean(ring2, axis=0)
    
    distance = calculate_distance(centroid1, centroid2)
    
    if distance < 5:
        plane1 = generate_plane(ring1)
        plane2 = generate_plane(ring2)
        
        angle = calculate_angle_between_planes(plane1, plane2)
        return angle
    else:
        return None


# if __name__ == "__main__":
def run(mol_file,pdb_file,ligand_id,ligand_id_mol2,dirct,structure):  
    # pdb_file = "C:/Users/ehsan/Desktop/pdf_interactions/revefenacin.pdb"
    data = rg.detect_rings_from_pdb(pdb_file,structure)
    prot_values = data.get('Amino Acid')
    lig_values = data.get('Ligand')
    aromatic= far.run(mol_file,pdb_file,ligand_id,ligand_id_mol2,dirct)
    enforce_chain = chain_guard_is_meaningful(prot_values, lig_values)
    tshape=[]
    ff=[]
    pl=[]
    for res in prot_values:
        if res.get('residue_name') == 'PRO':
            pass
        else:
            for cent in range(len(res['centroids'])):
                cent_p = res['centroids'][cent]
                for lig in lig_values:
                    if (not enforce_chain) or res['chain'].get_id() == lig['chain'].get_id():
                        for cn in range(len(lig['centroids'])):
                            cent_lig = lig['centroids'][cn]
    
                            # if len(lig_values) == 1:
                            lig_atoms = lig['atoms_in_ring']['atoms'][cn]
                            rings_pdb = lig['rings'][cn]
                            if (is_aromatic_ring(rings_pdb,aromatic)
                                    or far.ring_is_planar(lig['atom_coords'][cn])):
                                p_plane = generate_plane(res['atom_coords'][cent])
                                l_plane = generate_plane(lig['atom_coords'][cn])
                                angle = calculate_angle_between_planes(p_plane, l_plane)
                                distance = calculate_distance(cent_p, cent_lig)
                                offset = ring_centre_offset(cent_p, p_plane,
                                                            cent_lig, l_plane)
                                stacked = offset <= PISTACK_OFFSET_MAX
                                    
                                if stacked and angle <= 90 and angle >= 60 and distance <= 5.5:
                                   
                                    tshape.append((lig['residue_name'], lig['residue_id'][1],res['residue_name'], res['residue_id'][1],distance,angle))
                                    
            
                                    prot_cord = np.array(res['atom_coords'][cent] + [res['atom_coords'][cent][0]])
                                    lig_cord = np.array(lig['atom_coords'][cn] + [lig['atom_coords'][cn][0]])
                                   
                                elif stacked and angle <= 30 and angle >= 0.0 and distance <= 5.5:
                                   
                                    pl.append((lig['residue_name'], lig['residue_id'][1],res['residue_name'], res['residue_id'][1],distance,angle))
                                    
                                elif stacked and angle <= 60 and angle >= 30 and distance <= 4.0 :

                                    lig_atoms = lig['atoms_in_ring']['atoms'][cn]
                                    
                                    ff.append((lig['residue_name'], lig['residue_id'][1],res['residue_name'], res['residue_id'][1],distance,angle))
                        # else:
                        #     for lig in lig_values:
                        #         if res['chain'].get_id() == lig['chain'].get_id():
                        #             lig_atoms = lig['atoms_in_ring']['atoms'][cn] 
                        #             rings_pdb = lig_values[0]['rings'][cn]
                        #             if is_aromatic_ring(rings_pdb,aromatic):
                        #                 p_plane = generate_plane(res['atom_coords'][cent])
                        #                 l_plane = generate_plane(lig['atom_coords'][cn])
                        #                 angle = calculate_angle_between_planes(p_plane, l_plane)
                        #                 distance = calculate_distance(cent_p, cent_lig)
                                
                        #                 if angle <= 90 and angle >= 60 and distance <= 5.5  :
                                            
                        #                     tshape.append((lig['residue_name'], lig['residue_id'][1],res['residue_name'], res['residue_id'][1],distance,angle))
                                            
                        #                 elif angle <= 30 and angle >= 0.0 and distance <= 5.5:
                                            
                        #                     pl.append((lig['residue_name'], lig['residue_id'][1],res['residue_name'], res['residue_id'][1],distance,angle))
                        #                 elif angle <= 60 and angle >= 30 and distance <= 4.0 :
                                            
                        #                     ff.append((lig['residue_name'], lig['residue_id'][1],res['residue_name'], res['residue_id'][1],distance,angle))
    return tshape, pl, ff
