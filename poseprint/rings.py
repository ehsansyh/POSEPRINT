"""Ring perception in protein residues.

Carried over from the working copy of POSEPRINT with its detection logic unchanged.
See build_poseprint_repo.py in the CLINCH working tree for the exact edits that were
applied, which are limited to module naming, relative imports, dead imports and the
removal of one hard-coded scratch path.
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 15:58:49 2024

@author: ehsan
"""

import numpy as np
import networkx as nx
from Bio import PDB
import os


def residue_to_graph(residue):
    atoms = [atom for atom in residue.get_atoms() if atom.element != 'H']
    # atoms = {key: value for key, value in atoms.items() if 'H' not in key}    
    num_atoms = len(atoms)
    adjacency_matrix = np.zeros((num_atoms, num_atoms))

    for i in range(num_atoms):
        for j in range(i + 1, num_atoms):
            distance = atoms[i] - atoms[j]
            if distance <= 1.6 and distance >= 1.0:  # Typical C-C bond length is around 1.54 Å, so using 1.6 Å as a threshold
                adjacency_matrix[i, j] = 1
                adjacency_matrix[j, i] = 1

    G = nx.from_numpy_array(adjacency_matrix)
    return G

def calculate_ring_centroid(ring, atoms):
    centroid = np.mean([list(atoms.values())[i] for i in ring], axis=0)
    return centroid

def detect_rings_from_pdb(pdb_file,structure):

    rings_and_centroids_per_residue = {'Amino Acid': [], 'Ligand': []}

    for model in structure:
        for chain in model:
            for residue in chain:
                residue_id = residue.get_id()
                if residue.get_id()[0].strip() != 'W':
                    # Determine residue type (Amino Acid or Ligand)
                    residue_type = 'Amino Acid' if PDB.is_aa(residue) else 'Ligand'
                    # Convert residue to graph
                    G = residue_to_graph(residue)
    
                    # Detect rings
                    rings = list(nx.cycle_basis(G))

                    # Get atom coordinates
                    atoms = {atom.get_id(): atom.get_coord() for atom in residue.get_atoms()}
                    ring_atoms = []
                    r_atoms = []
                    atoms_r =[]
                    for ring_atom_ids in rings:
                        ring_atoms.append([atoms[atom_id] for atom_id in ring_atom_ids if atom_id in atoms])
    
                        #r_atoms.append([list(residue.get_atoms())[i] for i in len(list(residue.get_atoms()) if list(residue.get_atoms())[i].get_id()==atoms[atom_id])])
                    # Calculate centroids for each ring
                    centroids = [calculate_ring_centroid(ring, atoms) for ring in rings]
                    
                    if rings:  # Check if rings are detected
                        # Store ring atom coordinates
                        ring_atoms = [[list(atoms.keys())[i] for i in ring] for ring in rings]
                        
                        for d in range(len(ring_atoms)):
                            for at in range(len(list(residue.get_atoms()))):
                                if list(residue.get_atoms())[at].get_id() in ring_atoms[d]:
                                    r_atoms.append(list(residue.get_atoms())[at])
                            atoms_r.append(r_atoms)
                            r_atoms = []
                        atom_dict= {'ring': ring_atoms,
                                    'atoms': atoms_r}
                            
                        pdb_atoms = {list(residue.get_atoms())[i] for i in range(len(list(residue.get_atoms())))}
                        atom_coords = [[list(atoms.values())[i] for i in ring] for ring in rings]
                        rings_and_centroids_per_residue[residue_type].append({
                            'chain': chain,
                            'residue_name': residue.get_resname(),
                            'residue_id': residue_id,
                            'pdb_residue': residue,
                            'pdb_atoms': pdb_atoms,
                            'rings': rings,
                            'ring_atoms': ring_atoms,
                            'atoms_in_ring': atom_dict,
                            'atom_coords': atom_coords,
                            'centroids': centroids
                        })

    return rings_and_centroids_per_residue

# if __name__ == "__main__":
#     pdb_file = "C:/Users/ehsan/Desktop/pdf_interactions/3HF6.pdb"    
#     data_per_residue_type = detect_rings_from_pdb(pdb_file)
    
#     for residue_type, data_list in data_per_residue_type.items():
#         print(f"{residue_type}:")
        
#         for data in data_list:
#             chain = data['chain']
#             residue_name = data['residue_name']
#             residue_id = data['residue_id']
#             pdb_residue= data['pdb_residue']
#             pdb_atoms= data['pdb_atoms']
#             rings = data['rings']
#             ring_atoms = data['ring_atoms']
#             atoms_in_ring = data['atoms_in_ring']
#             atom_coords = data['atom_coords']
#             centroids = data['centroids']
            
#             print(f"  Residue {residue_name} {residue_id[1]}:")
            
#             # Print rings if they are detected
#             if rings:
#                 print("    Rings:")
#                 for i, ring in enumerate(rings):
#                     print(f"      Ring {i + 1}: {ring}")
#                     print("        Atom Coordinates:")
#                     for atom, coord in zip(ring_atoms[i], atom_coords[i]):
#                         print(f"          Atom {atom} coordinate is : {coord}")
                
#                 print("    Centroids:")
#                 for i, centroid in enumerate(centroids):
#                     print(f"      Centroid {i + 1}: {centroid}")
#             else:
#                 pass

