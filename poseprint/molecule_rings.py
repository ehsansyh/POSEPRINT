"""Ring perception helpers.

Carried over from the working copy of POSEPRINT with its detection logic unchanged.
See build_poseprint_repo.py in the CLINCH working tree for the exact edits that were
applied, which are limited to module naming, relative imports, dead imports and the
removal of one hard-coded scratch path.
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 18:23:11 2024

@author: ehsan
"""

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
import networkx as nx


def adjacency_matrix_to_graph(adjacency_matrix):
    G = nx.Graph()
    num_atoms = adjacency_matrix.shape[0]
    for i in range(num_atoms):
        G.add_node(i)
    for i in range(num_atoms):
        for j in range(i + 1, num_atoms):
            if adjacency_matrix[i, j] == 1:
                G.add_edge(i, j)
    return G

def detect_rings(mol):
    adjacency_matrix = Chem.GetAdjacencyMatrix(mol)
    G = adjacency_matrix_to_graph(adjacency_matrix)
    rings = list(nx.cycle_basis(G))
    return rings

def draw_molecule_with_rings(mol, rings):
    d = rdMolDraw2D.MolDraw2DCairo(400, 400)
    d.drawOptions().addAtomIndices = True
    
    # Draw molecule
    d.DrawMolecule(mol)
    d.FinishDrawing()
    
    # Get the image data as a binary string
    drawing_binary = d.GetDrawingText()
    
    # Convert binary string to numpy array for image display
    from io import BytesIO
    import matplotlib.image as mpimg
    
    imgdata = BytesIO()
    imgdata.write(drawing_binary)
    imgdata.seek(0)  # rewind the data
    
    # Read image data and display using matplotlib
    img = mpimg.imread(imgdata)
    plt.imshow(img)
    
    # Highlight rings
    for ring in rings:
        for i in range(len(ring)):
            j = (i + 1) % len(ring)
            atom1 = mol.GetAtomWithIdx(ring[i])
            atom2 = mol.GetAtomWithIdx(ring[j])
            plt.plot([atom1.GetIdx()], [atom2.GetIdx()], color='red', linewidth=2)
    
    plt.axis('off')
    plt.show()
# Example usage
if __name__ == "__main__":
    smiles = "c1ccccccc1CCO"  # Phenyl ethanol
    mol = Chem.MolFromSmiles(smiles)
    AllChem.Compute2DCoords(mol)
    rings = detect_rings(mol)
    draw_molecule_with_rings(mol, rings)

