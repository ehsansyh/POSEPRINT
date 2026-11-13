"""Salt bridges between charged ligand and protein groups.

Carried over from the working copy of POSEPRINT with its detection logic unchanged.
See build_poseprint_repo.py in the CLINCH working tree for the exact edits that were
applied, which are limited to module naming, relative imports, dead imports and the
removal of one hard-coded scratch path.
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 17:23:52 2024

@author: ehsan
"""

import numpy as np
import networkx as nx
from Bio import PDB
from Bio.PDB import *
from . import protein_anion as pa
from . import ligand_anion as la
from . import rings as rg
import math
from Bio.PDB import NeighborSearch
from . import aromatic as far
from . import ligand_cation as lc
from . import protein_cation as pc
import os


def calculate_distance(cent_p, cent_lig):
    """Calculate distance between two atoms."""
    distance = math.sqrt((cent_p[0]-(cent_lig[0]))**2+(cent_p[1]-(cent_lig[1]))**2+(cent_p[2]-(cent_lig[2]))**2)
    return distance

# if __name__ == "__main__":
def lig_cation_prot_anion(prot_values,lig_values,protein_anion,ligand_cation,structure): 
    
    sltbrd=[]
    sltbrd_id=[]
    if len(ligand_cation) == 0:
        pass
    else:
        for res in list(structure.get_residues()):
            if res.get_id()[0] == ' ':
            # print(res.get_resname())
            # cat_res=['HIS','ARG','LYS']
            # if res.get_resname() in cat_res:
                # print(res['residue_name'])
                for cnt in range(len(protein_anion)):
                    for atm in list(res.get_atoms()):
                        if atm.name == protein_anion[cnt][1].strip():
                            # print(atom.name,'------', protein_anion[cnt][1])
                            panion= atm.coord
                            if len(lig_values) == 1:
                                lig = lig_values.copy()
                                for cn in range(len(ligand_cation)):
                                    # cation = next((atom for atom in lig[0]['pdb_atoms'] if atom.name == ligand_cation[cn][1]))
                                    for atom in lig[0]['pdb_atoms']:
                                        if atom.name == ligand_cation[cn][1].strip():
                                            # print(atm,atom)

                                            cation=atom
                                            cat_coord = cation.coord
                                            distance = calculate_distance(panion, cat_coord)
                                            # print(distance)
                                            if distance <= 5.5:
                                                resid_no = (cation.get_parent().get_id()[1], atm.get_parent().get_id()[1])
                                                interaction = (cation, atm,distance)
                                                if resid_no not in sltbrd_id:
                                                    sltbrd_id.append(resid_no)
                                                    sltbrd.append(interaction)
                                                # result = print('Salt Bridge interactions between', cation.get_name(), 'of ', cation.get_parent().get_resname(),cation.get_parent().get_id()[1], cation.get_parent().get_parent().get_id(), 'and  ' , 
                                                #                atm.get_id(), atm.get_parent().get_resname(),atm.get_parent().get_id()[1], atm.get_parent().get_parent().get_id())
                                        else:
                                            pass
                            else:
                                for lig in lig_values:
                                    if res.get_parent().get_id() == lig['chain'].get_id():
                                        for cn in range(len(ligand_cation)):
                                            # cation = next((atom for atom in lig[0]['pdb_atoms'] if atom.name == ligand_cation[cn][1]))
                                            for atom in lig['pdb_atoms']:
                                                if atom.name == ligand_cation[cn][1].strip():
                                                    cation=atom
                                                    cat_coord = cation.coord
                                                    distance = calculate_distance(panion, cat_coord)
                                                    if distance <= 5.5:
                                                        resid_no = (cation.get_parent().get_id()[1], atm.get_parent().get_id()[1])
                                                        interaction = (cation, atm, distance)
                                                        if resid_no not in sltbrd_id:
                                                            sltbrd_id.append(resid_no)
                                                            sltbrd.append(interaction)
                                                        # result = print('Salt Bridge interactions between', cation.get_name(), 'of ', cation.get_parent().get_resname(),cation.get_parent().get_id()[1], cation.get_parent().get_parent().get_id(),'and  ' , 
                                                        #                atm.get_id(), atm.get_parent().get_resname(),atm.get_parent().get_id()[1], atm.get_parent().get_parent().get_id())
                                                else:
                                                    pass
    return sltbrd


def prot_cation_lig_anion(prot_values,lig_values,protein_cation,ligand_anion,structure):
    sltbrd_p=[]
    sltbrd_id_p=[]
    if len(ligand_anion) == 0:
        pass
    else:
        if len(lig_values) == 1:
            for cn in range(len(ligand_anion)):
                lig_atoms = ligand_anion[cn][1]
                lig = lig_values.copy()
                for atom in lig[0]['pdb_atoms']:
                    if atom.name == ligand_anion[cn][1].strip():
                        anion=atom
                        anion_coord = anion.coord
                        for res in list(structure.get_residues()):
                            # print(res.get_resname())
                            cat_res=['HIS','ARG','LYS']
                            if res.get_resname() in cat_res:
                                # print(res.get_resname())
                                for cnt in range(len(protein_cation)):
                                    for atmm in list(res.get_atoms()):
                                        if atmm.name == protein_cation[cnt][1].strip():
                                            pcat= atmm.coord
                                            distance = calculate_distance(anion_coord, pcat)
                                            # print(distance)
                                            if distance <= 5.5:
                                                resid_no = (anion.get_parent().get_id()[1], atmm.get_parent().get_id()[1])
                                                interaction = (anion, atmm, distance)
                                                if resid_no not in sltbrd_id_p:
                                                    sltbrd_id_p.append(resid_no)
                                                    sltbrd_p.append(interaction)
                                                # result=print('Salt Bridge interactions between', atom.get_name(), 'of ', atom.get_parent().get_resname(),atom.get_parent().get_id()[1], atom.get_parent().get_parent().get_id(),'and  ' , 
                                                #                atmm.get_id(), atmm.get_parent().get_resname(),atmm.get_parent().get_id()[1], atmm.get_parent().get_parent().get_id()) 
                    else:
                        pass
        else:
            for lig in lig_values:
                for cn in range(len(ligand_anion)):
                    for atom in lig['pdb_atoms']:
                        if atom.name == ligand_anion[cn][1].strip():
                            anion=atom
                            anion_coord = anion.coord
                            # print(anion)
                            for res in list(structure.get_residues()):
                                if res.get_parent().get_id() == lig['chain'].get_id():
                                    # print(res.get_resname())
                                    cat_res=['HIS','ARG','LYS']
                                    if res.get_resname() in cat_res:
                                        # print(res.get_resname())
                                        for cnt in range(len(protein_cation)):
                                            for atmm in list(res.get_atoms()):
                                                if atmm.name == protein_cation[cnt][1].strip():
                                                    pcat= atmm.coord
                                                    distance = calculate_distance(pcat,anion_coord)
                                                    # print(distance)
                                                    if distance <= 5.5:
                                                        resid_no = (anion.get_parent().get_id()[1], atmm.get_parent().get_id()[1])
                                                        interaction = (anion, atmm, distance)
                                                        if resid_no not in sltbrd_id_p:
                                                            sltbrd_id_p.append(resid_no)
                                                            sltbrd_p.append(interaction)
                                                    # result=print('Salt Bridge interactions between', atom.get_name(), 'of ', atom.get_parent().get_resname(),atom.get_parent().get_id()[1], atom.get_parent().get_parent().get_id(),'and  ' , 
                                                    #                atmm.get_id(), atmm.get_parent().get_resname(),atmm.get_parent().get_id()[1], atmm.get_parent().get_parent().get_id()) 
                        else:
                            pass           
    return sltbrd_p

def add_to_list(list1, list2):
    if not list1:  # If list1 is empty
        return list2
    if not list2:  # If list2 is empty
        return list1
    return list1 + list2

def run(pdb_file, mol_file, ligand_id,ligand_id_mol2, dirct,structure):

    data = rg.detect_rings_from_pdb(pdb_file,structure)
    prot_values = data.get('Amino Acid')
    lig_values = data.get('Ligand')

    # Run functions in parallel using concurrent.futures
    #import concurrent.futures
    #with concurrent.futures.ThreadPoolExecutor() as executor:
        #protein_cation = executor.submit(pc.run, pdb_file).result()
        #ligand_cation = executor.submit(lc.run, pdb_file, ligand_id, ligand_id_mol2,dirct).result()
        #protein_anion = executor.submit(pa.run, pdb_file).result()
        #ligand_anion = executor.submit(la.run, pdb_file, ligand_id,ligand_id_mol2, dirct).result()
    protein_cation = pc.run(pdb_file)
    ligand_cation = lc.run(pdb_file, ligand_id, ligand_id_mol2, dirct)
    protein_anion = pa.run(pdb_file)
    ligand_anion = la.run(pdb_file, ligand_id, ligand_id_mol2, dirct)

    lig_prot_cation = lig_cation_prot_anion(prot_values, lig_values, protein_anion, ligand_cation, structure)
    prot_lig_cation = prot_cation_lig_anion(prot_values, lig_values, protein_cation, ligand_anion, structure)

    # Use list comprehension to concatenate lists
    saltbridge = [item for sublist in [lig_prot_cation, prot_lig_cation] for item in sublist]

    return saltbridge
