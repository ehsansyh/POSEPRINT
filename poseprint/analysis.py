"""Run every detector over a set of complexes and tabulate the result.

Carried over from the working copy of POSEPRINT with its detection logic unchanged.
See build_poseprint_repo.py in the CLINCH working tree for the exact edits that were
applied, which are limited to module naming, relative imports, dead imports and the
removal of one hard-coded scratch path.
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Jul 26 13:52:34 2024

@author: ehsan
"""
import os
import tempfile
from . import hbond as interactions
from . import rings as detect_residue_ring_update
from . import pi_stacking as ring_compares
from . import pi_cation
from . import hydrophobic as hydrophobic_interactions
from . import salt_bridge as sb
from . import halogen as hl
from . import water_bridge as wb
import pandas as pd
from Bio import PDB
from Bio.PDB import *
import warnings


def check_is_lig(structure):
    for res in list(structure.get_residues()):
        if res.get_id()[0] != ' ' and res.get_id()[0] != 'W':
            print(res)
            break
def ext_lig(structure):
    ligid=[]
    ligid_mol2=[]
    for res in list(structure.get_residues()):
        if res.get_id()[0] != ' ' and res.get_id()[0] != 'W':
            # print(res.get_resname())
            ligand_id= f"{res.get_resname()}{res.get_id()[1]}"
            ligand_id_mol2 = f"{res.get_resname()} {res.get_id()[1]}"
            ligid.append(ligand_id)
            ligid_mol2.append(ligand_id_mol2)
    return ligid , ligid_mol2

def run_analyze(pdb_files, mol2_files, dirctf, save_dir=None, scratch=None,
                verbose=True):
    """Analyse each complex and return one table per ligand.

    scratch is the directory the detectors write their intermediate ligand files
    into. The original hard-coded one path for every call, which also made two
    parallel runs overwrite each other, so it is a parameter now and the caller
    is expected to give each run its own.

    save_dir keeps the original behaviour of writing a log and a CSV per ligand.
    Leave it as None to get the tables back without touching the disk.
    """
    results = []
    if scratch is None:
        scratch = tempfile.mkdtemp(prefix='poseprint_')
    dirct = os.path.join(scratch, '')
    os.makedirs(dirct, exist_ok=True)
    if save_dir is not None:
        save_dir = os.path.join(save_dir, '')
        os.makedirs(save_dir, exist_ok=True)
    for pdb, mol in zip(pdb_files, mol2_files):
        if verbose:
            print(pdb)
        pdb_file = os.path.join(dirctf, pdb)
        mol_file = os.path.join(dirctf, mol)
        parser = PDBParser()
        structure = parser.get_structure('complex', pdb_file)
        ligand, ligandm= ext_lig(structure)
        #print('\n==========================================================')
        #print(f"               {pdb.split('.pdb')[0]}  in process   ...    ")
        #print('==========================================================\n')
        warnings.filterwarnings("ignore", category=UserWarning)
        for f in range(len(ligand)):
            try :
                ligand_id = ligand[f]
                ligand_id_mol2= ligandm[f]
                if any(x in ligand_id for x in ["ZN", "CD", "FE", "MG", "CA", "NA", "NI", ]):
                    pass
                else:
                # if __name__ == "__main__":
                    #print(f"         Ligand name -->   {ligand_id_mol2}  \n")
                    hbonds = interactions.main(pdb_file,structure)
                    res_no = len(list(structure.get_residues()))
                    starter_res= list(structure.get_residues())[0].get_id()[1]
                    data_per_residue_type = detect_residue_ring_update.detect_rings_from_pdb(pdb_file,structure)
                    tshape, parallel, f_f = ring_compares.run(mol_file,pdb_file,ligand_id,ligand_id_mol2,dirct,structure)
                    lig_prot_cat, prot_lig_cat= pi_cation.run(pdb_file,mol_file,ligand_id,ligand_id_mol2,dirct,structure)
                    hydro = hydrophobic_interactions.hydrophobic(pdb_file, mol_file,ligand_id,ligand_id_mol2,dirct,structure)
                    salt_bridge = sb.run(pdb_file, mol_file,ligand_id,ligand_id_mol2,dirct ,structure)
                    halogen = hl.main(pdb_file,structure)
                    water_bridge= wb.main(pdb_file,structure)
                    model_Data= pd.DataFrame(columns=['Residues','H-bonds','Salt_bridge', 'Water_bridge', 'Pi-Stacking (T-shape)', 'Pi-Stacking(F-F)', 'Pi-Stacking(parallel)', 'Pi-cation', 'Hydrophobic interaction', 'Halogen bonds'],
                                             index=range(res_no))
                    
                    for i in range(res_no):
                        residue = list(structure.get_residues())[i]
                        res_name = residue.get_resname()
                        res_number = residue.get_id()[1]
                        model_Data.loc[i, 'Residues'] = f"{res_name} {res_number}"
                    
                    hbond_counts = {}
                    distance_counts = {}
                    angle_counts = {}
                    
                    for key in hbonds.keys():
                        for i in range(len(key) - 1):
                            res_id = '{} {}'.format(key[i], key[i+1])
                            hbond_counts[res_id] = hbond_counts.get(res_id, 0) + 1
                            for subelem in hbonds[key]:
                                distance_counts[res_id] = distance_counts.get(res_id, []) + [round(subelem[1], 2)]
                                angle_counts[res_id] = angle_counts.get(res_id, []) + [round(subelem[2], 2)]
            
                    for i, res_id in enumerate(model_Data['Residues']):
                        # resname, res_id = res_id.split()
                        # res_id = int(res_id)
                        model_Data.loc[i, 'H-bonds'] = hbond_counts.get(res_id, 0)
                        model_Data.loc[i, 'H-bond Distance'] = distance_counts.get(res_id, [0])[0] if distance_counts.get(res_id) else 0
                        model_Data.loc[i, 'H-bond Angle'] = angle_counts.get(res_id, [0])[0] if angle_counts.get(res_id) else 0
            
                    
                    salt_bridge_counts = {}
                    salt_bridge_distance = {}
                    try:
                        for key in salt_bridge:
                            for i in range(len(key)-1):
                                res_id= f"{key[i].get_parent().get_resname()} {key[i].get_parent().get_id()[1]}"
                                salt_bridge_counts[res_id] = salt_bridge_counts.get(res_id, 0) + 1
                                salt_bridge_distance[res_id] = salt_bridge_distance.get(res_id, []) + [round(subelem[2], 2)]
                    except Exception as e:
                        print(f"[Warning] Error in salt bridge analysis for {ligand_id}: {e}")
                        salt_bridge_counts = {}
                        salt_bridge_distance = {}
                    for i, res_id in enumerate(model_Data['Residues']):
                        model_Data.loc[i, 'Salt_bridge'] = salt_bridge_counts.get(res_id, 0)
                        model_Data.loc[i, 'Salt-Bridge Distance'] = salt_bridge_distance.get(res_id, [0])[0] if salt_bridge_distance.get(res_id) else 0
                    
                        
                    water_bridge_counts = {}
                    wd_distance_counts = {}
                    wd_angle_counts = {}
                    wa_distance_counts = {}
                    wa_angle_counts = {}
                    
                    for key in water_bridge.keys():
                        for i in range(len(key) - 1):
                            res_id = '{} {}'.format(key[i], key[i+1])
                            water_bridge_counts[res_id] = water_bridge_counts.get(res_id, 0) + 1
                            for subelem in water_bridge[key]:
                               wd_distance_counts[res_id] = wd_distance_counts.get(res_id, []) + [round(subelem[1], 2)]
                               wd_angle_counts[res_id] = wd_angle_counts.get(res_id, []) + [round(subelem[3], 2)]
                               wa_distance_counts[res_id] = wa_distance_counts.get(res_id, []) + [round(subelem[1], 2)]
                               wa_angle_counts[res_id] = wa_angle_counts.get(res_id, []) + [round(subelem[4], 2)]
            
                    for i, res_id in enumerate(model_Data['Residues']):
                        model_Data.loc[i, 'Water_bridge'] = water_bridge_counts.get(res_id, 0)
                        model_Data.loc[i, 'WD Distance'] = wd_distance_counts.get(res_id, [0])[0] if wd_distance_counts.get(res_id) else 0
                        model_Data.loc[i, 'WA Distance'] = wa_distance_counts.get(res_id, [0])[0] if wa_distance_counts.get(res_id) else 0
                        model_Data.loc[i, 'WD Angle'] = wd_angle_counts.get(res_id, [0])[0] if wd_angle_counts.get(res_id) else 0
                        model_Data.loc[i, 'WA Angle'] = wa_angle_counts.get(res_id, [0])[0] if wa_angle_counts.get(res_id) else 0
                            
                    tshape_counts = {}
                    parallel_counts = {}
                    ff_counts = {}
                    distance_values = {}
                    angle_values = {}
                    
                    for key in tshape:
                        for i in range(len(key) - 3):
                            res_id = '{} {}'.format(key[i], key[i+1])
                            tshape_counts[res_id] = tshape_counts.get(key[1], 0) + 1
                            distance_values[res_id] = round(key[4],2)
                            angle_values[res_id] = round(key[5],2)
                    
                    for key in parallel:
                        for i in range(len(key) - 3):
                            res_id = '{} {}'.format(key[i], key[i+1])
                            parallel_counts[res_id] = parallel_counts.get(key[1], 0) + 1
                            distance_values[res_id] = round(key[4],2)
                            angle_values[res_id] = round(key[5],2)
                    
                    for key in f_f:
                        for i in range(len(key) - 3):
                            res_id = '{} {}'.format(key[i], key[i+1])
                            ff_counts[res_id] = ff_counts.get(key[1], 0) + 1
                            distance_values[res_id] = round(key[4],2)
                            angle_values[res_id] = round(key[5],2)
                    
                    for i, res_id in enumerate(model_Data['Residues']):
            
                        model_Data.loc[i, 'Pi-Stacking (T-shape)'] = tshape_counts.get(res_id, 0)
                        model_Data.loc[i, 'Pi-Stacking(parallel)'] = parallel_counts.get(res_id, 0)
                        model_Data.loc[i, 'Pi-Stacking(F-F)'] = ff_counts.get(res_id, 0)
                        model_Data.loc[i, 'pi-pi Distance'] = distance_values.get(res_id, 0)
                        model_Data.loc[i, 'pi-pi Angle'] = angle_values.get(res_id, 0)
                    
                    
                    pi_cat_counts = {}
                    pi_cat_distance = {}
                    
                    pication = {}
                    for item in lig_prot_cat + prot_lig_cat:
                        key = (item[0], item[1], item[2], item[3])
                        if key not in pication or item[4] < pication[key][4]:
                            pication[key] = item
                    
                    pication = list(pication.values())
            
                    for key in pication:
                        for i in range(len(key) - 2):
                            res_id = '{} {}'.format(key[i], key[i+1])
                            pi_cat_counts[res_id] = pi_cat_counts.get(res_id, 0) + 1
                            pi_cat_distance[res_id] = round(key[4], 2)
                
                    for i, res_id in enumerate(model_Data['Residues']):
            
                        model_Data.loc[i, 'Pi-cation'] = pi_cat_counts.get(res_id, 0)
                        model_Data.loc[i, 'Pi-cation Distance'] = pi_cat_distance.get(res_id, 0)
                                
                    
                    result_hp = {}
                    for item in hydro:
                        key = (item[0], item[1], item[2], item[3])
                        if key not in result_hp or item[4] < result_hp[key][4]:
                            result_hp[key] = item
                    
                    result_hp = list(result_hp.values())
                    
                    hydro_counts = {}
                    hydro_distance = {}
                    for key in result_hp:
                        for i in range(len(key) - 2):
                            res_id = '{} {}'.format(key[i], key[i+1])
                            hydro_counts[res_id] = hydro_counts.get(res_id, 0) + 1
                            hydro_distance[res_id] = round(key[4], 2)
            
                    
                    for i, res_id in enumerate(model_Data['Residues']):
            
                        model_Data.loc[i, 'Hydrophobic interaction'] = hydro_counts.get(res_id, 0)
                        model_Data.loc[i, 'Hydrophobic Distance'] = hydro_distance.get(res_id, 0)
                        
                        
                    halogen_counts = {}
                    hl_distance_counts = {}
                    hl_angle_counts = {}
                    
                    for key in halogen.keys():
                        for i in range(len(key) - 1):
                            res_id = '{} {}'.format(key[i], key[i+1])
                            halogen_counts[res_id] = halogen_counts.get(res_id, 0) + 1
                            for subelem in halogen[key]:
                               hl_distance_counts[res_id] = hl_distance_counts.get(res_id, []) + [round(subelem[1], 2)]
                               hl_angle_counts[res_id] = hl_angle_counts.get(res_id, []) + [round(subelem[2], 2)]
            
                    for i, res_id in enumerate(model_Data['Residues']):
                        model_Data.loc[i, 'Halogen bonds'] = halogen_counts.get(res_id, 0)
                        model_Data.loc[i, 'Halogen Distance'] = hl_distance_counts.get(res_id, [0])[0] if hl_distance_counts.get(res_id) else 0
                        model_Data.loc[i, 'Halogen Angle'] = hl_angle_counts.get(res_id, [0])[0] if hl_angle_counts.get(res_id) else 0
                    
                        
                    _log = (save_dir + f"results_{pdb.split('.pdb')[0]}_{ligand_id}.log"
                            if save_dir else os.devnull)
                    with open(_log, 'w') as result:
                        result.write("Hydrogen bond interactions: \n========================HITMER========================\n")
                        for interaction_type, atoms in hbonds.items():
                            ligand_residue,ligand_id1, protein_residue, protein_id = interaction_type
                            for at in atoms:
                                result.write(f"Hydrogen bond interaction between {ligand_residue} {ligand_id1} and {protein_residue} {protein_id}: H-A Distance = {round(at[1],2)}, D-H-A Angle = {round(at[2],2)} --> Interaction between atoms {at[0]} \n")
                        result.write('\nSalt Bridge Interactions : \n=========================HITMER=========================\n')
                        for atom in salt_bridge:
                            result.write(f"Salt Bridge interactions between {atom[0].get_parent().get_resname()} {atom[0].get_parent().get_id()[1]} {atom[0].get_parent().get_parent().get_id()} and {atom[1].get_parent().get_resname()} {atom[1].get_parent().get_id()[1]} {atom[1].get_parent().get_parent().get_id()}: Distance = {round(atom[2],2)}  \n")
                        result.write("\nWater Bridge interactions: \n========================HITMER==========================\n")
                        for interaction_type, atoms in water_bridge.items():
                            ligand_residue,ligand_id1, water, water_id, protein_residue, protein_id = interaction_type
                            for dt in atoms:
                                 result.write(f"Watger Bridge interaction between {ligand_residue} {ligand_id1}, Water {water_id} and {protein_residue} {protein_id}: WD Distance = {round(dt[1],2)} , WA Distance = {round(dt[2],2)} , WD Angle = {round(dt[3],2)}, WA Angle = {round(dt[4],2)}\n")
                        result.write("\nPi-Pi interactions: \n========================HITMER==========================\n")
                        for t in range(len(tshape)):
                            result.write(f"T-shape Pi stacking interactions between {tshape[t][0]} {tshape[t][1]} and {tshape[t][2]} {tshape[t][3]}: center distance = {tshape[t][4]:.2f} , angle = {tshape[t][5]:.2f}\n")
                        for z in range(len(f_f)):
                            result.write(f"F-F or F-E Pi stacking interactions between {f_f[z][0]} {f_f[z][1]} and {f_f[z][2]} {f_f[z][3]}: center distance = {f_f[z][4]:.2f} , angle = {f_f[z][5]:.2f}\n")
                        for r in range(len(parallel)):
                            result.write(f"Parallel Pi stacking interactions between {parallel[r][0]} {parallel[r][1]} and {parallel[r][2]} {parallel[r][3]}: center distance = {parallel[r][4]:.2f} , angle = {parallel[r][5]:.2f}\n")
                        result.write("\nPi-Cation interactions: \n========================HITMER==========================\n")
                        for t in pication:
                            result.write(f"Pi-Cation interactions between {t[0]} {t[1]} and {t[2]} {t[3]}: Distance = {round(t[4],2)}\n")
                        result.write("\nHydrophobic interactions: \n========================HITMER==========================\n")
                        for t in range(len(result_hp)):
                            result.write(f"Hydrophobic interactions between {result_hp[t][0]} {result_hp[t][1]} and {result_hp[t][2]} {result_hp[t][3]}: Distance = {round(result_hp[t][4],2)}\n")
                        result.write("\nHalogen Bonds: \n========================HITMER==========================\n")
                        for interaction_type, atoms in halogen.items():
                            ligand_residue,ligand_id1, protein_residue, protein_id = interaction_type
                            for ad in atoms: 
                                result.write(f"Halogen bond interaction between {ligand_residue} {ligand_id1} and {protein_residue} {protein_id}: Distance = {round(ad[0][1], 2)}, Angle = {round(ad[0][2], 2)}\n")
                    
                    model_Data = model_Data.astype({col: 'float' if col != 'Residues' else 'object' for col in model_Data.columns})
                    results.append({"pose": pdb, "ligand": ligand_id,
                                    "table": model_Data})
                    if save_dir:
                        model_Data.to_csv(save_dir + f"model_data_{pdb.split('.pdb')[0]}_{ligand_id}.csv", index=False)
                        if verbose:
                            print(f"Data for {pdb.split('.pdb')[0]} saved in {save_dir}")

            except AttributeError:
                pass
            except ValueError as a:
                print(f"Skipping {ligand_id} of {pdb} pdb file due to error: {a} ")
                continue
            except IndexError as e:
                if str(e) == "Cannot write an AtomGroup with 0 atoms":
                    print(f"Skipping {ligand_id} of {pdb} pdb file due to error: {e}")
                    continue
                else:
                    raise
            except Exception as e:
                print(f"An error occurred: {e}")
                raise
        ligid=[]
        ligid_mol2=[]
    return results
