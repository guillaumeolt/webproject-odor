"""
Script to read and write molecules
"""
import os
from os import listdir, path
import argparse
import time, threading, queue, datetime
import sys
from io import StringIO
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D, MolToFile

from .conformers import generate_conformers


def load_multiple_molecule_smi(supplier, file):
    """
    Load Multiples molecules from a smile file
    """
    file_content = open(file).read()
    list_name = []
    res = []
    for mol in file_content.split("\n"):
        if len(mol) > 3:
            list_name.append(mol.split(" ")[1])
    k = 0
    for mol in supplier:
        if mol is not None:
            mol.SetProp("_Name", list_name[k])
            res.append(mol)
        else:
            print("ERROR LOAD : on Loading %s" % list_name[k])

        k += 1
    return res


def load_multiple_molecule_sdf(supplier, file):
    """
    Load Multiples molecules from a sdf file
    """
    file_content = open(file).read()
    list_name = []
    res = []
    for mol in file_content.split("$$$$\n"):
        list_name.append(mol.split("\n")[0])
    k = 0
    for mol in supplier:
        if mol is not None:
            mol.SetProp("_Name", list_name[k])
            res.append(mol)
        else:
            print("ERROR LOAD : on Loading %s" % list_name[k])
        k += 1
    return res


def load_multiple_mol2(file=None, sanitize=False):
    mols = []
    with open(file, 'r') as f:
        line = f.readline()
        while not f.tell() == os.fstat(f.fileno()).st_size:
            if line.startswith("@<TRIPOS>MOLECULE"):
                mol = []
                mol.append(line)
                line = f.readline()
                while not line.startswith("@<TRIPOS>MOLECULE"):
                    mol.append(line)
                    line = f.readline()
                    if f.tell() == os.fstat(f.fileno()).st_size:
                        mol.append(line)
                        break
                mol[-1] = mol[-1].rstrip()  # removes blank line at file end
                block = ",".join(mol).replace(',', '')
                m = Chem.MolFromMol2Block(block, sanitize=sanitize)
            mols.append(m)
    return (mols)

def load_all_mol_from_file(path_file_mol):
    """
    Load all mol from a specified path from sdf, smile and mol2 files
    """
    n_sdf = 0
    n_mol2 = 0
    n_smi = 0
    res = []
    file = path_file_mol
    print(file)
    name_mol = os.path.basename(path_file_mol)
    if ".sdf" in file:
        supplier = Chem.SDMolSupplier(file, sanitize=False, \
                                      removeHs=False, \
                                      strictParsing=True)
        if len(supplier) > 1:
            mols = load_multiple_molecule_sdf(supplier, file)
            n_sdf += len(mols)
            res += mols
        else:
            mol = supplier[0]
            if mol is not None:
                mol.SetProp("_Name", name_mol)
                res.append(mol)
                n_sdf += 1
            else:
                print("ERROR LOAD : on Loading %s" % name_mol)
    elif ".mol2" in file:
        mols = load_multiple_mol2(file)
        if mols is not None:
            n_mol2 += len(mols)
            res += mols
        else:
            print("ERROR LOAD : on Loading %s" % name_mol)
    elif ".smi" in file:
        supplier = Chem.SmilesMolSupplier(file, sanitize=False, \
                                          titleLine=False)
        if len(supplier) >= 1:
            mols = load_multiple_molecule_smi(supplier, file)
            n_smi += len(mols)
            res += mols
        else:
            mol = supplier[0]
            if mol is not None:
                mol.SetProp("_Name", name_mol)
                res.append(mol)
                n_smi += 1
            else:
                print("ERROR LOAD : on Loading %s" % name_mol)
    else:
        print("Format not recognized " + file)
    print("Loaded %d molecules, %d sdf %d mol2, %d smi" % (len(res), \
                                                           n_sdf, \
                                                           n_mol2, \
                                                           n_smi))
    return res
def load_all_mol_from(path_all_mol):
    """
    Load all mol from a specified path from sdf, smile and mol2 files
    """
    n_sdf = 0
    n_mol2 = 0
    n_smi = 0
    res = []
    for filename in listdir(path_all_mol):
        file = path.join(path_all_mol, filename)
        print(file)
        name_mol = filename.split('.')[0]
        if ".sdf" in file:
            supplier = Chem.SDMolSupplier(file, sanitize=False, \
                                          removeHs=False, \
                                          strictParsing=True)
            if len(supplier) > 1:
                mols = load_multiple_molecule_sdf(supplier, file)
                n_sdf += len(mols)
                res += mols
            else:
                mol = supplier[0]
                if mol is not None:
                    mol.SetProp("_Name", name_mol)
                    res.append(mol)
                    n_sdf += 1
                else:
                    print("ERROR LOAD : on Loading %s" % name_mol)
        elif ".mol2" in file:
            mols = load_multiple_mol2(file)
            if mols is not None:
                n_mol2 += len(mols)
                res += mols
            else:
                print("ERROR LOAD : on Loading %s" % name_mol)
        elif ".smi" in file:
            supplier = Chem.SmilesMolSupplier(file, sanitize=False, \
                                              titleLine=False)
            if len(supplier) >= 1:
                mols = load_multiple_molecule_smi(supplier, file)
                n_smi += len(mols)
                res += mols
            else:
                mol = supplier[0]
                if mol is not None:
                    mol.SetProp("_Name", name_mol)
                    res.append(mol)
                    n_smi += 1
                else:
                    print("ERROR LOAD : on Loading %s" % name_mol)
        else:
            print("Format not recognized " + file)
    #print("Loaded %d molecules, %d sdf %d mol2, %d smi" % (len(res), \
     #                                                      n_sdf, \
     #                                                      n_mol2, \
     #                                                      n_smi))
    return res

def save_2d_image_PNG_list(mol_list, path_output, name = "_Name"):
    """
    Save list of molecules to png files
    """
    mol_error=[]
    for mol in mol_list:
        mol.UpdatePropertyCache()
        Chem.rdDepictor.Compute2DCoords(mol)
        mol = Chem.RemoveHs(mol)
        path_out = os.path.join(path_output, mol.GetProp(name) + ".png")
        try:
            MolToFile(mol, path_out, size=(200, 200))
        except:
            mol_error.append(path_out)
            pass
    #print(mol_error)

def save_2d_image_SVG_list(mol_list, path_output, name = "_Name"):
    """
    Save list of molecules to png files
    """
    mol_error=[]
    for mol in mol_list:
        mol.UpdatePropertyCache()
        Chem.rdDepictor.Compute2DCoords(mol)
        mol = Chem.RemoveHs(mol)
        path_out = os.path.join(path_output, mol.GetProp(name) + ".svg")
        try:
            MolToFile(mol, path_out, size=(200, 200))
        except:
            mol_error.append(path_out)
            pass
    #print(mol_error)

def write_2d_pdb_list(mol_list, path_output, name="_Name"):
    """
    Write 2D PDB of list of molecules
    """
    i = 0
    for mol in mol_list:
        mol = Chem.AddHs(mol)
        Chem.rdDepictor.Compute2DCoords(mol)
        pdb_file = path.join(path_output, mol.GetProp("_Name") + ".pdb")
        writer_pdb = Chem.PDBWriter(pdb_file)
        writer_pdb.write(mol)

def write_3d_pdb_list(mol_list, path_output, name = "_Name"):
    """
    Write 3D PDB of list of molecules
    """
    i=0
    for mol in mol_list:
        i+=1
        mol, param = generate_conformers(mol,numberConf = 1)
        pdb_file = path.join(path_output, mol.GetProp("_Name") + ".pdb")
        writer_pdb = Chem.PDBWriter(pdb_file, flavor=4) # flavor 4 to Write CONECT records in both directions | flavor 8 Don't use multiple CONECTs to encode bond order
        writer_pdb.write(mol)