"""
Python script to generate descriptors from 3D molecules prepared with our
protocol. The script allow to generate RDKit descriptors or mordred descriptors.
This protocol take as assumption that you already standardize and compute isomers and conformers
of your molecules.
"""
import sys
import argparse
import time
import os
from io import StringIO
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem import Descriptors3D
from rdkit.ML.Descriptors.Descriptors import DescriptorCalculator
from rdkit.ML.Descriptors import MoleculeDescriptors
from mordred import Calculator, descriptors, is_missing
import re
# PATH

# FUNCTIONS
def get_molecules(path_dir):
    allMol = load_all_mol_from(path_dir)
    res_list = []
    mol_n=1
    for mol in allMol:
        mol_name = mol.GetProp("_Name")
        #print("{:*^50}".format("PROCESSING MOLECULE {} ({}/{})".format(mol_name,mol_n,len(allMol))))
        mol_n += 1
        #Standardize
        try:
            res_list.append(standardize(mol))
        except:
            #print("ERROR MOL: Standardize molecule " + mol_name)
            continue
    return res_list

def embed_molecule(mol_list):
    res = []
    name_list = []
    for mol in mol_list:
        #print(Chem.MolToSmiles(mol, isomericSmiles=True))
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        mol = Chem.RemoveHs(mol)
    return res
def get_decriptors_rdkit_fr(mol_list, list_fr = "All", compute_vec_desc = False):
    """
    Compute RDKit 2D  descriptors for a lsit of molecules and return a dataframe of the
    results. Only fragments descriptors
    """
    list_descriptors_mol = []
    list_names = []
    # 2D Descriptors
    #list_fr = []
    #for desc in Descriptors._descList:
    #    if "fr_" in desc[0]:
    #        list_fr.append(desc)
    if list_fr == "All":
        des_list_2D = [x[0] for x in Descriptors._descList[123:]]
    else:
        des_list_2D = []
        for desc in Descriptors._descList:
            if desc[0] in list_fr:
                des_list_2D.append(desc[0])

    des_list_all = des_list_2D
    calculator = MoleculeDescriptors.MolecularDescriptorCalculator(des_list_2D)
    n_mol = 0
    for mol in mol_list:
        # mol precessed
        n_mol += 1
        #print("MOL PROCESSED : {}/{}".format(n_mol,len(mol_list)))
        descriptors_mol = []
        # 2D Descriptors
        try:
            descriptors_2D = calculator.CalcDescriptors(mol)
        except:
            #print("ERROR MOL: 2D descripteurs RDKit. " + mol.GetProp("_Name"))
            continue
        if len(descriptors_2D) != len(des_list_2D):
            #print("ERROR MOL: 2D descriptors generated not same lengths as number available descriptors. "+ mol.GetProp("_Name"))
            continue
        descriptors_mol = list(descriptors_2D)
        list_descriptors_mol.append(descriptors_mol)
        list_names.append(mol.GetProp("_Name"))
    #print("FIRST DESCRIPTORS COMPUTED")
    df = pd.DataFrame(list_descriptors_mol, columns=des_list_all, index = list_names)

    return df

def get_decriptors_rdkit(mol_list, compute_vec_desc = False):
    """
    Compute RDKit 2D and 3D descriptors for a lsit of molecules and return a dataframe of the
    results
    """
    list_descriptors_mol = []
    list_names = []
    # 2D Descriptors
    des_list_2D = [x[0] for x in Descriptors._descList]
    des_list_3D = ["Asphericity",\
                   "Eccentricity",\
                   "InertialShapeFactor",\
                   "NPR1",\
                   "NPR2",\
                   "PMI1",\
                   "PMI2",\
                   "PMI3",\
                   "RadiusOfGyration",\
                   "SpherocityIndex"]
    des_list_all = des_list_2D + des_list_3D
    calculator = MoleculeDescriptors.MolecularDescriptorCalculator(des_list_2D)
    
    # Vector descriptors
    vector_descr_2d = {"AUTOCORR2D":'CalcAUTOCORR2D'}
    vector_descr_3d = {"AUTOCORR3D":'CalcAUTOCORR3D',\
                       "RDF":'CalcRDF',\
                       "MORSE":'CalcMORSE',\
                       "WHIM":'CalcWHIM',\
                       "GETAWAY":'CalcGETAWAY'}
    n_mol = 0
    for mol in mol_list:
        # mol precessed
        n_mol += 1
        #print("MOL PROCESSED : {}/{}".format(n_mol,len(mol_list)))
        descriptors_mol = []
        # 2D Descriptors
        try:
            descriptors_2D = calculator.CalcDescriptors(mol)
        except:
            #print("ERROR MOL: 2D descripteurs RDKit. " + mol.GetProp("_Name"))
            continue
        if len(descriptors_2D) != len(des_list_2D):
            #print("ERROR MOL: 2D descriptors generated not same lengths as number available descriptors. "+ mol.GetProp("_Name"))
            continue
        # 3D Descriptors
        descriptors_3D = []
        try:
            Asphericity = Descriptors3D.Asphericity(mol)
        except:
            Asphericity = "NA"
            #print("WARNING DESC: 3D descriptor Asphericity coulnd be computed. mol " + mol.GetProp("_Name"))
        descriptors_3D.append(Asphericity)
        try:
            Eccentricity = Descriptors3D.Eccentricity(mol)
        except:
            Eccentricity = "NA"
            #print("WARNING DESC: 3D descriptor Eccentricity coulnd be computed. mol " + mol.GetProp("_Name"))
        descriptors_3D.append(Eccentricity)
        try:
            InertialShapeFactor = Descriptors3D.InertialShapeFactor(mol)
        except:
            InertialShapeFactor = "NA"
            #print("WARNING DESC: 3D descriptor InertialShapeFactor coulnd be computed. mol " + mol.GetProp("_Name"))
        descriptors_3D.append(InertialShapeFactor)
        try:
            NPR1 = Descriptors3D.NPR1(mol)
        except:
            NPR1 = "NA"
            #print("WARNING DESC: 3D descriptor NPR1 coulnd be computed. mol " + mol.GetProp("_Name"))
        descriptors_3D.append(NPR1)
        try:
            NPR2 = Descriptors3D.NPR2(mol)
        except:
            NPR2 = "NA"
            #print("WARNING DESC: 3D descriptor NPR2 coulnd be computed. mol " + mol.GetProp("_Name"))
        descriptors_3D.append(NPR2)
        try:
            PMI1 = Descriptors3D.PMI1(mol)
        except:
            PMI1 = "NA"
            #print("WARNING DESC: 3D descriptor PMI1 coulnd be computed. mol " + mol.GetProp("_Name"))
        descriptors_3D.append(PMI1)
        try:
            PMI2 = Descriptors3D.PMI2(mol)
        except:
            PMI2 = "NA"
            #print("WARNING DESC: 3D descriptor PMI2 coulnd be computed. mol " + mol.GetProp("_Name"))
        descriptors_3D.append(PMI2)
        try:
            PMI3 = Descriptors3D.PMI3(mol)
        except:
            PMI3 = "NA"
            #print("WARNING DESC: 3D descriptor PMI3 coulnd be computed. mol " + mol.GetProp("_Name"))
        descriptors_3D.append(PMI3)
        try:
            RadiusOfGyration = Descriptors3D.RadiusOfGyration(mol)
        except:
            RadiusOfGyration = "NA"
            #print("WARNING DESC: 3D descriptor RadiusOfGyration coulnd be computed. mol " + mol.GetProp("_Name"))
        descriptors_3D.append(RadiusOfGyration)
        try:
            SpherocityIndex = Descriptors3D.SpherocityIndex(mol)
        except:
            SpherocityIndex = "NA"
            #print("WARNING DESC: 3D descriptor SpherocityIndex coulnd be computed. mol " + mol.GetProp("_Name"))
        descriptors_3D.append(SpherocityIndex)
        #
        descriptors_mol = list(descriptors_2D) + descriptors_3D
        if len(descriptors_mol) != 218:
            pass
            #print("ERROR MOL: 2D 3D descriptors generated not length available descriptors. "+ mol.GetProp("_Name"))
        else:
            list_descriptors_mol.append(descriptors_mol)
            list_names.append(mol.GetProp("_Name"))
    
    #print("FIRST DESCRIPTORS COMPUTED")
    df = pd.DataFrame(list_descriptors_mol, columns=des_list_all, index = list_names)

    if compute_vec_desc:
        n_mol = 0
        for mol in mol_list:
            # mol precessed
            n_mol += 1
            #print("MOL PROCESSED : {}/{}".format(n_mol,len(mol_list)))
            ## VECTORS 2D
            for name_d, v_d in vector_descr_2d.items():
                for n, num in enumerate(getattr(Chem.rdMolDescriptors, v_d)(mol)):
                    df.loc[mol.GetProp("_Name"), '{descr}_{k}'.format(descr=name_d, k=n)] = num
            ## VECTORS 3D
            for name_d, v_d in vector_descr_3d.items():
                for n, num in enumerate(getattr(Descriptors3D.rdMolDescriptors, v_d)(mol)):
                    df.loc[mol.GetProp("_Name"), '{descr}_{k}'.format(descr=name_d, k=n)] = num
        df = df.astype(float).fillna("NA")
  
    return df

def get_descriptors_mordred(mol_list):
    """
    Compute mordred 2D and 3D descriptors for a lsit of molecules and return a dataframe of the
    results
    """
    list_names=[mol.GetProp("_Name") for mol in mol_list]
    calc3D = Calculator(descriptors, ignore_3D=False)
    list_descriptors_mol = calc3D.pandas(mol_list)
    list_descriptors_mol = list_descriptors_mol.astype(float).fillna("NA")
    
    df = pd.DataFrame(list_descriptors_mol.values.tolist(), columns=list_descriptors_mol.columns, index = list_names)
    return df
def save_descriptors(mat_descriptors, path_output):
    mat_descriptors.to_csv(path_output)

#### MAIN ####
if __name__ == "__main__":
    ##Parser to deal with arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", nargs="?",\
           help="Path to the directory containing the molecules",\
           default=PATH_mol_sdf)
    parser.add_argument("-o",\
           help="Path to the output file to write the csv",\
           default=PATH_file_output)
    parser.add_argument("-name",\
           help="Path to the output file to write the csv",\
           default=NAME_file_output)
    args = parser.parse_args()
    # TODO: check 3D
    
    list_mol = get_molecules(args.d)
    n = len(list_mol)


    #list_mol = embed_molecule(list_mol)

    df_rdkit = get_decriptors_rdkit(list_mol)
    save_descriptors(df_rdkit,os.path.join(args.o,"mat_descriptors_"+args.name+"_rdkit_raw.csv"))
    
    #df_mordred = get_descriptors_mordred(list_mol)
    #save_descriptors(df_mordred,os.path.join(args.o,"mat_descriptors_"+args.name+"_mordred_raw.csv"))
    
   
