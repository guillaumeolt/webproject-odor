import numpy as np

import umap
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import DataStructs

import pickle
## UMAP ##
def get_umap_chem_odor(list_smile):
    fps = smi2fp(list_smile)
    X_fp = np.array(fps)
    mapper = umap.UMAP(random_state=2021).fit_transform(X_fp)
    #print(mapper)
    return(mapper)
def write_umap_chem_odor(mapper, path):
    with open(path, 'wb') as outp:
        pickle.dump(mapper, outp, pickle.HIGHEST_PROTOCOL)
def load_umap_chem_odor(path):
    with open(path, 'rb') as pickle_file:
        mapper = pickle.load(pickle_file)
    return mapper

def mol2svg(mol):
    d2d = rdMolDraw2D.MolDraw2DSVG(200, 100)
    try:
        d2d.DrawMolecule(mol)
        d2d.FinishDrawing()
        return d2d.GetDrawingText()
    except:
        return(d2d.GetDrawingText())

def mol2fp(mol, radi=2, nBits=1024):
    arr = np.zeros((1,))
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radi, nBits=nBits)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr

def smi2fp(Smiles_R) :
    mols = [Chem.MolFromSmiles(smi) for smi in Smiles_R]
    fpvect = [AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024) for mol in mols]
    fps_string = [j  for j in [fpvect[i].ToBitString() for i in range(len(fpvect))]]
    fps = []
    for i in fps_string :
        fps_int = []
        for j in i :
            fps_int.append(int(j))
        fps.append(fps_int)
    #X_fp = np.array(fps)
    return(fps)