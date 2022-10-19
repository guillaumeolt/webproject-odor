"""
Conformer generation.
"""

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

# Génération des conformeres et clusterisation

#Infos :
# https://rdkit.readthedocs.io/en/latest/Cookbook.html

def generate_conformers(molecule, numberConf, pruneThresh=1,\
                                              n_conf_thresh=50,\
                                              useRandomCoords=False):
    """
    # Parameter of conf =
    rdkit.org/docs/cppapi/structRDKit_1_1DGeomHelpers_1_1EmbedParameters.html
    Generate conformers with ETKDGv3 parameters
    Prune Generate conformers for those with an RMSD of 1
    Increase RMSD threshold of 0.1 until the final number of comformers
    is under 50 because it is not too high for the follow up of the processus.
    In case the molecule is too big the conformation geeneration is repeated
    by using random coordinates as a starting point instead of using a distance
    geometry embedding. A fix number of conformer is set to 50 but it can be
    a long process if the molecule is big (~1h for DB13928).
    """
    print("Conformer Generation {}".format(molecule.GetProp("_Name")))
    if useRandomCoords == False:
        molecule = Chem.AddHs(molecule)
        param = AllChem.ETKDGv3()
        param.useSmallRingTorsions = True
        param.randomSeed = 123
        param.pruneRmsThresh = float(pruneThresh)
        param.onlyHeavyAtomsForRMS = True
        param.numThreads = 0 #use all threads of the computer
        ids = AllChem.EmbedMultipleConfs(molecule,\
                                         numConfs=int(numberConf),\
                                         params=param)
        while len(ids) > n_conf_thresh:
            offset = 0.1 * (len(ids)//n_conf_thresh)
            param.pruneRmsThresh += offset
            print("{} Number of conformers > {} :{:^4}| Updating pruneRmsThresh to {}"\
            .format(molecule.GetProp("_Name"), n_conf_thresh, len(ids), param.pruneRmsThresh))
            ids = AllChem.EmbedMultipleConfs(molecule,\
                                             numConfs=numberConf,\
                                             params=param)
            if len(ids) > 50 and param.pruneRmsThresh >= 2:
                ids = AllChem.EmbedMultipleConfs(molecule,\
                                                 numConfs=n_conf_thresh,\
                                                 params=param)
                return molecule
    else:
        print("WARNING :Molecule may be too big using random coords as as a starting point instead of using a distance geometry embedding. Molecule " + molecule.GetProp("_Name"))
        molecule = Chem.AddHs(molecule)
        param = AllChem.ETKDGv3()
        param.useSmallRingTorsions = True
        param.randomSeed = 123
        param.pruneRmsThresh = float(pruneThresh)
        param.onlyHeavyAtomsForRMS = True
        param.numThreads = 0 #use all threads of the computer
        param.useRandomCoords=True
        ids = AllChem.EmbedMultipleConfs(molecule,\
                                         numConfs=1,\
                                         params=param)
    return molecule, param.pruneRmsThresh

def get_conformer_energies(mol, mmffVariant="MMFF94"):
    """
    Calculate conformer energies with MMFF94 forcefield (RDKit default).
    Parameters
    ----------
    mol : RDKit Mol
        Molecule.
    Returns
    -------
    energies : array_like
    """
    mmff_props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant=mmffVariant)
    energies = []
    for conf in mol.GetConformers():
        ff = AllChem.MMFFGetMoleculeForceField(mol,\
                                               mmff_props,\
                                               confId=conf.GetId())
        energy = ff.CalcEnergy()
        energies.append(energy)
        #print(confEs)
        if mol.HasProp("_ConfEnergies"):
            mol.SetProp("_ConfEnergies",\
                        "%s|%d:%.4f"%(mol.GetProp("_ConfEnergies"),\
                                      conf.GetId(),\
                                      energy))
        else:
            mol.SetProp("_ConfEnergies","%d:%.4f"%(conf.GetId(),energy))
    energies = np.asarray(energies, dtype=float)
    return energies

# region ConformerGen
# Return the list of mol and a list of list containing for each
# molecule the indices of the selected conformers
def generate_conformer_list(mol_list, n_conf, pruneRmsThresh, n_conf_thresh):
    """
    Generate conformers for a list of molecules
    """
    res = []
    listPruneRMS = []
    for mol in mol_list:
        mol, lastPruneRMSD = conformers.generate_conformers(mol, \
                                                            n_conf, \
                                                            pruneRmsThresh, \
                                                            n_conf_thresh)
        if mol.GetNumConformers() == 0:
            print("ERROR CONF : No conformers generated for molecule " + mol.GetProp('_Name') + \
                  "trying with random coords as starting point")
            mol, lastPruneRMSD = conformers.generate_conformers(mol, \
                                                                n_conf, \
                                                                pruneRmsThresh, \
                                                                n_conf_thresh, \
                                                                useRandomCoords=True)
            if mol.GetNumConformers() == 0:
                raise ValueError("ERROR CONF : No conformers generated for molecule " + mol.GetProp('_Name'))
        try:
            conformers.get_conformer_energies(mol)
        except AttributeError:
            print("ERROR CONF_NRG : Getting conformers energies molecule for " + mol.GetProp('_Name'))
        res.append(mol)
        listPruneRMS.append(lastPruneRMSD)
        print("Conformer Generated for %s , %d conf generated" % \
              (mol.GetProp("_Name"), mol.GetNumConformers()))
    return res, listPruneRMS

