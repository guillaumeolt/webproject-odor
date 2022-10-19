"""
Standardize molecules
"""
import os
from os import listdir, path
from collections import defaultdict, Counter
from rdkit import Chem
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, \
    StereoEnumerationOptions
from rdkit.Chem.MolStandardize.normalize import Normalization
from rdkit.Chem.MolStandardize.standardize import Standardizer, copy
from rdkit.Chem.MolStandardize import rdMolStandardize

NORMALIZATIONS_perso = (
    # Opposite of #2.1 in InChI technical manual? Covered by RDKit Sanitization.
    Normalization('Nitro to N+(O-)=O', \
                  '[N,P,As,Sb;X3:1](=[O,S,Se,Te:2])=[O,S,Se,Te:3]>>[*+1:1]([*-1:2])=[*:3]'),
    Normalization('Sulfone to S(=O)(=O)', \
                  '[S+2:1]([O-:2])([O-:3])>>[S+0:1](=[O-0:2])(=[O-0:3])'),
    Normalization('Pyridine oxide to n+O-', \
                  '[n:1]=[O:2]>>[n+:1][O-:2]'),
    Normalization('Azide to N=N+=N-', \
                  '[*,H:1][N:2]=[N:3]#[N:4]>>[*,H:1][N:2]=[N+:3]=[N-:4]'),
    Normalization('Diazo/azo to =N+=N-', \
                  '[*:1]=[N:2]#[N:3]>>[*:1]=[N+:2]=[N-:3]'),
    Normalization('Sulfoxide to -S+(O-)-', \
                  '[!O:1][S+0;X3:2](=[O:3])[!O:4]>>[*:1][S+1:2]([O-:3])[*:4]'),
    # Equivalent to #1.5 in InChI technical manual
    Normalization('Phosphate to P(O-)=O', \
                  '[O,S,Se,Te;-1:1][P+;D4:2][O,S,Se,Te;-1:3]>>[*+0:1]=[P+0;D5:2][*-1:3]'),
    # Equivalent to #1.8 in InChI technical manual
    Normalization('C/S+N to C/S=N+', \
                  '[C,S;X3+1:1]([NX3:2])[NX3!H0:3]>>[*+0:1]([N:2])=[N+:3]'),
    # Equivalent to #1.8 in InChI technical manual
    Normalization('P+N to P=N+', \
                  '[P;X4+1:1]([NX3:2])[NX3!H0:3]>>[*+0:1]([N:2])=[N+:3]'),
    Normalization('Normalize hydrazine-diazonium', \
                  '[CX4:1][NX3H:2]-[NX3H:3][CX4:4][NX2+:5]#[NX1:6]>>[CX4:1][NH0:2]=[NH+:3][C:4][N+0:5]=[NH:6]'),
    # Equivalent to #1.3 in InChI technical manual
    Normalization('Recombine 1,3-separated charges', \
                  '[N,P,As,Sb,O,S,Se,Te;-1:1]-[A+0:2]=[N,P,As,Sb,O,S,Se,Te;+1:3]>>[*-0:1]=[*:2]-[*+0:3]'),
    Normalization('Recombine 1,3-separated charges', \
                  '[n,o,p,s;-1:1]:[a:2]=[N,O,P,S;+1:3]>>[*-0:1]:[*:2]-[*+0:3]'),
    Normalization('Recombine 1,3-separated charges', \
                  '[N,O,P,S;-1:1]-[a:2]:[n,o,p,s;+1:3]>>[*-0:1]=[*:2]:[*+0:3]'),
    Normalization('Recombine 1,5-separated charges', \
                  '[N,P,As,Sb,O,S,Se,Te;-1:1]-[A+0:2]=[A:3]-[A:4]=[N,P,As,Sb,O,S,Se,Te;+1:5]' + \
                  '>>[*-0:1]=[*:2]-[*:3]=[*:4]-[*+0:5]'),
    Normalization('Recombine 1,5-separated charges', \
                  '[n,o,p,s;-1:1]:[a:2]:[a:3]:[c:4]=[N,O,P,S;+1:5]>>[*-0:1]:[*:2]:[*:3]:[c:4]-[*+0:5]'),
    Normalization('Recombine 1,5-separated charges', \
                  '[N,O,P,S;-1:1]-[c:2]:[a:3]:[a:4]:[n,o,p,s;+1:5]>>[*-0:1]=[c:2]:[*:3]:[*:4]:[*+0:5]'),
    # Conjugated cation rules taken from Francis Atkinson's standardiser.
    # Those that can reduce aromaticity aren't included
    Normalization('Normalize 1,3 conjugated cation', \
                  '[N,O;+0!H0:1]-[A:2]=[N!$(*[O-]),O;+1H0:3]>>[*+1:1]=[*:2]-[*+0:3]'),
    Normalization('Normalize 1,3 conjugated cation', \
                  '[n;+0!H0:1]:[c:2]=[N!$(*[O-]),O;+1H0:3]>>[*+1:1]:[*:2]-[*+0:3]'),
    # Normalization('Normalize 1,3 conjugated cation',
    # '[N,O;+0!H0:1]-[c:2]:[n!$(*[O-]),o;+1H0:3]>>[*+1:1]=[*:2]:[*+0:3]'),
    Normalization('Normalize 1,5 conjugated cation', \
                  '[N,O;+0!H0:1]-[A:2]=[A:3]-[A:4]=[N!$(*[O-]),O;+1H0:5]>>[*+1:1]=[*:2]-[*:3]=[*:4]-[*+0:5]'),
    Normalization('Normalize 1,5 conjugated cation', \
                  '[n;+0!H0:1]:[a:2]:[a:3]:[c:4]=[N!$(*[O-]),O;+1H0:5]>>[n+1:1]:[*:2]:[*:3]:[*:4]-[*+0:5]'),
    # Normalization('Normalize 1,5 conjugated cation',
    # '[N,O;+0!H0:1]-[c:2]:[a:3]:[a:4]:[n!$(*[O-]),o;+1H0:5]>>[*+1:1]=[c:2]:[*:3]:[*:4]:[*+0:5]'),
    # Normalization('Normalize 1,5 conjugated cation',
    # '[n;+0!H0:1]1:[a:2]:[a:3]:[a:4]:[n!$(*[O-]);+1H0:5]1>>[n+1:1]1:[*:2]:[*:3]:[*:4]:[n+0:5]1'),
    # Normalization('Normalize 1,5 conjugated cation',
    # '[n;+0!H0:1]:[a:2]:[a:3]:[a:4]:[n!$(*[O-]);+1H0:5]>>[n+1:1]:[*:2]:[*:3]:[*:4]:[n+0:5]'),
    # Equivalent to #1.6 in InChI technical manual. RDKit Sanitization handles this for perchlorate.
    Normalization('Charge normalization', \
                  '[F,Cl,Br,I,At;-1:1]=[O:2]>>[*-0:1][O-:2]'),
    Normalization('Charge recombination', \
                  '[N,P,As,Sb;-1:1]=[C+;v3:2]>>[*+0:1]#[C+0:2]'),
    # Normalization from previous protocole with ChemAxon 
    # (some may be redundant: 1,2,4,9,10,15,17)
    Normalization('1', '[#6][NH2+:1][O;X1-:2]>>[#6][N:1]=[O:2]'),
    Normalization('2', '[O:3]=[N:1]=[O:2]>>[#8-:2][N+:1]=[O:3]'),
    Normalization('3', '[#6]-[#7;X4:1]=[O:2]>>[#6][N+:1][O-:2]'),
    Normalization('4', '[n:1]=[O:2]>>[N+:1][#8-:2]'),
    Normalization('5', 'C=[N:1]=[O:2]>>[#8-:2][N+:1]=C'),
    Normalization('6', '[#6;X3+:1]-[#7;X3:2]>>[C:1]=[N+:2]'),
    Normalization('7', '[C;X2+:1]=[N;X2:2]>>[C:1]#[N+:2]'),
    Normalization('8', '[#6][N:1]=[N+:2]>>[#6][N+:1]#[N:2]'),
    Normalization('9', '[#6;X3-:1][N;X2+:2]#[N;X1:3]>>[C:1]=[N+:2]=[#7-:3]'),
    Normalization('10', '[N;X2-:1][N;X2+:2]#[N;X1:3]>>[N:1]=[N+:2]=[#7-:3]'),
    Normalization('11', '[N+:1][C-:2]=[O:3]>>[N:1]=[C:2]=[O:3]'),
    Normalization('12', '[#6]C([#6])=[S+:1][#8-:2]>>[#6]C([#6])=[S:1]=[O:2]'),
    Normalization('13', '[#6][S-:1]([#6])[C+:2]>>[#6][S:1]([#6])=[C:2]'),
    Normalization('14', '[#6][S;X3+:1]([#6])[#8-:2]>>[#6][S:1]([#6])=[O:2]'),
    Normalization('15', '[#6][S+:1]([#6])([#8-:2])=O>>[#6][S:1]([#6])(=[O:2])=O'),
    Normalization('16', '[#6][S+:1]([#6])([#8-:2])=C>>[#6][S:1]([#6])(=C)=[O:2]'),
    Normalization('17', '[#6][P+:1]([O;X2])([O;X2])[#8-:2]>>[#6][P:1](O)(O)=[O:2]'),
    Normalization('18', '[#6][P-:1]([#6])([#6])[C+:2]>>[#6][P:1]([#6])([#6])=[C:2]'),
    Normalization('19', '[O;X2][Se+:1]([O;X2])[#8-:2]>>O[Se:1](O)=[O:2]'),
    Normalization('20', '[O;X2][Si+:1]([O;X2])[#8-:2]>>O[Si:1](O)=[O:2]'),
    Normalization('21', \
                  '[*:1][N-:2][S:3]([*:4])(=[O:5])=[O:6]>>[*:1][N:2][S:3]([*:4])(=[O:5])=[O:6]'),
    Normalization('22', \
                  '[*:4][C:2](=[O:3])[N-:1][C:6]([*:5])=[O:7]>>[*:4][C:2](=[O:3])[N:1][C:6]([*:5])=[O:7]'),
    Normalization('23', \
                  '[#1A][N+:1]1=[C:2][C:3]=[C:4][C:5]=[C:6]1>>[C:4]1=[C:3][C:2]=[N:1][C:6]=[C:5]1')
)

def charge_nitrogen_atoms(mol):
    """
    Charge nitrogen atoms with a valence of 4 and a charge different of 1
    """
    pattern = Chem.MolFromSmarts("[#7]")
    at_matches = mol.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]
    if len(at_matches_list) > 0:
        for at_idx in at_matches_list:
            atom = mol.GetAtomWithIdx(at_idx)
            chg = atom.GetFormalCharge()
            if atom.GetExplicitValence() == 4 and chg != 1:
                atom.SetFormalCharge(1)
                print("WARNING NITROGEN: ({}) ".format(mol.GetProp('_Name')) + \
                      "Nitrogen explicit valence = 4 and charge != 1 " + \
                      "- Nitrogen charge have been changed")
                atom.UpdatePropertyCache()
    return mol


def standardize(mol):
    """
    Standardize a molecule by doing the following steps:
    - Charge nitrogen atoms with a valence of 4 without a charge of 1
    - Remove hydrogens (not essentials hydrogens)
    - Break covalent bonds between metals and organic atoms under certain
      conditions
    - Apply a series of Normalization transforms to correct functional groups
      and recombine charges
    - Enforce charges on certain atoms, then perform competitive reionization
    - Return the largest covalent unit
    - Neutralize molecule by adding/removing hydrogens. Attempts to preserve
      zwitterions
    - Assign chiral tags to atoms if the 3D structure of the molecule is known
    - Assign stereochemistry informations to atoms and bonds if not specified.
      Tag chiral atoms
    """
    print("Standardize {}".format(mol.GetProp('_Name')))
    name_mol = mol.GetProp('_Name')
    s = Standardizer(normalizations=NORMALIZATIONS_perso)
    # CORRECT N without Charge
    mol.UpdatePropertyCache(strict=False)
    mol = charge_nitrogen_atoms(mol)
    mol_props = mol.GetPropsAsDict()
    mol = copy.deepcopy(mol)
    Chem.SanitizeMol(mol)
    mol = Chem.RemoveHs(mol)
    mol = s.disconnect_metals(mol)
    mol = s.normalize(mol)
    mol = s.reionize(mol)
    mol = s.fragment_parent(mol)
    un = rdMolStandardize.Uncharger()
    mol = un.uncharge(mol)
    Chem.AssignAtomChiralTagsFromStructure(mol)  # Maybe not necessary here
    Chem.AssignStereochemistryFrom3D(mol)
    # should keep the predefine stereochemistry, keep chirality informations
    Chem.AssignStereochemistry(mol, flagPossibleStereoCenters=True)
    for k, v in mol_props.items():
        mol.SetProp(k, str(v))
    mol.SetProp("_Name", name_mol)
    return mol


def standardize_list(mol_list):
    """
    standardize a list of molecules
    """
    res_list = []
    num_mol = 0
    for mol in mol_list:
        try:
            res_list.append(standardize(mol))
            res_list[num_mol].SetProp("_Name", mol.GetProp('_Name'))
            num_mol += 1
        except ValueError:
            print("ERROR STD : in standardizing molecule " + mol.GetProp('_Name'))
    return res_list


# endregion
# Return the respective canonical tautomers of each molecule in the list
_normalization_transforms_phosphorus = """
P Valence7 to Valence5	[PH1:1](=[*:2])(=[O:3])[O:4]>>[PH0:1](=[*:2])([OH1:3])[*:4]
"""
_normalizer_params = rdMolStandardize.CleanupParameters()
_normalizer = rdMolStandardize.NormalizerFromData(_normalization_transforms_phosphorus,
                                                  _normalizer_params)


def canonical_tautomer(mol, enumerator, replace_phosphorus_grp=True):
    """
    Compute the canonical tautomere according a specific enumerator.
    Here the enumarator should preserve the stereochemistry
    
    Correct valency of 7 for phosphate caused by a valency of 7 accepted
    cause of Hexafluorophosphate. tranform fragment to one wiht a valence of 5
    for the phosphate
    """
    print("Canonicalize {}".format(mol.GetProp('_Name')))
    mol = enumerator.Canonicalize(mol)
    if replace_phosphorus_grp:
        pattern = Chem.MolFromSmarts("[PH](=O)(=O)O")
        at_matches = mol.GetSubstructMatches(pattern)
        at_matches_list = [y[0] for y in at_matches]
        if len(at_matches_list) > 0:
            print(at_matches_list)
            print("WARNING CANONICALIZE: Canonicalize has change phosphorus group," + \
                  " changing it back to a valence of 5 mol " + mol.GetProp('_Name'))
            name = mol.GetProp('_Name')
            print(Chem.MolToSmiles(mol, isomericSmiles=True))
            mol = Chem.RWMol(mol)
            mol = _normalizer.normalize(mol)
            Chem.SanitizeMol(mol)
            Chem.AssignAtomChiralTagsFromStructure(mol)
            Chem.AssignStereochemistryFrom3D(mol)
            # again to reassign chiral atoms
            Chem.AssignStereochemistry(mol, \
                                       force=True, \
                                       flagPossibleStereoCenters=True)
            mol.SetProp('_Name', name)
            print(Chem.MolToSmiles(mol, isomericSmiles=True))
            # mol = standardize(mol)
    #
    Chem.SanitizeMol(mol)
    Chem.AssignAtomChiralTagsFromStructure(mol)
    Chem.AssignStereochemistryFrom3D(mol)
    # again to reassign chiral atoms
    Chem.AssignStereochemistry(mol, \
                               force=True, \
                               flagPossibleStereoCenters=True)
    return mol


def canonical_tautomer_list(mol_list):
    """
    Return the canonical tautomer for a list of molecules and preserve its
    stereochemistry
    """
    params = rdMolStandardize.CleanupParameters()
    # Keep double bonds essential for the stereochemistry
    params.tautomerRemoveBondStereo = False
    # To keep the stereochemistry when its specified
    params.tautomerRemoveSp3Stereo = False
    enumerator = rdMolStandardize.TautomerEnumerator(params)
    res_list = []
    for mol in mol_list:
        res_list.append(canonical_tautomer(mol, enumerator))
    return res_list

