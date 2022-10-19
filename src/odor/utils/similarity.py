from rdkit.Chem import MACCSkeys
from rdkit import DataStructs

def get_dico_maccs(list_mol):
    """
    Compute maccs fingerprints for a list of molecules
    """
    dico_mol_maccs = {}
    for mol in list_mol:
        dico_mol_maccs[mol.GetProp("_Name")] = DataStructs.CreateFromFPSText(
                DataStructs.BitVectToFPSText(MACCSkeys.GenMACCSKeys(mol)))
    return dico_mol_maccs
    
    
def get_sorted_tanimoto_maccs(mol_maccs, db_maccs):
    '''
    mol_maccs: the maccs bitvector of the query molecule (list)
    db_maccs: all the maccs bitvector to compare with (list)
    return: all the tanimoto scores of the query molecule
    '''
    list_sim = []
    for name_mol in db_maccs.keys():
        sim = DataStructs.TanimotoSimilarity(list(mol_maccs.values())[0], db_maccs[name_mol])
        list_sim.append((name_mol, sim))
    list_sim.sort(reverse=True, key=lambda x: x[1])
    return list_sim
    
def get_dict_sim(mol, db_mol, similarity ="maccs"):
    '''
    Get a dictionary of tanimoto similarity for a molecule agaiinst a database
    '''
    dict_sim_mol_db = dict()
    if similarity == "maccs":
        mol_maccs = get_dico_maccs(mol) # into list when one molecule
        db_maccs = get_dico_maccs(db_mol)
        print(mol_maccs)
        #print(db_maccs)
        for sim in get_sorted_tanimoto_maccs(mol_maccs, db_maccs):
            dict_sim_mol_db[str(sim[0])] = sim[1]
    #### ADD OTHER SIMILARITY METRIC
    #### MORGAN etc ...
    return dict_sim_mol_db
