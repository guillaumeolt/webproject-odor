from django.core.files.storage import FileSystemStorage
import shutil
import re
import json
from django.utils.safestring import mark_safe

from .standardize import *
from .read_write import *
from .similarity import *
import pubchempy as pcp
from django.template.defaulttags import register
import subprocess
#from .prediction import *
# Check Files Type POST


@register.filter
def str_split_or(str_or):
    return str_or.split(";")
@register.filter
def get_item(dictionary, key):
    return dictionary.get(key)
@register.filter
def highlight_search(text, search):
    search = str(search)
    text = str(text)
    #pattern = re.compile(re.escape(search), re.IGNORECASE)
    rgx = re.compile(re.escape(search), re.IGNORECASE)
    #pattern = pattern.sub('<span class="highlight" style="background-color: #FFFF00">{}</span>'.format(search), text)
    pattern = rgx.sub(
        lambda m: '<b class="text text-danger">{}</b>'.format(m.group()),
        text
    )
    return mark_safe(pattern)
@register.filter
def get_cas(synonyms):
    try:
        for syn in synonyms:
            match = re.match('(\d{2,7}-\d\d-\d)', syn)
            if match:
                cas = match.group(1)
    except:
        return None
    return None


def get_path_svg_db(db_dict):
    for i in range(len(db_dict)):
        db_dict[i]["path_svg"] = os.path.join("media/db_mols_svg",str(db_dict[i]["idChemicals"])) + ".svg"
    return db_dict

# TRANSFORM DB INPUT
def tranform_db_dict(db_dict):
    """
    transform or dictionary to list
    """
    for i in range(len(db_dict)):
        if db_dict[i]["OlfRecept"] != None:
            db_dict[i]["OlfRecept"] = db_dict[i]["OlfRecept"].split(";")
    return(db_dict)
def tranform_db_dict_iduniprot(receptors_iduniprot):
    """
    transform or dictionary with or keys
    """
    receptors_iduniprot_dict = dict()
    for i in range(len(receptors_iduniprot)):
        receptors_iduniprot_dict[receptors_iduniprot[i]['GeneName']] = receptors_iduniprot[i]['idUniprot']
    return(receptors_iduniprot_dict)
# 1 Check Input:
def check_valid_input(request):
    """
    Process input from the form in the search path.
    return true/fals eand the molecule as rdkit
    """
    smile_post = request.POST.get('smilefield', "None")
    #print([smile_post],"-")
    if request.POST.get("btn_search") != None or request.POST.get("btn_predict") != None:
        if smile_post != "None" and smile_post != '':
            #print("here smile?")
            try:
                mol = Chem.MolFromSmiles(smile_post)
                mol.SetProp("_Name","ketcher_mol")
                return(True, mol, None)
            except:
                pass
            #print("here inchi?")
            try:
                mol = Chem.MolFromInchi(smile_post)
                mol.SetProp("_Name","ketcher_mol")
                return(True, mol, None)
            except:
                pass
            #InchiKey
            try:
                try:
                    results = pcp.get_compounds(smile_post,'inchikey')
                    compound = results[0]
                    try :
                        smile = compound.isomeric_smiles
                    except:
                        try:
                            smile = compound.canonical_smiles
                        except:
                            pass
                except:
                    pass
                mol = Chem.MolFromSmiles(smile)
                mol.SetProp("_Name","ketcher_mol")
                return(True, mol, compound)
            except:
                pass
    if request.POST.get("btn_match") != None:
        if smile_post != "None" and smile_post != '':
            #print("here?smart")
            try:
                mol = Chem.MolFromSmarts(smile_post)
                mol.SetProp("_Name","ketcher_mol")
                return(True, mol, None)
            except:
                pass
            try:
                mol = Chem.MolFromInchi(smile_post)
                mol.SetProp("_Name","ketcher_mol")
                return(True, mol, None)
            except:
                pass
            #InchiKey
            try:
                try:
                    results = pcp.get_compounds(smile_post,'inchikey')
                    compound = results[0]
                    try :
                        smile = compound.isomeric_smiles
                    except:
                        try:
                            smile = compound.canonical_smiles
                        except:
                            pass
                except:
                    pass
                mol = Chem.MolFromSmiles(smile)
                mol.SetProp("_Name","ketcher_mol")
                return(True, mol, compound)
            except:
                pass
    pubchem_name = request.POST.get('pubchem_name', "None")
    #print([pubchem_name],"-")
    if pubchem_name != "None" and pubchem_name != '':
        try:
            try:
                results = pcp.get_compounds(pubchem_name,'name')
                compound = results[0]
                try :
                    smile = compound.isomeric_smiles
                except:
                    try:
                        smile = compound.canonical_smiles
                    except:
                        pass
            except:
                pass
            mol = Chem.MolFromSmiles(smile)
            mol.SetProp("_Name","ketcher_mol")
            return(True, mol, compound)
        except:
            pass

    try:
        upload_file = request.FILES["document"]
        #print(request.FILES["document"],"----")
        fs = FileSystemStorage()
        fs.save("mol.sdf", upload_file)
        list_mol = load_all_mol_from_file("media/" + "mol.sdf" ) #upload_file.name)
        fs.delete("mol.sdf")
        mol = Chem.MolFromSmiles(Chem.MolToSmiles(list_mol[0]))
        mol.SetProp("_Name", upload_file.name)
        mol = standardize_list(list_mol)
        return (True, mol[0], None)
    except:
        #print("nop")
        pass
    return(False, "Error input files", None)
# 1 Check Input:
def get_pubchem_compound(mol):
    try:
        results = pcp.get_compounds(Chem.MolToSmiles(mol), 'smiles')
        compound = results[0]
        return(compound)
    except:
        return(None)

# 2 Check if mixture
def is_mixture(mol):
    """
    Check if a molecule is a mixture
    """
    if len(Chem.GetMolFrags(mol, asMols=False)) > 1:

        return True
    return False

# 4 Compute similarity
def  utils_get_similarity(ChemicalsOdors, mol_input, db_dict, sim_metric ="maccs"):
    """
    Compute similarity between the database and the query molecule
    """
    list_db = []
    for chem in ChemicalsOdors:
        mol = Chem.MolFromSmiles(chem.SMILE)
        mol.SetProp("_Name", str(chem.idChemicals))
        list_db.append(mol)

    if sim_metric == "maccs":
        dict_sim_macc = get_dict_sim([mol_input], list_db)
        for mol_info_db in db_dict:
            try:
                mol_info_db["Similarity"] = dict_sim_macc[str(mol_info_db["idChemicals"])]
            except:
                print("Error maccs")
    return db_dict


def run_2d_image_SVG_list_substructures(db_dict, patt, name="_Name"):

    db_dict_sub = []
    for mol_info_db in db_dict:

        mol = Chem.MolFromSmiles(mol_info_db["SMILE"])
        mol.SetProp("_Name",str(mol_info_db["idChemicals"]))

        hit_ats = list(mol.GetSubstructMatch(patt))
        if len(hit_ats) == 0:
            #print(len(hit_ats), Chem.MolToSmiles(patt), Chem.MolToSmiles(mol))
            continue
        else:
            pass

        hit_bonds = []
        for bond in patt.GetBonds():
            aid1 = hit_ats[bond.GetBeginAtomIdx()]
            aid2 = hit_ats[bond.GetEndAtomIdx()]
            hit_bonds.append(mol.GetBondBetweenAtoms(aid1, aid2).GetIdx())

        d = rdMolDraw2D.MolDraw2DSVG(200, 200)  # or MolDraw2DCairo to get PNGs
        rdMolDraw2D.PrepareAndDrawMolecule(d, mol, highlightAtoms=hit_ats, highlightBonds=hit_bonds)

        d.FinishDrawing()
        mol_info_db["path_svg"] = d.GetDrawingText()
        db_dict_sub.append(mol_info_db)
    return(db_dict_sub)

def run_2d_image_SVG_list_substructures_past(db_dict, patt, path_output, name="_Name"):
    # REMOVE ALL FILES
    folder = path_output
    for filename in os.listdir(folder):
        file_path = os.path.join(folder, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            print('Failed to delete %s. Reason: %s' % (file_path, e))
    db_dict_sub = []
    for mol_info_db in db_dict:
        #print(mol_info_db["SMILE"])
        mol = Chem.MolFromSmiles(mol_info_db["SMILE"])
        mol.SetProp("_Name",str(mol_info_db["idChemicals"]))
        #smi = 'c1cc(F)ccc1Cl'
        #mol = Chem.MolFromSmiles(smi)
        #mol.SetProp("_Name","test")
        #patt = Chem.MolFromSmarts('ClccccF')

        hit_ats = list(mol.GetSubstructMatch(patt))
        if len(hit_ats) == 0:
            #print(len(hit_ats), Chem.MolToSmiles(patt), Chem.MolToSmiles(mol))
            continue
        else:
            #print(Chem.MolToSmiles(patt), Chem.MolToSmiles(mol))
            pass

        hit_bonds = []
        for bond in patt.GetBonds():
            aid1 = hit_ats[bond.GetBeginAtomIdx()]
            aid2 = hit_ats[bond.GetEndAtomIdx()]
            hit_bonds.append(mol.GetBondBetweenAtoms(aid1, aid2).GetIdx())

        d = rdMolDraw2D.MolDraw2DSVG(500, 500)  # or MolDraw2DCairo to get PNGs
        rdMolDraw2D.PrepareAndDrawMolecule(d, mol, highlightAtoms=hit_ats, highlightBonds=hit_bonds)

        d.FinishDrawing()
        outs = open(path.join(path_output, mol.GetProp("_Name") + ".svg"), "w")
        print(d.GetDrawingText(), file=outs)
        #print(d.GetDrawingText())
        outs.close()
        mol_info_db["path_svg"] = path.join(path_output, mol.GetProp("_Name") + ".svg")
        db_dict_sub.append(mol_info_db)
    return(db_dict_sub)
## Past functions :
def utils_test(file_mol):
    list_db = load_all_mol_from("db")
    list_db = standardize_list(list_db)
    #list_mol = load_all_mol_from("media")
    list_mol = load_all_mol_from_file("media/"+file_mol)
    list_mol = standardize_list(list_mol)
    #Compute tanimoto


    list_mol[0].SetProp("_Name",os.path.splitext(file_mol)[0])
    #print(os.getcwd())
    dict_sim_macc = get_dict_sim(list_mol, list_db)

    save_2d_image_PNG_list(list_mol, "media", name = "_Name")
    return(True, list_mol[0].GetProp("_Name"))

def utils_get_maccs_similarity(ChemicalsOdors, db_dict, name_mol, smi_str=False):

    chemicals_odors = ChemicalsOdors
    list_db = []
    for chem in chemicals_odors:

        mol = Chem.MolFromSmiles(chem.SMILE)
        mol.SetProp("_Name", str(chem.idChemicals))
        list_db.append(mol)
    #list_db = load_all_mol_from("db")
    #list_db = standardize_list(list_db)
    #
    if smi_str!=False:
        try:
            list_mol = [Chem.MolFromSmiles(smi_str)]
            list_mol[0].SetProp("_Name","ketcher_mol")
            save_2d_image_PNG_list(list_mol, "media", name="_Name")
        except:
            pass
    else:
        try:

            list_mol = load_all_mol_from_file("media/" + name_mol)

        except:

            pass
    list_mol = standardize_list(list_mol)

    #Compute tanimoto

    dict_sim_macc = get_dict_sim([list_mol[0]], list_db)
    for mol_info_db in db_dict:
        try:
            mol_info_db["Similarity"] = dict_sim_macc[str(mol_info_db["idChemicals"])]
        except:
            #print(mol_info_db)
            #print(mol_info_db["idChemicals"], type(mol_info_db["idChemicals"]))
            #print(dict_sim_macc["6033"], dict_sim_macc[6033])
            #print("error")
            pass
    #print(db_dict,"--/--")
    return dict_sim_macc

### PREDICTION
import os

def utils_get_prediction_odor(BASE_DIR, query_smile = "CC=O", predict = "odor"):
    """
    Run subprocess usinng specific environnement to laucnch prediction_subprocess and predict a given smile odor
    """
    # Path to a Python interpreter that runs any Python script
    # under the virtualenv /path/to/virtualenv/
    #python_bin = "/home/guillaumeolt/miniconda3/envs/cpmli-predict-odor/bin/python"
    python_bin = str(BASE_DIR) + "/../.env-prediction/bin/python3.7"
    # Path to the script that must run under the virtualenv
    script_file = str(BASE_DIR) + "/odor/utils/prediction_subprocess.py"
    if predict == "odor":
        p = subprocess.Popen([python_bin, script_file,
                              "-smile", query_smile,
                              "-label", str(BASE_DIR) + "/odor/utils/infos_models/label_odors_names.txt",
                              "-PATH_GCC", str(BASE_DIR) + "/odor/utils/Trained_models/GNN_CODB_ckpt",
                              "-predict", "odor"], stdout=subprocess.PIPE) #, stderr=subprocess.PIPE)
    if predict == "or":
        p = subprocess.Popen([python_bin, script_file,
                              "-smile", query_smile,
                              "-label", str(BASE_DIR) + "/odor/utils/infos_models/label_or_names.txt",
                              "-PATH_GCC", str(BASE_DIR) + "/odor/utils/Trained_models/GNN_RCDB_HO_ckpt",
                              "-predict", "or"], stdout=subprocess.PIPE) # , stderr=subprocess.PIPE)
    out, err = p.communicate()
    #print([out.decode("utf-8") ],"----------------")
    #out = out.split("\n")
    """
    out = str(out)
    out = out.replace("b\'","")
    out = out.replace("\\n\'", "")
    out = out.rstrip().split(",")
    """
    out = out.decode("utf-8")
    out = out.rstrip()
    #print(out, type(out))
    out = out.replace("\\n", "")
    out = out.replace("\'", "\"")
    #print([out],"---------- - -")
    out = json.loads(out)
    #print(out)
    out = dict(sorted(out.items(), key=lambda x: x[1], reverse=True))
    return(out)
