import json
import random
import shutil
from django.core.serializers.json import DjangoJSONEncoder
from django.db.models import F, Q
from django.forms import JSONField
from django.http import HttpResponse
from django.shortcuts import render, redirect
from datetime import datetime
from django.core.files.storage import FileSystemStorage
from django.templatetags.static import static
from django.contrib.staticfiles.storage import staticfiles_storage
from odorwebsite.settings import BASE_DIR, STATIC_ROOT
from .utils.docking import get_url_dockign_seamdock
from .utils.get_descriptors import get_decriptors_rdkit
from .utils.test import *
from .utils.tools_bokeh import *
from .utils.tools_plotly import *
from .utils.tools_umap import *

from bokeh.plotting import figure
from bokeh.embed import components

import urllib.request

from .models import ChemicalsOdors,\
                          OlfactoryReceptors,\
                          Smell_Percepts,\
                          my_custom_sql,\
                          my_custom_sql_chem_id,\
                          my_custom_sql_chem_predict,\
                          my_custom_sql_chem_get_odor,\
                          my_custom_sql_odor_id,\
                          my_custom_sql_odor_get_chem,\
                          my_custom_sql_with_human_homologue


def OdorWebSite(request):
    # LOAD Database Infos : Chemicals info, Smells infos, Olfactory Receptors Infos
    db_dict = my_custom_sql()
    return render(request, "OdorWebSite_Home.html", context={'db': db_dict})

def OdorWebSite_design(request):
    # LOAD Database Infos : Chemicals info, Smells infos, Olfactory Receptors Infos
    db_dict = my_custom_sql()
    return render(request, "OdorWebSite_webesign.html", context={'db': db_dict})

def OdorWebSite_Home(request):
    # LOAD Database Infos : Chemicals info, Smells infos, Olfactory Receptors Infos
    db_dict = my_custom_sql()
    chemicals_odors = ChemicalsOdors.objects.all()
    olfactory_receptors = OlfactoryReceptors.objects.all()
    smell_percepts = Smell_Percepts.objects.all()
    return render(request, "OdorWebSite_Home.html", context={'db': db_dict,
                                                             "chemicals_odors": chemicals_odors,
                                                             "olfactory_receptors": olfactory_receptors,
                                                             "smell_percepts": smell_percepts})
def OdorWebSite_About(request):
    # LOAD Database Infos : Chemicals info, Smells infos, Olfactory Receptors Infos
    db_dict = my_custom_sql()
    return render(request, "OdorWebSite_About.html", context={'db': db_dict})


def OdorWebSite_Search(request):
    # LOAD Database Infos : Chemicals info, Smells infos, Olfactory Receptors Infos
    db_dict = my_custom_sql()
    receptors_iduniprot = OlfactoryReceptors.objects.values('GeneName', 'idUniprot')
    # Transform data
    db_dict = tranform_db_dict(db_dict)
    receptors_iduniprot_dict = tranform_db_dict_iduniprot(receptors_iduniprot)

    db_dict = get_path_svg_db(db_dict)
    if request.method == 'POST':
        # 1 Check Input
        is_valid, mol, pubchem_compound = check_valid_input(request)
        # 1 bis load pubchem infos
        if pubchem_compound == None and is_valid:
            pubchem_compound = get_pubchem_compound(mol)
        if not is_valid:
            return render(request, "OdorWebSite_Search-Structure.html", context={'db': db_dict,
                                                                       'dict_or_id': receptors_iduniprot_dict,
                                                                       'error_message': mol})
        # 2 Check if mixture
        if is_mixture(mol):
            return render(request, "OdorWebSite_Search-Structure.html", context={'db': db_dict,
                                                                       'dict_or_id': receptors_iduniprot_dict,
                                                                       'error_message': "Error input files mixtures"})
        # 3 Generate 2D image
        type_search_sub = False
        try:
            mol_svg = mol2svg(mol)
            #save_2d_image_PNG_list([mol], "media", name="_Name")
        except Exception as e:
            print("No image maybe smart input")
            pass
        if request.POST.get("btn_match") != None:
            #db_dict = run_2d_image_SVG_list_substructures(db_dict,mol,"media/result_substructure_match_png")
            db_dict = run_2d_image_SVG_list_substructures(db_dict, mol)
            type_search_sub = True
            pass
        # 4 Compute similarity
        if request.POST.get("btn_search") != None:
            chemicals_odors = ChemicalsOdors.objects.all()
            db_dict = utils_get_similarity(chemicals_odors, mol, db_dict, sim_metric="maccs")
        # 5 Sort similarity
        if request.POST.get("btn_search") != None:
            db_dict = sorted(db_dict, key=lambda i: i['Similarity'], reverse=True)
            #print(db_dict)

        # 2 Check if substructure retunr results
        if db_dict == []:
            if pubchem_compound is None:
                return render(request, "OdorWebSite_Search-Structure.html", context={'db': db_dict,
                                                                                     'dict_or_id': receptors_iduniprot_dict,
                                                                                     'error_message': "Error input substructure"})
            return render(request, "OdorWebSite_Search-Structure.html", context={'db': db_dict,
                                                                       'dict_or_id': receptors_iduniprot_dict,
                                                                       'error_message': "Error input substructure: no substructure match in the database"})
        #print(db_dict, pubchem_compound, mol)
        if pubchem_compound is None and Chem.MolToSmarts(mol) != None:
            return render(request, "OdorWebSite_Search-Structure.html",
                          context={"image": mol_svg,
                                   'db': db_dict,
                                   'dict_or_id': receptors_iduniprot_dict,
                                   'compound': pubchem_compound,
                                   "smart": Chem.MolToSmarts(mol) ,
                                   'type_search_sub': type_search_sub,
                                   'error_message': None})

        return render(request, "OdorWebSite_Search-Structure.html", context={"image": mol_svg,
                                                                   'db': db_dict,
                                                                   'dict_or_id': receptors_iduniprot_dict,
                                                                   'compound': pubchem_compound,
                                                                   'type_search_sub': type_search_sub,
                                                                   'smart': None,
                                                                   'error_message':None})


    return render(request, "OdorWebSite_Search-Structure.html", context={'db': db_dict,
                                                               'dict_or_id':receptors_iduniprot_dict,
                                                               'error_message':None,
                                                               'smart': None,
                                                               'compound': None})
def OdorWebSite_Predict(request):
    # LOAD Database Infos : Chemicals info, Smells infos, Olfactory Receptors Infos
    db_dict = my_custom_sql_with_human_homologue() #my_custom_sql_chem_predict() #my_custom_sql_with_human_homologue
    receptors_iduniprot = OlfactoryReceptors.objects.values('GeneName', 'idUniprot')
    # Transform data
    db_dict = tranform_db_dict(db_dict)
    receptors_iduniprot_dict = tranform_db_dict_iduniprot(receptors_iduniprot)

    # Add SVG
    db_dict = get_path_svg_db(db_dict)
    db_dict_all = db_dict
    if request.method == 'POST':
        # 1 Check Input
        is_valid, mol, pubchem_compound = check_valid_input(request)
        if not is_valid:
            return render(request, "OdorWebSite_Predict.html", context={'db': None,
                                                                       'dict_or_id': None,
                                                                       'error_message': "Error input files"})
        # 2 Check if mixture
        if is_mixture(mol):
            return render(request, "OdorWebSite_Predict.html", context={'db': None,
                                                                       'dict_or_id': None,
                                                                       'error_message': "Error input files mixtures"})
        # 3 Generate 2D image
        mol_svg = mol2svg(mol)
        #save_2d_image_PNG_list([mol], "media", name="_Name")

        dt_desc_mol = get_decriptors_rdkit([mol])
        dic_desc_mol = dict(dt_desc_mol.loc[mol.GetProp("_Name")])
        #print(dic_desc_mol)
        # 4 Compute similarity
        if request.POST.get("btn_predict") != None:
            chemicals_odors = ChemicalsOdors.objects.all()
            db_dict = utils_get_similarity(chemicals_odors, mol, db_dict, sim_metric="maccs")
        # 5 Sort similarity
        if request.POST.get("btn_predict") != None:
            db_dict = sorted(db_dict, key=lambda i: i['Similarity'], reverse=True)
        # 6 get 10 most similar
            db_dict = db_dict[:5]
        # 7 Prédict molécule
        #list = ["Odorless", "Odorant"]


        # 8 PREDICT / 9 Visualize prediction wit Bokeh
        mapper = load_umap_chem_odor("odor/static/media/umap/mapper.pkl")
        if request.POST.get("predict_model_field") == "Odor" or request.POST.get("predict_model_field") != "Olfactory Receptor (Human)":
            dict_pred = utils_get_prediction_odor(BASE_DIR, query_smile = Chem.MolToSmiles(mol), predict = "odor")
            list_pred = list(dict_pred.keys())
            script, div = get_bokeh_plot_odor_from_list_odors(mapper, db_dict_all, list_pred, path_svg_add="static/media/db_mols_svg/")

            # RADDAR PLOT PLOTLY
            tmp = list(dict_pred)
            tmp.append('All')
            #print(dict_pred,"------------------------",type(dict_pred))
            df = get_data_desc_plotly_list_odor(db_dict_all, tmp)
            div_radar_plot = get_radar_plot_from_list_odor(df)

        if request.POST.get("predict_model_field") == "Olfactory Receptor (Human)":
            dict_pred = utils_get_prediction_odor(BASE_DIR, query_smile = Chem.MolToSmiles(mol), predict="or")
            
            #list_pred = list(dict_pred.keys())
            #script, div = get_bokeh_plot_odor_from_list_or(mapper, db_dict_all, list_pred, path_svg_add="static/media/db_mols_svg/")


            # RADDAR PLOT PLOTLY
            #tmp = list(list_pred)
            #tmp.append('All')
            #print(list_pred,"------------------------",type(list_pred))
            #df = get_data_desc_plotly_list_or(db_dict_all, tmp)
            #div_radar_plot = get_radar_plot_from_list_or(df)

            #div_radar_plot = None
        #res_mol_predict = random.choice(list)



        return render(request, "OdorWebSite_Predict.html", context={"image": mol_svg,
                                                                   'db': db_dict,
                                                                   'dict_or_id': receptors_iduniprot_dict,
                                                                   'dic_desc_mol':dic_desc_mol,
                                                                   'prediction':dict_pred,
                                                                   'error_message':None,
                                                                    'script': script, 'div': div,
                                                                    'div_radar_plot':div_radar_plot})


    return render(request, "OdorWebSite_Predict.html", context={'db': None,
                                                               'dict_or_id':None,
                                                               'error_message':None})

def OdorWebSite_Contact(request):
    # LOAD Database Infos : Chemicals info, Smells infos, Olfactory Receptors Infos
    db_dict = my_custom_sql()
    return render(request, "OdorWebSite_Contact.html", context={'db': db_dict})


def index(request):
    date = datetime.today()
    chemicals_odors = ChemicalsOdors.objects.all()
    db_dict = my_custom_sql()
    #db_dict = json.dumps(db_dict)
    #print(date, len(db_dict))
    if request.method == 'POST':
        test=False
        try:
            smile_post = request.POST.get('smilefield',"---what--")
            #print(smile_post,"---------------------------")
            #test, mol_name = utils_test("ketcher_mol")
            dict_sim_macc = utils_get_maccs_similarity(chemicals_odors, db_dict, "ketcher_mol",smile_post)
            js_data = json.dumps(dict_sim_macc, sort_keys=False)
            test = True
            mol_name = "ketcher_mol"
        except:
            pass
        try:
            upload_file = request.FILES["document"]
            fs = FileSystemStorage()
            fs.save(upload_file.name, upload_file)
            test, mol_name = utils_test(upload_file.name)
            dict_sim_macc = utils_get_maccs_similarity(chemicals_odors, db_dict, upload_file.name)
            js_data = json.dumps(dict_sim_macc, sort_keys=False)
        except:
            print("nop")
            pass
        smile_post = request.POST.get('smilefield', "---what--")
        if test:
            #print(db_dict,"--//")
            db_dict = sorted(db_dict, key=lambda i: i['Similarity'],reverse=True)
            return render(request, "index.html", context={"prenom": "Guillaume_OAUIP",
                                                  "date": date,
                                                  "image":"media/"+mol_name+".png",
                                                  "show_id":"True",
                                                  'chemicals_odors':chemicals_odors,
                                                  'dict_sim_macc':None,
                                                  'db': db_dict})

    return render(request, "index.html", context={"prenom": "Guillaume",
                                                  "date": date,
                                                  'chemicals_odors': chemicals_odors,
                                                  'db': db_dict})

def OdorWebSite_OlfactoryReceptor_template(request, GeneName=None):
    GeneName_or = OlfactoryReceptors.objects.get(GeneName=GeneName)
    db_dict_all = my_custom_sql()
    db_dict_all = tranform_db_dict(db_dict_all)

    # get chem with known or
    print(GeneName_or.idOlfactoryReceptors , "-------------")
    dic_or_chem = my_custom_sql_odor_get_chem(GeneName_or.idOlfactoryReceptors)[0]
    #docking list chemicals
    chemicals_odors = ChemicalsOdors.objects.all()


    # plotly
    path_odor_plotly = os.path.join(STATIC_ROOT,'media/plotly_or/' + GeneName + ".html")
    with open(path_odor_plotly, 'r') as file:
    	div_radar_plot = file.read()

    # bokeh
    path_odor_bokeh_div = os.path.join(STATIC_ROOT,'media/bokeh_or/' +GeneName + ".div")
    path_odor_bokeh_script = os.path.join(STATIC_ROOT,'media/bokeh_or/' + GeneName + ".script")
    with open(path_odor_bokeh_div, 'r') as file:
    	div = file.read()
    with open(path_odor_bokeh_script, 'r') as file:
    	script = file.read()
    return render(request, "OdorWebSite_OR.html", context={"GeneName":GeneName_or,
                                                           'script': script, 'div': div,
                                                           'div_radar_plot': div_radar_plot,
                                                           'chemicals_odors':chemicals_odors,
                                                           'dic_or_chem': dic_or_chem})
def OdorWebSite_Chemical_template(request, chem_id=None):
    #print(chem_id, "-------------------")
    olfactory_receptors = OlfactoryReceptors.objects.all()
    chem = ChemicalsOdors.objects.get(idChemicals=chem_id)
    db_dict = my_custom_sql_chem_id(chem_id)[0]
    dic_chem_odor = my_custom_sql_chem_get_odor(chem_id)[0]

    db_dict = {str(key): str(value) for key, value in db_dict.items()}

    return render(request, "OdorWebSite_Chemical.html", context={"chem":chem,
                                                                 "db":db_dict,
                                                                 'olfactory_receptors' : olfactory_receptors,
                                                                 "dic_chem_odor":dic_chem_odor})


def test(request):
    # create a plot
    plot = figure(plot_width=400, plot_height=400)

    # add a circle renderer with a size, color, and alpha

    plot.circle([1, 2, 3, 4, 5], [6, 7, 2, 4, 5], size=20, color="navy", alpha=0.5)

    script, div = components(plot)

    return render(request, 'test.html', {'script': script, 'div': div})

def OdorWebSite_Odor_template(request, odor_id=None, db_dict_all=None):
    #print(chem_id, "-------------------")
    odor = Smell_Percepts.objects.get(idSmell_Percepts=odor_id)
    db_dict = my_custom_sql_odor_id(odor_id)[0]
    dic_odor_chem = my_custom_sql_odor_get_chem(odor_id)[0]

    db_dict_all = my_custom_sql()
    db_dict_all = tranform_db_dict(db_dict_all)

    # Add SVG
    db_dict_all = get_path_svg_db(db_dict_all)
    d = dict()
    for key, value in dic_odor_chem.items():
        d[key] = {str(k): str(v) for k, v in value.items()}
    dic_odor_chem = d
    db_dict = {str(key): str(value) for key, value in db_dict.items()}

    # plotly
    path_odor_plotly = os.path.join(STATIC_ROOT,'media/plotly_odors/' + odor.Odor + ".html")
    with open(path_odor_plotly, 'r') as file:
    	div_radar_plot = file.read()

    # bokeh
    path_odor_bokeh_div = os.path.join(STATIC_ROOT,'media/bokeh_odors/' + odor.Odor + ".div")
    path_odor_bokeh_script = os.path.join(STATIC_ROOT,'media/bokeh_odors/' + odor.Odor + ".script")
    with open(path_odor_bokeh_div, 'r') as file:
    	div = file.read()
    with open(path_odor_bokeh_script, 'r') as file:
    	script = file.read()

    return render(request, "OdorWebSite_Odor.html", context={"odor":odor,
                                                             "db":db_dict,
                                                             "dic_odor_chem":dic_odor_chem,
                                                             'script': script, 'div': div,
                                                             'div_radar_plot': div_radar_plot})
def OdorWebSite_search_chem_or_odor(request):
    if request.method == 'POST':
        search_input_chemical = request.POST.get("search_chemical","None")
        search_input_OR = request.POST.get("search_or","None")
        search_input_odor = request.POST.get("search_odor", "None")

        if search_input_chemical != "":
            db_dict = ChemicalsOdors.objects.filter(Q(Name__icontains=search_input_chemical) |
                                                    Q(idChemicals__icontains=search_input_chemical) |
                                                    Q(IUPAC_name__icontains=search_input_chemical) |
                                                    Q(Pubchem_CID__icontains=search_input_chemical) |
                                                    Q(SMILE__icontains=search_input_chemical) |
                                                    Q(CAS__icontains=search_input_chemical))
            return render(request, "OdorWebSite_Search_chemical.html", context={"search_type":"Chemicals",
                                                                               "search_value":search_input_chemical,
                                                                               "db_dict":db_dict})
        if search_input_OR != "":
            db_dict = OlfactoryReceptors.objects.filter(Q(GeneName__icontains=search_input_OR) |
                                                        Q(idUniprot__icontains=search_input_OR)|
                                                        Q(Species__icontains=search_input_OR)  |
                                                        Q(Synonym__icontains=search_input_OR))
            return render(request, "OdorWebSite_Search_or.html", context={"search_type":"Olfactory Receptors",
                                                                               "search_value":search_input_OR,
                                                                               "db_dict":db_dict})
        if search_input_odor != "":
            db_dict = Smell_Percepts.objects.filter(Q(Odor__icontains=search_input_odor))
            return render(request, "OdorWebSite_Search_odor.html", context={"search_type":"Odor",
                                                                            "search_value":search_input_odor,
                                                                            "db_dict":db_dict})
    return render(request, "OdorWebSite_Search.html", context={"search_type": "Search"})

def OdorWebSite_docking_chem_or(request):
    olfactory_receptors = OlfactoryReceptors.objects.all()
    chemicals_odors = ChemicalsOdors.objects.all()
    params = True
    if request.method == 'POST':
        # Search Input molecule structure
        if request.POST.get("search_chemical", "None") == "":
            # 1 Check Input
            is_valid, mol, pubchem_compound = check_valid_input(request)

            if not is_valid:
                return render(request, "OdorWebSite_Docking.html", context={"olfactory_receptors": olfactory_receptors,
                                                                            "chemicals_odors": chemicals_odors,
                                                                            'error_message': "Error inputs"})
            # 2 Check if mixture
            if is_mixture(mol):
                return render(request, "OdorWebSite_Docking.html", context={"olfactory_receptors": olfactory_receptors,
                                                                            "chemicals_odors": chemicals_odors,
                                                                            'error_message': "Error inputs"})
            # 3 Generate 2D image and 3d pdb

            mol.SetProp("_Name", "lig")
            write_3d_pdb_list([mol], str(BASE_DIR) + '/media/docking/', name="_Name")
            save_2d_image_PNG_list([mol], "media/docking", name="_Name")
            docking_input_chemical = str(BASE_DIR) + '/media/docking/lig.pdb'
        else:
            docking_input_chemical = request.POST.get("search_chemical", "None")

        # Search Input protein structure
        print(request.POST.get("search_or", "None"))
        print(request.POST.get("search_or", "None")," sdferf")
        if request.POST.get("search_or", "None") == "":
            print(request.POST.get('pdb_name', "None"), "---------------------<>")
            try:
                pdb_name = request.POST.get('pdb_name', "None")
                print(pdb_name, "---------------------<>")
                urllib.request.urlretrieve('http://files.rcsb.org/download/' + pdb_name + '.pdb', str(BASE_DIR) + '/media/docking/prot.pdb')

                docking_input_OR = str(BASE_DIR) + "/media/docking/prot.pdb"
                params = False
            except:
                path_prot = None
                return render(request, "OdorWebSite_Docking.html", context={"olfactory_receptors": olfactory_receptors,
                                                                            "chemicals_odors": chemicals_odors,
                                                                            'error_message': "Error inputs"})
        else:
            try:
                docking_input_OR = request.POST.get("search_or", "None")
                print("odor/static/media/db_prots/"+docking_input_OR, str(BASE_DIR) + '/media/docking/prot.pdb')
                shutil.copyfile("odor/static/media/db_prots/"+docking_input_OR, str(BASE_DIR) + '/media/docking/prot.pdb')
            except:
                path_prot = None
                return render(request, "OdorWebSite_Docking.html", context={"olfactory_receptors": olfactory_receptors,
                                                                            "chemicals_odors": chemicals_odors,
                                                                            'error_message': "Error inputs"})


        url = get_url_dockign_seamdock(str(BASE_DIR) + "/media/docking/lig.pdb", \
                                       str(BASE_DIR) + "/media/docking/prot.pdb" , website_param = params)
        #print(url)s
        return redirect(url)
    return render(request, "OdorWebSite_Docking.html", context={"olfactory_receptors": olfactory_receptors,
                                                                "chemicals_odors": chemicals_odors,
                                                                'error_message': None})

def ketcher(request):
    return render(request, "ketcher-2.4.0/index.html")