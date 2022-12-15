from django.apps import AppConfig

class OdorConfig(AppConfig):
    default_auto_field = 'django.db.models.BigAutoField'
    name = 'odor'
    def ready(self):
        """
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
        from os import listdir, path
        import os
        import argparse
        from io import StringIO
        from rdkit import Chem
        from rdkit.Chem.Draw import rdMolDraw2D, MolToFile
        from .utils.read_write import write_2d_pdb_list, write_3d_pdb_list, save_2d_image_SVG_list
        from .utils.test import tranform_db_dict, get_path_svg_db
        from .utils.tools_plotly import get_data_desc_plotly_list_odor, get_radar_plot_from_list_odor,\
                                        get_data_desc_plotly_list_or, get_radar_plot_from_list_or
        from .utils.tools_bokeh import get_bokeh_plot_odor_from_list_odors, get_bokeh_plot_odor_from_list_odors,\
                                       get_bokeh_plot_odor_from_list_or, get_bokeh_plot_odor_from_list_chem

        from .utils.tools_umap import get_umap_chem_odor, write_umap_chem_odor, load_umap_chem_odor
        from .utils.test import utils_get_phy_tree, tranform_db_dict_iduniprot_bis
        from odorwebsite.settings import BASE_DIR, STATIC_ROOT,STATIC_URL
        from sys import path

        print("launching server")
         #  COMPUTE CHEMICALS IMAGES
        #print(path)
        #return(0)
        # Add SVG

        chemicals_odors = ChemicalsOdors.objects.all()
        list_db = []
        for chem in chemicals_odors:
            print("--",type(chem.SMILE), chem.SMILE)
            try:
                mol = Chem.MolFromSmiles(chem.SMILE)
                mol.SetProp("_Name", str(chem.idChemicals))
                list_db.append(mol)
            except:
                pass

        save_2d_image_SVG_list(list_db,"odor/static/media/db_mols_svg")
        write_2d_pdb_list(list_db,"odor/static/media/db_mols_2d")
        write_3d_pdb_list(list_db,"odor/static/media/db_mols_3d")


        ## MAPPER

        list_db = []
        list_smi = []
        list_name = []
        list_path_svg = []
        list_color = []
        list_legend = []
        list_odors = []
        list_chem_id = []
        chemicals_odors = ChemicalsOdors.objects.all()
        for chem in chemicals_odors:
            print("--",type(chem.SMILE), chem.SMILE)

            mol = Chem.MolFromSmiles(chem.SMILE)
            mol.SetProp("_Name", str(chem.idChemicals))
            list_db.append(mol)
            list_smi.append(chem.SMILE)
            list_name.append(chem.Name)
            list_path_svg.append("odor/static/media/db_mols_svg/"+str(chem.idChemicals)+".svg")
            list_odors.append("")
            list_chem_id.append(chem.idChemicals)


        mapper = get_umap_chem_odor(list_smi)
        write_umap_chem_odor(mapper,"odor/static/media/umap/mapper.pkl")
        
        mapper = load_umap_chem_odor("odor/static/media/umap/mapper.pkl")
        ## COMPUTE ODORS FIGURES
        db_dict_all = my_custom_sql()
        db_dict_all = tranform_db_dict(db_dict_all)

        # Add SVG
        db_dict_all = get_path_svg_db(db_dict_all)

        smell_percepts = Smell_Percepts.objects.all()
        for odor in smell_percepts:
            print(odor.Odor)
            # plotly
            df = get_data_desc_plotly_list_odor(db_dict_all, ["All", odor.Odor])#odor.Odor])
            div_radar_plot = get_radar_plot_from_list_odor(df)

            #open text file
            text_file = open("odor/static/media/plotly_odors/"+odor.Odor+".html", "w")
            #write string to file
            text_file.write(div_radar_plot)
            #close file
            text_file.close()
            # bokeh

            path_svg_add = '../static/media/db_mols_svg/' #os.path.join(STATIC_URL,'media/db_mols_svg/') 
            #path_svg_add = os.path.join("static", "media/db_mols_svg/")
            print(path_svg_add, "-----------")
            script, div = get_bokeh_plot_odor_from_list_odors(mapper, db_dict_all, [odor.Odor], path_svg_add =path_svg_add)
            print(script, div)
            #open text file
            text_file = open("odor/static/media/bokeh_odors/"+odor.Odor+".div", "w")
            #write string to file
            text_file.write(div)
            #close file
            text_file.close()

            #open text file
            text_file = open("odor/static/media/bokeh_odors/"+odor.Odor+".script", "w")
            #write string to file
            text_file.write(script)
            #close file
            text_file.close()
            print(odor.Odor, path_svg_add)

        db_dict_all = my_custom_sql()
        db_dict_all = tranform_db_dict(db_dict_all)

        # Add SVG
        db_dict_all = get_path_svg_db(db_dict_all)
        ## COMPUTE OR FIGURES
        receptors_idOR = OlfactoryReceptors.objects.values('GeneName', 'idOlfactoryReceptors')#'idUniprot')
        # Transform data
        receptors_idOR_dict = tranform_db_dict_iduniprot_bis(receptors_idOR)

        olfactory_receptors = OlfactoryReceptors.objects.all()
        for or_i in olfactory_receptors:
            print(or_i.GeneName)
            # plotly
            df = get_data_desc_plotly_list_or(db_dict_all, ["All", str(or_i.idOlfactoryReceptors)], receptors_idOR_dict)#odor.Odor])
            div_radar_plot = get_radar_plot_from_list_or(df)

            #open text file
            text_file = open("odor/static/media/plotly_or/"+str(or_i.idOlfactoryReceptors)+".html", "w")
            #write string to file
            text_file.write(div_radar_plot)
            #close file
            text_file.close()

            path_svg_add = '../static/media/db_mols_svg/'#os.path.join(STATIC_ROOT,'media/db_mols_svg/')
            script, div = get_bokeh_plot_odor_from_list_or(mapper, db_dict_all, [str(or_i.idOlfactoryReceptors)], receptors_idOR_dict, path_svg_add =path_svg_add)
            print(script, div)
            #open text file
            text_file = open("odor/static/media/bokeh_or/"+str(or_i.idOlfactoryReceptors)+".div", "w")
            #write string to file
            text_file.write(div)
            #close file
            text_file.close()

            #open text file
            text_file = open("odor/static/media/bokeh_or/"+str(or_i.idOlfactoryReceptors)+".script", "w")
            #write string to file
            text_file.write(script)
            #close file
            text_file.close()


        ## COMPUTE PHYLOGENIQUE TREE mouse
        # OR mouse
        olfactory_receptors = OlfactoryReceptors.objects.all()
        for or_i in olfactory_receptors:
            print(or_i.GeneName)
            utils_get_phy_tree(BASE_DIR, query_or=or_i.GeneName,
                                         path_tree="odor/static/media/phylo_tree_PhyML_Or10ad1.tree",
                                         path_output="odor/static/media/phylogenic_tree_OR/"+str(or_i.idOlfactoryReceptors)+".svg",
                                         path_homology="odor/static/media/dt_homology_human_mouse.csv",
                                         species="mouse")

        # CHEM mouse
        db_dict = my_custom_sql()
        # Transform data
        db_dict = tranform_db_dict(db_dict)
        for k in db_dict:
            if k["OlfRecept"] is not None:
                utils_get_phy_tree(BASE_DIR, query_or=";".join(k["OlfRecept"]),
                                             path_tree="odor/static/media/phylo_tree_PhyML_Or10ad1.tree",
                                             path_output="odor/static/media/phylogenic_tree_chem/"+str(k["idChemicals"])+".svg",
                                             path_homology="odor/static/media/dt_homology_human_mouse.csv",
                                             species="mouse")
                print(k["idChemicals"], ";".join(k["OlfRecept"]))

        ## COMPUTE PHYLOGENIQUE TREE human
        # OR mouse
        olfactory_receptors = OlfactoryReceptors.objects.all()
        for or_i in olfactory_receptors:
            print(or_i.GeneName)
            utils_get_phy_tree(BASE_DIR, query_or=or_i.GeneName,
                                         path_tree="odor/static/media/phylo_tree_PhyML_human.tree",
                                         path_output="odor/static/media/phylogenic_tree_human_OR/"+str(or_i.idOlfactoryReceptors)+".svg",
                                         path_homology="odor/static/media/dt_homology_mouse_human.csv",
                                         species="human")

        # CHEM mouse
        db_dict = my_custom_sql()
        # Transform data
        db_dict = tranform_db_dict(db_dict)
        for k in db_dict:
            if k["OlfRecept"] is not None:
                utils_get_phy_tree(BASE_DIR, query_or=";".join(k["OlfRecept"]),
                                             path_tree="odor/static/media/phylo_tree_PhyML_human.tree",
                                             path_output="odor/static/media/phylogenic_tree_human_chem/"+str(k["idChemicals"])+".svg",
                                             path_homology="odor/static/media/dt_homology_mouse_human.csv",
                                             species="human")
                print(k["idChemicals"], ";".join(k["OlfRecept"]))

        ## COMPUTE CHEM UMAP FIGURES
        db_dict_all = my_custom_sql()
        db_dict_all = tranform_db_dict(db_dict_all)

        # Add SVG
        db_dict_all = get_path_svg_db(db_dict_all)
        mapper = load_umap_chem_odor("odor/static/media/umap/mapper.pkl")
        chemicals_odors = ChemicalsOdors.objects.all()
        for chem_i in chemicals_odors:
            path_svg_add = '../static/media/db_mols_svg/'#os.path.join(STATIC_ROOT,'media/db_mols_svg/')
            script, div = get_bokeh_plot_odor_from_list_chem(mapper, db_dict_all, [str(chem_i.idChemicals)], path_svg_add =path_svg_add)
            #print(script, div)
            #open text file
            text_file = open("odor/static/media/bokeh_chem/"+str(chem_i.idChemicals)+".div", "w")
            #write string to file
            text_file.write(div)
            #close file
            text_file.close()

            #open text file
            text_file = open("odor/static/media/bokeh_chem/"+str(chem_i.idChemicals)+".script", "w")
            #write string to file
            text_file.write(script)
            #close file
            text_file.close()
        quit()
        """