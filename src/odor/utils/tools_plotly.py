import pandas as pd
#import plotly
from plotly.offline import plot
import plotly.express as px
import plotly.graph_objs as go
from rdkit import Chem

from .get_descriptors import get_decriptors_rdkit_fr
DIC_FRAGMENTS = {
# Oxygens
'fr_C_O' : "Number of carbonyl O",
'fr_Al_OH' : "Number of aliphatic hydroxyl groups",
'fr_Ar_OH' : "Number of aromatic hydroxyl groups",
'fr_methoxy' : "Number of methoxy groups -OCH3",
'fr_ester' : "Number of esters",
'fr_Al_COO' : "Number of aliphatic carboxylic acids",
'fr_Ar_COO' : "Number of Aromatic carboxylic acide",
'fr_COO' : "Number of carboxylic acids",
'fr_ketone' : "Number of ketones",
'fr_ether' : "Number of ether oxygens (including phenoxy)",
'fr_phenol' : "Number of phenols",
'fr_aldehyde' : "Number of aldehydes",
# Nitrogens
'fr_NH2' : "Number of carbonyl O",
'fr_NH1' : "Number of carbonyl O",
'fr_NH0' : "Number of carbonyl O",
'fr_Ar_N' : "Number of carbonyl O",
'fr_Ar_NH' : "Number of carbonyl O",
'fr_aniline' : "Number of carbonyl O",
'fr_Imine' : "Number of carbonyl O",
'fr_nitrile' : "Number of carbonyl O",
'fr_nitro' : "Number of carbonyl O",
'fr_amide' : "Number of carbonyl O",
# Halogens
'fr_halogen' : "Number of halogens",
'fr_alkyl_halide' : "Number of alkyl halides",
# Sulfurs
'fr_sulfide' : "Number of thioether",
'fr_SH' : "Number of thiol groups",
# Miscellaneous Functional Groups
'fr_furan' : "Number of furan rings",
'fr_pyridine' : "Number of pyridine rings",
'fr_piperdine' : "Number of piperdine rings",
'fr_morpholine' : "Number of morpholine rings",
'fr_lactone' : "Number of cyclic esters (lactones)",
'fr_epoxide' : "Number of epoxide rings",
'fr_unbrch_alkane' : "Number of unbranched alkanes  of at least 4 members (excludes halogenated alkanes)",
'fr_bicyclic' : "Bicyclic",
'fr_benzene' : "Number of benzene rings",
# Topliss Metabolism
#'fr_nitro_arom_nonortho' : "Number of non-ortho nitro benzene ring substituents",
#'fr_phenol_noOrthoHbond' : "Number of phenolic OH excluding ortho intramolecular Hbond substituents",
#'fr_Al_OH_noTert' : "Number of aliphatic hydroxyl groups excluding tert-OH",
#'fr_para_hydroxylation' : "Number of para-hydroxylation sites",
#'fr_allylic_oxid' : "Number of allylic oxidation sites excluding steroid dienone",
#'fr_aryl_methyl' : "Number of aryl methyl sites for hydroxylation",
#'fr_Ndealkylation2' : "Number of tert-alicyclic amines (no heteroatoms, not quinine-like bridged N)",
#'fr_ketone_Topliss' : "Number of ketones excluding diaryl, a,b-unsat. dienones, heteroatom on Calpha",
#'fr_ArN' : "Number of N functional groups attached to aromatics"
}


def get_radar_plot_from_list_odor(df):
    """
    fig = go.Figure(data=go.Scatterpolar(
        r=[1, 5, 2, 2, 3],
        theta=['processing cost', 'mechanical properties', 'chemical stability', 'thermal stability',
               'device integration'],
        fill='toself'
    ))

    fig.update_layout(
        polar=dict(
            radialaxis=dict(
                visible=True
            ),
        ),
        showlegend=False
    )


    df = px.data.wind()
    print(df)
    fig = px.line_polar(df, r="frequency", theta="direction", color="strength", line_close=True)
    """
    svg_test = '<svg xmlns="http://www.w3.org/2000/svg"><g><rect x="0" y="0" width="100" height="100" fill="red"></rect><text x="0" y="50" font-family="Verdana" font-size="35" fill="blue">Hello</text></g></svg>'

    fig = px.line_polar(df, r="odor_values", theta="Fragment", color="odor",line_close=True, custom_data=['Fragment', 'Fragment_desc', 'odor_count_occ'])
                        #range_r = [-1,1],)

                        #hover_data={"Fragment": True,
                        #            "Fragment_desc":"Fragment description:%s",
                        #            "odor": True,
                        #            "odor_values":"Occurence:.2f"
                        #            })

    #fig.update_traces(hovertemplate='test: %{r}'+svg_test)
    fig.update_traces(    hovertemplate="<br>".join([
                                        "Fragment: %{customdata[0]}",
                                        "Fragment description: %{customdata[1]}",
                                        "Occurrence: %{customdata[2]:.2f} %"
    ]))
    #plotly.offline.plot(fig, filename='media/test.html')
    fig.write_image("/home/guillaumeolt/CMPLI/Projet_odeurs/Results/BDD_odors/plotly_radarplot_odor.svg")
    return plot(fig, output_type="div")

def get_data_desc_plotly_list_odor(db_dict,list_odors):
    df = None
    #idx = 0
    for odor in list_odors:
        list_mol = []
        for chem in db_dict:
            try:
                if odor == "All" or odor in chem["smell"].split(";"):
                    mol = Chem.MolFromSmiles(chem['SMILE'])
                    mol.SetProp("_Name",chem["Name"])
                    mol.SetProp("_idChemicals", str(chem["idChemicals"]))
                    try:
                        mol.SetProp("_smell", chem["smell"])
                    except:
                        pass
                    list_mol.append(mol)
            except:
                pass
        data_desc = get_decriptors_rdkit_fr(list_mol, list_fr = DIC_FRAGMENTS.keys())
        data_desc[data_desc != 0] = 1
        if odor == "All":
            data_desc_all = data_desc.mean(axis=0)
            dict_data_all = data_desc_all.to_dict()
        if df is None:
            data_desc_mean = data_desc.mean(axis=0)
            list_fr_desc = [DIC_FRAGMENTS[i] for i in data_desc_mean.keys()]
            df = pd.DataFrame({'Fragment': data_desc_mean.index,
                               "Fragment_desc": list_fr_desc,
                                 'odor_values': data_desc_mean.values,
                                 'odor_count': data_desc.sum(axis=0).values,
                               'odor_count_occ': data_desc.sum(axis=0).values / len(list_mol) * 100,
                                 'odor_count_all_odor': [len(list_mol) for i in range(len(data_desc_mean.values))],
                                 'odor': [odor for i in range(len(data_desc_mean.values))]})
                              #index = range(idx, idx+len(data_desc.values)))
            #idx = idx + len(data_desc.values)
        else:
            data_desc_mean = data_desc.mean(axis=0)
            list_fr_desc = [DIC_FRAGMENTS[i] for i in data_desc_mean.keys()]
            df_bis = pd.DataFrame({'Fragment': data_desc_mean.index,
                                   "Fragment_desc":list_fr_desc,
                                 'odor_values': data_desc_mean.values,
                                 'odor_count': data_desc.sum(axis=0).values,
                                 'odor_count_occ': data_desc.sum(axis=0).values/len(list_mol)*100,
                                 'odor_count_all_odor': [len(list_mol) for i in range(len(data_desc_mean.values))],
                                 'odor': [odor for i in range(len(data_desc_mean.values))]})
                                  #index = range(idx, idx+len(data_desc.values)))
            #idx = idx + len(data_desc.values)
            df = pd.concat([df,df_bis])


    # REMOVE ROWS WHERE ALL SCORE = 0

    #df = df.drop(df[(df['odor_values'] == 0) & (df['odor'] == "All")].index)
    #df = df.drop(df[(df['odor_count'] < 20) & (df['odor'] == "All")].index)
    # Change values to differences against all dataset

    df.reset_index(drop=True ,inplace=True)
    for index, row in df.iterrows():


        #print(row['descriptors'],"_", dict_data_all[row['descriptors']], index)
        #print(row['odor_values'], row['descriptors'], dict_data_all[row['descriptors']] )

        df.loc[index,'odor_values'] = row['odor_values'] - dict_data_all[row['Fragment']]



    return(df)
    #get_decriptors_rdkit
def get_radar_plot_from_list_or(df):

    fig = px.line_polar(df, r="or_values", theta="Fragment", color="or",line_close=True, custom_data=['Fragment', 'Fragment_desc', 'or_count_occ'])
                        #range_r = [-1,1],
                        #hover_data=["Fragment", "or_values", "or_count", "or_count_all_or"])
    fig.update_traces(    hovertemplate="<br>".join([
                                        "Fragment: %{customdata[0]}",
                                        "Fragment description: %{customdata[1]}",
                                        "Occurrence: %{customdata[2]:.2f} %"
    ]))
    fig.write_image("/home/guillaumeolt/CMPLI/Projet_odeurs/Results/BDD_odors/plotly_radarplot_or.svg")

    """
        fig = px.line_polar(df, r="odor_values", theta="Fragment", color="odor",line_close=True, custom_data=['Fragment', 'Fragment_desc', 'odor_count_occ'])
                        #range_r = [-1,1],)

                        #hover_data={"Fragment": True,
                        #            "Fragment_desc":"Fragment description:%s",
                        #            "odor": True,
                        #            "odor_values":"Occurence:.2f"
                        #            })

    #fig.update_traces(hovertemplate='test: %{r}'+svg_test)
    fig.update_traces(    hovertemplate="<br>".join([
                                        "Fragment: %{customdata[0]}",
                                        "Fragment description: %{customdata[1]}",
                                        "Occurrence: %{customdata[2]:.2f} %"
    ]))
    """


    return plot(fig, output_type="div")
def get_data_desc_plotly_list_or(db_dict,list_or):
    df = None
    #idx = 0
    n=0
    for OR in list_or:
        list_mol = []
        for chem in db_dict:
            try:
                if OR == "All" or OR in chem["OlfRecept"]:
                    n += 1
                    mol = Chem.MolFromSmiles(chem['SMILE'])
                    mol.SetProp("_Name",chem["Name"])
                    mol.SetProp("_idChemicals", str(chem["idChemicals"]))
                    try:
                        mol.SetProp("_OlfRecept", ";".join(chem["OlfRecept"]))
                    except:
                        pass
                    list_mol.append(mol)
            except:
                pass
        data_desc = get_decriptors_rdkit_fr(list_mol, list_fr = DIC_FRAGMENTS.keys())
        data_desc[data_desc != 0] = 1
        if OR == "All":
            data_desc_all = data_desc.mean(axis=0)
            dict_data_all = data_desc_all.to_dict()
        if df is None:
            data_desc_mean = data_desc.mean(axis=0)
            list_fr_desc = [DIC_FRAGMENTS[i] for i in data_desc_mean.keys()]
            df = pd.DataFrame({'Fragment': data_desc_mean.index,
                               "Fragment_desc": list_fr_desc,
                                 'or_values': data_desc_mean.values,
                                 'or_count': data_desc.sum(axis=0).values,
                                 'or_count_all_or': [len(list_mol) for i in range(len(data_desc_mean.values))],
                               'or_count_occ': data_desc.sum(axis=0).values/len(list_mol)*100,
                                 'or': [OR for i in range(len(data_desc_mean.values))]})
                              #index = range(idx, idx+len(data_desc.values)))
            #idx = idx + len(data_desc.values)
        else:
            data_desc_mean = data_desc.mean(axis=0)
            list_fr_desc = [DIC_FRAGMENTS[i] for i in data_desc_mean.keys()]
            df_bis = pd.DataFrame({'Fragment': data_desc_mean.index,
                                   "Fragment_desc": list_fr_desc,
                                 'or_values': data_desc_mean.values,
                                 'or_count': data_desc.sum(axis=0).values,
                                 'or_count_all_or': [len(list_mol) for i in range(len(data_desc_mean.values))],
                                   'or_count_occ': data_desc.sum(axis=0).values / len(list_mol) * 100,
                                 'or': [OR for i in range(len(data_desc_mean.values))]})
                                  #index = range(idx, idx+len(data_desc.values)))
            #idx = idx + len(data_desc.values)
            df = pd.concat([df,df_bis])

    #print(df)
    # REMOVE ROWS WHERE ALL SCORE = 0
    #print(df['or_values'])

    # df = df.drop(df[(df['or_values'] == 0) & (df['or'] == "All")].index)

    #df = df.drop(df[(df['or_count'] < 20) & (df['or'] == "All")].index)
    # Change values to differences against all dataset
    #print(df)
    df.reset_index(drop=True ,inplace=True)
    for index, row in df.iterrows():
        #print(index, df.loc[index,:])

        #print(row['descriptors'],"_", dict_data_all[row['descriptors']], index)
        #print(row['odor_values'], row['descriptors'], dict_data_all[row['descriptors']] )

        df.loc[index,'or_values'] = row['or_values'] - dict_data_all[row['Fragment']]

    #print(dict_data_all)

    return(df)
    #get_decriptors_rdkit