import argparse
import json
from prediction import *
import pandas as pd
import itertools

PATH_label_odors_names = '/home/guillaumeolt/CMPLI/Projet_odeurs/web-projects-django/src/ServerOdors/utils/infos_models/label_odors_names.txt'
PATH_GCC_CODB = "/home/guillaumeolt/CMPLI/Projet_odeurs/web-projects-django/src/ServerOdors/utils/Trained_models/GNN_CODB_ckpt"
PATH_GCC_CODB = "/home/guillaumeolt/CMPLI/Projet_odeurs/stage_M2_rayane/Trained_models/GNN_CODB_ckpt"
def utils_get_prediction_odor():
    dbClean = pd.read_excel("/home/guillaumeolt/CMPLI/Projet_odeurs/stage_M2_rayane/Databases/CODB_updated.xlsx")
    L_Odor_Cat = 'ODOR1,ODOR2,ODOR3,ODOR4,ODOR5,ODOR6,ODOR7,ODOR8,ODOR9,ODOR10'.split(',')

    # ### Encodage binaire des odeurs pour chacune des molécules
    y_multilabel, label_odors_names = CODB_Binarizer_full(L_Odor_Cat, dbClean)
    #print(label_odors_names)
    #print(y_multilabel)
    #print(len(label_odors_names))
    tasks = label_odors_names

    # SMILE
    use_chirality = True  # False pour ne pas prendre en compte
    random_state = 2020  # None pour aléatoire
    smiles = dbClean.SMILES.tolist()
    x = Smiles2GraphObj(smiles, use_chirality)
    naf = len(x[0].atom_features[0])
    # Restore model
    model = GNNPerso(len(tasks), batch_size=32, graph_conv_layers=[64, 64], mode='classification', dropout=0.40,
                     tensorboard=True, model_dir="/home/guillaumeolt/CMPLI/Projet_odeurs/stage_M2_rayane/Trained_models/GNN_CODB_ckpt", number_atom_features=naf)
    model.restore()

def predict_odor(smile, label_odors_names, PATH_GCC_CODB):
    """
    For a given smile predict its odors using the GNN_CODB model
    """
    tasks = label_odors_names
    # SMILE
    use_chirality = True  # False pour ne pas prendre en compte
    random_state = 2020  # None pour aléatoire
    #print([smile],"infunction")
    smiles = [smile]
    x = Smiles2GraphObj(smiles, use_chirality)
    naf = len(x[0].atom_features[0])

    # Restore model
    if tf.executing_eagerly():
        tf.compat.v1.disable_eager_execution()
    model = GNNPerso(len(tasks), batch_size=32, graph_conv_layers=[64, 64], mode='classification', dropout=0.40,
                     tensorboard=True,
                     model_dir=PATH_GCC_CODB,
                     number_atom_features=naf)
    model.restore()
    predicted_Ext_Testset = PredictFromSmiles([smile], model, tasks, use_chirality=True)
    dict_pred = get_json_predict(predicted_Ext_Testset, tasks)
    #test_odors_all, t = getOdorExtSet(predicted_Ext_Testset, tasks, pourcentage=True)
    #print(predicted_Ext_Testset, t)
    return(dict_pred)
def predict_OR(smile, label_or_names, PATH_GCC_RCDB):
    """
    For a given smile predict its OR using the GNN_CODB model
    """
    tasks = label_or_names
    # SMILE
    use_chirality = True  # False pour ne pas prendre en compte
    random_state = 2020  # None pour aléatoire
    #print([smile],"infunction")
    smiles = [smile]
    x = Smiles2GraphObj(smiles, use_chirality)
    naf = len(x[0].atom_features[0])

    # Restore model
    if tf.executing_eagerly():
        tf.compat.v1.disable_eager_execution()
    model = GNNPerso(len(tasks), batch_size=32, graph_conv_layers=[64, 64], mode='classification', dropout=0.40,
                     tensorboard=True,
                     model_dir=PATH_GCC_RCDB,
                     number_atom_features=naf)
    model.restore()
    predicted_Ext_Testset = PredictFromSmiles([smile], model, tasks, use_chirality=True)
    dict_pred = get_json_predict(predicted_Ext_Testset, tasks)
    #test_or_all = getOdorExtSet(predicted_Ext_Testset, tasks)
    return(dict_pred)
def get_label_names(path_label_names):
    """
    read file containing  odors or or of the model. return list
    """
    with open(path_label_names) as f:
        list_label = f.read().splitlines()
    return(list_label)

######

if __name__ == "__main__":
    ##Parser to deal with arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-label",\
           help="Path to the file containing the odors to predict",\
           default=PATH_label_odors_names)
    parser.add_argument("-smile",\
           help="String of the smile to predict",\
           default=None)
    parser.add_argument("-PATH_GCC",\
           help="PATH to the model to use",\
           default=PATH_GCC_CODB)
    parser.add_argument("-predict",\
           help="which model to use",\
           default="odor")
    args = parser.parse_args()

    #test_odors_all = utils_get_prediction_odor()
    if args.predict == "odor":
        if args.smile != None:
            label_odors_names = get_label_names(args.label)
            predicted_odors = predict_odor(args.smile, label_odors_names, args.PATH_GCC)
            #print(predicted_odors,','.join(predicted_odors))
            #print(f','.join(predicted_odors[0]))#, ','.join(list(itertools.chain.from_iterable(predicted_odors))))
            #print(json.dumps(str(predicted_odors)))
            print(predicted_odors)
            #print(str(predicted_odors))
            quit()
    if args.predict == "or":
        if args.smile != None:
            label_or_names = get_label_names(args.label)
            predicted_or = predict_OR(args.smile, label_or_names, args.PATH_GCC)
            #print(predicted_odors,','.join(predicted_odors))
            #print(f','.join(predicted_or[0]))#, ','.join(list(itertools.chain.from_iterable(predicted_odors))))
            #print(json.dumps(str(predicted_or)))
            print(predicted_or)
            #print(str(predicted_or))
            quit()