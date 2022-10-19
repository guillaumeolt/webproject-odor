import collections
import copy
import sklearn
import numpy as np
import pandas as pd 
import deepchem as dc
import tensorflow as tf
#import matplotlib.pyplot as mlp
#
from deepchem.models import GraphConvModel
from deepchem.data.datasets import NumpyDataset

from deepchem.metrics import to_one_hot
from deepchem.feat.graph_features import ConvMolFeaturizer
from deepchem.feat.mol_graphs import ConvMol
from deepchem.models import KerasModel, layers
from deepchem.models.losses import L2Loss, SoftmaxCrossEntropy
from tensorflow.keras.layers import Input, Dense, Reshape, Softmax, Dropout, Activation, BatchNormalization
#
from rdkit import Chem
#
from sklearn.preprocessing import StandardScaler, MultiLabelBinarizer, Normalizer
from sklearn.model_selection import train_test_split
from sklearn.feature_selection import VarianceThreshold



class TrimGraphOutput(tf.keras.layers.Layer):
    
    """Trim the output to the correct number of samples.
    GraphGather always outputs fixed size batches.  This layer trims the output
    to the number of samples that were in the actual input tensors.
    """

    def __init__(self, **kwargs):
        super(TrimGraphOutput, self).__init__(**kwargs)

    def call(self, inputs):
        n_samples = tf.squeeze(inputs[1])
        return inputs[0][0:n_samples]


class GNNPerso(KerasModel):

    def __init__(self,
                 n_tasks,
                 graph_conv_layers=[64, 64],
                 dense_layer_size=128,
                 dropout=0.0,mode="classification",
                 number_atom_features=75,
                 n_classes=2,
                 uncertainty=False,
                 batch_size=100,
                 **kwargs):
        """
        Parameters
        ----------
        n_tasks: int
          Number of tasks
        graph_conv_layers: list of int
          Width of channels for the Graph Convolution Layers
        dense_layer_size: int
          Width of channels for Atom Level Dense Layer before GraphPool
        dropout: list or float
          the dropout probablity to use for each layer.  The length of this list should equal
          len(graph_conv_layers)+1 (one value for each convolution layer, and one for the
          dense layer).  Alternatively this may be a single value instead of a list, in which
          case the same value is used for every layer.
        mode: str
          Either "classification" or "regression"
        number_atom_features: int
            75 is the default number of atom features created, but
            this can vary if various options are passed to the
            function atom_features in graph_features
        n_classes: int
          the number of classes to predict (only used in classification mode)
        uncertainty: bool
          if True, include extra outputs and loss terms to enable the uncertainty
          in outputs to be predicted
        """
        if mode not in ['classification', 'regression']:
            raise ValueError("mode must be either 'classification' or 'regression'")
        self.n_tasks = n_tasks
        self.mode = mode
        self.dense_layer_size = dense_layer_size
        self.graph_conv_layers = graph_conv_layers
        self.number_atom_features = number_atom_features
        self.n_classes = n_classes
        self.uncertainty = uncertainty
        if not isinstance(dropout, collections.Sequence):
            dropout = [dropout] * (len(graph_conv_layers) + 1)
        if len(dropout) != len(graph_conv_layers) + 1:
            raise ValueError('Wrong number of dropout probabilities provided')
        self.dropout = dropout
        if uncertainty:
            if mode != "regression":
                raise ValueError("Uncertainty is only supported in regression mode")
            if any(d == 0.0 for d in dropout):
                raise ValueError(
            'Dropout must be included in every layer to predict uncertainty')

        # Build the model.

        atom_features = Input(shape=(self.number_atom_features,))
        degree_slice = Input(shape=(2,), dtype=tf.int32)
        membership = Input(shape=tuple(), dtype=tf.int32)
        n_samples = Input(shape=tuple(), dtype=tf.int32)
        dropout_switch = tf.keras.Input(shape=tuple())

        self.deg_adjs = []
        for i in range(0, 10 + 1):
            deg_adj = Input(shape=(i + 1,), dtype=tf.int32)
            self.deg_adjs.append(deg_adj)
        in_layer = atom_features
        for layer_size, dropout in zip(self.graph_conv_layers, self.dropout):
            gc1_in = [in_layer, degree_slice, membership] + self.deg_adjs
            gc1 = layers.GraphConv(layer_size, activation_fn=tf.nn.relu)(gc1_in)
            batch_norm1 = BatchNormalization(fused=False)(gc1)
            if dropout > 0.0:
                batch_norm1 = layers.SwitchedDropout(rate=dropout)(
                    [batch_norm1, dropout_switch])
            gp_in = [batch_norm1, degree_slice, membership] + self.deg_adjs
            in_layer = layers.GraphPool()(gp_in)
        dense = Dense(self.dense_layer_size, activation=tf.nn.relu)(in_layer)
        batch_norm3 = BatchNormalization(fused=False)(dense)
        if self.dropout[-1] > 0.0:
            batch_norm3 = layers.SwitchedDropout(rate=self.dropout[-1])(
                [batch_norm3, dropout_switch])
        self.neural_fingerprint = layers.GraphGather(
            batch_size=batch_size,
            activation_fn=tf.nn.tanh)([batch_norm3, degree_slice, membership] +
                                      self.deg_adjs)

        n_tasks = self.n_tasks
        if self.mode == 'classification':
            n_classes = self.n_classes
            logits = Reshape((n_tasks, n_classes))(Dense(n_tasks * n_classes)(
                self.neural_fingerprint))
            logits = TrimGraphOutput()([logits, n_samples])
            output = Softmax()(logits)
            outputs = [output, logits]
            output_types = ['prediction', 'loss']
            loss = SoftmaxCrossEntropy()
        else:
            output = Dense(n_tasks)(self.neural_fingerprint)
            output = TrimGraphOutput()([output, n_samples])
            if self.uncertainty:
                log_var = Dense(n_tasks)(self.neural_fingerprint)
                log_var = TrimGraphOutput()([log_var, n_samples])
                var = Activation(tf.exp)(log_var)
                outputs = [output, var, output, log_var]
                output_types = ['prediction', 'variance', 'loss', 'loss']

                def loss(outputs, labels, weights):
                    diff = labels[0] - outputs[0]
                    return tf.reduce_mean(diff * diff / tf.exp(outputs[1]) + outputs[1])
            else:
                outputs = [output]
                output_types = ['prediction']
                loss = L2Loss()
        model = tf.keras.Model(
            inputs=[
                atom_features, degree_slice, membership, n_samples, dropout_switch
            ] + self.deg_adjs,
            outputs=outputs)
        super(GNNPerso, self).__init__(
            model, loss, output_types=output_types, batch_size=batch_size, **kwargs)

    def default_generator(self,
                        dataset,
                        epochs=1,
                        mode='fit',
                        deterministic=True,
                        pad_batches=True):
        for epoch in range(epochs):
            for (X_b, y_b, w_b, ids_b) in dataset.iterbatches(
                batch_size=self.batch_size,
                deterministic=deterministic,
                pad_batches=pad_batches):
                if self.mode == 'classification':
                    y_b = to_one_hot(y_b.flatten(), self.n_classes).reshape(
                        -1, self.n_tasks, self.n_classes)
                multiConvMol = ConvMol.agglomerate_mols(X_b)
                n_samples = np.array(X_b.shape[0])
                if mode == 'predict':
                    dropout = np.array(0.0)
                else:
                    dropout = np.array(1.0)
                inputs = [
                    multiConvMol.get_atom_features(), multiConvMol.deg_slice,
                    np.array(multiConvMol.membership), n_samples, dropout
                ]
                for i in range(1, len(multiConvMol.get_deg_adjacency_lists())):
                    inputs.append(multiConvMol.get_deg_adjacency_lists()[i])
                yield (inputs, [y_b], [w_b])

########################################################"""
def CODB_Binarizer_full(L_Odor_Cat,dbClean):
    categories = L_Odor_Cat
    #dataset_df_python[categories]
    L_vectors = []
    for i in range(len(dbClean)) :

        vector  = ()
        for cat in categories :
            odor = dbClean[cat].iloc[i]
            if str(odor) != 'nan' :
                vector += (odor,)
        L_vectors.append(vector)
    ##
    ##


    mlb = MultiLabelBinarizer()
    y_multilabel = mlb.fit_transform(L_vectors)
    label_odors_names = mlb.classes_
    return(y_multilabel,label_odors_names)
    
    
#####################

def Smiles2GraphObj(smiles, use_chirality = False):
    """Génère une list d'objet de type graph moléculaire à partir d'une liste de smiles"""

    # Conversion du format de molécule 'Smiles' au format 'Mol'
    mols = [Chem.MolFromSmiles(s) for s in smiles]

    #Featurizer de type ConvMol
    featurizer = dc.feat.ConvMolFeaturizer(use_chirality=use_chirality)
    x = featurizer.featurize(mols)
    return(x)
    
    
###################


def SmiFromOneOdor(odortest = 'popcorn'):
    
    L_Odor_Cat = 'ODOR1,ODOR2,ODOR3,ODOR4,ODOR5,ODOR6,ODOR7,ODOR8,ODOR9,ODOR10'.split(',')
    smilestotest = []
    
    for i in range(len(dbClean)):
        
        
        if odortest in dbClean.iloc[i].tolist() :
            vector = ()
            stti = dbClean.iloc[i]['SMILES']
            #print(stti)
            smilestotest.append(stti)
            
            for cat in L_Odor_Cat :
                odor = dbClean[cat].iloc[i]
                if str(odor) != 'nan' :
                    vector += (odor,)
            L_vectors.append(vector)
            
    mlb = MultiLabelBinarizer()
    y = mlb.fit_transform(L_vectors)
    odor_names = mlb.classes_        
    return smilestotest, y, odor_names
    #print('END'))
    
    
####### 
def train_GNN(model, train_dataset, num_epochs, plotlosses=True) :
    losses = []
    for i in range(num_epochs):
        loss = model.fit(train_dataset, nb_epoch=1)
        print("Epoch %d loss: %f" % (i, loss))
        losses.append(loss)
    if plotlosses == True :
        mlp.plot(range(num_epochs), losses)
#######
def test_GNN(model, test_dataset,  train_dataset=None) :
    metrics = [dc.metrics.Metric(dc.metrics.roc_auc_score, np.mean, mode="classification"),
               dc.metrics.Metric(dc.metrics.prc_auc_score, np.mean, mode="classification")]
               #dc.metrics.Metric(dc.metrics.balanced_accuracy_score, np.mean, mode="classification")]
    
    if train_dataset :
        train_scores = model.evaluate(train_dataset, metrics, per_task_metrics=True, )
        test_scores = model.evaluate(test_dataset, metrics, per_task_metrics=True)
        return(train_scores,test_scores)
    
    else :
        test_scores = model.evaluate(test_dataset, metrics, per_task_metrics=True)
        return(test_scores)
########



def PredictOneOdorMols(odor, model, use_chirality=False):
    stt, y_stt, odor_names_stt = SmiFromOneOdor(odor)

    molstotest = [Chem.MolFromSmiles(smi) for smi in stt]
    featurizer = dc.feat.ConvMolFeaturizer(use_chirality=use_chirality)
    xtt = featurizer.featurize(molstotest)

    xtt_dataset = NumpyDataset(xtt, np.array(len(xtt)*[len(tasks)*[0]]))
    xtt_dataset_val = NumpyDataset(xtt, y_stt )

    predicted_odors = model.predict(xtt_dataset)
    #predicted_odors_val = model.predict(xtt_dataset_val)
    return(predicted_odors)

###########

def PredictFromSmiles(smiles, model, tasks, use_chirality=False):

    molstotest = [Chem.MolFromSmiles(smi) for smi in smiles]
    featurizer = dc.feat.ConvMolFeaturizer(use_chirality=use_chirality)
    xtt = featurizer.featurize(molstotest)

    xtt_dataset = NumpyDataset(xtt, np.array(len(xtt)*[len(tasks)*[0]]))
    #xtt_dataset_val = NumpyDataset(xtt, y_stt )

    predicted_odors = model.predict(xtt_dataset)
    #predicted_odors_val = model.predict(xtt_dataset_val)
    return(predicted_odors)


#############

def get_json_predict(predicted_Ext_Testset, tasks):
    dict_pred=dict()
    all_pred_odors=[]
    for i, molecule in enumerate(predicted_Ext_Testset):
        L_pred_odors = []
        for j,odpred in enumerate(molecule) :
            for k,binval in enumerate(odpred) :
                if binval >= 0.5 and k == 1 :
                    dict_pred[tasks[j]] = round(binval,2)
                    L_pred_odors.append([tasks[j], round(binval,2)])
        all_pred_odors.append(L_pred_odors)
    return(dict_pred)

def getOdorExtSet(predicted_Ext_Testset, tasks,  pourcentage=False):
    all_pred_odors=[]
    for i, molecule in enumerate(predicted_Ext_Testset):
        L_pred_odors = []
        for j,odpred in enumerate(molecule) :
            for k,binval in enumerate(odpred) :
                if binval >= 0.5 and k == 1 :
                    #print(k,binval)
                    #print(j,label_odors_names[j])
                    #if label_odors_names[j] == 'popcorn' :
                    #    print('OK')
                    L_pred_odors.append([tasks[j], round(binval,2)])
        all_pred_odors.append(L_pred_odors)
        #print( L_pred_odors,'\n')
        #stt[i],


    test_odors_all = []
    for i in all_pred_odors:
        test_odors_i=[]
        for j in i :
            #print(j[0])
            test_odors_i.append(j[0])
        test_odors_all.append(test_odors_i)

    if pourcentage == True :
        return(all_pred_odors, test_odors_all)
    else :
        return(test_odors_all)
########



def create_dico_RCDB(Train_RCDB_DB) :
    L_smiles = []
    dico_train_RCDB = {}
    full_recep_list = Train_RCDB_DB['Recepteur olfactif'].tolist()+Train_RCDB_DB['Récepteur homologue humain'].tolist()
    for i in range(len(Train_RCDB_DB)):
        smiles_i = Train_RCDB_DB.iloc[i]['SMILES']
        rec_hmlg = Train_RCDB_DB.iloc[i]['Récepteur homologue humain']
        rec  = Train_RCDB_DB['Recepteur olfactif'].iloc[i]
        compound = Train_RCDB_DB['Nom'].iloc[i]

        #print(rec_hmlg, rec,compound )
        if compound in dico_train_RCDB :
            if str(rec_hmlg) != 'nan'  and  full_recep_list.count(rec_hmlg)>=1 :
                if rec_hmlg not in dico_train_RCDB[compound] :
                    dico_train_RCDB[compound].append(rec_hmlg)
                if smiles_i not in L_smiles : 
                        L_smiles.append(smiles_i)
            else :
                if full_recep_list.count(rec)>=1 :
                    if rec not in dico_train_RCDB[compound] :
                        dico_train_RCDB[compound].append(rec)
                    if smiles_i not in L_smiles : 
                        L_smiles.append(smiles_i)
        else :
            if str(rec_hmlg) != 'nan' and  full_recep_list.count(rec_hmlg)>=1 :
                dico_train_RCDB[compound] = [rec_hmlg]
                if smiles_i not in L_smiles : 
                        L_smiles.append(smiles_i)
            else :
                if full_recep_list.count(rec)>=1 : 
                    dico_train_RCDB[compound] =  [rec]
                    if smiles_i not in L_smiles : 
                        L_smiles.append(smiles_i)
    return(dico_train_RCDB, L_smiles)
##############
def create_dico_RCDB_v2(RCDB) :
    dico = {}
    L_smiles = []
    for i in range(len(Full_RCDB)) :
        name = Full_RCDB['Nom'].iloc[i]
        recep = Full_RCDB['Recepteur olfactif'].iloc[i]
        recep_hum = Full_RCDB['Synonyme Recepteur olfactif'].iloc[i]
        smile = Full_RCDB['SMILES'].iloc[i]
        if name in dico :
        
            if type(recep_hum) == type('str') :
                dico[name] += [recep_hum]
            else :
                dico[name] += [recep]
        else :
            L_smiles.append(smile)
            if type(recep_hum) == type('str') :
                dico[name] = [recep_hum]
            else :
                dico[name] = [recep]
    return(dico, L_smiles)
##############
def create_dico_RCDB_human(Train_RCDB_DB) :
    """Create dict for human + orthologs"""
    L_smiles = []
    dico_train_RCDB = {}
    full_recep_list = Train_RCDB_DB[Train_RCDB_DB.Espèce == 'humain']['Recepteur olfactif'].tolist()+Train_RCDB_DB['Récepteur homologue humain'].tolist()
    for i in range(len(Train_RCDB_DB)):
        smiles_i = Train_RCDB_DB.iloc[i]['SMILES']
        rec_hmlg = Train_RCDB_DB.iloc[i]['Récepteur homologue humain']
        rec  = Train_RCDB_DB['Recepteur olfactif'].iloc[i]
        compound = Train_RCDB_DB['Nom'].iloc[i]

        #print(rec_hmlg, rec,compound )
        if compound in dico_train_RCDB :
            if str(rec_hmlg) != 'nan'  and  full_recep_list.count(rec_hmlg)>=1 :
                if rec_hmlg not in dico_train_RCDB[compound] :
                    dico_train_RCDB[compound].append(rec_hmlg)
                if smiles_i not in L_smiles : 
                        L_smiles.append(smiles_i)
            else :
                if full_recep_list.count(rec)>=1 :
                    if rec not in dico_train_RCDB[compound] :
                        dico_train_RCDB[compound].append(rec)
                    if smiles_i not in L_smiles : 
                        L_smiles.append(smiles_i)
        else :
            if str(rec_hmlg) != 'nan' and  full_recep_list.count(rec_hmlg)>=1 :
                dico_train_RCDB[compound] = [rec_hmlg]
                if smiles_i not in L_smiles : 
                        L_smiles.append(smiles_i)
            else :
                if full_recep_list.count(rec)>=1 : 
                    dico_train_RCDB[compound] =  [rec]
                    if smiles_i not in L_smiles : 
                        L_smiles.append(smiles_i)
    return(dico_train_RCDB, L_smiles)
##############

##############
def create_dico_RCDB_human_only(Train_RCDB_DB) :
    L_smiles = []
    dico_train_RCDB = {}
    RCDB_df = Train_RCDB_DB[Train_RCDB_DB.Espèce == 'humain']
    full_recep_list = RCDB_df['Recepteur olfactif'].tolist()
    for i in range(len(RCDB_df)):
        smiles_i = RCDB_df.iloc[i]['SMILES']
        rec  = RCDB_df['Recepteur olfactif'].iloc[i]
        compound = RCDB_df['Nom'].iloc[i]

        #print(rec_hmlg, rec,compound )
        if compound in dico_train_RCDB :
            if full_recep_list.count(rec)>=1 :
                if rec not in dico_train_RCDB[compound] :
                    dico_train_RCDB[compound].append(rec)
                if smiles_i not in L_smiles : 
                    L_smiles.append(smiles_i)
        else :
            if full_recep_list.count(rec)>=1 : 
                dico_train_RCDB[compound] =  [rec]
                if smiles_i not in L_smiles : 
                    L_smiles.append(smiles_i)
    return(dico_train_RCDB, L_smiles)
##############
def RCDB_Binarizer(dico):
    L_vectors = []

    for i in dico:
        L_vectors.append(dico[i])
        
    mlb = MultiLabelBinarizer()
    y_multilabel = mlb.fit_transform(L_vectors)
    label_recepteurs_names = mlb.classes_    
    return(y_multilabel, label_recepteurs_names)
############
def create_train_Dataset_fromRCDB(Train_RCDB_DB, use_chirality=False):
    """Cette fonction prend en entrée un jeu de données RCDB et renvoi un Dataset qui servira à entrainer le modèle GNN.
    On utilisera cette fonction pour entrainer le modèle sur l'ensemble du jeu de données"""
    
    dico_full_RCDB, SMILES_Train = create_dico_RCDB(Train_RCDB_DB)
    
    print('dico', len(dico_full_RCDB))
    
    # Création de graphes à partir des smiles
    x = Smiles2GraphObj(SMILES_Train, use_chirality=use_chirality)
    print('x', len(x))
    '''try :
        mols = []
        for i,s in enumerate(SMILES_Train):
            mols.append(Chem.MolFromSmiles(s))
    except :
        print('Erreur pour le SMILES  : ',i, s, 'Vérifier le jeu de données')

    #Featurizer de type ConvMol
    featurizer = dc.feat.ConvMolFeaturizer()
    x = featurizer.featurize(mols)'''    
    
    # Encodage binaire des labels
    '''L_vectors = []
    for i in dico_full_RCDB:
        L_vectors.append(dico_full_RCDB[i])

    mlba = MultiLabelBinarizer()
    y_multilabel = mlba.fit_transform(L_vectors)
    label_recepteurs_names = mlba.classes_'''
    y_multilabel, label_recepteurs_names = RCDB_Binarizer(dico_full_RCDB)
    
    return(x,y_multilabel,label_recepteurs_names )
#############
def create_train_Dataset_fromRCDB_human(Train_RCDB_DB, use_chirality=False):
    """Cette fonction prend en entrée un jeu de données RCDB et renvoi un Dataset qui servira à entrainer le modèle GNN.
    On utilisera cette fonction pour entrainer le modèle sur l'ensemble du jeu de données"""
    
    dico_full_RCDB, SMILES_Train = create_dico_RCDB_human(Train_RCDB_DB)  #create_dico_RCDB_human
    
    print('dico', len(dico_full_RCDB))
    
    # Création de graphes à partir des smiles
    x = Smiles2GraphObj(SMILES_Train, use_chirality=use_chirality)
    print('x', len(x))
    
    y_multilabel, label_recepteurs_names = RCDB_Binarizer(dico_full_RCDB)    
    return(x,y_multilabel,label_recepteurs_names )
################

#############
def create_train_Dataset_fromRCDB_human_only(Train_RCDB_DB, use_chirality=False):
    """Cette fonction prend en entrée un jeu de données RCDB et renvoi un Dataset qui servira à entrainer le modèle GNN.
    On utilisera cette fonction pour entrainer le modèle sur l'ensemble du jeu de données"""
    
    dico_full_RCDB, SMILES_Train = create_dico_RCDB_human_only(Train_RCDB_DB)  #create_dico_RCDB_human
    
    print('dico', len(dico_full_RCDB))
    
    # Création de graphes à partir des smiles
    x = Smiles2GraphObj(SMILES_Train, use_chirality=use_chirality)
    print('x', len(x))
    
    y_multilabel, label_recepteurs_names = RCDB_Binarizer(dico_full_RCDB)    
    return(x,y_multilabel,label_recepteurs_names )
################
################


def create_dataset_fromdic(dico_RCDB, Full_RCDB, use_chirality=False):
    L_smiles_dico=[]
    for comp in dico_RCDB:
        for i in range(len(Full_RCDB)):
            if comp == Full_RCDB.iloc[i]['Nom']:
                L_smiles_dico.append(Full_RCDB.iloc[i]['SMILES'])
                print(comp,Full_RCDB.iloc[i]['SMILES'])
                break
    # Création de graphes à partir des smiles
    try :
        mols = []
        for i,s in enumerate(L_smiles_dico):
            mols.append(Chem.MolFromSmiles(s))
    except :
        print('Erreur pour le SMILES  : ',i, s, 'Vérifier le jeu de données')

    #Featurizer de type ConvMol
    featurizer = dc.feat.ConvMolFeaturizer(use_chirality=use_chirality)
    x = featurizer.featurize(mols)    

    # Encodage binaire des labels
    L_vectors = []
    for i in dico_RCDB:
        L_vectors.append(dico_RCDB[i])

    mlb = MultiLabelBinarizer()
    y_multilabel = mlb.fit_transform(L_vectors)
    label_recepteurs_names = mlb.classes_


    return(x,y_multilabel,label_recepteurs_names )
##########
##########

def CV_GNN(Full_Dataset, kfolds, epochs, tasks, DB, naf) :
    splitter = dc.splits.RandomSplitter()
    all_sets = splitter.k_fold_split(Full_Dataset, kfolds)
    scores_test = []
    for i in range(len(all_sets)):
        print('\nFold n°', i+1, ' :  ')
        
        if DB == 'CODB':
            model  = GNNPerso(len(tasks), batch_size=32, mode='classification', dropout=0.47,tensorboard=True, number_atom_features = naf)
        elif DB == 'RCDB' :
            model  = GNNPerso(len(tasks), batch_size=32, mode='classification', dropout=0.25,tensorboard=True, number_atom_features = naf)

        for j in range(len(all_sets)):
            if i == j :
                n_totest = i
            else :
                train_GNN(model, all_sets[j][1],epochs, plotlosses=False)
        test_score_i = test_GNN(model, all_sets[n_totest][1], train_dataset=None )
        scores_test.append(test_score_i)
    
    return(scores_test)
#############

def CV_GNN_2(Full_Dataset, kfolds, epochs, tasks, DB) :
    splitter = dc.splits.RandomStratifiedSplitter()
    all_sets = splitter.k_fold_split(Full_Dataset, kfolds)
    scores_test = []
    for i in range(len(all_sets)):
        print('\nFold n°', i+1, ' :  ')
        
        if DB == 'CODB':
            model  = GNNPerso(len(tasks), batch_size=32, mode='classification', dropout=0.47,tensorboard=True, model_dir='models', number_atom_features = 75)
        elif DB == 'RCDB' :
            model  = GNNPerso(len(tasks), batch_size=32, mode='classification', dropout=0.25,tensorboard=True, model_dir='models', number_atom_features = 78)

        for j in range(len(all_sets)):
            if i == j :
                n_totest = i
            else :
                train_GNN(model, all_sets[j][0],epochs, plotlosses=False)
        test_score_i = test_GNN(model, all_sets[n_totest][1], train_dataset=None )
        scores_test.append(test_score_i)
    
    return(scores_test)

##############
def CV_perso_GNN(preDatasets, epochs) :
    all_sets = [NumpyDataset(datasets[0], datasets[1])for datasets in preDatasets]
    scores_test = []
    dico_index={}
    dico_test={}
    for i in range(4):
        for j in range(4):
            if i == j :
                dico_test[i]=[j]
                pass
            else :
                if i not in dico_index:
                    dico_index[i] = [j]

                else :
                    dico_index[i].append(j)

    for i in dico_index:
        print('\nFold n°', i+1, ' :  ')

        common_recep = []
        for i_dico in dico_index[i]:
            print(i_dico)
            for recep in preDatasets[i_dico][2] :
                if recep not in common_recep :
                    common_recep.append(recep)        
        tasks = common_recep
        model  = GNNPerso(len(tasks), batch_size=32, mode='classification', dropout=0.25,tensorboard=True, model_dir='models', number_atom_features = 75)

        for j in range(len(all_sets)):
            if i == j :
                n_totest = i
            else : 
                train_GNN(model, all_sets[j],epochs, plotlosses=False)
        test_score_i = test_GNN(model, all_sets[n_totest] )
        scores_test.append(test_score_i)
        
        
##########


def split_4fold_CV(dico_full_RCDB) :
    """A partir du dictionnaire regroupant l'ensemble des données, cette fonction créée 4 dictionnaires/ datasets,\
    où chaque récepteur a au moins 2 composés associés, de manière à pouvoir effectuer une Cross-Validation correcte de type 4fold.\
    Cette fonction renvoi une liste de dictionnaires."""
    import random
    all_sets = [{},{},{},{}]
    set1 = {}
    set2 = {}
    set3 = {}
    set4 = {}
    L_recepteurs_totaux = []
    L_recepteurs_done = []

    for comp in dico_full_RCDB : 
        for rec in dico_full_RCDB[comp]:
            L_recepteurs_totaux.append(rec)

    for rec in L_recepteurs_totaux :

        if rec not in L_recepteurs_done :
            L_rand_i = []
        # 8 ou Plus  -> 4sets
            if L_recepteurs_totaux.count(rec) >=8 : 
                c1=0
                #print('\n',L_recepteurs_totaux.count(rec), rec, '    CAS >=8 \n')
                for comp in dico_full_RCDB:
                    if c1 == 4:
                        c1=0                
                    if rec in dico_full_RCDB[comp]:
                        #print(rec, comp)
                        if comp not in all_sets[c1] :
                            all_sets[c1][comp] = [rec]
                            c1+=1
                        else :
                            all_sets[c1][comp].append(rec)
                            c1+=1                    

            # 6 et 7 -> 3 sets 
            elif L_recepteurs_totaux.count(rec) in range(6,8) :


                while len(L_rand_i)!=3 :
                    rand_i = random.randint(0,3)
                    if rand_i not in L_rand_i :
                        L_rand_i.append(rand_i)
                c2 = 0        
                #print('\n',L_recepteurs_totaux.count(rec), rec, '    CAS 6,7 \n')
                for comp in dico_full_RCDB:
                    if c2 == 3:
                        c2=0 
                    if rec in dico_full_RCDB[comp]:
                        #print(rec, comp)
                        #print(L_rand_i)
                        if comp not in all_sets[L_rand_i[c2]] :
                            all_sets[L_rand_i[c2]][comp] = [rec]
                            c2+=1
                        else :
                            all_sets[L_rand_i[c2]][comp].append(rec)
                            c2+=1


            # 4 à 5 -> 2 sets
            elif L_recepteurs_totaux.count(rec) in range(4,6) : 


                while len(L_rand_i)!=2 :
                    rand_i = random.randint(0,3)
                    if rand_i not in L_rand_i :
                        L_rand_i.append(rand_i)
                c3 = 0 

                #print('\n',L_recepteurs_totaux.count(rec), rec, '    CAS 4,5 \n')
                for comp in dico_full_RCDB:
                    if c3 == 2:
                        c3=0 
                    if rec in dico_full_RCDB[comp]:
                        #print(rec, comp)
                        if comp not in all_sets[L_rand_i[c3]] :
                            all_sets[L_rand_i[c3]][comp] = [rec]
                            c3+=1
                        else :
                            all_sets[L_rand_i[c3]][comp].append(rec)
                            c3+=1

            # 3 -> set            
            elif L_recepteurs_totaux.count(rec) <4 : 
                L_rand_i = [random.randint(0,3)]

                #print('\n',L_recepteurs_totaux.count(rec), rec,'    CAS <=3 \n')
                for comp in dico_full_RCDB:

                    if rec in dico_full_RCDB[comp]:
                        #print(rec, comp)
                        if comp not in all_sets[L_rand_i[0]] :
                            all_sets[L_rand_i[0]][comp] = [rec]
                        else :
                            all_sets[L_rand_i[0]][comp].append(rec)

            L_recepteurs_done.append(rec)
            #print(L_rand_i)
    return(all_sets)

##########
    
