#!F:\Projects\20202507_dNNDR\src\ python
# -*- coding: utf-8 -*-

# Import essential libraries

import numpy as np
import pandas as pd
import pickle
import urllib.parse
import urllib.request
#import wget
import os
import seaborn as sns
import matplotlib.pyplot as plt
from collections import Counter
from sklearn.model_selection import train_test_split
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score
from sklearn.metrics import roc_curve, auc
# from rdkit import Chem
# from rdkit.ML.Descriptors.MoleculeDescriptors import MolecularDescriptorCalculator as md
# from rdkit import DataStructs
#from chembl_webresource_client.new_client import new_client
from tqdm import tqdm

class PROCESS:
    def __init__(self, data):
        self.data = data

    # Extaract terget data for 52k drugs from CHEMBL
    def extractDrugTargets(self, drugList):
        main_dict = dict.fromkeys(drugList, [])
        for drug in main_dict:
            main_dict[drug] = dict.fromkeys(['target_organism', 'target_chembl_id', 'type', 'units', 'value'],[])

        for drug in tqdm(main_dict): 
            records = new_client.activity.filter(molecule_chembl_id=drug).only(['target_organism', 'target_chembl_id', 'type', 'units', 'value'])

            target_chembl_id, target_organism, type, units, value = [], [], [], [], []
            for a in records:
                target_chembl_id.append(a['target_chembl_id'])
                target_organism.append(a['target_organism'])
                type.append(a['type'])
                units.append(a['units'])
                value.append(a['value'])
                
            main_dict[drug]['target_organism'] = target_chembl_id
            main_dict[drug]['target_chembl_id'] = target_chembl_id
            main_dict[drug]['type'] = type
            main_dict[drug]['units'] = units
            main_dict[drug]['value'] = value
        return main_dict
    
    # Extract Drug-Target pairs for Homo-Sapiens with IC50 values
    def getDTI(self, drug_target):
        DTI = pd.DataFrame(columns=['target_organism', 'drug', 'target', 'IC50', 'unit'])
        for drug in tqdm(drug_target):
            for i in range(len(drug_target[drug]['target_chembl_id'])):
                if drug_target[drug]['type'][i]=='IC50' and drug_target[drug]['target_organism'][i]=='Homo sapiens':
                    dict = {'target_organism':drug_target[drug]['target_organism'][i], 
                            'drug':drug,
                            'target':drug_target[drug]['target_chembl_id'][i],
                            'IC50':drug_target[drug]['value'][i],
                            'unit':drug_target[drug]['units'][i]}
                    DTI = DTI.append(dict, True)
        # Save to file
        DTI.to_csv('data/DTI2.csv')
        return DTI
    
    # Screen Data based on sequence availiability
    def screenDTI(self, DTI, pathToMapping, pathToOutput):
        chembl2uniprot = pd.read_csv(pathToMapping, sep='\t', header=None)  #Import CHEMBL_ID to uniprot_ID mapping
        units=['nM','uM','pM','mM']     # Units to be selected
        DTI_screened_units = DTI.loc[DTI['unit'].isin(units)]   # Extract datapoints with required units
        DTI_screened_units = DTI_screened_units[DTI_screened_units['IC50'].notnull()]
        DTI_screened_mapping = DTI_screened_units[DTI_screened_units['target'].isin(chembl2uniprot[0].tolist())]    # Extract datapoints for which uniprot_ID ia available
        DTI_screened_mapping.to_csv(pathToOutput)

        # Summary statistics
        print('Total targets and datapoints acquired   : %s | %s'%(len(DTI['target'].unique()), len(DTI)))
        print('Targets and datapoints after IC50 screen: %s | %s'%(len(DTI_screened_units['target'].unique()), len(DTI_screened_units)))
        print('Removing NaNs and unannotated CHEMBL_ID : %s | %s'%(len(DTI_screened_mapping['target'].unique()), len(DTI_screened_mapping)))
        return DTI_screened_mapping

    # Save dict 
    def save_obj(self, obj, name):
        with open(name, 'wb') as f:
            pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

    # Load dict
    def load_obj(self, name):
        with open(name, 'rb') as f:
            return pickle.load(f)

    # Extract drug descriptors from smile string
    def mol2des(self, chembl_id, attribs):
        try: 
            m = Chem.MolFromSmiles(new_client.activity.filter(molecule_chembl_id=chembl_id).only(['canonical_smiles'])[0]['canonical_smiles'])
            calc = md(attribs)
            des = pd.DataFrame(calc.CalcDescriptors(m)).T
            des.columns = attribs
            des['Drug'] = str(chembl_id)
            status = True
        except:
            des = []
            status = False
        return des, status
    
    # Map CHEMBL traget IDs to uniprot ID
    def chembl2uniprot(self, pathToData, pathToOutput):
        url = 'https://www.uniprot.org/uploadlists/'
        targets = pd.read_csv(pathToData)
        
        for i in tqdm(range(len(targets['target'].unique()))):
            params = {
            'from': 'CHEMBL_ID',
            'to': 'ACC',
            'format': 'tab',
            'query': targets['target'].unique()[i]
            }

            data = urllib.parse.urlencode(params)
            data = data.encode('utf-8')
            req = urllib.request.Request(url, data)
            with urllib.request.urlopen(req) as f:
                response = f.read()
            a = response.decode('utf-8')

            with open(pathToOutput, "a") as file: 
                file.write('\n')
                file.write(a.split('\n')[1])
        return None
    
    # Extract protein descriptors for retrieved targets
    def extractProteinDes(self, des_type, pathToData, pathToOutput):
        fset = pd.DataFrame()
        chembl2uniprot = pd.read_csv(pathToData, sep='\t', header=None)
        des_type = des_type
        for CHEMBL_ID, uniprot_ID in zip(chembl2uniprot[0], chembl2uniprot[1]):
            #print(CHEMBL_ID, uniprot_ID)
            try:
                push ='python /home/dell15/KING/Work_ubuntu/Projects/20200703_drugTarget2/src2/dNNDR/iFeature-master/iFeature.py --file %s --type %s'%('/home/dell15/KING/Work_ubuntu/Projects/20200703_drugTarget2/src2/dNNDR/data/fasta_968/'+uniprot_ID+'.fasta', des_type)
            except:
                print('Fasta not found')
            os.system(push)
            enc = pd.read_csv('encoding.tsv', delimiter='\t', encoding='utf-8')
            enc['uniprot_ID'] = uniprot_ID
            enc['CHEMBL_ID'] = CHEMBL_ID
            fset = fset.append(enc)
        fset.to_csv(pathToOutput+'/'+ des_type + '.csv')
        return fset

class ANALYSIS:
    def __init__(self, EXP):
        self.EXP = EXP

    def plot_seq_count(self, df, data_name):
        sns.distplot(df['seq_char_count'].values)
        plt.title(f'Sequence char count: {data_name}')
        plt.grid(True)

    def get_code_freq(self, df, data_name):
        df = df.apply(lambda x: " ".join(x))
        codes = []
        for i in df: # concatination of all codes
            codes.extend(i)
        codes_dict= Counter(codes)
        codes_dict.pop(' ') # removing white space

        print(f'Codes: {data_name}')
        print(f'Total unique codes: {len(codes_dict.keys())}')

        df = pd.DataFrame({'Code': list(codes_dict.keys()), 'Freq': list(codes_dict.values())})
        return df.sort_values('Freq', ascending=False).reset_index()[['Code', 'Freq']]

    def plot_code_freq(self, df, data_name):
        plt.title(f'Code frequency: {data_name}')
        sns.barplot(x='Code', y='Freq', data=df)
        plt.grid(True)

    def create_dict(self, codes):
        char_dict = {}
        for index, val in enumerate(codes):
            char_dict[val] = index+1
        return char_dict

    def integer_encoding(self, data, dict):
        """
        - Encodes code sequence to integer values.
        - 20 common amino acids are taken into consideration
        and rest 4 are categorized as 0.
        """
        encode_list = []
        for row in data['seq'].values:
            row_encode = []
            for code in row:
                row_encode.append(dict.get(code, 0))
            encode_list.append(np.array(row_encode))
        return encode_list


    def roc(self, model, label, **kwargs):
        data = []
        for key in kwargs:
            data.append(kwargs[key])
        y_pred_keras = model.predict(data)

        # Compute ROC/AUC
        fpr = dict()
        tpr = dict()
        roc_auc = dict()

        for i in range(3):
            fpr[i], tpr[i], _ = roc_curve(label[:, i], y_pred_keras[:, i])
            roc_auc[i] = auc(fpr[i], tpr[i])

        return roc_auc, fpr, tpr

    def aupr(self, model, label, **kwargs):
        data = []
        for key in kwargs:
            data.append(kwargs[key])
        y_pred_keras = model.predict(data)

        # compute PR/AUPR
        precision, recall, average_precision = dict(), dict(), dict()

        for i in range(3):
            precision[i], recall[i], _ = precision_recall_curve(label[:, i], y_pred_keras[:, i])
            average_precision[i] = average_precision_score(label[:, i], y_pred_keras[:, i])

        # A "micro-average": quantifying score on all classes jointly
        precision["micro"], recall["micro"], _ = precision_recall_curve(label.ravel(), y_pred_keras.ravel())
        average_precision["micro"] = average_precision_score(label, y_pred_keras, average="micro")
        #print('Average precision score , micro-averaged over all classes: {0:0.2f}'
        #    .format(average_precision["micro"]))
        return  precision, recall, average_precision

    def plotROC_PR(self, fpr, tpr, roc_auc, precision, recall, average_precision, save):
        fig, ax = plt.subplots(1,2,figsize=(12,5))

        for i, activity in zip(range(3),['Active', 'Intermediate','Inactive']):
            ax[0].plot(fpr[i], tpr[i], label=activity+' (AUC: %0.2f)' % roc_auc[i], alpha=1)
        ax[0].plot([0, 1], [0, 1], 'k--')
        ax[0].set_ylim([0.0, 1.01])
        ax[0].set_xlim([0.0, 1.01])
        ax[0].set_yticks(np.arange(0, 1.1, 0.1))
        ax[0].set_xticks(np.arange(0, 1.1, 0.1))
        ax[0].set_title('ROC curves for all three classes', fontsize=14)
        ax[0].set_xlabel('False Positive Rate', fontsize=14)
        ax[0].set_ylabel('True Positive Rate', fontsize=14)
        #ax.set_title('Receiver operating characteristic ('+train_test+')', fontsize=14)
        ax[0].tick_params(axis="x", labelsize=12)
        ax[0].tick_params(axis="y", labelsize=12) 
        ax[0].grid(linestyle='-.', linewidth=0.7)
        ax[0].legend(fontsize=12)

        ax[1].step(recall['micro'], precision['micro'], where='post')
        ax[1].set_xlabel('Recall', fontsize=14)
        ax[1].set_ylabel('Precision', fontsize=14)
        ax[1].plot([0, 1], [1, 0], 'k--')
        ax[1].set_ylim([0.0, 1.01])
        ax[1].set_xlim([0.0, 1.00])
        ax[1].set_yticks(np.arange(0, 1.1, 0.1))
        ax[1].set_xticks(np.arange(0, 1.1, 0.1))
        ax[1].set_title(
            'Micro-averaged AUPR for classes: AP={0:0.2f}'
            .format(average_precision["micro"]), fontsize=14)
        ax[1].tick_params(axis="x", labelsize=12)
        ax[1].tick_params(axis="y", labelsize=12) 
        ax[1].grid(linestyle='-.', linewidth=0.7)

        if save:
            plt.savefig('plots/'+self.EXP+'_ROC_PR.png', dpi=500, format = 'png', bbox_inches='tight')
        return None

    def plotTrainingPerf(self, history, save):
        #  Training performance
        fig, ax = plt.subplots(1,2,figsize=(12,5))

        acc = history.history['accuracy']
        val_acc = history.history['val_accuracy']
        loss = history.history['loss']
        val_loss = history.history['val_loss']
        x = range(1, len(acc) + 1)

        ax[0].plot(x, acc, 'b', label='Training acc')
        ax[0].plot(x, val_acc, 'r', label='Validation acc')
        ax[0].set_title('Training and validation accuracy')
        ax[0].set_xlabel('Epoch', fontsize=14)
        ax[0].set_ylabel('Accuracy', fontsize=14)
        ax[0].tick_params(axis="x", labelsize=12)
        ax[0].tick_params(axis="y", labelsize=12)
        ax[0].set_ylim([0,1])
        ax[0].set_xticks(np.arange(0, 350, 50))
        #ax[0].set_yticks(np.arange(0, 1, 0.1))
        ax[0].legend(fontsize=12)
        ax[0].grid(linestyle='-.', linewidth=0.7)

        ax[1].plot(x, loss, 'b', label='Training loss')
        ax[1].plot(x, val_loss, 'r', label='Validation loss')
        ax[1].set_title('Training and validation loss')
        ax[1].set_xlabel('Epoch', fontsize=14)
        ax[1].set_ylabel('Loss (categorical crossentropy)', fontsize=14)
        ax[1].tick_params(axis="x", labelsize=12)
        ax[1].tick_params(axis="y", labelsize=12)
        ax[1].set_ylim([0,2])
        ax[1].set_xticks(np.arange(0, 350, 50))
        ax[1].legend(fontsize=12)
        ax[1].grid(linestyle='-.', linewidth=0.7)
        if save:
            plt.savefig('plots/'+self.EXP+'_training_perf.png', dpi=500, format = 'png', bbox_inches='tight')
        return None
