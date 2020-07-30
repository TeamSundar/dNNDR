#!E:\Projects\202004_MultiLabelClass python
# -*- coding: utf-8 -*-

# Import essential libraries

import numpy as np
import pandas as pd
import pickle
import urllib.parse
import urllib.request
import wget
import os
from rdkit import Chem
from rdkit.ML.Descriptors.MoleculeDescriptors import MolecularDescriptorCalculator as md
from rdkit import DataStructs
from chembl_webresource_client.new_client import new_client
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
    def extractProteinDes(self, des_type, pathToData):
        fset = pd.DataFrame()
        chembl2uniprot = pd.read_csv(pathToData, sep='\t', header=None)
        des_type = des_type
        for CHEMBL_ID, uniprot_ID in zip(chembl2uniprot[0], chembl2uniprot[1]):
            print(CHEMBL_ID, uniprot_ID)
            push ='python F:/Projects/202004_DrugTargetNN/iFeature/iFeature.py --file %s --type %s'%('data/fasta/'+uniprot_ID+'.fasta', des_type)
            os.system(push)
            enc = pd.read_csv('encoding.tsv', delimiter='\t', encoding='utf-8')
            enc['uniprot_ID'] = uniprot_ID
            enc['CHEMBL_ID'] = CHEMBL_ID
            fset = fset.append(enc)
        fset.to_csv('data/target_'+ des_type + '.csv')
        return fset


