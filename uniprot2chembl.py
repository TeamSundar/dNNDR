import urllib.parse
import urllib.request
import pandas as pd
import wget
import os
from tqdm import tqdm

url = 'https://www.uniprot.org/uploadlists/'
targets = pd.read_csv('data/DTI.csv')
fset = pd.DataFrame()

for CHEMBL_ID in tqdm(targets['target'].unique()):
    params = {
    'from': 'CHEMBL_ID',
    'to': 'ACC',
    'format': 'tab',
    'query': CHEMBL_ID
    }

    try:
        data = urllib.parse.urlencode(params)
        data = data.encode('utf-8')
        req = urllib.request.Request(url, data)
        with urllib.request.urlopen(req) as f:
            response = f.read()
        a = response.decode('utf-8')
        #print(a)
        #print(a.split('\n')[1].split('\t'))

        #desc = ['AAC', 'CKSAAP', 'CTDC', 'CTDD', 'CTDT', 'DDE', 'DPC', 'GAAC', 'GDPC', 'GTPC', 'TPC', 'QSOrder', 'KSCTriad', 'CTriad']
        des_type = 'AAC'
        #for uniprotID in tqdm(protein_list):
        wget.download('https://www.uniprot.org/uniprot/%s.fasta'%(a.split('\n')[1].split('\t')[1]), out='fasta')
        push ='python F:/Projects/202004_DrugTargetNN/iFeature/iFeature.py --file %s --type %s'%('fasta/'+(a.split('\n')[1].split('\t')[1])+'.fasta', des_type)
        os.system(push)
        enc = pd.read_csv('encoding.tsv', delimiter='\t', encoding='utf-8')
        enc['#'] = a.split('\n')[1].split('\t')[1]
        fset = fset.append(enc)
    except:
        pass
        #print('Uniprot ID for %s not found.'%CHEMBL_ID)

fset.to_csv('target_'+ des_type + '.csv')

# with open("sample.txt", "a") as a_file:
#     file1 = open("chembl2uniprot.txt","w") 
#     file1.write('\n')
#     file1.write(a)
#     file1.close()

# df = pd.read_csv('myfile.txt', sep='\t')
# print(df['To'])