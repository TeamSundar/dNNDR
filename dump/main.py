import pandas as pd
# Python modules used for API access...
from chembl_webresource_client.new_client import new_client

#df = pd.read_excel (r'data/ci400709d_si_002.xlsx', sheet_name='KIBA')
#print(df.head)

# chembl_id = "CHEMBL5554"
# activities = new_client.mechanism.filter(target_chembl_id=chembl_id)
# compound_ids = [x['molecule_chembl_id'] for x in activities]
# approved_drugs = new_client.molecule.filter(molecule_chembl_id__in=compound_ids).filter(max_phase=4)

# for record in approved_drugs:
#     print("{:10s} : {}".format(record["molecule_chembl_id"], record["pref_name"]))

#Get compound record using client...
chembl_id = 'CHEMBL10'
records = new_client.activity.filter(molecule_chembl_id=chembl_id)
target_df = pd.DataFrame(columns=['target_organism', 'target_chembl_id', 'type', 'units', 'value'])
# for i in range(len(records)):
#     target = {'target_organism': records[i]['target_organism'], 'target_chembl_id': records[i]['target_chembl_id'], 
#             'type': records[i]['type'], 'units':records[i]['units'],
#             'value':records[i]['value']}

#     target_df = target_df.append(target, ignore_index=True)

# print(target_df.head())
print(len(records))


