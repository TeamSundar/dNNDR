import urllib.parse
import urllib.request

url = 'https://www.uniprot.org/uploadlists/'

params = {
'from': 'ACC+ID',
'to': 'CHEMBL_ID',
'format': 'tab',
'query': 'O00141 O00142 O00238 O00311 O00329 O00418 O00443 O00444 O00506 O00750'
}

data = urllib.parse.urlencode(params)
data = data.encode('utf-8')
req = urllib.request.Request(url, data)
with urllib.request.urlopen(req) as f:
   response = f.read()
a = response.decode('utf-8')
file1 = open("myfile.txt","w") 
file1.write(a)
file1.close()
print(a)