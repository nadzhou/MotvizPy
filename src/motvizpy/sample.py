import urllib.parse
import urllib.request

example_list = ["XP_004053513", "XP_008146952", "XP_027445634"]
example_str = " ".join(example_list)


url = 'https://www.uniprot.org/uploadlists/'

params = {
'from': 'P_REFSEQ_AC',
'to': 'ID',
'format': 'tab',
'query': example_str
}

data = urllib.parse.urlencode(params)
data = data.encode('utf-8')
req = urllib.request.Request(url, data)
with urllib.request.urlopen(req) as f:
   response = f.read()
print(response.decode('utf-8'))