import bacdive


def start_client(info_file: str) -> bacdive.BacdiveClient:
    with open(info_file, "r") as f:
        tokens = f.readlines()
    email = tokens[0].strip()
    pw = tokens[1].strip()
    client = bacdive.BacdiveClient(email, pw)
    return client


client = start_client(".bacdive_info")

# help(client.search)
# help(client.retrieve)

# the search method fetches all BacDive-IDs matching your query
# and returns the number of IDs found
print("Taxonomy search 1")
query = "Acetanaerobacterium_elongatum_gca_900103835.IMG-taxon_2667527408_annotated_assembly.49"
count = client.search(taxonomy='"SELECT *')

query = {"genome": "GCA_900103835"} # genome sequence
query = {"genome": "GCA_002005445"} # genome sequence

# run query
client.search(**query)

print("Taxonomy search 2")
count = client.search(taxonomy='Acetobacter aceti')


print("Taxonomy search 3")
client.search(taxonomy='Acetanaerobacterium elongatum')

print(count, 'strains found.')
filter=["phylum","range","metabolite",'utilization activity','oxygen tolerance',"enzymes"]

bac_list_attr = list()
filter=["phylum"]
result = client.retrieve(filter)
bac_dict = {k:v for x in result for k,v in x.items()}
for i in bac_dict.values():
    for j in i:
        print("Phylum:",j["phylum"])
        bac_list_attr.append(j["phylum"])

filter=["metabolite utilization"]
result = client.retrieve(filter)
bac_dict = {k:v for x in result for k,v in x.items()}
for i in bac_dict.values():
    for j in i[0].values():
        for z in j:
            print("Metabolite_utilisation:", z["metabolite"]+"_"+z["utilization activity"])
            bac_list_attr.append(z["metabolite"]+"_"+z["utilization activity"])

# the retrieve method lets you iterate over all strains
# and returns the full entry as dict
# Entries can be further filtered using a list of keys (e.g. ['keywords'])
strains = client.retrieve()
print(next(strains).keys())
for strain in strains:
    print(strain, "\n")

