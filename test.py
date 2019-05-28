import json

with open('cache_compounds_new.json') as cache:
    ter = json.load(cache)


to_check = []
for system in ter:
    for compound in list(system.keys()):
        comp = system[compound][0]
        if len(comp) == 3:
            to_check.append([compound, system])


print(len(to_check))
#with open('to_check.json','w') as toc:
#    json.dump(to_check, toc, indent=4, separators=(',',': '))
