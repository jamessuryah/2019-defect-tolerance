import json
from matplotlib import pyplot as plt
import numpy as np


with open('binary_ranges.json', 'r') as oq:
    q = json.load(oq)


q = [[x[0], abs(x[1][1]-x[1][0])] for x in q]
q.sort(key=lambda x: x[1], reverse=True)
for item in q:
    print(item)
items = []
newitems = []
for item in q:
    if item[0] in items:
        pass
        #print(item[0])
        #print(item)
    else:
        items.append(item[0])
        newitems.append(item)


elements = {}
for items in newitems:
    if items[1][-1] not in list(elements.keys()):
        elements[items[1][-1]] = 1
    else:
        elements[items[1][-1]] += 1


#with open('output_q_fixed.json', 'w') as oqf:
#    json.dump(newitems, oqf)

q.sort(key=lambda x: x[2], reverse=True)
n = 0
for item in q:
    if n < 200 and '32' not in item[0]:
        print([item[0], item[-1]])