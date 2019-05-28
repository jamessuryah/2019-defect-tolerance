import json
import numpy as np
import matplotlib.pyplot as plt

with open('output_file_mesh.json') as ofm:
    meshes = json.load(ofm)


dic_of_thicc = {}
f = 0
c = 0
for item in list(meshes.keys()):
    ap = []
    mesh = meshes[item]
    meshx = mesh['x']
    meshy = mesh['y']
    anion = mesh['anion']
    xlength = []
    ylength = []
    for xpoint in list(meshx.keys()):
        line = meshx[xpoint]
        if len(line)<2:
            f+= 1
        else:
            line_length = abs(line[0][1]-line[1][1])
            xlength.append(line_length)
    for ypoint in list(meshy.keys()):
        line = meshy[ypoint]
        if len(line)<2:
            f +=1
        else:
            line_length=abs(line[0][0]-line[1][0])
            ylength.append(line_length)
    if len(xlength)>3 and len(ylength)>3:
        c += 1
        if np.mean(xlength)<np.mean(ylength):
            dic_of_thicc[item] = [max(xlength), anion]
        else:
            dic_of_thicc[item] = [max(ylength), anion]

overall = {}
for i in dic_of_thicc:
    item = dic_of_thicc[i]
    if item[1] not in list(overall.keys()):
        print(item[1])
        overall[item[1]] = [item[0]]
    else:
        overall[item[1]].append(item[0])

total = overall['O']+overall['S']+overall['Se']+overall['Te']
xlog = [np.log10(x) for x in total]
bins = np.linspace(min(xlog), max(xlog), 40)
for an in ['O', 'S', 'Se', 'Te']:
    stuff = overall[an]
    plt.hist([np.log10(x) for x in stuff], bins, alpha=0.5, label=an)
plt.xlabel('Log of limiting potential for vacancy, eV')
plt.ylabel('Number of ions')
plt.legend(loc='upper right')
plt.show()

