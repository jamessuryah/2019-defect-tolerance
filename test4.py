import json
from matplotlib import pyplot as plt
import numpy as np


with open('output_file_tern_c.json', 'r') as ofm:
    m = json.load(ofm)


hress=[]
lress=[]
areas = []
for mat in list(m.keys()):
    test = True
    mdic = m[mat]
    area = mdic['area']
    rad = np.sqrt(area/(4*np.pi))
    hr = mdic['hi_res']
    lr = mdic['low_res']
    hrx = hr['x']
    hry = hr['y']
    lrx = lr['x']
    lry = lr['y']
    xlines = []
    ylines = []
    lranges = []
    for xmes in list(hrx.keys()):
        range = hrx[xmes]
        if len(range) > 1:
            xlines.append(abs(range[1][1]-range[0][1]))
    for ymes in list(hry.keys()):
        range = hry[ymes]
        if len(range) > 1:
            ylines.append(abs(range[1][0]-range[0][0]))
    try:
        if sum(xlines)/len(xlines) > sum(ylines)/len(ylines):
            hres = max(ylines)
        else:
            hres = max(xlines)
    except ZeroDivisionError:
        test = False
    try:
        lranges.append(abs(lrx[list(lrx.keys())[0]][1][1]-lrx[list(lrx.keys())[0]][0][1]))
    except IndexError:
        print(lrx)
    try:
        lranges.append(abs(lry[list(lry.keys())[0]][1][0]-lry[list(lry.keys())[0]][0][0]))
    except IndexError:
        print(lry)
    if len(lranges) > 0 and test:
        lres = min(lranges)
        hress.append(hres)
        lress.append(lres)
        areas.append(area)

plt.plot(hress, lress,'k.', hress, areas, 'b.')
plt.xlabel('From many meshes, eV')
plt.ylabel('From centre point, eV')
plt.show()