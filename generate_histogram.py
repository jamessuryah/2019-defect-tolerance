import json
from matplotlib import pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


def func(x, a, b):
    return a * np.exp(-b*x)


#Ternary stuff
with open('output_file_tern_c.json', 'r') as oftc:
    t_data = json.load(oftc)


with open('output.json') as output:
    data_t = json.load(output)


with open('output_q_fixed.json') as output:
    data_q = json.load(output)


with open('output_b.json') as biput:
    data_b = json.load(biput)


sorted_b = {}
for datum in data_b:
    anion = [x for x in datum[1] if x in ['O', 'S', 'Se', 'Te']][0]
    if anion not in list(sorted_b.keys()):
        sorted_b[anion] = [datum]
    else:
        sorted_b[anion].append(datum)
sorted_t = {}
sorted_t_hr = {}
for datum in data_t:
    anion = [x for x in datum[1] if x in ['O', 'S', 'Se', 'Te']][0]
    if anion not in list(sorted_t.keys()):
        sorted_t[anion] = [datum]
    else:
        sorted_t[anion].append(datum)
    if datum[0] in list(t_data.keys()):
        tdic = t_data[datum[0]]
        try:
            hrdx = tdic["hi_res"]["x"]
            hrdy = tdic["hi_res"]["y"]
            ly = [abs(hrdx[x][1][1] - hrdx[x][0][1]) for x in list(hrdx.keys())]
            lx = [abs(hrdy[x][1][0] - hrdy[x][0][0]) for x in list(hrdy.keys())]
            l = min([max(lx), max(ly)])
            if anion not in list(sorted_t_hr.keys()):
                sorted_t_hr[anion] = [l]
            else:
                sorted_t_hr[anion].append(l)
        except IndexError:
            pass

sorted_q = {}
for datum in data_q:
    anion = [x for x in datum[1] if x in ['O', 'S', 'Se', 'Te']][0]
    if anion not in list(sorted_q.keys()):
        sorted_q[anion] = [datum]
    else:
        sorted_q[anion].append(datum)


for an in list(sorted_q.keys()):
    print(an)
    print(np.mean([x[-1] for x in sorted_q[an]]))
    print(np.power(10, np.mean([np.log10(x[-1]) for x in sorted_q[an]])))

ob = sorted_b['O']
ot = sorted_t['O']
oq = sorted_q['O']
sizesb = [x[-1]/2 for x in ob]
sizest = [np.sqrt(x[-1]/np.pi) for x in ot]
sizesq = [np.cbrt(3*x[-1]/(4*np.pi)) for x in oq]
sizesall = sizesb+sizest+sizesq
linrange = np.linspace(min(sizesall), 4, 40)
logrange = np.linspace(np.log10(min(sizesall)), np.log10(max(sizesall)), 40)
logb= [np.log10(x) for x in sizesb]
logt= [np.log10(x) for x in sizest]
logq = [np.log10(x) for x in sizesq]

print(np.mean(sizesb))
print(np.power(10, np.mean(logb)))
print(np.mean(sizest))
print(np.power(10, np.mean(logt)))
print(np.mean(sizesq))
print(np.power(10, np.mean(logq)))

# hist_b = plt.hist(logb, logrange, alpha=0.3, label="Binary")
# hist_t = plt.hist(logt, logrange, alpha=0.3, label="Ternary")
hist_q = plt.hist(sizesq, linrange, alpha=1)
yhq = list(hist_q[0])
xhq_r = list(hist_q[1])
xhq = [(xhq_r[a]+xhq_r[a+1])/2 for a in range(0, len(xhq_r)-1)]
max_index = yhq.index(max(yhq))
popt, pcov = curve_fit(func, xhq[max_index:], yhq[max_index:])
fit = []
for i in xhq:
    if xhq.index(i) < max_index:
        print(func(i, *popt)-yhq[xhq.index(i)])
    fit.append(func(i, *popt))
plt.plot(xhq, fit, 'k-')
#
plt.legend(loc="upper right")
plt.xlabel('Stability metric, eV')
plt.ylabel('Number of phases')
# plt.show()
#bv = []
#for anion in list(sorted_b.keys()):
#    bv += [x[-1] for x in sorted_b[anion]]
#for anion in list(sorted_b.keys()):
#    vol = [x[-1] for x in sorted_b[anion]]
#    linbin = np.linspace(min(bv), max(bv), 40)
#    logvol = [np.log10(v) for v in vol]
#    lbins = np.linspace(-7, 3, 40)
#    y, binEdges = np.histogram(vol,bins=40)
#    bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
#    plt.plot(bincenters,y,'-',label=anion)
#    print(np.mean(vol))
#    print(np.power(10, np.mean([np.log10(x) for x in vol])))
#    print("{0} arithmetic average = {1}".format(anion, np.mean(vol)))
#    print("{0} geometric average = {1}".format(anion, np.power(10, np.mean(logvol))))
    #plt.axvline(np.log10(np.mean(vol)), color='k', linestyle='dashed', linewidth=1)
    #plt.axvline(np.mean(logvol), color='k', linestyle='-', linewidth=1)
#plt.legend(loc='upper right')
#plt.rc('text', usetex=True)
#plt.title('Distribution for binary stability range')
#plt.ylabel('Number of phases')
#plt.xlabel('Stability range, eV')
print('----')
print(popt)
print(pcov)
perr = np.sqrt(np.diag(pcov))
print(perr)
plt.show()

