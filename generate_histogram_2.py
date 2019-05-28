import json
from matplotlib import pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import scipy.stats as sstats
from scipy.stats import norm, sem
maxwell = sstats.maxwell

with open('output_q.json') as output:
    data_q = json.load(output)


with open('output.json') as output:
    data_t = json.load(output)

sorted_t = []
sorted_q  = []
for datum in data_t:
    anion = [x for x in datum[1] if x in ['O', 'S', 'Se', 'Te']][0]
    if anion == 'O':
        sorted_t.append(datum)
    if datum in list(t_data.keys()):


for datum in data_q:
    anion = [x for x in datum[1] if x in ['O', 'S', 'Se', 'Te']][0]
    if anion == 'O':
        sorted_q.append(datum)
sorted_t_calc = []
sorted_q_calc = []
for item in sorted_t:
    itvol = item[-1]
    itadj = (itvol/np.pi) ** (1. / 2)
    nitem = item[:-1]+[itadj]
    sorted_t_calc.append(nitem)
for item in sorted_q:
    itvol = item[-1]
    itadj = (3*itvol/(4*np.pi)) ** (1. / 3)
    nitem = item[:-1]+[itadj]
    sorted_q_calc.append(nitem)
print(sorted_t_calc[201])
print(sorted_t[201])
print(sorted_q_calc[201])
print(sorted_q[201])
sorteds = {'Ternary': sorted_t_calc, 'Quartenary': sorted_q_calc}
bins = np.linspace(-3, 1, 80)
for names in list(sorteds.keys()):
    sorts = sorteds[names]
    vol = [x[-1] for x in sorts]
    print(len(vol))
    logvol = [np.log10(v) for v in vol]
    mu, std = norm.fit(logvol)
    mfit = maxwell.fit(vol, floc=0)
    data_sem = sem(logvol)
    norm_p = norm.pdf(bins, mu, std)
    dnorm_p = [x*len(vol)/sum(norm_p) for x in norm_p]
    maxw_p = maxwell.pdf(np.linspace(0,5,100), *mfit)
    dmaxw_p = [x*len(vol)/sum(maxw_p) for x in maxw_p]
    plt.hist(logvol,bins, alpha=0.5, label=names)
    plt.plot(bins, dnorm_p, 'k')
    #plt.axvline((np.mean(vol)), color='k', linestyle='dashed', linewidth=1)
    plt.axvline(np.mean(logvol), color='k', linestyle='-', linewidth=1)
    print(names)
    print(np.mean(vol))
    print(mu)
    print(data_sem)
plt.legend(loc='upper right')
plt.xlim(-3, 1)
plt.ylabel('Number of phases')
plt.xlabel('Log stability vector, eV')
plt.show()

# Quantify the axes q
