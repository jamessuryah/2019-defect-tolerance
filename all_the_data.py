import json
import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import expon
from scipy.optimize import curve_fit
# Quartenary volumes
with open('output_q.json', 'r') as oq:
    q_volumes = json.load(oq)


#Ternary stuff
with open('output_file_tern_c.json', 'r') as oftc:
    t_data = json.load(oftc)


#Binary stuff
with open('binary_ranges.json','r') as br:
    b_ranges = json.load(br)

#Others
with open('int_points_ttt_fixed.json', 'r') as ip:
    m = json.load(ip)


def linfunc(x, a):
    return a*x

sum_dic = {}
binrange = []
sum_dic['b'] = {}
for bitem in b_ranges:
    name = bitem[0]
    dist = abs(bitem[1][1]-bitem[1][0])
    sum_dic['b'][name] = dist
    binrange.append(dist)

sum_dic['t'] = {}


sum_dic['q'] = {}
for material in q_volumes:
    vol = material[-1]
    anion = material[1][-1]
    name = material[0]
    rad = np.cbrt(3*vol/(4*np.pi))
    sum_dic['q'][name] = {'anion': anion, "volume": rad}


c = 0
for n in list(m.keys()):
    if n not in ["failure", "failures"]:
        mat = m[n]
        axes = ['x', 'y', 'z']
        diss = []
        anion = mat['anion']
        for ax in axes:
            idx = axes.index(ax)
            try:
                mat_ax = mat[ax]["0"]
                #print(mat_ax)
                #print(idx)
                ord = [x[idx] for x in mat_ax]
                if len(ord) > 1:
                    dis = max(ord)-min(ord)
                    diss.append(dis)
            except TypeError:
                print(n)
        if len(diss) != 0:
            res = min(diss)
            if res > 0.00000000001:
                c += 1
                if n in list(sum_dic['q'].keys()):
                    sum_dic['q'][n]['center'] = [res, anion]
                else:
                    sum_dic['q'][n] = {'center': [res, anion]}


radiis = []
aradiis = []
t_hires = []
t_lowres = []
t_hires_a = {}
for t_name in list(t_data.keys()):
    tdic = t_data[t_name]
    if "area" in list(tdic.keys()):
        aradiis.append(np.sqrt(tdic["area"]/np.pi))
    if "area" in list(tdic.keys()) and "hi_res" in list(tdic.keys()) and "low_res" in list(tdic.keys()):
        try:
            anion = tdic['anion']
            area = tdic["area"]
            radii = np.sqrt(area/np.pi)
            hrdx = tdic["hi_res"]["x"]
            hrdy = tdic["hi_res"]["y"]
            ly = [abs(hrdx[x][1][1]-hrdx[x][0][1]) for x in list(hrdx.keys())]
            lx = [abs(hrdy[x][1][0]-hrdy[x][0][0]) for x in list(hrdy.keys())]
            l = min([max(lx), max(ly)])
            lrdx = tdic["low_res"]["x"]
            lrdy = tdic["low_res"]["y"]
            lrly = [abs(lrdx[x][1][1]-lrdx[x][0][1]) for x in list(lrdx.keys())]
            lrlx = [abs(lrdy[x][1][0]-lrdy[x][0][0]) for x in list(lrdy.keys())]
            ll = min([max(lrly), max(lrlx)])
            radiis.append(radii)
            t_hires.append(l)
            t_lowres.append(ll)
        except Exception:
            pass


logradiis = [np.log10(x) for x in radiis]
logt_hires = [np.log10(x) for x in t_hires]
logt_lowres = [np.log10(x) for x in t_lowres]
tmin = min(radiis+t_hires+t_lowres)
tmax = max(radiis+t_hires+t_lowres)
tlinbins = np.linspace(tmin, tmax, 40)
tlogbins = np.linspace(np.log10(tmin), np.log10(tmax), 40)
#plt.hist(logradiis, tlogbins, alpha=0.5, label="Radii")
#plt.hist(logt_hires, tlogbins, alpha=0.5, label="Hi-Res")
#plt.hist(logt_lowres, tlogbins, alpha=0.5, label="Low-Res")
#plt.legend(loc='upper right')
#plt.plot(radiis, t_hires, 'b.', label = "High-res")
#plt.plot(radiis, t_lowres, 'r.', label = "Single")
#plt.plot(radiis, radiis, 'k-', label="Equals radii")
#plt.xlabel("Radii metric, eV")
#plt.ylabel("Metric, eV")
#plt.legend(loc='upper left')
#plt.show()


vols = []
vv = []
cen = []
sdq = sum_dic['q']
for name in list(sdq.keys()):
    sdqn = sdq[name]
    if 'volume' in list(sdqn.keys()):
        vv.append(sdqn['volume'])
    if 'center' in list(sdqn.keys()) and 'volume' in list(sdqn.keys()):
        cen.append(sdqn['center'])
        vols.append(sdqn['volume'])

print('Quaternary means')
print(np.mean(vv))
print(np.power(10, np.mean([np.log10(x) for x in vv])))
#min_x = min((min(cen), min(vols)))
#max_x = max((max(cen), max(vols)))
#linbins = np.linspace(min_x, max_x, 40)
#logbins = np.linspace(np.log10(min_x), np.log10(max_x), 40)
#logvols = [np.log10(x) for x in vv]

#logbinrange = [np.log10(x) for x in binrange]
#allrange = binrange+aradiis+vv
#alllin = np.linspace(min(allrange), 10, 50)
#lllog = np.linspace(np.log10(min(allrange)), np.log10(max(allrange)), 25)
#aradiis = [x*2 for x in aradiis]
#vv = [x*2 for x in vv]
#bin_fit = expon.fit(binrange, floc=0)
#ter_fit = expon.fit(aradiis, floc=0)
#qua_fit = expon.fit(vv, floc=0)
#print(qua_fit)
#for var in [t_hires, t_lowres, cen]:
#    print('A')
#    print(np.mean(var))
#    print(np.power(10, np.mean([np.log10(x) for x in var])))


#plt.hist(vv, alllin, alpha=0.4, label="Quartenary")
#plt.hist(binrange, alllin, normed=True, alpha=0.4, label="Binary")
#plt.plot(alllin, bin_reg, 'k-')
#plt.show()
#plt.hist([x*2 for x in aradiis], alllin, alpha=0.4, label="Ternary")

#plt.legend(loc="upper right")
#print(np.median(binrange))
#print(np.median(aradiis)*2)
#print(np.median(vv)*2)
#plt.show()

# Fit from peak to end with exponential

#print(len(cen))
#print(len(vols))
#plt.hist(radiis, tlinbins, alpha=0.3, label="Radius", color="red")
#plt.hist(t_hires, tlinbins, alpha=0.3, label="Mesh", color="blue")
#plt.hist(t_lowres, tlinbins, alpha=0.3, label="Center-point", color="yellow")
popt1, pcov1 = curve_fit(linfunc, vols, [x[0] for x in cen])
perr1 = np.sqrt(np.diag(pcov1))
#popt2, pcov2 = curve_fit(linfunc, radiis, t_lowres)
#perr2 = np.sqrt(np.diag(pcov2))
plt.plot(vols, vols, 'k-', label="Equals radii")
plt.plot(vols, [x[0] for x in cen], 'b.', label="Center-point")
plt.plot(vols, [x*popt1[0] for x in vols], 'b--', label="Center fit")
#plt.plot(radiis, t_lowres, 'b.', label="Center-point")
#plt.plot(radiis, [x*popt2[0] for x in radiis], 'b--', label="Center fit")
#plt.hist(cen, linbins, alpha=0.5, label="Center-point")
#plt.hist(vols, linbins, alpha=0.5, label="Radius")
#plt.hist([np.log10(x) for x in cen], logbins, alpha=0.5, label="Center-point")
#plt.hist([np.log10(x) for x in vols], logbins, alpha=0.5, label="Radius")
key_q = {}
all = []
for item in cen:
    anion = item[-1]
    all.append(item[0])
    if anion not in list(key_q.keys()):
        key_q[anion] = [item[0]]
    else:
        key_q[anion].append(item[0])

bins = np.linspace(min(all), max(all), 40)
#for anion in list(key_q.keys()):
#    plt.hist(key_q[anion], bins, alpha=0.4, label=anion)
plt.title('Distribution of the metric for quaternary materials')
plt.legend(loc='upper left')
#plt.ylabel('Number of materials')
plt.ylabel('Metric, eV')
plt.xlabel('Radius, eV')
print(popt1)
print(perr1)
#print(popt2)
#print(perr2)
plt.show()