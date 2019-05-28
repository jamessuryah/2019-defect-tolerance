import time
start_time = time.time()
import json
from pymatgen import MPRester
import random
API_key = 'TnKFxKZTef8PsbvDG5C9'
mpr = MPRester(API_key)
# criteria = {"elements": {'$in': ['O', 'S', 'Se', 'Te']}}
# properties = ['pretty_formula', 'formula', 'formation_energy_per_atom', 'final_energy']
# cur = mpr.query(criteria, properties)
# with open('temp_cache_q.json', 'w') as tc:
#     json.dump(cur, tc)

with open('temp_cache_q.json', 'r') as tc:
    cache = json.load(tc)


cations = ['Li', 'Be', 'Na', 'Mg', 'Al', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Cs', 'Ba', 'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi']
print(len(cations))
anions = ['O', 'S', 'Se', 'Te']


def convertformula(f):
    nf = []
    for e in list(f.keys()):
        if e not in ['O', 'S', 'Se', 'Te']:
            nf.append([e, f[e]])
    for e in list(f.keys()):
        if e in ['O', 'S', 'Se', 'Te']:
            nf.append([e, f[e]])
    return nf

anion_count = [0, 0, 0, 0]
# filt = []
# for item in cache:
#      # Filters out all items that contains stuff outside above
#      comp = item['formula']
#      elements = list(comp.keys())
#      cmatch = 0
#      amatch = 0
#      for k in elements:
#          if k in cations:
#              cmatch += 1
#          elif k in anions:
#              amatch += 1
#      if cmatch + amatch == len(elements) and amatch == 1 and len(elements) > 1:
#          filt.append(item)
#
#
# with open('prefiltered_cache.json', 'w') as pfc:
#     json.dump(filt, pfc, indent=4)


with open('prefiltered_cache.json', 'r') as pfc:
    filtered = json.load(pfc)


q_phases = []
for c in filtered:
    cf = c['formula']
    if len(list(cf.keys())) == 4:
        q_phases.append(c)


final_form = []
for k in range(len(q_phases)):
    qp = q_phases[k]
    qpf = qp['formula']
    elements = set(qpf.keys())
    compete = {}
    for obj in filtered:
        of = obj['formula']
        if set(of.keys()).issubset(elements):
            compete[obj['pretty_formula']] = [convertformula(of), obj['formation_energy_per_atom']*sum(of.values()), obj['final_energy']]
    final_form.append([qp['pretty_formula'], compete])
    if k%250 == 0:
        print(k)


print(len(final_form))
with open('to_check_q.json', 'w') as tcq:
    json.dump(final_form, tcq)
# Filters out all items that contains stuff outside above
print(time.time()-start_time)


