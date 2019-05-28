import numpy as np
from pymatgen import MPRester
from matplotlib import pyplot as plt
m = MPRester("TnKFxKZTef8PsbvDG5C9")
import json


def get_compounds(eoi, compulsory = None):
    #print('Starting...')
    # For a query we set criteria that have to match
    # and properties we want to get from the database
    compounds = []
    if compulsory is None:
        criteria = {"elements": {"$all": eoi}}
        properties = ['pretty_formula', 'formula', 'formation_energy_per_atom','final_energy']
        compounds += m.query(criteria, properties)
    else:
        # Usually the anion would be the compulsory
        # This function gives the combinatorics so we can get say BaO and SnO in a Ba/Sn/O system
        criteria = {"elements": {"$in": [x for x in eoi if x != compulsory], '$all': [compulsory]}}
        properties = ['pretty_formula', 'formula', 'formation_energy_per_atom', 'final_energy']
        cur = m.query(criteria, properties)
        compounds = compounds + cur
    # Remove duplicates
    singular_c = []
    names = []
    # Filter out multiple phases of the same material - pick the one with the lowest formation energy
    # By definition these materials will always be more stable at any potential anyway
    for c in compounds:
        if c['formula'] not in names:
            names.append(c['formula'])
            singular_c.append(c)
        else:
            incumbent = [k for k in singular_c if k['formula'] == c['formula']][0]
            if c['formation_energy_per_atom'] < incumbent['formation_energy_per_atom']:
                singular_c[singular_c.index(incumbent)] = c
    compounds = singular_c


    #Filter out stuff that contains other elements
    p_compounds = []
    for comp in compounds:
        form = list(comp['formula'].keys())
        check = True
        for el in form:
            if el not in eoi:
                check = False
        if check:
            p_compounds.append(comp)
    compounds = p_compounds


    output_dic = {}
    for dic in compounds:
        comp = []
        name = dic['pretty_formula']
        compo = dic['formula']
        elements = []
        natom = 0
        # Reorganize data
        for el in eoi:
            if el in list(compo.keys()):
                elements.append([el, compo[el]])
                natom += compo[el]
        comp.append(elements)
    # Convert "per atom"  to "per molecule"
        hf = dic['formation_energy_per_atom']
        hf_fixed = hf*natom
        comp.append(hf_fixed)
        comp.append(dic['final_energy'])
        output_dic[name] = comp
    return(output_dic)


def solve_binary(comp, dic):
    tardic = dic[comp]
    tcc = tardic[0][0][1]
    tac= tardic[0][1][1]
    te = tardic[1]
    lt_i = [0]
    gt_i = [te/tcc]
    for mat in list(dic.keys()):
        if mat != comp:
            matdic = dic[mat]
            formula = matdic[0]
            energy = matdic[1]
            anion_coeff = formula[1][1]
            cat_coeff = formula[0][1]
            mult_e = te*anion_coeff-energy*tac
            mult_c = tcc*anion_coeff-cat_coeff*tac
            if mult_c == 0:
                if energy/cat_coeff < te/tac:
                    return None
                else:
                    pass
            elif mult_c < 0:
                lt_i.append(mult_e/mult_c)
            else:
                gt_i.append(mult_e/mult_c)
    high = min(lt_i)
    low = max(gt_i)
    if low < high:
        return [low, high]
    else:
        return None
with open('comb.txt', 'r') as comb:
    a = comb.readlines()

anions = ['O', 'S', 'Se', 'Te']
combs = []
elements = []
for line in a:
    k = list(eval(line))
    for item in k:
        if item not in elements and item not in anions:
            elements.append(item)
            for anion in anions:
                combs.append([item, anion])



#cache_list = []
#for comb_i in range(len(combs)):
#    comb = combs[comb_i]
#    if comb_i % 10 == 0:
#        print('{0}/{1}'.format(comb_i, len(combs)))
#    c = get_compounds(comb, comb[-1])
#    cache_list.append([comb, c])


#print(cache_list[2])
with open('binarycache.json', 'r') as bc:
    cache_list = json.load(bc)
print(cache_list[2])
count = 0
ranges = []
for s in range(len(cache_list)):
    sets = cache_list[s]
    eoi = sets[0]
    anion = eoi[-1]
    data = sets[1]
    for comp in list(data.keys()):
        count += 1
        r = solve_binary(comp, data)
        if r is not None:
            ranges.append((comp, r, anion))
print(count)
print(len(ranges))
sizes = {}
for item in ranges:
    if item[-1] not in list(sizes.keys()):
        sizes[item[-1]] = 1
    else:
        sizes[item[-1]] += 1
print(sizes)
sizes = {}
total = []
for mat_range in ranges:
    size = abs(mat_range[1][1]-mat_range[1][0])
    if mat_range[-1] not in list(sizes.keys()):
        sizes[mat_range[-1]] = [size]
    else:
        sizes[mat_range[-1]].append(size)
    total.append(size)

histbins = np.linspace(0, max(total), 40)
lhistbins = np.linspace(np.log10(min(total)), np.log10(max(total)), 50)
plt.xlabel('Log stability range, eV')
plt.ylabel('Number of materials')
plt.hist([np.log10(x) for x in total], lhistbins)
for anions in list(sizes.keys()):
    print(anions)
    print(np.power(10, np.mean([np.log10(x) for x in sizes[anions]])))
    print(np.mean(sizes[anions]))
plt.show()

