from pymatgen import MPRester
import pymatgen as pmg
from itertools import combinations
import re
# I have an enviornment variable set up for my API key
# But you can get your API key from your materials project dashboard
# and put it directly in MPRester() as a string
m = MPRester("TnKFxKZTef8PsbvDG5C9")


# Just a convenience function
def lcm(x, y):
   if x > y:
       greater = x
   else:
       greater = y
   while(True):
       if((greater % x == 0) and (greater % y == 0)):
           lcm = greater
           break
       greater += 1

   return lcm


# Get total energies of compounds of interest
def get_compounds(eoi, compulsory = None, raw=False):
    # For a query we set criteria that have to match
    # and properties we want to get from the database
    if compulsory is None:
        criteria = {"elements": {"$all": eoi}}
        properties = ['pretty_formula', 'material_id', 'final_energy']
        compounds = m.query(criteria,properties)
    else:
    # Usually the anion would be the compulsory
    # This function gives the combinatorics so we can get say BaO and SnO in a Ba/Sn/O system
        opt = []
        for k in range(1, len(eoi)):
            opt = opt + [list(a) for a in list(combinations([x for x in eoi if x != compulsory], k))]
        compounds = []
        for com in opt:
            com.append(compulsory)
            criteria = {"elements": {"$all": com}}
            properties = ['pretty_formula', 'material_id', 'final_energy']
            cur = m.query(criteria, properties)
            compounds = compounds + cur

    # Remove duplicates
    singular_c = []
    names = []
    for c in compounds:
        if c['pretty_formula'] not in names:
            names.append(c['pretty_formula'])
            singular_c.append(c)
    compounds = singular_c


    # Results and filter out any additions
    filtered_c = []
    for compound in compounds:
        comp = pmg.Composition(compound['pretty_formula'])
        n = 0
        for element in list(comp.keys()):
            if str(element) not in eoi:
                n += 1
        if n == 0:
            filtered_c.append(compound)
    compounds = filtered_c
    if raw:
        return compounds
    output_dic = {}
    for dic in compounds:
    # String manipulation of the formula to split them to atoms
    # Probably faster than using PMG API
        comp = []
        name = dic['pretty_formula']
        splitname = re.findall('[A-Z][^A-Z]*', name)
        elements = []
        for part in splitname:
            el = ''.join([i for i in part if not i.isdigit()])
            dig = ''.join([i for i in part if i.isdigit()])
            if dig == '':
                dig = '1'
            elements.append([el, dig])
        comp.append(elements)
        comp.append(dic['final_energy'])
        output_dic[name] = comp
    return(output_dic)



# Function generates the inequalities that are used in other uses
# The function iterates over all compounds in the list and
# generates all possible inequalities of them with other phases
# in form of a dictionary
def create_inequalities(compounds, eoi, control_element = None):
    # Assumes that the anion is placed last
    inequalities = {}
    formulas = {}
    if control_element is None:
        control_element = eoi[-1]
    for compound in compounds:
        comp = pmg.Composition(compound['pretty_formula'])
        dict = comp.get_el_amt_dict()
        lhs = []
        for k in list(dict.keys()):
            lhs.append([k, dict[k]])
        formulas[compound['pretty_formula']] = lhs
    for compound in compounds:
        per_comp = []
        c_formula = formulas[compound['pretty_formula']]
        # Basic inequalities
        per_comp.append([c_formula, '<', compound['final_energy']])
        for items in c_formula:
            if items[0] != control_element:
                per_comp.append([[items], '>', compound['final_energy']])
                per_comp.append([[items], '<', 0])
        # Generates inequalities outside the first one
        for o_comp in compounds:
            if o_comp != compound:
                o_formula = formulas[o_comp['pretty_formula']]
                compul = [x[1] for x in c_formula if x[0] == control_element][0]
                o_compul = [x[1] for x in o_formula if x[0] == control_element][0]
                f = lcm(compul, o_compul)
                constants = [[x[0], x[1]*f/compul] for x in c_formula]
                o_constants = [[x[0], x[1]*f/o_compul] for x in o_formula]
                f_constants = []
                for elm in [x for x in eoi if x!= control_element]:
                    el = [x[1] for x in constants if x[0] == elm]
                    sel = [x[1] for x in o_constants if x[0] == elm]
                    if len(el) == 0:
                        el = [0]
                    if len(sel) == 0:
                        sel = [0]
                    f_constants.append([elm, el[0]-sel[0]])
                per_comp.append([o_comp['pretty_formula'], f_constants, '>',
                                f*(compound['final_energy']/compul-o_comp['final_energy']/o_compul)])
        inequalities[compound['pretty_formula']] = per_comp
    return inequalities


def print_inequalities(ineq_dic):
    for k in list(ineq_dic.keys()):
        print(k)
        ineq_list = ineq_dic[k]
        for ineq_raw in ineq_list:
            if len(ineq_raw) == 4:
                ineq_raw = ineq_raw[1::]
            ineq_str = ''
            ineq_els = ineq_raw[0]
            for ineq_el in ineq_els:
                if ineq_el[1] != 0:
                    if len(ineq_str) == 0 or ineq_el[1]<0:
                        ineq_str += '{0}{1}'.format(int(ineq_el[1]), ineq_el[0])
                    else:
                        ineq_str += '+{0}{1}'.format(int(ineq_el[1]), ineq_el[0])
            ineq_str += ' {0} {1:.3f}'.format(ineq_raw[1], ineq_raw[2])
            print(ineq_str)
        print('\n')


# Create an output file that follows the format of input.dat
def create_output_file(compound_of_interest, dic, dependent, filename = 'input.dat'):
    if compound_of_interest not in list(dic.keys()):
        print('Compound not listed!')
        return None
    lines = []
# Setup the part for the compound of interest
    coi_data = dic[compound_of_interest]
    elements_coi = coi_data[0]
    energy_coi = coi_data[-1]
    lines.append(str(len(elements_coi)))
    line = []
    for element_coi in elements_coi:
        line.append(str(int(element_coi[1])))
        line.append(element_coi[0])
    line.append(str(energy_coi))
    lines.append(' '.join(line))
# Dependent element and number of competing phases
    lines.append(dependent)
    lines.append(str(len(list(dic.keys()))-1))
# Individual competing phases
    for comp in list(dic.keys()):
        if comp != compound_of_interest:
            line = []
            comp_data = dic[comp]
            elements_comp = comp_data[0]
            energy_comp = comp_data[-1]
            # Number of elements in compound
            lines.append(str(len(elements_comp)))
            for element_comp in elements_comp:
                line.append(str(int(element_comp[1])))
                line.append(element_comp[0])
            line.append(str(energy_comp))
            lines.append(' '.join(line))


    # Output
    with open(filename, 'w') as file:
        for fline in lines:
            file.write(fline)
            file.write('\n')


elements_of_interest = ['Cu', 'Sn', 'S']
c = get_compounds(elements_of_interest, 'S', raw=True)
print(str(c).replace('_', '$_$'))
create_output_file('ZnCu2SnS4', c, 'S', 'input1.dat')
