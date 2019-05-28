from pymatgen import MPRester
import pymatgen as pmg
from itertools import combinations
import re
import numpy as np
import matplotlib.pyplot as plt
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


def find_2d_intersect(lineA, lineB, axesnames = ['x', 'y']):
    lineAC, lineAE = lineA[0], lineA[1]
    lineBC, lineBE = lineB[0], lineB[1]
    Adic = {}
    Bdic = {}
    for axis in axesnames:
        if axis in [x[0] for x in lineAC]:
            Adic[axis] = [x[1] for x in lineAC if x[0] == axis][0]
        else:
            Adic[axis] = 0
        if axis in [x[0] for x in lineBC]:
            Bdic[axis] = [x[1] for x in lineBC if x[0] == axis][0]
        else:
            Bdic[axis] = 0
    axis = axesnames
    # Normal case
    if 0 not in [Adic[k] for k in list(Adic.keys())] and 0 not in [Bdic[k] for k in list(Bdic.keys())]:
        if Adic[axis[1]] / Adic[axis[0]] == Bdic[axis[1]] / Bdic[axis[0]]:
            return None
        else:
            axis_one_i = (Adic[axis[1]]*lineBE-Bdic[axis[1]]*lineAE)/(Adic[axis[1]]*Bdic[axis[0]]-Bdic[axis[1]]*Adic[axis[0]])
            axis_two_i = (lineAE - axis_one_i*Adic[axis[0]])/Adic[axis[1]]
    # If there is a zero
    else:
        p = [[k, Adic[k]] for k in list(Adic.keys()) if Adic[k] == 0]+[[k, Bdic[k]] for k in list(Bdic.keys()) if Bdic[k] == 0]
        if len(p) == 1:
            if [[k, Adic[k]] for k in list(Adic.keys()) if Adic[k] == 0] != [[]]:
                stat = Bdic
                E = lineBE
            else:
                stat = Adic
                E = lineAE
            if p[0][0] == axis[0]:
                axis_one_i = 0
                axis_two_i = E/(stat[[axis[1]]])
            else:
                axis_two_i = 0
                axis_one_i = E / (stat[[axis[0]]])
        else:
            return {axis[0]: [x[1] for x in p if x[0] == axis[0]][0], axis[1]:[x[1] for x in p if x[0] ==axis[1]][0]}
    return {axis[0]: axis_one_i, axis[1]: axis_two_i}


# Get total energies of compounds of interest
def get_compounds(eoi, compulsory = None):
    print('Starting...')
    # For a query we set criteria that have to match
    # and properties we want to get from the database
    compounds = []
    if compulsory is None:
        criteria = {"elements": {"$all": eoi}}
        properties = ['pretty_formula', 'material_id', 'formation_energy_per_atom']
        compounds += m.query(criteria, properties)
    else:
        # Usually the anion would be the compulsory
        # This function gives the combinatorics so we can get say BaO and SnO in a Ba/Sn/O system
        opt = []
        for k in range(1, len(eoi)+1):
            opt = opt + [list(a) for a in list(combinations([x for x in eoi], k))]
        for com in opt:
            # com.append(compulsory)
            criteria = {"elements": {"$all": com}, "nelements": len(com)}
            properties = ['unit_cell_formula', 'final_energy', 'pretty_formula', 'material_id', 'formation_energy_per_atom']
            cur = m.query(criteria, properties)
            compounds = compounds + cur

    # Remove duplicates
    singular_c = []
    names = []
    for c in compounds:
        if c['pretty_formula'] not in names:
            names.append(c['pretty_formula'])
            singular_c.append(c)
        else:
            incumbent = [k for k in singular_c if k['pretty_formula'] == c['pretty_formula']][0]
            if c['final_energy'] < incumbent['final_energy']:
                singular_c[singular_c.index(incumbent)] = c
    compounds = singular_c


    output_dic = {}
    for dic in compounds:
    # String manipulation of the formula to split them to atoms
    # Probably faster than using PMG API
        comp = []
        name = dic['pretty_formula']
        comp.append(dic['unit_cell_formula'])
        hf = dic['final_energy']
        comp.append(hf)
        # if natom < 3*len(eoi):
        output_dic[name] = comp
    return(output_dic)



# Function generates the inequalities that are used in other uses
# The function iterates over all compounds in the list and
# generates all possible inequalities of them with other phases
# in form of a dictionary
def create_inequality(target, compounds, eoi, control_element = None):
    # Assumes that the anion is placed last
    # For now, 2D
    if len(eoi) == 3 or len(eoi) == 2:
        inequalities = {}
        if control_element is None:
            control_element = eoi[-1]
        target_comp = compounds[target]

        inequalities[target] = [[x for x in target_comp[0].items() if x[0] != control_element], target_comp[1]]

        # Number of control atoms (dep var) in the target compound
        control_q = [x[1] for x in target_comp[0].items() if x[0] == control_element][0]

        for element in eoi:
            if element != control_element:
                inequalities[element] = [[[element, 1]], 0]
        for compname in [x for x in list(compounds.keys()) if x != target]:
            compound = compounds[compname]
            form_e = compound[-1]
            elements = compound[0]
            print(elements)
            # Solve such that the control element is 0
            if control_element not in [x[0] for x in elements]:
                # ? why?
                pass
            else:
                eqn_num = []
                # comp_q = [x[1] for x in elements.items() if x[0] == control_element][0]
                # lcm_q = lcm(comp_q, control_q)
                # comp_f, cont_f = lcm_q/comp_q, lcm_q/control_q
                # for el in [x for x in eoi if x != control_element]:
        #             cont_el_cur = [x[1] for x in target_comp[0] if x[0] == el]
        #             comp_el_cur = [x[1] for x in elements if x[0] == el]
        #             if len(cont_el_cur) == 0:
        #                 cont_el_cur = [0]
        #             if len(comp_el_cur) == 0:
        #                 comp_el_cur = [0]
        #             eqn_num.append([el, cont_el_cur[0]*cont_f-comp_el_cur[0]*comp_f])
                # result = [eqn_num, target_comp[1]*cont_f-form_e*comp_f]
                # inequalities[compname] = result

        # return inequalities


def plot_inequalities(i, target, axes, ppoints = None):
    if len(axes) == 2:
        f_eq = i[target]
        xlim = [f_eq[1]/f_eq[0][0][1], 0]
        ylim = [f_eq[1]/f_eq[0][1][1], 0]
        xrange = np.linspace(xlim[0], xlim[1], 20)
        yvals = [(f_eq[1]-x*f_eq[0][0][1])/f_eq[0][1][1] for x in xrange]
        eqns = {target: [xrange, yvals]}
        for ineq in list(i.keys()):
            if len(i[ineq][0])!=1:
                eqn = i[ineq]
                yvals = [(eqn[1]-x*eqn[0][0][1])/eqn[0][1][1] for x in xrange]
                eqns[ineq] = [xrange, yvals]
        for n in list(eqns.keys()):
            yvals = eqns[n][1]
            plt.plot(xrange, yvals, label=n)
        plt.xlim(xlim[0], xlim[1])
        plt.ylim(ylim[0], ylim[1])
        plt.xlabel('Potential of {0}'.format(axes[0]))
        plt.ylabel('Potential of {0}'.format(axes[1]))
        plt.legend()
        if ppoints is not None:
            ppoints_x = [x[axes[0]] for x in ppoints]
            ppoints_y = [x[axes[1]] for x in ppoints]
            plt.fill(ppoints_x, ppoints_y)
        plt.show()
    elif len(axes) == 1:
        f_eq = i[target]
        xlim = [f_eq[1] / f_eq[0][0][1], 0]
        plt.figure()
        for x in list(i.keys()):
            comp = i[x]
            plt.plot([comp[1]/comp[0][0][1]],[0], '.', label = x)
        plt.xlabel('Potential of {0}'.format(axes[0]))
        plt.xlim(xlim[0]-1, xlim[1])
        plt.ylim(-0.5, 0.5)
        plt.legend()
        plt.show()

    else:
        pass


def get_polygon_points(i, target, axes):
    polygon_points = []
    competing_phases = []
    f_eq = i[target]
    xlim = [f_eq[1] / f_eq[0][0][1], 0]
    ylim = [f_eq[1] / f_eq[0][1][1], 0]
    ipoints = []
    for k in list(i.keys()):
        if k != target:
            if len(i[k][0]) > 1:
                ipoint = find_2d_intersect(f_eq, i[k], axes)
                ipoints.append([k, ipoint])
    # Sort by x-intercept
    ipoints.sort(key=lambda x: x[1][axes[0]])
    stats = []
    test = True
    for stuff in ipoints:
        ineq = i[stuff[0]]
        if len(ineq[0]) >1:
        # Until the "turning point"
            if ineq[0][1][1] > 0 and len(stats) == 0:
                stats.append(1)
                before = ipoints[ipoints.index(stuff)-1]
                after = stuff
                polygon_points.append(after[1])
                polygon_points.append(before[1])
                competing_phases.append(after[0])
                competing_phases.append(before[0])
            elif ineq[0][1][1] < 0 and len(stats) != 0:
                stats.append(-1)
                if ineq[0][1][1]*ineq[0][0][1] > 0:
                    pass
                else:
                    print(stuff)
                    return None
    # Focus on the "before" line


	#   Loops around the Polygon
    before_name = before[0]
    before_line = i[before[0]]
    before_int = []
    a_number = 1
    axiscode = 0
    prev = target
    while True:
        for k in list(i.keys()):
            if k != before_name:
                if len(i[k][0]) > 1:
                    bpoint = find_2d_intersect(before_line, i[k], axes)
                    if bpoint is not None:
                        before_int.append([k, bpoint])
        before_int.sort(key=lambda x: x[1][axes[axiscode]])
        # Find the index number of target - then the next one on the X-axis
        p = before_int.index([x for x in before_int if x[0] == prev][0])
        next_item = before_int[p+a_number]
        if next_item[0] != after[0]:
            polygon_points.append(next_item[1])
            competing_phases.append(next_item[0])
            prev = before_name
            before_name = next_item[0]
            before_line = i[next_item[0]]
            before_int = []
            c = i[next_item[0]]
            if len(c[0]) > 1:
                if c[0][1][1] > 0:
                    axiscode = 0
                    a_number = -1
                else:
                    axiscode = 0
                    a_number = 1
            # X-intercept
            elif c[0][0][0] == axes[0]:
                axiscode = 0
                a_number = 1
            # Y-intercept
            else:
                axiscode = 1
                a_number = -1
        else:
            polygon_points.append(next_item[1])
            break
    return [polygon_points, competing_phases]


# elements_of_interest = ['Zn', 'Sn','S']
elements_of_interest = ['Ba', 'Sn','O']
# c = get_compounds(elements_of_interest, 'S')
c = get_compounds(elements_of_interest, 'O')
# print(c)
# i = create_inequality('ZnSnS3', c, elements_of_interest, 'S')
i = create_inequality('BaSnO3', c, elements_of_interest, 'O')
# print(i)
# # p_points = get_polygon_points(i, 'ZnSnS3', ['Zn', 'Sn'])
# p_points = get_polygon_points(i, 'BaSnO3', ['Ba', 'Sn'])
# if p_points is not None:
#     plot_inequalities(i, 'BaSnO3', ['Ba', 'Sn'], p_points[0])
# else:
#     plot_inequalities(i, 'BaSnO3', ['Ba', 'Sn'])
## print(i)

