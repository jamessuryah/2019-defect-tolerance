from pymatgen import MPRester
import pymatgen as pmg
from itertools import combinations
import re
from operator import add, sub, mul
from sympy import Plane, Point3D, Line3D
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Path3DCollection, Line3DCollection
from scipy.spatial import ConvexHull as CHull
from polyhedron import Vrep, Hrep
import time
import json
import scipy
np.set_printoptions(precision=4)
# I have an enviornment variable set up for my API key
# But you can get your API key from your materials project dashboard
# and put it directly in MPRester() as a string
m = MPRester("TnKFxKZTef8PsbvDG5C9")


def ncross(a, b):
    c = [a[1]*b[2] - a[2]*b[1],
         a[2]*b[0] - a[0]*b[2],
         a[0]*b[1] - a[1]*b[0]]
    size = np.sqrt(c[0]**2+c[1]**2+c[2]**2)
    c = [x/size for x in c]
    return c


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
        inequalities[target] = [[x for x in target_comp[0] if x[0] != control_element], target_comp[1]]
        # Number of control atoms (dep var) in the target compound
        control_q = [x[1] for x in target_comp[0] if x[0] == control_element][0]
        for element in eoi:
            if element != control_element:
                inequalities[element] = [[[element, 1]], 0]
        for compname in [x for x in list(compounds.keys()) if x != target]:
            compound = compounds[compname]
            form_e = compound[-1]
            elements = compound[0]
        # Solve such that the control element is 0
            if control_element not in [x[0] for x in elements]:
                pass
            else:
                eqn_num = []
                comp_q = [x[1] for x in elements if x[0] == control_element][0]
                lcm_q = lcm(comp_q, control_q)
                comp_f, cont_f = lcm_q/comp_q, lcm_q/control_q
                for el in [x for x in eoi if x != control_element]:
                    cont_el_cur = [x[1] for x in target_comp[0] if x[0] == el]
                    comp_el_cur = [x[1] for x in elements if x[0] == el]
                    if len(cont_el_cur) == 0:
                        cont_el_cur = [0]
                    if len(comp_el_cur) == 0:
                        comp_el_cur = [0]
                    eqn_num.append([el, cont_el_cur[0]*cont_f-comp_el_cur[0]*comp_f])
                result = [eqn_num, target_comp[1]*cont_f-form_e*comp_f]
                inequalities[compname] = result
            # Basic inequalities
            #per_comp.append([c_formula, '<', compound['formation_energy_per_atom']])
            #for items in c_formula:
            #    if items[0] != control_element:
            #        per_comp.append([items, '>', compound['formation_energy_per_atom']])
            #        per_comp.append([items, '<', 0])
            # Generates inequalities outside the first one
            #for o_comp in compounds:
            #    if o_comp != compound:
            #        o_formula = formulas[o_comp['pretty_formula']]
            #        compul = [x[1] for x in c_formula if x[0] == control_element][0]
            #        o_compul = [x[1] for x in o_formula if x[0] == control_element][0]
            #        f = lcm(compul, o_compul)
            #        constants = [[x[0], x[1]*f/compul] for x in c_formula]
            #        o_constants = [[x[0], x[1]*f/o_compul] for x in o_formula]
            #        f_constants = []
            #        for elm in [x for x in eoi if x!= control_element]:
            #            el = [x[1] for x in constants if x[0] == elm]
            #            sel = [x[1] for x in o_constants if x[0] == elm]
            #            if len(el) == 0:
            #                el = [0]
            #            if len(sel) == 0:
            #                sel = [0]
            #            f_constants.append([elm, el[0]-sel[0]])
            #        per_comp.append([o_comp['pretty_formula'], f_constants, '>',
            #                        f*(compound['formation_energy_per_atom']/compul-o_comp['formation_energy_per_atom']/o_compul)])
            #inequalities[compound['pretty_formula']] = per_comp
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


# Create an output file that follows the input format of cplapy
def create_output_file(compound_of_interest, dic, dependent, filename = 'input.dat'):
    if compound_of_interest not in list(dic.keys()):
        print('Compound not listed!')
        return None
    lines = []
# Setup the part for the compound of interest
    coi_data = dic[compound_of_interest]
    elements_coi = coi_data[0]
    energy_coi = coi_data[-2]
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
            energy_comp = comp_data[-2]
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


def get_mu_dep(A0, b0, mu, indx_dep):
    a_dep = A0[indx_dep]
    A0 = np.delete(A0, indx_dep)
    mu_dep = (b0 - np.dot(A0, mu)) / a_dep
    return mu_dep[0]


def read_cplap_input(file):
    '''
    read cplap input file
    and return A, b, and lists of secondary phase and elements
    '''
    A = []
    b = []
    # elements
    sec_phases = []

    with open(file) as f:
        lines = [line for line in f.readlines() if not line[0] == '#']
        n_elements = int(lines.pop(0))

        # read product
        line = lines.pop(0).split()
        b_i = float(line.pop())
        elements = line[1::2]
        A_i = np.array(line[0::2], dtype=int)
        dep_elem = lines.pop(0).strip()
        indx_dep = elements.index(dep_elem)
        dep_coeff = float(A_i[indx_dep])

        # product
        A.append(-A_i)
        b.append(-b_i)
        sec_phases.append(''.join(['{}{}'.format(i, j) for i, j in zip(elements, A_i)]))
        # elemental each
        for i in np.delete(range(len(elements)), [indx_dep]):
            A_ie = np.zeros(len(elements))
            A_ie[i] = A_i[i]
            A.append(A_ie)
            b.append(0)
            sec_phases.append(''.join(['{}{}'.format(i, j) for i, j in zip(elements, A_ie)]))
            A.append(-A_ie)
            b.append(-b_i)
            sec_phases.append(''.join(['{}{}'.format(i, j) for i, j in zip(elements, A_ie)]))
        # secondary phases
        no_sec = int(lines.pop(0))
        for i_sec in range(no_sec):
            n_elem = int(lines.pop(0))  # ignore

            line = lines.pop(0).split()
            b_i = float(line.pop())

            elem_i = line[1::2]
            indx_elem = [elements.index(elem) for elem in elem_i]
            A_it = np.array(line[0::2], dtype=int)
            A_i = np.array([A_it[indx_elem.index(i)]
                            if i in indx_elem else 0
                            for i in range(len(elements))])

            A.append(A_i + A[0] / dep_coeff * A_i[indx_dep])
            b.append(b_i + b[0] / dep_coeff * A_i[indx_dep])
            sec_phases.append(''.join(['{}{}'.format(i, j) for i, j in zip(elem_i, A_it)]))
    # remove dependent column
    A0 = A[0]  # to calculate mu_dependent
    A = np.delete(A, indx_dep, axis=1)
    b = np.array(b).reshape((len(b), 1))
    elements.pop(indx_dep)

    mu_dep = lambda mu: get_mu_dep(A0, b[0], mu, indx_dep)
    return A, b, elements, sec_phases, mu_dep, dep_elem


def sort_vert(ininc_i, adj):
    ininc_sorted = []
    ininc_i = list(ininc_i)
    while len(ininc_i) > 0:
        v = ininc_i.pop()
        ininc_sorted.append(v)
        # find adj
        adj_i = adj[v]
        ininc_i = sorted(ininc_i, reverse=True,
                         key=lambda x: np.where(np.concatenate([adj_i, np.arange(1000)]) == x)[0][0])
    return ininc_sorted


def draw_plane(ax, verts, ininc, adj, sec_phases=[]):
    # find number of polygons with verts
    if len([half for half in ininc if len(half) > 0]) > 8:
        cmap = plt.get_cmap("tab10_r")
    else:
        cmap = plt.get_cmap("Set1")

    # draw plane
    color_counter = 0
    for i, ininc_i in enumerate(ininc):
        if len(ininc_i) < 3:
            continue

        ininc_i = sort_vert(ininc_i, adj)
        x = []
        y = []
        z = []
        for v in ininc_i:
            x.append(verts[v][0])
            y.append(verts[v][1])
            z.append(verts[v][2])
        x.append(verts[ininc_i[0]][0])
        y.append(verts[ininc_i[0]][1])
        z.append(verts[ininc_i[0]][2])
        coord = [list(zip(x, y, z))]

        label = sec_phases[i]
        polygon = Poly3DCollection(coord, alpha=0.9, label=label, closed=True)
        polygon.set_facecolor(cmap(color_counter))
        color_counter += 1
        polygon._facecolors2d = polygon._facecolors3d
        polygon._edgecolors2d = polygon._edgecolors3d

        ax.add_collection3d(polygon)
        path = Line3DCollection(coord, lw=2, color='k')
        ax.add_collection3d(path)
    ax.legend(loc='center left')
    ax.legend(loc='upper center', ncol=5)


def write_vert_output(generators, inc, elements, sec_phases, mu_dep, dep_elem, file='output_vert.dat'):
    print('gen', len(elements))
    with open(file, 'w') as f:
        [f.write('# {0: <7}'.format(elem)) for elem in elements + [dep_elem] + [': secondary phases']]
        f.write('\n')
        for i, vert in enumerate(generators):
            # write chemical potential
            [f.write('{0:0.4f}  '.format(v)) for v in np.hstack((vert, [mu_dep(vert)]))]
            f.write('  :')
            # write secondary phase
            [f.write(' {0: <15}'.format(sec_phases[inc_i])) for inc_i in inc[i]]
            f.write('\n')
        f.write('\n')


def write_half_output(generators, ininc, elements, sec_phases, mu_dep, dep_elem, file='output_half.dat'):
    print("no", len(ininc))
    print(ininc)
    with open(file, 'w') as f:
        f.write('# {0: <13}:'.format('Sec phase'))
        f.write(' Vertices \n')
        for i, half in enumerate(ininc):
            if len(half) == 0:
                continue
            f.write('{0: <15}:'.format(sec_phases[i]))
            for v_i in half:
                [f.write('  {0:0.4f}'.format(v)) for v in np.hstack((generators[v_i], [mu_dep(generators[v_i])]))]
                f.write(', ')
            f.write('\n')


def draw_pd(ax, file):
    """ tool_tip_missing
    """
    A, b, elements, sec_phases, mu_dep, dep_elem = read_cplap_input(file)

    p = Hrep(A, b)
    #print(p)
    #print(len(A))
    #print(len(p.inc))
    #print(len(p.ininc))
    #print(p.ininc)
    #print(p.inadj)
    #print(p.is_vertex)

    draw_plane(ax, p.generators, p.ininc, p.adj, sec_phases)
    write_vert_output(p.generators, p.inc, elements, sec_phases, mu_dep, dep_elem)
    write_half_output(p.generators, p.ininc, elements, sec_phases, mu_dep, dep_elem)

    set_axis(ax, A, b, elements)
    return p.generators


# draw_pd() withouu actually drawing
def calc_pd(file):
    A, b, elements, sec_phases, mu_dep, dep_elem = read_cplap_input(file)
    p = Hrep(A, b)
    return p.generators


def calc_pd_c(file):
    A, b, elements, sec_phases, mu_dep, dep_elem = read_cplap_input(file)
    p = Hrep(A, b)
    return p

def set_axis(ax, A, b, elements):
    # calc elemental chemical potential and set them 0
    buffer = 0.2  # eV

    mu_elements = []
    for i, ele in enumerate(elements):
        mu_i = []
        for k, a in enumerate(A):
            if b[k] == 0.: continue
            l_elemental = True
            for j, coeff in enumerate(a):
                if i != j and coeff != 0:
                    l_elemental = False
            if l_elemental:
                mu_i.append(b[k] / A[k, i])
        mu_elements.append(np.max(mu_i))
    print('mu_elements', elements, np.ravel(mu_elements))

    ax.set_xlim([-1 - buffer + mu_elements[0], 0 + buffer + mu_elements[0]])
    ax.set_ylim([-2 - buffer + mu_elements[1], -1 + buffer + mu_elements[1]])
    # ax.set_zlim([-2 -buffer + mu_elements[2], 0 + buffer + mu_elements[2]])

    ax.set_xticks(-np.arange(2) + mu_elements[0])
    ax.set_yticks(-np.arange(2) + mu_elements[1] - 1)
    # ax.set_zticks(-np.arange(3) + mu_elements[2])

    ax.set_xticklabels(-1 * np.arange(2))
    ax.set_yticklabels(-1 * np.arange(2) - 1)
    # ax.set_zticklabels(-1*np.arange(3))

    ax.set_xlabel(r'$\mu_{{{}}}$ $(eV)$'.format(elements[0]))
    ax.set_ylabel(r'$\mu_{{{}}}$ $(eV)$'.format(elements[1]))
    # ax.set_zlabel(r'$\mu_{{{}}}$ $(eV)$'.format(elements[2]))


def line_intersection(line1, line2, line1t='n', line2t='n'):
    xdiff = (line1[0][0] - line1[1][0], line2[0][0] - line2[1][0])
    ydiff = (line1[0][1] - line1[1][1], line2[0][1] - line2[1][1]) #Typo was here

    def det(a, b):
        return a[0] * b[1] - a[1] * b[0]

    div = det(xdiff, ydiff)
    if div == 0:
       return None

    d = (det(*line1), det(*line2))
    x = det(d, xdiff) / div
    y = det(d, ydiff) / div
    # Check if the points are in range
    if line1t=='v':
        testr = (line1[0][1], line1[1][1])
        test1 = min(testr)<=y<=max(testr)
    else:
        testr = (line1[0][0], line1[1][0])
        test1 = min(testr)<=x<=max(testr)
    if line1t != 'v':
        testr = (line2[0][1], line2[1][1])
        test2 = min(testr) <= y <= max(testr)
    else:
        testr = (line2[0][0], line2[1][0])
        test2 = min(testr) <= x <= max(testr)
    if test1 and test2:
        return x, y
    else:
        return None


def lineplane_intersection(line, plane, la=0):
    bbox = ((min([x[0] for x in plane]), min([x[1] for x in plane]), min([x[2] for x in plane])),
                  (max([x[0] for x in plane]), max([x[1] for x in plane]), max([x[2] for x in plane])))
    if bbox[0][la-1]<=line[0][la-1]<=bbox[1][la-1] and bbox[0][la-2]<=line[0][la-2]<=bbox[1][la-2]:
        s_line = Line3D(Point3D(line[0]), Point3D(line[1]))
        s_plane = Plane(Point3D(plane[0]), Point3D(plane[1]), Point3D(plane[2]))
        itp = s_plane.intersection(s_line)
        in_range = False
        if len(itp) > 0:
            intpoint = itp[0]
            if min([x[la] for x in line])<=intpoint[la]<=max([x[la] for x in line]):
                if bbox[0][0]<=intpoint[0]<=bbox[1][0]:
                    if bbox[0][1] <= intpoint[1] <= bbox[1][1]:
                        if bbox[0][2] <= intpoint[2] <= bbox[1][2]:
                            in_range = True
        if in_range:
            return [float(x) for x in intpoint]
        else:
            return None
    else:
        return None

def ppdist(coord, normal, ppoint):
    d = 0
    for p in range(len(coord)):
        d += normal[p]*(coord[p]-ppoint[p])
    return d

def calc_width_2d(bbox, vert, res = 20):
    xrange = bbox[0]
    yrange = bbox[1]
    xpoints = np.linspace(xrange[0], xrange[1], num=res+1, endpoint=False)[1:]
    ypoints = np.linspace(yrange[0], yrange[1], num=res+1, endpoint=False)[1:]
    intpoints = {'x':{}, 'y':{}}
    vhull = CHull(vert)
    baselines = []
    for simplex in vhull.simplices:
        line = (vert[simplex, 0], vert[simplex, 1])
        or_line = ((line[0][0], line[1][0]), (line[0][1], line[1][1]))
        baselines.append(or_line)
    for xpoint in xpoints:
        pone = (xpoint, yrange[0])
        ptwo = (xpoint, yrange[1])
        line1 = (pone, ptwo)
        for line2 in baselines:
            w = line_intersection(line1, line2, line1t='v')
            if w is not None:
                if xpoint not in list(intpoints['x'].keys()):
                    intpoints['x'][xpoint] = [w]
                else:
                    intpoints['x'][xpoint] = intpoints['x'][xpoint]+[w]
    for ypoint in ypoints:
        pone = (xrange[0], ypoint)
        ptwo = (xrange[1], ypoint)
        line1 = (pone, ptwo)
        for line2 in baselines:
            w = line_intersection(line1, line2)
            if w is not None:
                if ypoint not in list(intpoints['y'].keys()):
                    intpoints['y'][ypoint] = [w]
                else:
                    intpoints['y'][ypoint] = intpoints['y'][ypoint] + [w]
    return intpoints


def linepoly(line, n, vers, adjc, axis=0):
    # First define the points and normal
    ppoint = line[0]
    n1 = n[0]
    n2 = n[1]
    dsts = []
    lt1 = []
    gt1 = []
    for p in range(len(vers)):
        ver = vers[p]
        d = ppdist(ver, n1, ppoint)
        if d < 0:
            lt1.append(p)
        else:
            gt1.append(p)
        dsts.append(d)
    if len(lt1) == 0 or len(gt1) == 0:
        return None
    newpoly = []
    newpd = []
    newps = []
    for gt in gt1:
        adjc_points = adjc[p]
        for adjc_point in adjc_points:
            if adjc_point in lt1:
                vi = vers[gt]
                vj = vers[adjc_point]
                di = dsts[gt]
                dj = dsts[adjc_point]
                divisor = abs(di)+abs(dj)
                via = [abs(dj)*x for x in vi]
                vja = [abs(di)*x for x in vj]
                vs = [x/divisor for x in list(map(add, via, vja))]
                newpoly.append(vs)
                dtwo = ppdist(vs, n2, ppoint)
                newpd.append(dtwo)
                if dtwo < 0:
                    newps.append(-1)
                else:
                    newps.append(1)
    if 1 not in newps or -1 not in newps:
        return None
    dleft = []
    dright = []
    uleft = []
    uright = []
    pmin = newpoly[newpd.index(min(newpd))]
    pmax = newpoly[newpd.index(max(newpd))]
    new_normal = ncross(list(map(sub, pmax, pmin)), n2)
    for point in newpoly:
        idx = newpoly.index(point)
        if idx == newpoly.index(pmin):
            dleft.append([idx, 0])
            uleft.append([idx, 0])
        elif idx == newpoly.index(pmax):
            dright.append([idx, 0])
            uright.append([idx, 0])
        else:
            inplane_d = ppdist(point, new_normal, pmin)
            if inplane_d < 0:
                if newpd[idx] < 0:
                    dleft.append([idx, inplane_d])
                else:
                    dright.append([idx, inplane_d])
            else:
                if newpd[idx] < 0:
                    uleft.append([idx, inplane_d])
                else:
                    uright.append([idx, inplane_d])
    dlp=newpoly[[x for x in dleft if x[0]==max([x[0] for x in dleft])][0][0]]
    drp=newpoly[[x for x in dright if x[0]==max([x[0] for x in dright])][0][0]]
    ulp=newpoly[[x for x in uleft if x[0] == max([x[0] for x in uleft])][0][0]]
    urp=newpoly[[x for x in uright if x[0]==max([x[0] for x in uright])][0][0]]
    target = line[0][axis-1]
    intone = [0, 0, 0]
    inttwo = [0, 0, 0]
    intone[axis-1]=target
    inttwo[axis-1]=target
    intone[axis-2]=line[0][axis-2]
    inttwo[axis-2]=line[0][axis-2]
    dperc = (target-dlp[axis-1])/(drp[axis-1]-dlp[axis-1])
    uperc = (target-ulp[axis-1])/(urp[axis-1]-ulp[axis-1])
    d_axis = dlp[axis]+dperc*(drp[axis]-dlp[axis])
    u_axis = ulp[axis]+uperc*(urp[axis]-ulp[axis])
    intone[axis]=d_axis
    inttwo[axis]=u_axis
    return [intone, inttwo]


def calc_width_3d(v, res = 5):
    vers = v.generators
    adjc = v.adj
    if len(vers) > 3:
        xrange = (min(vers[:, 0]), max(vers[:, 0]))
        yrange = (min(vers[:, 1]), max(vers[:, 1]))
        zrange = (min(vers[:, 2]), max(vers[:, 2]))
        xpoints = np.linspace(xrange[0], xrange[1], num=res+1, endpoint=False)[1:]
        ypoints = np.linspace(yrange[0], yrange[1], num=res+1, endpoint=False)[1:]
        zpoints = np.linspace(zrange[0], zrange[1], num=res+1, endpoint=False)[1:]
        xyline = []
        xzline = []
        yzline = []
        for xpoint in xpoints:
            for ypoint in ypoints:
                xyline.append(((xpoint, ypoint, zrange[0]), (xpoint, ypoint, zrange[1])))
            for zpoint in zpoints:
                xzline.append(((xpoint, yrange[0], zpoint), (xpoint, yrange[1], zpoint)))
        for ypoint in ypoints:
            for zpoint in zpoints:
                yzline.append(((xrange[0], ypoint, zpoint), (xrange[1], ypoint, zpoint)))
        int_points = {'x': {}, 'y': {}, 'z':{}}
        ipx = {}
        ipy = {}
        ipz = {}
        nfaces = []
        for veridx in range(len(vers)):
            f = []
            adjs = adjc[veridx]
            pos = [list(x)+[veridx] for x in list(combinations(adjs, 2))]
            nfaces += pos
        faces = []
        for face in nfaces:
            a = list(set(face))
            a.sort()
            if a not in faces:
                faces.append(a)
        for l in range(len(xyline)):
            line = xyline[l]
            int_line = []
            for face in faces:
                plane = [vers[x] for x in face]
                ip = lineplane_intersection(line, plane, la=2)
                if ip is not None and len(ip)!=0:
                    int_line.append(ip)
            if len(int_line) > 2:
                int_line.sort(key=lambda x: x[2])
                int_line = [int_line[0], int_line[-1]]
            ipz[l] = int_line
        for l in range(len(xzline)):
            line = xzline[l]
            int_line = []
            for face in faces:
                plane = [vers[x] for x in face]
                ip = lineplane_intersection(line, plane, la=1)
                if ip is not None and len(ip)!=0:
                    int_line.append(ip)
            if len(int_line) > 2:
                int_line.sort(key=lambda x: x[1])
                int_line = [int_line[0], int_line[-1]]
            ipy[l] = int_line
        for l in range(len(yzline)):
            line =yzline[l]
            int_line = []
            for face in faces:
                plane = [vers[x] for x in face]
                ip = lineplane_intersection(line, plane, la=0)
                if ip is not None and len(ip)!=0:
                    int_line.append(ip)
            if len(int_line) > 2:
                int_line.sort(key=lambda x: x[0])
                int_line = [int_line[0], int_line[-1]]
            ipx[l]=int_line
        int_points['x'] = ipx
        int_points['y'] = ipy
        int_points['z'] = ipz
        return int_points
    else:
        return None


with open('to_check_q.json', 'r') as tcq:
     d = json.load(tcq)
#
#
#with open('int_points.json', 'r') as ip:
#    int_points = json.load(ip)

# -1 because failures
#

t = time.time()
failures = []
int_points = {'failures': []}
print(len(d))
for k in range(len(d)):
    print(k)
    material = d[k]
    name = material[0]
    if name not in list(int_points.keys()) and name not in int_points['failures']:
        try:
            mat_dic = material[1]
            self_dic = mat_dic[name]
            elements_of_interest = [x[0] for x in self_dic[0]]
            anion = [x for x in elements_of_interest if x in ['O', 'S', 'Se', 'Te']][0]
            create_output_file(name, mat_dic, anion, 'inputcache.dat')
            v = calc_pd_c('inputcache.dat')
            if len(v.generators)>3:
                a = calc_width_3d(v, res=1)
                int_points[name] = a
            if k%5 == 0:
                print('Writing temp...')
                print('There are {0} entries'.format(len(int_points)-1))
                with open('int_points_ttt.json', 'w') as ip:
                    json.dump(int_points, ip)
                print(time.time()-t)
        except Exception:
            failures.append([name, k])
            int_points['failure'] = failures
    else:
        print(name)
#  for item in rvertices:
#    #print(item)
#      for ad in ad_list:
#            np = vertices[ad]
#            x = [item[0], np[0]]
#            y = [item[1], np[1]]
#            z = [item[2], np[2]]
            #print(vertices[ad])
#            ax.plot(x, y, z, 'b-')
    #print(vars(v))
    #vertices = scipy.array([[0, 0, 0], [0, 1, 0], [1, 0.5, 0], [0.5, 0.5, 1]])
    #if len(vertices)>3:
    #    ax.plot([x[0] for x in vertices], [x[1] for x in vertices], [x[2] for x in vertices], 'k.')
    #    vhull = CHull(vertices, qhull_options="Qc")
    #    for simplex in vhull.simplices:
    #        ax.plot(vertices[simplex, 0], vertices[simplex, 1], vertices[simplex, 2], 'b-')
    #    print(vhull.coplanar)
#plt.show()
read_from_set = 'x'
# Single case
if read_from_set == 0:
    elements_of_interest = ['Cu', 'Sn','S']
    c = get_compounds(elements_of_interest, elements_of_interest[-1])
    create_output_file('CuSn2S3', c, elements_of_interest[-1], 'input_BaTiO3.dat')
    fig = plt.figure()
    file = "input_BaTiO3.dat"
    vertices = calc_pd(file)
    #print(vertices)
    vhull = CHull(vertices)
    #vertsort = []
    for simplex in vhull.simplices:
        plt.plot(vertices[simplex, 0], vertices[simplex, 1], 'k-')
    #vertsort = [x for _,x in sorted(zip(vhull,vertices))]
    #print(vertsort)
    xbound = ((min([x[0] for x in vertices]), max([x[0] for x in vertices])))
    ybound = ((min([x[1] for x in vertices]), max([x[1] for x in vertices])))
    bound = (xbound, ybound)
    #plt.plot([x[0] for x in vertsort], [y[1] for y in vertsort])
    #plt.show()
    ip = calc_width_2d(bound, vertices)
    for item in list(ip['x'].keys()):
        plt.plot([x[0] for x in ip['x'][item]], [x[1] for x in ip['x'][item]], 'r-')
    for item in list(ip['y'].keys()):
        plt.plot([x[0] for x in ip['y'][item]], [x[1] for x in ip['y'][item]], 'b-')
    plt.xlabel('Potential of {0}, eV'.format(elements_of_interest[0]))
    plt.ylabel('Potential of {0}, eV'.format(elements_of_interest[1]))
    plt.show()
    #for m in list(ip['y'].keys()):
    #    if len(ip['y'][m]) != 2:
    #        print(ip['y'][m])
# Read from a file
elif read_from_set == 1:
    start_time = time.time()
    failures = []
    with open('to_check.json') as to_chek:
        to_check = json.load(to_chek)
    vertlim = {}
    results = {}
    for k in range(len(to_check)):
        if k%5 == 0 and k != 0:
            if k%25 == 0:
                print(k)

        material = to_check[k]
        name = material[0]
        mat_dic = material[1]
        self_dic = mat_dic[name]
        elements_of_interest = [x[0] for x in self_dic[0]]
        anion = [x for x in elements_of_interest if x in ['O', 'S', 'Se', 'Te']][0]
        create_output_file(name, mat_dic, anion,'inputcache.dat')
        vertices = calc_pd('inputcache.dat')
        if len(vertices) > 2:
            # Bounding box of the polygon
            # All of the polygon is included in the rectangle
            vhull = CHull(vertices)
            area = vhull.volume
            xbound=((min([x[0] for x in vertices]), max([x[0] for x in vertices])))
            ybound=((min([x[1] for x in vertices]), max([x[1] for x in vertices])))
            bound = (xbound, ybound)
            ip = calc_width_2d(bound, vertices)
            ipalt = calc_width_2d(bound, vertices, res=1)
            results[name] = {}
            results[name]['area'] = area
            results[name]['hi_res'] = ip
            results[name]['low_res'] = ipalt
            results[name]['anion'] = anion
            #if len(vertices) > 2:
            #    area = CHull(vertices).volume
            #    areas.append((name, elements_of_interest, area))
    #with open('output_q.json', 'w') as output:
    #    json.dump(areas, output, indent=4, separators=(',',': '))
    #print(len(areas))
    with open('output_file_tern_c.json','w') as ofj:
        json.dump(results, ofj)
    print(time.time()-start_time)
    print(failures)
elif read_from_set == 'binary':
    pass