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


# Create an output file that follows the format of input.dat
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
                # Takes only the point with the highest z-value and the lowest
                # Same for ipx and ipy but x and y values instead (i.e. x[0] and x[1] sort)
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


# This defines the number of interception lines calculated per axis
# For 3D, the total number of lines is 3R^2
resolution = 5
failures = []
int_points = {'failures': []}
for k in range(len(d)):
    print(len(d))
    if k%5 == 0:
        print(k)
    material = d[k]
    name = material[0]
    # Checks if the material has been used before
    if name not in list(int_points.keys()) and name not in int_points['failures']:
        try:
            mat_dic = material[1]
            self_dic = mat_dic[name]
            elements_of_interest = [x[0] for x in self_dic[0]]
            anion = [x for x in elements_of_interest if x in ['O', 'S', 'Se', 'Te']][0]
            # Generates the vertices
            create_output_file(name, mat_dic, anion, 'inputcache.dat')
            v = calc_pd_c('inputcache.dat')
            # Vertices and adjacents
            if len(v.generators)>3:
                a = calc_width_3d(v, res=resolution)
                # variable a is a dic with keys x y and z
                # Which represents interception lines along that axis
                # Inside each item (e.g. a['x'] there are keys numbered (e.g. a['x'][0])
                # Which are intersection points of a line with the polyhedron
                # There should be only 2 points
                # But because some planes aren't perfect triangles sometimes there are more than 1
                # The easiest way to sort out duplicates is min() and max()
                int_points[name] = a
            if k%5 == 0:
                print('Writing temp...')
                with open('int_points.json', 'w') as ip:
                    json.dump(int_points, ip)
        except Exception:
            failures.append([name, k])
            int_points['failure'] = failures
    else:
        # Indicates duplication in database
        print(name)
# Print all failed ones
print([x[0] for x in failures])


with open('int_points_ttt.json', 'r') as source:
    int_points = json.load(source)

# int_points can then be plotted

