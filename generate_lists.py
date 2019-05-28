from pymatgen import MPRester
import pymatgen as pmg
m = MPRester("TnKFxKZTef8PsbvDG5C9")


with open('comb.txt') as file:
    lines = file.readlines()


total = '['+''.join(lines).replace('\n',',')+']'
total_list = eval(total)
for k in range(0,10):
    element_set = total_list[k]



