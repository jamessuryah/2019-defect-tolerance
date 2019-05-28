import json


with open('int_points_ttt.json') as ttt:
    t = json.load(ttt)


with open('output_q_fixed.json') as oqf:
    q = json.load(oqf)


nt = t
for name in list(t.keys()):
    t = [x for x in q if x[0]==name]
    if len(t) != 1:
        print(name)
    else:
        anion = t[0][1][-1]
        nt[name]['anion'] = anion

with open('int_points_ttt_fixed.json', 'w') as file:
    json.dump(nt, file)