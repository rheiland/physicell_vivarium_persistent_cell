import graphviz
import math

f = graphviz.Digraph('openvt-persistent-cell', filename='rwh.gv', engine='neato')

fs1 = '16'
f.node('Vivarium', shape='rectangle', style='filled', fillcolor='orange', fontsize='24', pos='0,0.!')

all_fw = ['Artistoo','Chaste','CompuCell3D','Morpheus','PhysiCell','Tissue Forge','TST']
R=3
theta = math.pi 
theta_del = math.pi / 6
# theta = math.pi - theta_del 
for fw in all_fw:
    mypos = f'{R*math.cos(theta)},{R*math.sin(theta)}!'
    # print("pos=",mypos)
    f.node(fw, shape='ellipse', style='filled', fillcolor='lightblue', fontsize=fs1, pos=mypos)
    # f.node(fw, shape='ellipse', style='filled', fillcolor='lightblue', fontsize=fs1)
    theta -= theta_del
    f.edge('Vivarium', fw)
    f.edge(fw, 'Vivarium')
# f.node('TST', shape='ellipse', style='filled', fillcolor='lightblue', fontsize=fs1, pos='2.5,1.5!')

f.node('EFECT', shape='ellipse', style='filled', fillcolor='yellow', fontsize='24', pos='0,-1.5!')

f.edge('Vivarium', 'EFECT')
f.edge( 'EFECT', 'Vivarium')

f.view()