# -*- coding: utf-8 -*-
"""
Created on Sat Dec 29 15:08:34 2018

@author: ldx
"""
#!/usr/bin/python
from gurobipy import GRB
from gurobipy import tuplelist
from gurobipy import tupledict
from gurobipy import multidict
from gurobipy import Model
from gurobipy import quicksum
import matplotlib.pyplot as plt     
import networkx as nx
import os
import time
import xlrd

path_of_file = os.path.abspath('')
name_of_file = "\\default_date"
source_type = "xlsx"

if_paint_solution = True
if_print_solution = True
if_print_variable = False
if_view_topology = True
if_write_module = False
if_timing = True
size_of_batch = 100
route_type_set = ["least_leap","fastest_reach"]
route_type = 0
medium_rate = 300
programming_type = GRB.INTEGER
δ = 2.5 

def draw_topology(topology_graph): 
    plt.rcParams['figure.dpi'] = 100                          
    nx.draw(topology_graph,
            pos = nx.spring_layout(topology_graph),
            arrows = True,
            with_labels = True,
            nodelist = topology_graph.nodes(),
            style = 'dashed',
            edge_color = 'b',
            width = 2,
            node_color = 'y',
            alpha = 0.5)
    #保存绘制图像
    plt.savefig("result.png")
    plt.show()   
    
def view_topology(path_of_file):
    topology_graph = nx.Graph()
    with xlrd.open_workbook(os.path.join(path_of_file)) as topology_book:
        topology_sheet = topology_book.sheet_by_name("E")
        i = 1
        while True:
            try:
                topology_graph.add_edge(topology_sheet.cell_value(i, 0),
                                        topology_sheet.cell_value(i, 1))
                i= i + 1
            except IndexError:
                break
    draw_topology(topology_graph)   

def lowest_common_multiple(x, y):
    if x > y:
        greater = x
    else:
        greater = y
    while(True):
        if((greater % x == 0) and (greater % y == 0)):
            lowest_common_multiple = greater
            break
        greater += 1
    return lowest_common_multiple 
    
def print_variable():
    for v in md.getVars():
        print('%s %g' % (v.varName, v.x))
            
def print_solution():
    if md.status == GRB.Status.OPTIMAL:
        print('\nCost: %g' % md.objVal)
        print('\noffset:')
        for m in f:
            print('%s %g' % (m, φ[m].x))
        print('\nchannel:')
        for k in F:
            print('%s %g' % (k, ρ[k].x))
    else:
        print('No solution')
        
def paint_solution():
    from matplotlib import colors
    color_index = list(colors.CSS4_COLORS.keys())
    lowest_T = 1
    for m in old_f:
        lowest_T = lowest_common_multiple(old_T[m],lowest_T)
    fig, ax = plt.subplots(1)
    i = 0
    for va,vb in E:
        target_f = old_f.select(va,vb,'*','*')
        for m in target_f:
            mark_m = int(m[2]) % len(color_index)
            for α in range(round(lowest_T/old_T[m])):
                rect = plt.Rectangle((old_φ[m] + α*old_T[m],i-0.5),
                                     old_e2e[m],
                                     1,
                                     color = color_index[mark_m])
                ax.add_patch(rect)
                ax.annotate(m[2]+'.'+str(round(old_ρ[m[0],m[1],m[2]])), 
                            (old_φ[m] + α*old_T[m] + old_e2e[m]/2, i), 
                            color='black', 
                            weight='bold', 
                            fontsize=12, 
                            ha='center', 
                            va='center')
        i += 1 
    plt.xlim((0, lowest_T))    
    plt.ylim((-0.5, len(E)-0.5)) 
    plt.yticks(range(0, len(E)), E)
    plt.show()  

if source_type == "xlsx":         
    source_book = xlrd.open_workbook(os.path.join('cyber.xlsx'))
    if if_view_topology:
        view_topology("cyber.xlsx")
    #边集，存储所有双向边
    source_sheet = source_book.sheet_by_name("E")
    E = {}
    i = 1
    while True:
        try:
            start_node = source_sheet.cell_value(i, 0)
            end_node = source_sheet.cell_value(i, 1)
            arc_rate = source_sheet.cell_value(i, 2)
            arc_length = source_sheet.cell_value(i, 3)
            E[(start_node,end_node)] = [arc_rate,arc_length]
            E[(end_node,start_node)] = [arc_rate,arc_length]
            i= i + 1
        except IndexError:
            break  

    #流集，存储流与流属性        
    source_sheet = source_book.sheet_by_name("S")
    full_S = {}
    i = 1
    while True:
        try:
            task_id = str(int(source_sheet.cell_value(i, 0)))
            task_period = source_sheet.cell_value(i, 1)
            task_deadline = source_sheet.cell_value(i, 2)
            task_length = source_sheet.cell_value(i, 3)
            full_S[task_id] = [task_period,task_deadline,task_length]
            i = i + 1
        except IndexError:
            break
           
    #有序存储流实例
    source_sheet = source_book.sheet_by_name("F")
    full_F = tuplelist([])
    i = 1
    while True:
        try:
            task_id = str(int(source_sheet.cell_value(i, 0)))
            start_node = source_sheet.cell_value(i, 1)
            end_node = source_sheet.cell_value(i, 2)
            full_F.append((start_node, 
                           end_node,
                           task_id))
            i = i + 1
        except IndexError:
            break 
elif source_type == "txt":
    with open(path_of_file+name_of_file,'r') as input_file:
        input_content = input_file.readlines()
        num_of_node = int(input_content[0])
        num_of_arc = int(input_content[1])
        num_of_task = int(input_content[2]) 
        
        topology_graph=nx.Graph()
        E = {}
        for i in range(num_of_arc):
            arc_i = input_content[3+i].split()
            start_node = arc_i[0]
            end_node = arc_i[1]
            arc_rate = int(arc_i[2])
            arc_length = int(arc_i[3])
            E[(start_node,end_node)] = [arc_rate,arc_length]
            E[(end_node,start_node)] = [arc_rate,arc_length]
            topology_graph.add_edge(start_node,end_node)
        if if_view_topology:
            draw_topology(topology_graph)
        full_S = {}
        full_F = tuplelist([])
        for i in range(num_of_task):
            task_i = input_content[3+num_of_arc+i].split()
            start_node = task_i[0]
            end_node = task_i[1]
            task_period = float(task_i[2])
            task_deadline = float(task_i[3])
            task_length = int(task_i[4])
            full_S[str(i)] = [task_period,task_deadline,task_length]
            route = nx.shortest_path(topology_graph,
                                     source=start_node,
                                     target=end_node)
            for t in range(len(route)-1):
                full_F.append((route[t],route[t+1],str(i)))
        
E,R,LE = multidict(E)
old_κ = tupledict({})
for va,vb in E:
    old_κ[va,vb] = -1
old_φ = tupledict({})
old_ρ = tupledict({})
old_F = tuplelist([])
old_T = tupledict({})
old_e2e = tupledict({})
old_f = tuplelist([])
    
batch = 0
while len(full_S) > batch*size_of_batch:    
    #有序存储帧实例及属性
    S = dict(list(full_S.items())[size_of_batch*batch:size_of_batch*(batch+1)])
    S,task_period,D,task_size = multidict(S)
    F = tuplelist([])
    for s in S:
        F += full_F.select('*','*',s)
    f = {}
    for k_start,k_end,k_id in F:
        i = 0
        size_k = int(task_size[k_id])
        while size_k > 1500:
            f[(k_start,k_end,k_id,i)] = [task_period[k_id],1500]
            size_k -= 1500
            i += 1
        f[(k_start,k_end,k_id,i)] = [task_period[k_id],1500]
    f,T,L =  multidict(f)
    e2e = tupledict({})
    for m in f:
        m_start,m_end,m_id,m_i = m
        e2e[m] = L[m]*8/R[m_start,m_end] + LE[m_start,m_end]/medium_rate
    md = Model('cyber')
    
    φ = md.addVars(f, name='phi',vtype = programming_type)
    ρ = md.addVars(F, name='rho',vtype = programming_type)
    κ = md.addVars(E, name='chi',vtype = programming_type)
    
    md.setObjective(quicksum(κ[va,vb] for va,vb in E), GRB.MINIMIZE)
    
    for va,vb in E:
        for k in F.select(va,vb,'*'):
            md.addConstr((ρ[k] <= κ[va,vb]), 'chi Constaints')
        md.addConstr((old_κ[va,vb] <= κ[va,vb]), 'oldchi Constaints')
            
    for m in f:
        md.addConstr((φ[m] <= T[m] - e2e[m]), 'Frame Constaints')
    
    for va,vb in E:
        target_f = f.select(va,vb,'*','*')
        for m in target_f[:-1]:
            for n in target_f[target_f.index(m)+1:]:
                for α in range(round(lowest_common_multiple(T[m],T[n])/T[m])):
                    for β in range(round(lowest_common_multiple(T[m],T[n])/T[n])):
                        if β*T[n] > (α+1)*T[m]:
                            break
                        if β*T[n] < (α-1)*T[m]:
                            continue
                        σ = md.addVar()
                        md.addGenConstrIndicator(σ, 
                                                 True, 
                                                 φ[m] - φ[n], 
                                                 GRB.LESS_EQUAL,
                                                 - e2e[m] - α*T[m] + β*T[n], 
                                                 'Link Constaints 1') 
                        md.addGenConstrIndicator(σ, 
                                                 False, 
                                                 φ[n] - φ[m], 
                                                 GRB.LESS_EQUAL,
                                                  - e2e[n] - β*T[n] + α*T[m], 
                                                 'Link Constaints 2') 
        for m in target_f:
            for n in old_f.select(va,vb,'*','*'):
                for α in range(round(lowest_common_multiple(T[m],old_T[n])/T[m])):
                    for β in range(round(lowest_common_multiple(T[m],old_T[n])/old_T[n])):
                        if β*old_T[n] > (α+1)*T[m]:
                            break
                        if β*old_T[n] < (α-1)*T[m]:
                            continue
                        σ = md.addVar()
                        md.addGenConstrIndicator(σ, 
                                                 True,
                                                 φ[m],
                                                 GRB.LESS_EQUAL, 
                                                 - e2e[m] - α*T[m] + β*old_T[n] + old_φ[n], 
                                                 'Link Constaints 1') 
                        md.addGenConstrIndicator(σ, 
                                                 False,
                                                 -φ[m], 
                                                 GRB.LESS_EQUAL, 
                                                 - old_e2e[n] - β*old_T[n] - old_φ[n] + α*T[m], 
                                                 'Link Constaints 2')
    
    
    for m in f:
        for n in f.select(m[1],'*',m[2],m[3]):
            md.addConstr((φ[m] + e2e[m] + δ <= φ[n]), 'Flow Transmission Constaints')
    for s in S:
        target_f = f.select('*','*',s,'*')
        md.addConstr((φ[target_f[-1]] + e2e[target_f[-1]] - φ[target_f[0]] <= D[s]), 'End-to-End Constaints')
                                
    for va,vb in E:
        target_F = F.select(va,vb,'*')
        for k in target_F[:-1]:
            for l in target_F[target_F.index(k)+1:]: 
                εlk = md.addVar()
                md.addGenConstrIndicator(εlk,
                                         True,
                                         ρ[l] - ρ[k], 
                                         GRB.GREATER_EQUAL, 
                                         1, 
                                         'rho lk Constaints 1') 
                md.addGenConstrIndicator(εlk,
                                         False, 
                                         ρ[l] - ρ[k], 
                                         GRB.LESS_EQUAL, 
                                         0, 
                                         'rho lk Constaints 2')
                εkl = md.addVar()
                md.addGenConstrIndicator(εkl,
                                         True,
                                         ρ[k] - ρ[l], 
                                         GRB.GREATER_EQUAL,
                                         1, 
                                         'rho kl Constaints 1')
                md.addGenConstrIndicator(εlk,
                                         False, 
                                         ρ[k] - ρ[l], 
                                         GRB.LESS_EQUAL, 
                                         0, 
                                         'rho kl Constaints 2')
                for m in f.select(k[0],k[1],k[2],'*'):
                    for n in f.select(l[0],l[1],l[2],'*'):
                        for α in range(round(lowest_common_multiple(T[m],T[n])/T[m])):
                            for β in range(round(lowest_common_multiple(T[m],T[n])/T[n])):
                                if β*T[n] > (α+1)*T[m]:
                                    break
                                if β*T[n] < (α-1)*T[m]:
                                    continue                        
                                ω = md.addVar()
                                ωε = md.addVar()
                                md.addGenConstrIndicator(ωε,
                                                         True, 
                                                         ω+εlk+εkl, 
                                                         GRB.GREATER_EQUAL, 
                                                         1, 
                                                         'support Constaints 1')
                                εω = md.addVar()
                                md.addGenConstrIndicator(εω,
                                                         True, 
                                                         -ω+εlk+εkl, 
                                                         GRB.GREATER_EQUAL, 
                                                         0, 
                                                         'support Constaints 2')
                                for mp in f.select('*',m[0],m[2],m[3]):
                                    md.addGenConstrIndicator(ωε,
                                                             False, 
                                                             φ[n] - φ[mp], 
                                                             GRB.LESS_EQUAL, 
                                                             - δ - β*T[n] + α*T[m], 
                                                             'Frame Isolation Constraints 1')
                                for np in f.select('*',n[0],n[2],n[3]):
                                    md.addGenConstrIndicator(εω,
                                                             False, 
                                                             φ[m] - φ[np], 
                                                             GRB.LESS_EQUAL, 
                                                             - δ - α*T[m] + β*T[n], 
                                                             'Frame Isolation Constraints 2')
        for k in target_F:
            for l in old_F.select(va,vb,'*'): 
                εlk = md.addVar()
                md.addGenConstrIndicator(εlk,
                                         True,  
                                         - ρ[k], 
                                         GRB.GREATER_EQUAL,
                                         1-old_ρ[l], 
                                         'rho lk Constaints 1') 
                md.addGenConstrIndicator(εlk,
                                         False,  
                                         - ρ[k], 
                                         GRB.LESS_EQUAL, 
                                         -old_ρ[l], 
                                         'rho lk Constaints 2')
                εkl = md.addVar()
                md.addGenConstrIndicator(εkl,
                                         True, 
                                         ρ[k], 
                                         GRB.GREATER_EQUAL, 
                                         1+old_ρ[l], 
                                         'rho kl Constaints 1')
                md.addGenConstrIndicator(εlk,
                                         False, 
                                         ρ[k], 
                                         GRB.LESS_EQUAL, 
                                         old_ρ[l], 
                                         'rho kl Constaints 2')
                for m in f.select(k[0],k[1],k[2],'*'):
                    for n in old_f.select(l[0],l[1],l[2],'*'):
                        for α in range(round(lowest_common_multiple(T[m],old_T[n])/T[m])):
                            for β in range(round(lowest_common_multiple(T[m],old_T[n])/old_T[n])):
                                if β*old_T[n] > (α+1)*T[m]:
                                    break
                                if β*old_T[n] < (α-1)*T[m]:
                                    continue                        
                                ω = md.addVar()
                                ωε = md.addVar()
                                md.addGenConstrIndicator(ωε,
                                                         True, 
                                                         ω+εlk+εkl, 
                                                         GRB.GREATER_EQUAL, 
                                                         1, 
                                                         'support Constaints 1')
                                εω = md.addVar()
                                md.addGenConstrIndicator(εω,
                                                         True, 
                                                         -ω+εlk+εkl, 
                                                         GRB.GREATER_EQUAL, 
                                                         0, 
                                                         'support Constaints 2')
                                for mp in f.select('*',m[0],m[2],m[3]):
                                    md.addGenConstrIndicator(ωε,
                                                             False, 
                                                             - φ[mp], 
                                                             GRB.LESS_EQUAL, 
                                                             - δ - β*old_T[n] + α*T[m] - old_φ[n], 
                                                             'Frame Isolation Constraints 1')
                                for np in old_f.select('*',n[0],n[2],n[3]):
                                    md.addGenConstrIndicator(εω,
                                                             False, 
                                                             φ[m], 
                                                             GRB.LESS_EQUAL, 
                                                             - δ - α*T[m] + β*old_T[n] + old_φ[np], 
                                                             'Frame Isolation Constraints 2')
    if if_timing:
        start_point = time.clock()                                
        md.optimize()
        elapsed = (time.clock() - start_point)
        print("Time used: ",elapsed)
    else:
        md.optimize()
    if md.status != GRB.Status.OPTIMAL:
        print("no solution")
        break
    for va,vb in E:
        old_κ[va,vb] = κ[va,vb].x 
    for m in f:
        old_φ[m]=φ[m].x
    for k in F:
        old_ρ[k]=ρ[k].x
        
    old_f += f
    old_T.update(T)   
    old_e2e.update(e2e)
    old_F += F
    batch += 1
    if print_variable:
        print_variable()
    if print_solution:
        print_solution()
if if_paint_solution:
    paint_solution()
if if_write_module:
    md.write('cyber.lp')