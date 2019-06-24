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

batch = 100
δ = 2.5 
source = "txt"
programming_type = GRB.CONTINUOUS

tur = 0
fpath = os.path.abspath('..')

def draw(G): 
    plt.rcParams['figure.dpi'] = 100                          
    pos=nx.spring_layout(G)
    nx.draw(G,pos,arrows=True,with_labels=True,nodelist=G.nodes(),style='dashed',
            edge_color='b',width=2,node_color='y',alpha=0.5)
    #保存绘制图像
    plt.savefig("result.png")
    plt.show()   
    
def viewe(path):
    G = nx.Graph()
    book = xlrd.open_workbook(os.path.join(path))
    sh = book.sheet_by_name("E")
    i = 1
    while True:
        try:
            G.add_edge(sh.cell_value(i, 0),sh.cell_value(i, 1))
            i= i + 1
        except IndexError:
            break
    draw(G)   

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
    
def printSol():
    for v in md.getVars():
        print('%s %g' % (v.varName, v.x))
            
def printSolution():
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
        
def painte():
    from matplotlib import colors
    colorindex = list(colors.CSS4_COLORS.keys())
    lcmT = 1
    for m in oldf:
        lcmT = lcm(oldT[m],lcmT)
    fig, ax = plt.subplots(1)
    i=0
    for va,vb in E:
        fs = oldf.select(va,vb,'*','*')
        for m in fs:
            markm = int(m[2]) % len(colorindex)
            for α in range(round(lcmT/oldT[m])):
                rect = plt.Rectangle((oldφ[m] + α*oldT[m],i-0.5),oldL[m],1,
                                     color = colorindex[markm])
                ax.add_patch(rect)
                ax.annotate(m[2]+'.'+str(round(oldρ[m[0],m[1],m[2]])), 
                            (oldφ[m] + α*oldT[m] + oldL[m]/2, i), color='black', 
                            weight='bold', fontsize=12, ha='center', va='center')
        i += 1 
    plt.xlim((0, lcmT))    
    plt.ylim((-0.5, len(E)-0.5)) 
    plt.yticks(range(0, len(E)), E)
    plt.show()  

if source == "xlsx":         
    book = xlrd.open_workbook(os.path.join('cyber.xlsx'))
    viewe("cyber.xlsx")
    #边集，存储所有双向边
    sh = book.sheet_by_name("E")
    E = tuplelist([])
    i = 1
    while True:
        try:
            E.append((str(int(sh.cell_value(i, 0))),str(int(sh.cell_value(i, 1)))))
            E.append((str(int(sh.cell_value(i, 1))),str(int(sh.cell_value(i, 0)))))
            i= i + 1
        except IndexError:
            break  
    
    #流集，存储流与流属性        
    sh = book.sheet_by_name("S")
    fullS = {}
    i = 1
    while True:
        try:    
            fullS[str(int(sh.cell_value(i, 0)))] = [sh.cell_value(i, 1),sh.cell_value(i, 2)]
            i = i + 1
        except IndexError:
            break
           
    #有序存储流实例
    sh = book.sheet_by_name("F")
    fullF = tuplelist([])
    i = 1
    while True:
        try:
            fullF.append((str(int(sh.cell_value(i, 1))),str(int(sh.cell_value(i, 2))),str(int(sh.cell_value(i, 0)))))
            i = i + 1
        except IndexError:
            break 
elif source == "txt":
    with open(fpath+"\\grid\\grid7",'r') as mnet:
        grid = mnet.readlines()
        numofnode = int(grid[0])
        numofarc = int(grid[1])
        numoftask = int(grid[2]) 
        
        G=nx.Graph()    
        E = tuplelist([])
        superT = 0
        for i in range(numofarc):
            arci = grid[3+i].split()
            E.append((arci[0],arci[1]))
            G.add_edge(arci[0],arci[1])
            superT = max(superT,int(arci[3]))
        draw(G)
        fullS = {}
        fullF = tuplelist([])
        for i in range(numoftask):
            taski = grid[3+numofarc+i].split()
            fullS[str(i)] = [superT,int(taski[2])]
            route = nx.shortest_path(G,source=taski[0],target=taski[1])
            for t in range(len(route)-1):
                fullF.append((route[t],route[t+1],str(i)))
        
oldκ = tupledict({})
for va,vb in E:
    oldκ[va,vb] = -1
oldφ = tupledict({})
oldρ = tupledict({})
oldF = tuplelist([])
oldT = tupledict({})
oldL = tupledict({})
oldf = tuplelist([])
    
while len(fullS) > tur*batch: 
    
    #有序存储帧实例及属性
    S = dict(list(fullS.items())[batch*tur:batch*(tur+1)])
    S,e2e,size = multidict(S)
    F = tuplelist([])
    for s in S:
        F += fullF.select('*','*',s)
    f = {}
    for k in F:
        i = 0
        sizek = size[k[2]]
        while sizek > 15:
            f[(k[0],k[1],k[2],i)] = [e2e[k[2]],15.0]
            sizek -= 15
            i += 1
        f[(k[0],k[1],k[2],i)] = [e2e[k[2]],sizek]
    f,T,L =  multidict(f)
     
    md = Model('cyber')
    
    φ = md.addVars(f, name='phi',vtype = programming_type)
    ρ = md.addVars(F, name='rho',vtype = programming_type)
    κ = md.addVars(E, name='chi',vtype = programming_type)
    
    md.setObjective(quicksum(κ[va,vb] for va,vb in E), GRB.MINIMIZE)
    
    for va,vb in E:
        for k in F.select(va,vb,'*'):
            md.addConstr((ρ[k] <= κ[va,vb]), 'chi Constaints')
        md.addConstr((oldκ[va,vb] <= κ[va,vb]), 'oldchi Constaints')
            
    for m in f:
        md.addConstr((φ[m] <= T[m] - L[m]), 'Frame Constaints')
    
    for va,vb in E:
        fs = f.select(va,vb,'*','*')
        for m in fs[:-1]:
            for n in fs[fs.index(m)+1:]:
                for α in range(round(lcm(T[m],T[n])/T[m])):
                    for β in range(round(lcm(T[m],T[n])/T[n])):
                        if β*T[n] > (α+1)*T[m]:
                            break
                        if β*T[n] < (α-1)*T[m]:
                            continue
                        σ = md.addVar()
                        md.addGenConstrIndicator(σ, True, φ[m] - φ[n]  , 
                                                 GRB.LESS_EQUAL , - L[m] - α*T[m] + β*T[n], 
                                                 'Link Constaints 1') 
                        md.addGenConstrIndicator(σ, False, φ[n] - φ[m]  , 
                                                 GRB.LESS_EQUAL , - L[n] - β*T[n] + α*T[m], 
                                                 'Link Constaints 2') 
        for m in fs:
            for n in oldf.select(va,vb,'*','*'):
                for α in range(round(lcm(T[m],oldT[n])/T[m])):
                    for β in range(round(lcm(T[m],oldT[n])/oldT[n])):
                        if β*oldT[n] > (α+1)*T[m]:
                            break
                        if β*oldT[n] < (α-1)*T[m]:
                            continue
                        σ = md.addVar()
                        md.addGenConstrIndicator(σ, True, φ[m], GRB.LESS_EQUAL, 
                                                 - L[m] - α*T[m] + β*oldT[n] + oldφ[n], 
                                                 'Link Constaints 1') 
                        md.addGenConstrIndicator(σ, False, -φ[m], GRB.LESS_EQUAL, 
                                                 - oldL[n] - β*oldT[n] - oldφ[n] + α*T[m], 
                                                 'Link Constaints 2')
    
    
    for m in f:
        for n in f.select(m[1],'*',m[2],m[3]):
            md.addConstr((φ[m] + L[m] + δ <= φ[n]), 'Flow Transmission Constaints')
    for s in S:
        fs = f.select('*','*',s,'*')
        print(fs)
        md.addConstr((φ[fs[-1]] + L[fs[-1]] - φ[fs[0]] <= e2e[s]), 'End-to-End Constaints')
                                
    for va,vb in E:
        Fs = F.select(va,vb,'*')
        for k in Fs[:-1]:
            for l in Fs[Fs.index(k)+1:]: 
                εlk = md.addVar()
                md.addGenConstrIndicator(εlk, True, ρ[l] - ρ[k] , 
                                         GRB.GREATER_EQUAL , 1, 'rho lk Constaints 1') 
                md.addGenConstrIndicator(εlk, False, ρ[l] - ρ[k] , 
                                         GRB.LESS_EQUAL , 0, 'rho lk Constaints 2')
                εkl = md.addVar()
                md.addGenConstrIndicator(εkl, True, ρ[k] - ρ[l] , 
                                         GRB.GREATER_EQUAL , 1, 'rho kl Constaints 1')
                md.addGenConstrIndicator(εlk, False, ρ[k] - ρ[l] , 
                                         GRB.LESS_EQUAL , 0, 'rho kl Constaints 2')
                for m in f.select(k[0],k[1],k[2],'*'):
                    for n in f.select(l[0],l[1],l[2],'*'):
                        for α in range(round(lcm(T[m],T[n])/T[m])):
                            for β in range(round(lcm(T[m],T[n])/T[n])):
                                if β*T[n] > (α+1)*T[m]:
                                    break
                                if β*T[n] < (α-1)*T[m]:
                                    continue                        
                                ω = md.addVar()
                                ωε = md.addVar()
                                md.addGenConstrIndicator(ωε, True, ω+εlk+εkl , 
                                                         GRB.GREATER_EQUAL , 1, 
                                                         'support Constaints 1')
                                εω = md.addVar()
                                md.addGenConstrIndicator(εω, True, -ω+εlk+εkl , 
                                                         GRB.GREATER_EQUAL , 0, 
                                                         'support Constaints 2')
                                for mp in f.select('*',m[0],m[2],m[3]):
                                    md.addGenConstrIndicator(ωε, False, 
                                                             φ[n] - φ[mp], 
                                                             GRB.LESS_EQUAL, 
                                                             - δ - β*T[n] + α*T[m], 
                                                             'Frame Isolation Constraints 1')
                                for np in f.select('*',n[0],n[2],n[3]):
                                    md.addGenConstrIndicator(εω, False, 
                                                             φ[m] - φ[np], 
                                                             GRB.LESS_EQUAL, 
                                                             - δ - α*T[m] + β*T[n], 
                                                             'Frame Isolation Constraints 2')
        for k in Fs:
            for l in oldF.select(va,vb,'*'): 
                εlk = md.addVar()
                md.addGenConstrIndicator(εlk, True,  - ρ[k] , GRB.GREATER_EQUAL ,
                                         1-oldρ[l], 'rho lk Constaints 1') 
                md.addGenConstrIndicator(εlk, False,  - ρ[k] , GRB.LESS_EQUAL , 
                                         -oldρ[l], 'rho lk Constaints 2')
                εkl = md.addVar()
                md.addGenConstrIndicator(εkl, True, ρ[k]  , GRB.GREATER_EQUAL , 
                                         1+oldρ[l], 'rho kl Constaints 1')
                md.addGenConstrIndicator(εlk, False, ρ[k]  , GRB.LESS_EQUAL , 
                                         oldρ[l], 'rho kl Constaints 2')
                for m in f.select(k[0],k[1],k[2],'*'):
                    for n in oldf.select(l[0],l[1],l[2],'*'):
                        for α in range(round(lcm(T[m],oldT[n])/T[m])):
                            for β in range(round(lcm(T[m],oldT[n])/oldT[n])):
                                if β*oldT[n] > (α+1)*T[m]:
                                    break
                                if β*oldT[n] < (α-1)*T[m]:
                                    continue                        
                                ω = md.addVar()
                                ωε = md.addVar()
                                md.addGenConstrIndicator(ωε, True, ω+εlk+εkl , 
                                                         GRB.GREATER_EQUAL , 1, 
                                                         'support Constaints 1')
                                εω = md.addVar()
                                md.addGenConstrIndicator(εω, True, -ω+εlk+εkl , 
                                                         GRB.GREATER_EQUAL , 0, 
                                                         'support Constaints 2')
                                for mp in f.select('*',m[0],m[2],m[3]):
                                    md.addGenConstrIndicator(ωε, False, - φ[mp], 
                                                             GRB.LESS_EQUAL, 
                                                             - δ - β*oldT[n] + α*T[m] - oldφ[n], 
                                                             'Frame Isolation Constraints 1')
                                for np in oldf.select('*',n[0],n[2],n[3]):
                                    md.addGenConstrIndicator(εω, False, φ[m], 
                                                             GRB.LESS_EQUAL, 
                                                             - δ - α*T[m] + β*oldT[n] + oldφ[np], 
                                                             'Frame Isolation Constraints 2')
    start_point = time.clock()                                
    md.optimize()
    elapsed = (time.clock() - start_point)
    print("Time used: ",elapsed)
    if md.status != GRB.Status.OPTIMAL:
        print("no solution")
        break
    for va,vb in E:
        oldκ[va,vb] = κ[va,vb].x 
    for m in f:
        oldφ[m]=φ[m].x
    for k in F:
        oldρ[k]=ρ[k].x
        
    oldf += f
    oldT.update(T)   
    oldL.update(L)
    oldF += F
    tur += 1
    #printSolution()
#painte()
#md.write('cyber.lp')