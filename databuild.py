# -*- coding: utf-8 -*-
"""
Created on Mon May 27 16:11:34 2019

@author: ldx
"""
def view_topology(path):
    import networkx as nx
    import matplotlib.pyplot as plt
    topology_grape = nx.Graph()
    topology_file = open(path,'r')
    for topology_line in topology_file.readlines():
        node1,node2 = topology_line.split()
        topology_grape.add_edge(node1,node2)    
    nx.draw(topology_grape,with_labels=True,font_weight='bold')
    #保存绘制图像
    #plt.savefig("topology.png")
    plt.show()
def star_topology_build(count,size = None):
    import random
    if not size:
        size = random.randint(2,20)
    star_topology = open("star_topology"+str(count),'w')
    for i in range(size):
        star_topology.write("SW0 ES"+str(i)+"\n")
    star_topology.flush()
    star_topology.close()
def tree_topology_build(count,size = None):
    import random
    if not size:
        size = random.randint(2,20)
    tree_topology = open("tree_topology"+str(count),'w')
    ES_size = -1
    SW_size = 0
    for i in range(size-1):
        node_is_SW = random.randint(0,1)
        if not node_is_SW:
            ES_size += 1
            tree_topology.write("SW"+str(random.randint(0,SW_size))
            +" ES"+str(ES_size)+"\n")
        else:
            SW_size += 1
            tree_topology.write("SW"+str(random.ranint(0,SW_size-1))
            +" SW"+str(SW_size)+"\n")
    tree_topology.flush()
    tree_topology.close()
def snowflake_topology_build(count,size = None):
    import random
    if not size:
        size = random.randint(2,20)
    snowflake_topology = open("snowflake_topology"+str(count),'w')
    ES_size = -1
    SW_size = 1
    for i in range(size-1):
        node_is_SW = random.randint(0,1)
        if not node_is_SW:
            ES_size += 1
            snowflake_topology.write("SW"+str(random.randint(1,SW_size))
            +" ES"+str(ES_size)+"\n")
        else:
            SW_size += 1
            snowflake_topology.write("SW0 SW"+str(SW_size)+"\n")
    snowflake_topology.flush()
    snowflake_topology.close()       
star_topology_build(1)
view_topology("star_topology1")