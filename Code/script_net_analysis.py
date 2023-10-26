# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 15:38:10 2023

@author: arnau
"""
import pandas as pd
import networkx as nx
from networkx.algorithms import bipartite
import json

import util_net_analysis as util


# Data loading
# Task1
human_df = pd.read_csv('../datasets/HmPtC1_armsrace.csv')
Chimp_df = pd.read_csv('../datasets/HmPtC1_armsrace_pt.csv') 
# Task 2
human_alu_df = pd.read_csv('../datasets/HmC1_Alu.csv') 
human_csv_df = pd.read_csv('../datasets/HmC1_SVA.csv') 
# Task 3 
lost_AD_dfn = pd.read_csv('../datasets/tcx_lostAD_link.csv') 


'''Default script'''
dataframe = lost_AD_dfn
# Constructing the bipartite network and loading the ZNF and TE lists.
net, ZNFnodes, TEnodes = util.bipartite_network_constr(dataframe)
# Computing the node property dataframe.
node_met = util.nodes_metrics(net, ZNFnodes, TEnodes)
# Saving the dataframe.
node_met.to_excel('../results/node_met.xlsx')
# Computing the bipartite density.
dens = bipartite.density(net, ZNFnodes)
print("The bipartite density is : ")
print(dens)
# Ploting node properties distributions.
util.violin_met(node_met)
# Computing the list of hubs based on the threshold.
util.hub_list(node_met,0.05)
# Converting the list of hubs into a dictionnary and saving into an excel file.
util.hub_dict(node_met, 0.05, save = True, save_name = 'lost_ad_hubdf')
# Computing the assortativity analysis
util.assortativity(net, dataframe)
# Importing the community data.
module_dict = json.load(open('../results/node_condor_cluster.json'))
comu_list = list(pd.Series(module_dict['com']).unique())
module_dict = dict(zip( module_dict['id'], module_dict['com']))
# Setting the community node attribute.
nx.set_node_attributes(net, module_dict, name = 'comu' )
#Preparing the assortativity submodule dataframes.
assortativity_df = pd.DataFrame(columns = ['comu', 'assortativity'])
assortativity_df.loc[0, 'comu'] = 'all'
assortativity_df.loc[0, 'assortativity'] = util.assortativity(net, dataframe)
n = 0
# Subseting the largest connected component.
largest_cc = max(nx.connected_components(net), key=len)
largest_cc = net.subgraph(largest_cc)
# Filling the dataframe.
for i in comu_list:
    n +=1
    com_nodes = [x for x,y in largest_cc.nodes(data=True) if y['comu']==i]
    sub_net = largest_cc.subgraph(com_nodes)
    assort = nx.attribute_assortativity_coefficient(sub_net, 'age', nodes=None)
    assortativity_df.loc[n, 'comu'] = i
    assortativity_df.loc[n, 'assortativity'] = assort
# Saving the dataframe.
assortativity_df.to_excel('../results/assort.xlsx')

