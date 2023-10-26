# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 14:20:52 2023

@author: arnau
"""

''' Questions : Full network with all the possible nodes for density?

Suggestions : 
    
Network plotting :
net_plot(net)

Skewness analysis : 
net_metrics(node_met, ZNFnodes, TEnodes)


#problème de community detection

# modularity
default_comu_dict = comu_dict_big_cluster = commu_det(net, 1)
comu_dict_big_cluster = commu_det(net, 0.2)

save_dict_json(default_comu_dict, name = 'default_comu')
save_dict_json(comu_dict_big_cluster, name = 'big_cluster')

netsave(net, name = 'lost_ad_gephi', weight_name =  'coef')

#fonctiion dataframe:
node_met.loc[['SVA_E','SVA_F']]
#terminer les foncitons et ensuite faire toutes les analyses de a à z. 
proj_net = bipartite.weighted_projected_graph(net, TEnodes, ratio=False)


te_node_attr = dataframe.loc[:, ['teName','repFamily']]
te_node_attr.set_index('teName', inplace =True)
te_node_attr = te_node_attr.to_dict()
te_node_attr = te_node_attr['repFamily']
nx.set_node_attributes(proj_net, te_node_attr, name = 'repFamily' )
nx.attribute_assortativity_coefficient(proj_net, 'repFamily', nodes=None)
nx.get_node_attributes(proj_net, 'repFamily')

save_dict_json(te_node_attr, name = 'te_node_attr')

netsave(net, name = 'proj_net_gephi', weight_name =  'weight')

test = multiresol_modul_plot(net)


'''