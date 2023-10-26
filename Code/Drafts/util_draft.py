# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 14:21:08 2023

@author: arnau
"""


'''Suggestions : 
    
network ploting :
    def net_plot(net):
        for n, d in net.nodes(data=True):
            net.nodes[n]["degree"] = net.degree(n)
        edge_graphic_dict = nx.get_edge_attributes(net,"coef")
        edge_graphic_dict = {k: v/100 for k, v in edge_graphic_dict.items()}
        nx.set_edge_attributes(net, edge_graphic_dict, name="edge_graphic")
        fig, ax = plt.subplots(figsize=(7, 7))
        nv.circos(net, sort_by="degree", group_by="bipartite", node_color_by="bipartite", node_enc_kwargs={"size_scale": 1}, edge_alpha_by="edge_graphic")

Skewness analysis :
    def net_metrics(metricsdf, ZNFnodes, TEnodes):
        
        partition_prop = pd.DataFrame(data=None, index=['Size'], columns=['ZNF', 'TE'])
        partition_prop.loc['Size'] = [len(ZNFnodes),len(TEnodes)]
        partition_propZNF = pd.DataFrame(data = metricsdf.loc[metricsdf['type']=='ZNF'].iloc[:,0:5].apply(stats.skew), columns = ['ZNF'])
        partition_propTe = pd.DataFrame(metricsdf.loc[metricsdf['type']=='TE'].iloc[:,0:5].apply(stats.skew), columns = ['TE'])
        partition_propZNF = partition_propZNF.merge(partition_propTe, left_index =True, right_index=True)
        partition_prop = pd.concat([partition_prop , partition_propZNF])  
        return partition_prop
    
def commu_det(net, resolution = 1):
    commu = nx.community.greedy_modularity_communities(
        net, weight='coef', resolution=resolution)
    n_subgroup = len(commu)
    modu = nx.community.modularity(
        net, communities=commu, weight='coef', resolution=resolution)
    commu_dict = {}
    commu_list = [list(x) for x in commu]
    key = 0
    for i in commu_list :
        commu_dict[key] = i
        key += 1
    print('Modularity is')
    print(modu)
    return commu_dict

def n_submodul(net, resolution = 4):
    commu = nx.community.greedy_modularity_communities(
        net, weight='coef', resolution=resolution)
    n_subgroup = len(commu)
    
    return n_subgroup
def multiresol_modul_plot(net):
    values_df = pd.DataFrame(columns=['resolution', 'number of subgroups'])
    n = 0
    for res in np.arange(0,1,0.1):
        n_sub = n_submodul(net, res)
        values_df.loc[n, 'resolution'] = res
        values_df.loc[n, 'number of subgroups'] = n_sub
        n +=1
    values_df.plot(x = "resolution", y = "number of subgroups")
    return values_df

def save_dict_json(dict_value, name = ''):
    name = '../datasets/'+ name +'.json'
    with open(name, 'w') as fp:
        json.dump(dict_value, fp)
        return None
    
    
def netsave(net, name ='', weight_name = ''):
    adj = nx.to_pandas_adjacency(net, nodelist=None, dtype=None, weight=weight_name)
    adj = adj.abs()
    adj.to_excel('../datasets/'+name +'.xlsx')
    return None
    
    
'''