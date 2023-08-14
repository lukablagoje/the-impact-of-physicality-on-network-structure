import pickle
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from os import listdir
from os.path import isfile, join
from itertools import combinations
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['axes.titlesize'] = 16
import itertools
from neuprint import fetch_synapses, NeuronCriteria as NC, SynapseCriteria as SC
from neuprint import Client
from neuprint import fetch_adjacencies 
from neuprint import fetch_neurons
from neuprint import merge_neuron_properties
from neuprint import fetch_synapse_connections
from neuprint import Client
TOKEN = "eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJlbWFpbCI6Imx1a2FibGFnb2pldmljMTk5NUBnbWFpbC5jb20iLCJsZXZlbCI6Im5vYXV0aCIsImltYWdlLXVybCI6Imh0dHBzOi8vbGgzLmdvb2dsZXVzZXJjb250ZW50LmNvbS9hLS9BT2gxNEdqdDZpdFFGR2xTSWZUTElNTjRmcEt1QzZ3QmE2Rlp0WU1XYmpKV1ZBPXM5Ni1jP3N6PTUwP3N6PTUwIiwiZXhwIjoxODA5MTE1OTEyfQ.0h6CJp8xfQEpkW8a2_gqJUBrEA5GyBiZkNvDjRpoXoY" # <--- Paste your token here
           # (or define NEUPRINT_APPLICATION CREDENTIALS in your environment)

    
c = Client('neuprint.janelia.org', 'hemibrain:v1.2.1', TOKEN)

step  = 0.01
end = 1.00
radius_list = list(np.arange(0,end  ,step))
radius_list_below = [np.round(i,2) for i in radius_list]



radius_list = radius_list_below  + [1, 2, 3, 4, 5, 10, 15, 20, 25, 30, 40, 45, 50, 60, 70, 80, 90, 100, 150,200]
metagraph_radius_spherical_dict = {}
connectome_radius_spherical_dict = {}

number_of_neurons_spherical_dict = {}
number_of_metagraph_links_spherical_dict = {}
number_of_connectome_links_spherical_dict = {}
metagraph_links_spherical_dict = {}
connectome_links_spherical_dict = {}
infile = open('ME(R)_edge_dict_spherical_'+str(200)+'.pkl','rb')
metagraph_dictionary = pickle.load(infile)

all_neurons = []
for key in metagraph_dictionary.keys():
    for value in metagraph_dictionary[key]:
        all_neurons.append(key)
        all_neurons.append(value)
    all_neurons = list(set(all_neurons))    
    
criteria = NC(bodyId=all_neurons)
neuron_df,connectome_edgelist = fetch_adjacencies(criteria,criteria,properties=['size'])
connectome_pair_list = connectome_edgelist[['bodyId_pre','bodyId_post']].values.tolist()
for radius in radius_list:
    infile = open('ME(R)_edge_dict_spherical_'+str(radius)+'.pkl','rb')
    
    metagraph_dictionary = pickle.load(infile)
    
    ##ADD ALL NEURONS AT ONCE, CONNECTOME SHOULD BE FIXED, GRAPH SHOULD HAVE ISOLATED 
    #FILTER IF THEY ARE IN CONNECTOME, DON"T USE THEM USE THEM FOR METAGRAPH
    #CONNECTOME BETWEENEES VS METAGRAPH DEGREE, in csv file, or text so you can add it to plot
    #PART OF THE APPENDIX -HOW THE NEURONS ARE CHOSEN, HOW IS THE GRAPH CONSTRUCTED,  WHAT IS THE APPROXIMATION METHOD,
    #PICK DIFFERENT BRAIN REGIONS not only M(E)R, just to check 
    # CHECK for whole skeleton
    #Write about end to end distances and about radius lengths, sanity check for approximations
    
    #SHOULD CONNECTOME AND METAGRAPH LINKS SHOULD COUNT IN FOR THE NETWORK OR NOT?
    metagraph_pair_list = []
    metagraph_neurons = []
    for key in metagraph_dictionary.keys():
        for value in metagraph_dictionary[key]:
            if (([key,value] in connectome_pair_list ) or ([value,key] in connectome_pair_list )):
                if not([key, value] in metagraph_pair_list) and not([value,key] in metagraph_pair_list):
                    metagraph_pair_list.append([key,value])
                    metagraph_neurons.append(key)
                    metagraph_neurons.append(value)
    metagraph_radius_spherical_dict[radius] = metagraph_pair_list
    connectome_radius_spherical_dict[radius] =  connectome_pair_list
    number_of_metagraph_links_spherical_dict[radius] = len(metagraph_pair_list)

    
    print('Radius',radius)
    print("Number of neurons",len(metagraph_neurons))
    print("Number of links in the metagraph: " + str(len(metagraph_pair_list)))
    print("Number of links in the connectome: " + str(len(connectome_pair_list )))
    
metagraph_links_spherical_dictionary =  metagraph_radius_spherical_dict
connectome_links_spherical_dictionary    =  connectome_radius_spherical_dict


with open('metagraph_links_spherical_dictionary.pickle', 'wb') as handle:
    pickle.dump(metagraph_links_spherical_dictionary , handle, protocol=pickle.HIGHEST_PROTOCOL)
with open('connectome_links_spherical_dictionary.pickle', 'wb') as handle:
    pickle.dump(connectome_links_spherical_dictionary , handle, protocol=pickle.HIGHEST_PROTOCOL)

largest_threshold = max(metagraph_links_spherical_dictionary.keys())
link_and_threshold_dict = {}
for link in metagraph_links_spherical_dictionary[largest_threshold]:
    link_and_threshold_dict[(link[0],link[1])] = [largest_threshold]

for minimum_distance in sorted(radius_list,reverse=True):
    print(minimum_distance)
    for link in metagraph_links_spherical_dictionary[minimum_distance]:
        link_and_threshold_dict[(link[0],link[1])].append(minimum_distance)
        
        
link_and_threshold_list = []
for link in link_and_threshold_dict.keys():
    min_threshold = min(link_and_threshold_dict[link])
    link_and_threshold_list.append([link[0],link[1],min_threshold])
    
    
    
spherical_csv = pd.DataFrame(link_and_threshold_list,columns=['bodyId_1','bodyId_2','minimum_distance_threshold'])
spherical_csv.to_csv('ME(R)_synapse_only_metagraph_edgelist_min_dist_threshold_weight.csv')