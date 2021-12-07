# -*- coding: utf-8 -*-
"""
Created on Tue Nov 30 21:49:03 2021

@author: lukab
"""

from Important_functions_optimized import *
from vedo import Mesh, merge, show
import numpy as np
import pickle
import pandas as pd
import timeit

from os import listdir
from os.path import isfile, join

mypath = 'ME(R)/CSV/'
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
volumes = {}
areas = {}
bounding_boxes = {}
times_loading = {}
segments = {}
segmentwise_lengths = {}
shortest_lengths = {}
count = 0
segmentwise_areas = {}
segmentwise_volumes = {}

def kd_tree_from_point_list(point_list_1):
    X = np.array( point_list_1)
    tree = KDTree(X,compact_nodes=True,balanced_tree=True,leafsize = 25)
            # with open('tree_1', 'wb') as f:
            #     pickle.dump(tree_1, f)                

            # with open('tree_1', 'rb') as f:
            #     tree_1 = pickle.load(f)
    return tree

# Potential speed-up - load everything in the RAM memory
# Potential speed-up - paralelize the process
# Sci spatial

for bodyid in onlyfiles[0:1]:
    print('Count: ' + str(count) + ' of ' + str(len(onlyfiles)))
    neuron_properties = load_neuron_data(mypath+bodyid)
    point_list = create_neuron_point_cloud(neuron_properties,120,1)
    # tree = kd_tree_from_point_list(point_list)
    # start = timeit.default_timer()
    # with open('tree', 'wb') as f:
    #         pickle.dump(tree, f)                
    # end = timeit.default_timer()
    
    # print(end - start)
    # start = timeit.default_timer()
    # with open('tree', 'rb') as f:
    #        tree_= pickle.load(f)
    # end = timeit.default_timer()
    
    # print(end - start)   
point_cloud = Points(point_list)
show(point_cloud).close()
#point_cloud_1 = Points(point_list_1).c('black')
#for bodyid in onlyfiles[1:2]:
#     print('Count: ' + str(count+1) + ' of ' + str(len(onlyfiles)))
#     neuron_properties = load_neuron_data(mypath+bodyid)
#     point_list_2 = create_neuron_point_cloud(neuron_properties,6,1)
#     point_list_2 = point_list_2
#     point_list_label_2 = []
#     for point in point_list_2:
#         point_list_label_2.append((point[0],point[1],point[2],int(bodyid[:-4])))
#     point_cloud_2 = Points(point_list_2).c('green')


import timeit

import numpy as np
from scipy.spatial import distance
from scipy.spatial import KDTree
import matplotlib.pyplot as plt


def metagraph_edge_between_two_kd_trees(tree_1,tree_2,threshold_distance,find_true_minimum_distance=False):
    neighbors = np.array(tree_1.query_ball_tree(tree_2,threshold_distance))
    if neighbors.size == 0:
        return False
    else:
        if find_true_minimum_distance == True:
                    #This can be speed up, not needed for metagraph specifically
                    neighbor_pairs_distances = []
                    for index, x in np.ndenumerate(neighbors):
                        if len(x)>0:
                            for i in x:
                                neighbor_pairs_distances.append([index[0], i, distance.euclidean(X[index[0]],Y[i])])
                    minimum_distance = np.min([item[2] for item in neighbor_pairs_distances])
                    indices = np.argwhere(neighbor_pairs_distances == minimum_distance) 
                    coord_1 = X[neighbor_pairs_distances[int(indices[0][0])][0]]
                    coord_2 = Y[neighbor_pairs_distances[int(indices[0][0])][1]]
                    print('Minimum distance:', minimum_distance)
                    print('It is between the coordinates: ' + str(coord_1)+ ' and ' + str(coord_2))
        return True




